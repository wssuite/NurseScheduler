/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "solvers/mp/sp/RotationSP.h"

#include <memory>
#include <utility>

#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"

#ifdef BOOST
#include "solvers/mp/sp/rcspp/boost/RotationSP.h"
#endif

RotationSP::RotationSP(PScenario scenario,
                       int firstDayId,
                       int nbDays,
                       PLiveNurse nurse,
                       std::vector<PResource> pResources,
                       const SubProblemParam &param) :
    RCSPPSubProblem(std::move(scenario),
                    firstDayId,
                    nbDays,
                    std::move(nurse),
                    std::move(pResources),
                    param) {
}


void RotationSP::createNodes(const PRCGraph &pRCGraph) {
  // create a source for every day except last one (from -1 to nDays_-2)
  int initT = pLiveNurse_->pStateIni_->pShift_->type;
  const auto &pASSource = initT < pScenario_->nShiftTypes() ?
      pScenario_->pAnyTypeShift(initT) :
      pLiveNurse_->pStateIni_->pShift_;
  for (int d = -1; d < nDays_ - 1; d++)
    pRCGraph->addSingleNode(
        SOURCE_NODE, pScenario_->firstDay_.addAndGet(d),
        d >= 0 ? pScenario_->shiftsFactory().pAnyRestShift() : pASSource);

  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d <= nDays_ - 1; ++d)
    for (int st=0; st < pScenario_->nShiftTypes(); st++) {
      const auto &pAS = pScenario_->pAnyTypeShift(st);
      // create a sink node if rest or last day
      pNodesPerDay_[d][st] = pRCGraph->addSingleNode(
          (pAS->isRest() || d == nDays_ - 1) ? SINK_NODE : PRINCIPAL_NETWORK,
          pScenario_->firstDay_.addAndGet(d),
          pAS);
    }
}


void RotationSP::createArcs(const PRCGraph &pRCGraph) {
  // arcs from sources
  const PShift &pRestShift = pScenario_->pRestShift();
  const PShift &pEndShift = pScenario_->shiftsFactory().pEndShift();
  // previous shift is rest except for the first day
  PShift prevS = pLiveNurse_->pStateIni_->pShift_;
  int day = 0;
  for (const PRCNode &pSource : pRCGraph_->pSources()) {
    for (auto shiftId : prevS->successors) {
      if (shiftId == 0) continue;  // no resting arc from the source, must work
      if (pLiveNurse_->isShiftNotAvail(shiftId)) continue;
      const PShift &pS = pScenario_->pShift(shiftId);
      PRCNode pN = pNodesPerDay_[day][pS->type];
      addSingleArc(pRCGraph, pSource, pN, pS, day);
    }
    // update previous shift to rest (allow all successors) and day
    prevS = pScenario_->pRestShift();
    day++;
  }

  // arcs from the previous day to current day;
  // Can only go from work to work or work to rest
  // those from nDays-2 to nDays-1 are the arcs to the sinks
  for (int d = 1; d < nDays_; ++d)
    for (const PShift &pS : pScenario_->pShifts()) {
      PRCNode pOrigin = pNodesPerDay_[d-1][pS->type];
      if (pOrigin->type ==  SINK_NODE) continue;  // rest nodes are the sinks
      for (int succId : pS->successors) {
        if (pLiveNurse_->isShiftNotAvail(succId)) continue;
        const PShift &pS = pScenario_->pShift(succId);
        PRCNode pTarget = pNodesPerDay_[d][pS->type];
        PRCArc pArc = addSingleArc(pRCGraph, pOrigin, pTarget,
                                   pS->isRest() ? pEndShift : pS, d);
      }
    }
}

void RotationSP::createInitialLabels() {
  std::vector<PRCLabel> pLabels;
  pLabels.reserve(pRCGraph_->pSources().size());
  for (const PRCNode &pS : pRCGraph_->pSources()) {
    PRCLabel pL;
    // if initial day, use initial state
    if (pS->dayId == -1) {
      pL = std::make_shared<RCLabel>(
          pRCGraph_->pResources(), *pLiveNurse_->pStateIni_);
    } else {
      // state corresponding to the min rest shift done if any
      int nCons = pLiveNurse_->minConsDaysOff();
      State state(pScenario_->pRestShift(),
                  pS->dayId, 0, 0, 0, 0, nCons, nCons);
      pL = std::make_shared<RCLabel>(pRCGraph_->pResources(), state);
    }
    pL->setNode(pS);
    pLabels.push_back(pL);  // increment before as d starts at -1
  }
  pRcsppSolver_->setSourceLabels(pLabels);
}

bool RotationSP::postprocess() {
  const Shift &endShift = *pScenario_->shiftsFactory().pEndShift();
  for (RCSolution &sol : theSolutions_) {
    // if none shift
    if (sol.pShifts().back()->equals(endShift))
      sol.popBack();
    else if (sol.lastDayId() + 1 < nDays_)
      Tools::throwError("Rotation SP should produce solutions that end with "
                        "an end shift when not ending on the last");
  }
  return true;
}

void RotationSP::computeCost(MasterProblem *pMaster, RCSolution *rcSol) const {
  /************************************************
   * Compute all the costs of a roster:
   ************************************************/
#ifdef NS_DEBUG
  double cost = rcSol->cost();
  PRCLabel pL = rcSol->pLabel_;
#endif
  /*
  * Compute resources costs
  */
  // if previous day is the initial state of the nurse
  State state;
  if (rcSol->firstDayId() == 0) {
    state = *pLiveNurse_->pStateIni_;
  } else {
    // create a fake initial state
    int nCons = pLiveNurse_->minConsDaysOff();
    state = State(pScenario_->pRestShift(), rcSol->firstDayId() - 1,
                  0, 0, 0, 0, nCons, nCons);
  }

  bool addEndShift = rcSol->pShifts().back()->isWork() &&
      rcSol->lastDayId() < pMaster->nDays() - 1;
  // add a rest shift (if needed) at the end of the stretch to ensure
  // that all resources are priced (only if not the last day)
  if (addEndShift)
    rcSol->pushBack(pScenario_->shiftsFactory().pEndShift());
  // price resources
  computeResourcesCosts(state, rcSol);
  // remove rest shift if added
  if (addEndShift) rcSol->popBack();

#ifdef NS_DEBUG
  if (cost < DBL_MAX-1 && std::abs(cost - rcSol->cost()) > EPSILON) {
    std::cerr << "# " << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: rcspp " << cost
              << " != recomputed " << rcSol->cost()
              << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "#   | Base cost     : + " << rcSol->cost() << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    if (addEndShift)
      rcSol->pushBack(pScenario_->shiftsFactory().pEndShift());
    std::pair<vector<PResource>, PRCGraph> pActiveResources2 =
        computeResourcesCosts(state, rcSol);
    if (rcSol->pLabel_ && pL) {
      std::cerr << "Associated rcspp label: " << std::endl;
      std::cerr << pL->toStringRecursive(pResources_) << std::endl;
      std::cerr << "Associated new label: " << std::endl;
      std::cerr << rcSol->pLabel_->toStringRecursive(pResources_);
    }
    Tools::throwError("RotationSP::computeCost does not get the same cost.");

#ifdef BOOST
    // check with boost if default resources
    if (pScenario_->isWeightsDefined()) {
      boostRCSPP::RotationSP sp(pScenario_, nDays(), pLiveNurse_, param_);
      sp.computeCost(nullptr, rcSol);
    }
#endif
  }

  // restore label
  rcSol->pLabel_ = pL;
#endif
}
