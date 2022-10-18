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
#include "solvers/mp/sp/rcspp/boost/RotationSP.h"

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
                    std::move(param)) {
}


void RotationSP::createNodes(const PRCGraph &pRCGraph) {
  // create a source for every day except last one (from -1 to nDays_-2)
  const PAbstractShift &pRestShift = pScenario_->pRestShift();
  for (int d = -1; d < nDays_ - 1; d++)
    pRCGraph->addSingleNode(
        SOURCE_NODE, pScenario_->firstDay_.addAndGet(d),
        d >= 0 ? pRestShift : pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d < nDays_ - 1; ++d)
    for (const auto& pShift : pScenario_->pShifts()) {
      PRCNode &pN = pNodesPerDay_[d][pShift->id];
      // if rest, create a sink node
      if (pShift->isRest())
        pN = pRCGraph->addSingleNode(SINK_NODE,
                                     pScenario_->firstDay_.addAndGet(d),
                                     pShift);
      else
        // else, create a work arc
        pN = pRCGraph->addSingleNode(PRINCIPAL_NETWORK,
                                     pScenario_->firstDay_.addAndGet(d),
                                     pShift);
    }

  // every shift on the last day is a sink
  for (const auto& pShift : pScenario_->pShifts())
    pNodesPerDay_[nDays_ - 1][pShift->id] =
        pRCGraph->addSingleNode(SINK_NODE,
                                pScenario_->firstDay_.addAndGet(nDays_-1),
                                pShift);
}


void RotationSP::createArcs(const PRCGraph &pRCGraph) {
  // arcs from sources
  const PShift &pRestShift = pScenario_->pRestShift();
  // previous shift is rest except for the first day
  PShift prevS = pLiveNurse_->pStateIni_->pShift_;
  int day = 0;
  for (const PRCNode &pSource : pRCGraph_->pSources()) {
    for (auto shiftId : prevS->successors) {
      if (shiftId == 0) continue;  // no resting arc from the source, must work
      if (pLiveNurse_->isShiftNotAvailNorAlt(shiftId)) continue;
      PRCNode pN = pNodesPerDay_[day][shiftId];
      addSingleArc(pRCGraph, pSource, pN, pScenario_->pShift(shiftId), day);
    }
    // update previous shift and day
    prevS = pRestShift;
    day++;
  }

  // arcs from the previous day to current day;
  // Can only go from work to work or work to rest
  // those from nDays-2 to nDays-1 are the arcs to the sinks
  for (int d = 1; d < nDays_; ++d)
    for (const PShift &pS : pScenario_->pShifts()) {
      PRCNode pOrigin = pNodesPerDay_[d-1][pS->id];
      if (pOrigin->type ==  SINK_NODE) continue;  // rest nodes are the sinks
      for (int succId : pS->successors) {
        if (pLiveNurse_->isShiftNotAvailNorAlt(succId)) continue;
        PRCNode pTarget = pNodesPerDay_[d][succId];
        addSingleArc(pRCGraph, pOrigin, pTarget,
                     pScenario_->pShift(succId), d);
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
      State state(pS->dayId, 0, 0, 0, nCons, nCons, pScenario_->pRestShift());
      pL = std::make_shared<RCLabel>(pRCGraph_->pResources(), state);
    }
    pL->setNode(pS);
    pLabels.push_back(pL);  // increment before as d starts at -1
  }
  pRcsppSolver_->setSourceLabels(pLabels);
}

bool RotationSP::postprocess() {
  for (RCSolution &sol : theSolutions_) {
    if (sol.pShifts().back()->isRest())
      sol.popBack();
    else if (sol.lastDayId() + 1 < nDays_)
      Tools::throwError("Rotation SP should produce solutions that end with "
                        "a rest shift when not ending on the last");
  }
  return true;
}

void RotationSP::computeCost(MasterProblem *pMaster, RCSolution *rcSol) const {
  /************************************************
   * Compute all the costs of a roster:
   ************************************************/
#ifdef DBG
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
    state = State(
        rcSol->firstDayId() - 1, 0, 0, 0,
        nCons, nCons, pScenario_->pRestShift());
  }
  // add a rest shift at the end of the stretch to ensure
  // that all resources are priced (only if not the last day)
  bool restShiftAdded  = false;
  if (rcSol->pShifts().back()->isWork() &&
      rcSol->lastDayId() < pMaster->nDays() - 1) {
    rcSol->pushBack(pScenario_->pRestShift());
    restShiftAdded = true;
  }
  computeResourcesCosts(state, rcSol);
  if (restShiftAdded) rcSol->popBack();

#ifdef DBG
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
    if (rcSol->pLabel_ && pL) {
      std::cerr << "Associated rcspp label: " << std::endl;
      std::cerr << pL->toStringRecursive(pResources_) << std::endl;
      std::cerr << "Associated new label: " << std::endl;
      std::cerr << rcSol->pLabel_->toStringRecursive(pResources_);
    }
    Tools::throwError("RotationSP::computeCost does not get the same cost.");

    // check with boost if default resources
    if (pMaster->useDefaultResources() && pMaster->isINRC2()) {
      boostRCSPP::RotationSP sp(pScenario_, nDays(), pLiveNurse_, param_);
      sp.computeCost(nullptr, rcSol);
    }
  }

  // restore label
  rcSol->pLabel_ = pL;
#endif
}
