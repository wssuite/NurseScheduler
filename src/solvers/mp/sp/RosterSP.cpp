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

#include "solvers/mp/sp/RosterSP.h"

#include <list>
#include <memory>
#include <utility>

#include "solvers/mp/sp/rcspp/boost/RosterSP.h"

RosterSP::RosterSP(PScenario pScenario,
                   int firstDayId,
                   int nDays,
                   PLiveNurse nurse,
                   std::vector<PResource> pResources,
                   const SubProblemParam& param) :
    RCSPPSubProblem(std::move(pScenario),
                    firstDayId,
                    nDays,
                    std::move(nurse),
                    std::move(pResources),
                    param) {}


void RosterSP::createNodes(const PRCGraph &pRCGraph) {
  auto pDay = pScenario_->firstDay_.previous();
  pRCGraph->addSingleNode(SOURCE_NODE, pDay,
                          pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d < nDays_ - 1 ; ++d) {
    pDay = pDay->next();
    for (const auto &pShift : pScenario_->pShifts())
      pNodesPerDay_[d][pShift->id] =
          pRCGraph->addSingleNode(PRINCIPAL_NETWORK, pDay, pShift);
  }

  // every shift on the last day is a sink
  pDay = pDay->next();
  for (const auto &pShift : pScenario_->pShifts())
    pNodesPerDay_[nDays_ - 1][pShift->id] =
        pRCGraph->addSingleNode(SINK_NODE, pDay, pShift);
}


void RosterSP::createArcs(const PRCGraph &pRCGraph) {
  if (!pLiveNurse_->pStateIni_->pShift_)
    Tools::throwError("Nurse %s does not have any initial shift",
                      pLiveNurse_->name_.c_str());

  // arcs from source to first day
  for (auto shiftId : pLiveNurse_->pStateIni_->pShift_->successors) {
    if (pLiveNurse_->isShiftNotAvailNorAlt(shiftId)) continue;
    const PRCNode &pN = pNodesPerDay_[0][shiftId];
    addSingleArc(
        pRCGraph, pRCGraph->pSource(0), pN, pScenario_->pShift(shiftId), 0);
  }

  // arcs from the previous day to current day;
  // those from nDays-2 to nDays-1 are the arcs to the sinks
  for (int d = 1; d < nDays_; ++d) {
    for (const PShift &pS : pScenario_->pShifts()) {
      const PRCNode &pOrigin = pNodesPerDay_[d-1][pS->id];
      for (int succId : pS->successors) {
        if (pLiveNurse_->isShiftNotAvailNorAlt(succId)) continue;
        const PRCNode &pTarget = pNodesPerDay_[d][succId];
        addSingleArc(pRCGraph, pOrigin, pTarget, pScenario_->pShift(succId), d);
      }
    }
  }
}

void RosterSP::createInitialLabels() {
  PRCLabel pL = std::make_shared<RCLabel>(pRCGraph_->pResources(),
                                          *pLiveNurse_->pStateIni_);
  pL->setNode(pRCGraph_->pSource(0));
  pRcsppSolver_->setSourceLabels({pL});
}

void RosterSP::computeCost(MasterProblem *pMaster, RCSolution *rcSol) const {
  /************************************************
   * Compute all the costs of a roster:
   ************************************************/
#ifdef DBG
  double cost = rcSol->cost();
  PRCLabel pL = rcSol->pLabel_;
#endif
  /*
  * Compute resources costs: this includes all the individual soft costs of the
   * nurse, including preferences, complete weekends, forbidden patterns, etc.
  */
  vector<PResource> pActiveResources =
      computeResourcesCosts(*pLiveNurse_->pStateIni_, rcSol);

#ifdef DBG
  if (cost < DBL_MAX - 1 && std::abs(cost - rcSol->cost()) > EPSILON) {
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: rcspp " << cost
              << " != recomputed " << rcSol->cost()
              << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    vector<PResource> pActiveResources2 =
        computeResourcesCosts(*pLiveNurse_->pStateIni_, rcSol);
    if (rcSol->pLabel_ && pL) {
      std::cerr << "Associated rcspp label: " << std::endl;
      std::cerr << pL->toStringRecursive(pActiveResources) << std::endl;
      std::cerr << "Associated new label: " << std::endl;
      std::cerr << rcSol->pLabel_->toStringRecursive(pActiveResources2);
    }
    Tools::throwError("RosterSP::computeCost does not get the same cost.");
  }

  // check with boost if default resources
  if (pMaster->useDefaultResources() && pMaster->isINRC2() &&
      pMaster->getDynamicWeights().version() > 0) {
    boostRCSPP::RosterSP sp(pScenario_, nDays(), pLiveNurse_, param_);
    sp.computeCost(nullptr, rcSol);
  }

  // restore label
  rcSol->pLabel_ = pL;
#endif
}
