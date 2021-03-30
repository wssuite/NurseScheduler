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
                   int nDays,
                   PLiveNurse nurse,
                   std::vector<PResource> pResources,
                   SubProblemParam param) :
    RCSPPSubProblem(std::move(pScenario),
                    nDays,
                    std::move(nurse),
                    std::move(pResources),
                    std::move(param)) {}


void RosterSP::createNodes(const PRCGraph &pRCGraph) {
  pRCGraph->addSingleNode(SOURCE_NODE, -1, pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d < nDays_ - 1 ; ++d)
    for (const auto& pShift : pScenario_->pShifts())
      pNodesPerDay_[d][pShift->id] =
          pRCGraph->addSingleNode(PRINCIPAL_NETWORK, d, pShift);

  // every shift on the last day is a sink
  for (const auto &pShift : pScenario_->pShifts())
    pNodesPerDay_[nDays_ - 1][pShift->id] =
        pRCGraph->addSingleNode(SINK_NODE, nDays_ - 1, pShift);
}


void RosterSP::createArcs(const PRCGraph &pRCGraph) {
  // arcs from source to first day
  PShift pShiftIni = pScenario_->pShift(pLiveNurse_->pStateIni_->shift_);
  for (auto shiftId : pShiftIni->successors) {
    if (!pLiveNurse_->isShiftAvailable(shiftId)) continue;
    PRCNode pN = pNodesPerDay_[0][shiftId];
    addSingleArc(
        pRCGraph, pRCGraph->pSource(0), pN, pScenario_->pShift(shiftId), 0);
  }

  // arcs from the previous day to current day;
  // those from nDays-2 to nDays-1 are the arcs to the sinks
  for (int d = 1; d < nDays_; ++d) {
    for (const PShift &pS : pScenario_->pShifts()) {
      PRCNode pOrigin = pNodesPerDay_[d-1][pS->id];
      for (int succId : pS->successors) {
        if (!pLiveNurse_->isShiftAvailable(succId)) continue;
        PRCNode pTarget = pNodesPerDay_[d][succId];
        addSingleArc(pRCGraph, pOrigin, pTarget, pScenario_->pShift(succId), d);
      }
    }
  }
}

double RosterSP::dualCost(const PRCArc &pArc) {
  double dualCost = 0;
  // if start, add constant
  if (pArc->origin->type == SOURCE_NODE)
    dualCost -= pCosts_->constant();
  // iterate through the shift to update the cost
  int curDay = pArc->stretch.firstDay();
  for (const auto& pS : pArc->stretch.pShifts()) {
    if (pS->isWork())
      dualCost -= pCosts_->workedDayShiftCost(curDay, pS->id);
    curDay++;
  }
  return dualCost;
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
#endif
  /*
  * Compute resources costs
  */
  computeResourcesCosts(*pLiveNurse_->pStateIni_, pMaster, rcSol);

  /*
   * Compute complete weekend
   */
  if (pLiveNurse_->needCompleteWeekends()) {
    int k = 0;
    bool rest = false;
    for (const PShift &pS : rcSol->pShifts()) {
      // on sunday, if complete weekend, it's either:
      // work on saturday (rest=false) and sunday
      // rest on saturday (rest=true) and sunday
      if (Tools::isSunday(k) && (rest ^ pS->isRest()))
        rcSol->addCost(pScenario_->weights().WEIGHT_COMPLETE_WEEKEND,
                       COMPLETE_WEEKEND_COST);
      rest = pS->isRest();
      k++;
    }
  }

  computePreferencesCost(rcSol);

#ifdef DBG
  if (std::abs(cost - rcSol->cost()) > 1e-3) {
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: " << rcSol->cost() << " != " << cost
              << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    Tools::throwError("RosterSP::computeCost does not get the same cost.");
  }

  // check with boost if default resources
  if (pMaster->useDefaultResources()) {
    boostRCSPP::RosterSP sp(pScenario_, nDays(), pLiveNurse_, param_);
    sp.computeCost(nullptr, rcSol);
  }
#endif
}
