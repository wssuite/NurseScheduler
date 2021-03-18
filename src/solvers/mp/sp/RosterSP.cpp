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

RosterSP::RosterSP(PScenario pScenario,
                   int nDays,
                   PLiveNurse nurse,
                   std::vector<PResource> pResources,
                   SubproblemParam param) :
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
