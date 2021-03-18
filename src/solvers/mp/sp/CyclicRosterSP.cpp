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

#include "CyclicRosterSP.h"

#include <map>
#include <memory>
#include <set>
#include <utility>


OffsetRosterSP::OffsetRosterSP(PScenario pScenario,
                               int firstDay,
                               int nDays,
                               PLiveNurse nurse,
                               std::vector<PResource> pResources,
                               SubproblemParam param) :
    RosterSP(std::move(pScenario),
             nDays,
             std::move(nurse),
             std::move(pResources),
             std::move(param)),
    firstDay_(firstDay) {}

void OffsetRosterSP::build() {
  // set first day of the resources
  std::map<PResource, int> firstDays;
  for (const auto &pR : pResources_) {
    firstDays[pR] = pR->firstDay();
    pR->firstDay(firstDay_);
  }
  // build
  RosterSP::build();
  // reset first days
  for (const auto &p : firstDays)
    p.first->firstDay(p.second);
}

void OffsetRosterSP::createNodes(const PRCGraph &pRCGraph) {
  pRCGraph->addSingleNode(
      SOURCE_NODE, firstDay_ - 1, pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-1
  for (int d = 0; d < nDays_ ; ++d)
    for (const auto& pShift : pScenario_->pShifts()) {
      // the first day is always a rest shift
      if (d == 0 && pShift->isWork()) continue;
      // the real last day is always a work shift
      if (d == nDays_ - 1 && pShift->isRest()) continue;
      pNodesPerDay_[d][pShift->id] =
          pRCGraph->addSingleNode(PRINCIPAL_NETWORK,
                                  (d+firstDay_) % nDays_,
                                  pShift);
    }

  // create a last rest sink at the end
  pRCGraph->addSingleNode(SINK_NODE, firstDay_, pScenario_->pShift(0));
}


void OffsetRosterSP::createArcs(const PRCGraph &pRCGraph) {
  // arcs from source to first day
  PShift pShiftIni = pScenario_->pShift(pLiveNurse_->pStateIni_->shift_);
  for (auto shiftId : pShiftIni->successors) {
    if (!pLiveNurse_->isShiftAvailable(shiftId)) continue;
    PRCNode pN = pNodesPerDay_[0][shiftId];
    if (pN == nullptr)
      continue;
    addSingleArc(pRCGraph, pRCGraph->pSource(0), pN,
                 pScenario_->pShift(shiftId), pN->day);
  }

  // arcs from the previous day to current day;
  for (int d = 1; d < nDays_; ++d) {
    for (const PShift &pS : pScenario_->pShifts()) {
      PRCNode pOrigin = pNodesPerDay_[d-1][pS->id];
      // The node can be a nullptr if:
      //   - d = 1 and a work shift or
      //   - d = nDays_ - 1 and a rest shift
      if (pOrigin == nullptr)
        continue;
      for (int succId : pS->successors) {
        if (!pLiveNurse_->isShiftAvailable(succId)) continue;
        PRCNode pTarget = pNodesPerDay_[d][succId];
        if (pTarget == nullptr)
          continue;
        addSingleArc(pRCGraph, pOrigin, pTarget,
                     pScenario_->pShift(succId), pTarget->day);
      }
    }
  }

  // create the arcs from rest shifts to the sink
  for (const PShift &pS : pScenario_->pShifts()) {
    PRCNode pOrigin = pNodesPerDay_[nDays_-1][pS->id];
    if (pOrigin == nullptr)
      continue;
    addSingleArc(pRCGraph, pOrigin, pRCGraph->pSink(0),
                 pScenario_->pShift(0), firstDay_);
  }
}

CyclicRosterSP::CyclicRosterSP(PScenario pScenario,
                   int nDays,
                   PLiveNurse nurse,
                   std::vector<PResource> pResources,
                   SubproblemParam param) :
    RosterSP(std::move(pScenario),
             nDays,
             std::move(nurse),
             std::move(pResources),
             std::move(param)),
    maxOffset_(pLiveNurse_->maxConsDaysWork() + 1) {
  // check if cyclic
  if (!pScenario_->isCyclic())
    Tools::throwError("As the scenario is not cyclic, "
                      "it cannot use CyclicRosterSP.");
  // check if the maxOffset is not too big
  if (2*maxOffset_ >= nDays_)
    std::cout << "WARNING: there are a lot of subproblems (" << maxOffset_
              << ") for the cyclic roster. You should set a smaller value "
                 "for the upper bound on the maximum number of consecutive "
                 "worked days." << std::endl;

  // create the subgraphs
  int firstDay = pLiveNurse_->num_ % maxOffset_;
  for (int  i = 0; i < maxOffset_; i++, firstDay++) {
    if (firstDay >= maxOffset_) firstDay -= maxOffset_;
    offsetSPs_.emplace_back(std::make_shared<OffsetRosterSP>(
        pScenario_, firstDay, pRCGraph_->nDays(),
        pLiveNurse_, pResources_, param_));
  }
}

void CyclicRosterSP::build() {
  for (const auto &pSP : offsetSPs_)
    pSP->build();
}

// Solve : Returns TRUE if negative reduced costs path were found;
// FALSE otherwise.
bool CyclicRosterSP::solve(
    PLiveNurse nurse,
    const PDualCosts &costs,
    const SubproblemParam &param,
    const std::set<std::pair<int, int>> &forbiddenDayShifts,
    double redCostBound) {
  // resolve SPs while ont of them produces at least solution or
  // they have all been solved
  POffsetRosterSP pSP = popOffsetFront();
  std::list<POffsetRosterSP> offsetSolved = {pSP};
  timerSolve_.start();
  while (!pSP->solve(nurse, costs, param, forbiddenDayShifts, redCostBound) &&
      !offsetSPs_.empty()) {
    // update offsets
    pSP = popOffsetFront();
    offsetSolved.push_back(pSP);
  }
  timerSolve_.stop();

  // update offsets
  addOffsetSPs(offsetSolved);

  // update solution
  bestReducedCost_ = pSP->bestReducedCost();
  theSolutions_ = pSP->getSolutions();

  return !theSolutions_.empty();
}

POffsetRosterSP CyclicRosterSP::popOffsetFront() {
  POffsetRosterSP pSP = offsetSPs_.front();
  offsetSPs_.pop_front();
  return pSP;
}
