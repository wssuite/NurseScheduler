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
                               SubProblemParam param) :
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
  pRCGraph->addSingleNode(SINK_NODE, firstDay_, pScenario_->pRestShift());
}


void OffsetRosterSP::createArcs(const PRCGraph &pRCGraph) {
  // arcs from source to first day
  PShift pShiftIni = pLiveNurse_->pStateIni_->pShift_;
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
                 pScenario_->pRestShift(), firstDay_);
  }
}

bool OffsetRosterSP::postprocess() {
  for (RCSolution &sol : theSolutions_) {
#ifdef DBG
    if (sol.firstDay() != firstDay_)
      Tools::throwError("The RC solution and the offset grapb do not start "
                        "on the same day.");
    if (sol.pShifts().back()->isWork())
      Tools::throwError(
          "A cyclic roster cannot work on the last day before being erased.");
#endif
    sol.popBack();  // erase last rest shift
#ifdef DBG
    if (sol.pShifts().back()->isRest())
      Tools::throwError(
          "A cyclic roster cannot rest on the last day before being rotated "
          "back to a normal roster.");
#endif
    // put the last n shifts at the start of the roster
    sol.rotate(sol.firstDay());

#ifdef DBG
    if (sol.firstDay() != 0)
      Tools::throwError("A roster cannot have a first day different of 0");
    if (nDays() != sol.nDays())
      Tools::throwError("A roster cannot have a length "
                        "that is different than the number of days.");
#endif
  }
  return true;
}

CyclicRosterSP::CyclicRosterSP(PScenario pScenario,
                   int nDays,
                   PLiveNurse nurse,
                   std::vector<PResource> pResources,
                   SubProblemParam param) :
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
    const PDualCosts &costs,
    const std::set<std::pair<int, int>> &forbiddenDayShifts,
    double redCostBound) {
  // resolve SPs while ont of them produces at least solution or
  // they have all been solved
  POffsetRosterSP pSP = popOffsetFront();
  std::list<POffsetRosterSP> offsetSolved = {pSP};
  timerSolve_.start();
  while (!pSP->solve(costs, forbiddenDayShifts, redCostBound) &&
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

void CyclicRosterSP::computeCost(
    MasterProblem *pMaster, RCSolution *rcSol) const {
  // check if need to rotate
  bool rotate = rcSol->firstDay() == 0 && rcSol->pShifts().front()->isWork();
#ifdef DBG
  if (rcSol->pShifts().front()->isWork())
    Tools::throwError(
        "A cyclic roster cannot work on the first day before being rotated.");
  if (rcSol->pShifts().back()->isWork())
    Tools::throwError(
        "A cyclic roster cannot work on the last day before being erased.");
#endif

  if (rotate) {
    // Find the first rest day to rotate and put the first day on a rest.
    int fDay = 0;
    for (const PShift &pS : rcSol->pShifts()) {
      if (pS->isRest()) break;
      ++fDay;
    }
    // rotate to set the first day to fDay
    rcSol->rotate(-fDay-1);
    // add a rest shift at the end
    rcSol->pushBack(Stretch(rcSol->lastDay(), pScenario_->pRestShift()));
  }

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
    int k = rcSol->firstDay();
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

  /*
   * Compute preferencesCost
   */
  for (int k = rcSol->firstDay(); k <= rcSol->lastDay(); ++k) {
    int level =
        pLiveNurse_->wishesOffLevel(k % rcSol->nDays(), rcSol->shift(k));
    if (level != -1)
      rcSol->addCost(pScenario_->weights().WEIGHT_PREFERENCES_OFF[level],
                     PREFERENCE_COST);
    level = pLiveNurse_->wishesOnLevel(k % rcSol->nDays(), rcSol->shift(k));
    if (level != -1)
      rcSol->addCost(pScenario_->weights().WEIGHT_PREFERENCES_ON[level],
                     PREFERENCE_COST);
  }
#ifdef DBG
  if (std::abs(cost - rcSol->cost()) > EPSILON) {
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: " << rcSol->cost() << " != " << cost
              << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    Tools::throwError("RosterSP::computeCost does not get the same cost.");
  }
#endif

  if (rotate) {
#ifdef DBG
    if (rcSol->pShifts().back()->isWork())
      Tools::throwError(
          "A cyclic roster cannot work on the last day before being erased.");
#endif
    rcSol->popBack();  // erase last rest shift
#ifdef DBG
    if (rcSol->pShifts().back()->isRest())
      Tools::throwError(
          "A cyclic roster cannot rest on the last day before being rotated "
          "back to a normal roster.");
#endif
    // put the last n shifts at the start of the roster
    rcSol->rotate(rcSol->firstDay());

#ifdef DBG
    if (rcSol->firstDay() != 0)
      Tools::throwError("A roster cannot have a first day different of 0");
    if (nDays() != pMaster->nDays())
      Tools::throwError("A roster cannot have a length "
                        "that is different than the number of days.");
#endif
  }
}
