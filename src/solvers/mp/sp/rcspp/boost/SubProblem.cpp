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

#include "SubProblem.h"

#include <algorithm>
#include <map>

namespace boostRCSPP {

void SubProblem::build() {
  g_ = RCGraph(nDays_);

  //  initShortSuccessions();

  createNodes();
  createArcs();

  // Set all status to authorized
  Tools::initVector2D(&dayShiftStatus_, nDays_, pScenario_->nbShifts_, true);
  nPathsMin_ = 0;
}

BoostRCSPPSolver *SubProblem::initRCSSPSolver() {
  return new BoostRCSPPSolver(&g_,
                              maxReducedCostBound_,
                              param_.verbose_,
                              param_.epsilon_,
                              param_.search_strategy_,
                              param_.rcsppMinNegativeLabels_,
                              nullptr);
}

Penalties SubProblem::initPenalties() const {
  Penalties penalties;
  // Add each label
  // 0: CONS_DAYS
  penalties.addLabel(pContract_->minConsDaysWork_,
                     pContract_->maxConsDaysWork_,
                     pScenario_->weights().WEIGHT_CONS_DAYS_WORK);

  // if a live nurse is defined (true when solving)
  if (pLiveNurse_) {
    // 1: DAYS
    penalties.addLabel(pLiveNurse_->minTotalShifts(),
                       pLiveNurse_->maxTotalShifts(),
                       pScenario_->weights().WEIGHT_TOTAL_SHIFTS);
    // 2: WEEKEND
    penalties.addLabel(0, pLiveNurse_->maxTotalWeekends(),
                       pScenario_->weights().WEIGHT_TOTAL_WEEKENDS);
  }

  return penalties;
}

// shortest path problem is to be solved
bool SubProblem::solveRCGraph() {
  // solve the RC SPP
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> sinks = g_.sinks();

  if (param_.oneSinkNodePerLastDay_ && sinks.size() > 1) {
    sinks.resize(sinks.size() - 1);  // remove last sink (it's the main one)
    for (int a : arcsTosink_)
      g_.forbidArc(a);
  } else {
    // keep just the main one
    sinks = {sinks.back()};
  }

  // 0 - Remove all forbidden edges
  std::map<int, Arc_Properties > arcs_removed =
      g_.removeForbiddenArcsFromBoost();

  // 1 - solve the resource constraints shortest path problem
  if (sinks.empty()) sinks = g_.sinks();
  Penalties penalties = initPenalties();
  BoostRCSPPSolver *solver = initRCSSPSolver();

  std::vector<RCSolution>
      solutions = solver->solve(labels_, penalties, sinks);
  delete solver;

  // 2 - Add back all forbidden edges
  g_.restoreForbiddenArcsToBoost(arcs_removed);

  // Extract the best reduced cost
  for (const RCSolution &sol : solutions) {
    theSolutions_.push_back(sol);
    nPaths_++;
    nFound_++;
    if (bestReducedCost_ > sol.cost)
      bestReducedCost_ = sol.cost;
  }

  return !solutions.empty();
}

//--------------------------------------------
//
// Functions for the ARCS of the rcspp
//
//--------------------------------------------

// Function that creates the arcs of the network
void SubProblem::createArcs() {
  // Initialization
  Tools::initVector4D(&arcsFromSource_,
                      pScenario_->nbShiftsType_,
                      nDays_,
                      0,
                      0,
                      -1);
  Tools::initVector3D(&principalToPrincipal_,
                      pScenario_->nbShiftsType_,
                      pScenario_->nbShiftsType_,
                      nDays_,
                      -1);

  // create arcs
  createArcsSourceToPrincipal();
  createArcsPrincipalToPrincipal();
  createArcsPrincipalToSink();
}

//--------------------------------------------
//
// Functions for the costs
//
//--------------------------------------------
void SubProblem::updateArcDualCosts() {
  // A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
  for (PrincipalGraph &pg : principalGraphs_)
    for (unsigned int k = minConsDays_ - 1;
         k < arcsFromSource_[pg.shiftType()].size(); k++)
      for (int n = 0; n <= pg.maxCons(); ++n)
        for (int a : arcsFromSource_[pg.shiftType()][k][n]) {
          const Arc_Properties &arc_prop = g_.arc(a);
          if (!arc_prop.forbidden && canSuccStartHere(arc_prop) &&
              pg.checkFeasibilityEntranceArc(arc_prop, n)) {
            double c = 0;
            // if first day -> add the historical costs
            if (k == minConsDays_ - 1)
              c += historicalCost(arc_prop.shifts.front());
            // if rest shift, just add shift cost
            if (pg.shiftType() == 0)
              c += shiftCost(a);
            else
              // otherwise, call startWorkCost method
              c += startWorkCost(a);
            g_.updateCost(a, c);
            // For an arc that starts on the first day, must update the
            // consumption based on the historical state
            if (k == minConsDays_ - 1)
              g_.updateConsumptions(a, startConsumption(k, arc_prop.shifts));
          } else {
            g_.forbidArc(a);
          }
        }

  // B. ARCS : PRINCIPAL GRAPH
  for (PrincipalGraph &pg : principalGraphs_)
    pg.updateArcCosts();
}

// take into account historical state depending on current shift
// (should be called wisely)
double SubProblem::historicalCost(int currentShift) const {
  double cost = 0;
  int shiftTypeIni = pLiveNurse_->pStateIni_->shiftType_;
  int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
  int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

  // if resting
  if (pScenario_->isRestShift(currentShift)) {
    // 1. The nurse was resting
    if (shiftTypeIni == 0) {
      // if the nurse has already exceeded its max amount of rest,
      // add one penalty as current shift is already over the max
      if (pLiveNurse_->pStateIni_->consDaysOff_
          >= pLiveNurse_->maxConsDaysOff())
        cost += pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
    } else {
      // 2. The nurse was working
      // pay just penalty for min
      int diff = pLiveNurse_->minConsDaysWork() - nConsWorkIni;
      cost += std::max(.0, diff * pScenario_->weights().WEIGHT_CONS_DAYS_WORK);

      int diff2 = pScenario_->minConsShiftsOf(shiftTypeIni) - nConsShiftIni;
      cost += std::max(.0, diff2 * pScenario_->weights().WEIGHT_CONS_SHIFTS);
    }
  } else {
    // otherwise, currently working
    // 1. The nurse was resting: pay more only if the rest is too short
    if (shiftTypeIni == 0) {
      int diffRest =
          pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      cost +=
          std::max(.0, diffRest * pScenario_->weights().WEIGHT_CONS_DAYS_OFF);
    } else {
      // 2. The nurse was working
      // a. If the number of consecutive days worked has already exceeded the
      // max, subtract now the cost that will be added later
      int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
      cost -=
          std::max(.0, diffWork * pScenario_->weights().WEIGHT_CONS_DAYS_WORK);

      // b. The nurse was working on a different shift: if too short,
      // add the corresponding cost
      int shiftType = pScenario_->shiftIDToShiftTypeID_[currentShift];
      if (shiftTypeIni != shiftType) {
        int diff = pScenario_->minConsShiftsOf(shiftTypeIni) - nConsShiftIni;
        cost += std::max(.0, diff * pScenario_->weights().WEIGHT_CONS_SHIFTS);
      } else if (nConsShiftIni >= pScenario_->maxConsShiftsOf(shiftTypeIni)) {
        // c. If working on the same shift type, need to update the
        // consecutive shift cost just if exceeding the max
        cost += pScenario_->weights().WEIGHT_CONS_SHIFTS;
      }
    }
  }

  return cost;
}

std::vector<int> SubProblem::startConsumption(
    int day, std::vector<int> shifts) const {
  if (pScenario_->isRestShift(shifts.back()))
    return {0, 0, 0};

  int timeDuration = 0, size = 0;
  for (int s : shifts) {
    if (pScenario_->isRestShift(s)) {
      timeDuration = 0;
      size = 0;
    } else {
      timeDuration += pScenario_->timeDurationToWork_[s];
      ++size;
    }
  }

  // set the consumption
  std::vector<int> c = {
      size,
      timeDuration,
      Tools::nWeekendsInInterval(day - size + 1, day)
  };

  // if need to take the historical state
  if (day == size - 1 && pLiveNurse_)
    c[CONS_DAYS] += pLiveNurse_->pStateIni_->consDaysWorked_;
  return c;
}

//--------------------------------------------
//
// Functions to update the maximum length of a rotation
//
//--------------------------------------------
void SubProblem::updatedMaxRotationLengthOnNodes(int maxRotationLength) {
  if (maxRotationLength != maxRotationLength_) {
    maxRotationLength_ = maxRotationLength;
    for (int v = 0; v < g_.nodesSize(); v++) {
      std::vector<int> ubs = g_.nodeUBs(v);
      ubs[CONS_DAYS] = maxRotationLength;
      g_.updateUBs(v, ubs);
    }
  }
}

//--------------------------------------------
//
// Functions to forbid / authorize arcs and nodes
//
//--------------------------------------------


// Returns true if the succession succ starting on day k does not violate
// any forbidden day-shift
bool SubProblem::canSuccStartHere(int a) const {
  return canSuccStartHere(g_.arc(a));
}

bool SubProblem::canSuccStartHere(const Arc_Properties &arc_prop) const {
  return canSuccStartHere(arc_prop.day, arc_prop.shifts);
}

bool SubProblem::canSuccStartHere(int k, const std::vector<int> &shifts) const {
  // If the succession with the previous shift (day -1) is not allowed
  if (k == 0 &&
      pScenario_->isForbiddenSuccessorShift_Shift(
          shifts.front(), pLiveNurse_->pStateIni_->shift_))
    return false;
  // If some day-shift is forbidden
  for (int s : shifts)
    if (!dayShiftStatus_[k++][s])
      return false;
  return true;
}

// Forbids a day-shift couple
void SubProblem::forbidDayShift(int k, int s) {
  SP::forbidDayShift(k, s);

  int sh = pScenario_->shiftIDToShiftTypeID_[s];
  principalGraphs_[sh].forbidDayShift(k, s);
}

// (re)Authorizes the day-shift couple
void SubProblem::authorizeDayShift(int k, int s) {
  SP::authorizeDayShift(k, s);

  int sh = pScenario_->shiftIDToShiftTypeID_[s];
  principalGraphs_[sh].authorizeDayShift(k, s);
}

// Reset all authorizations to true
void SubProblem::resetAuthorizations() {
  SP::resetAuthorizations();

  // reset authorizations for all arcs and nodes
  g_.resetAuthorizations();
}

// forbid any arc that authorizes the violation of a consecutive constraint
void SubProblem::forbidViolationConsecutiveConstraints() {
  for (PrincipalGraph &pg : principalGraphs_)
    pg.forbidViolationConsecutiveConstraints();
}

}  // namespace boostRCSPP
