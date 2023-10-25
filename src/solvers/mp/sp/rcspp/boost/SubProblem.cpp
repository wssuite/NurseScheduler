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
#include <memory>
#include <utility>

#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

namespace boostRCSPP {

const int SubProblem::maxSubproblemStrategyLevel_ = 3;

void SubProblem::build() {
  g_ = RCGraph(nDays_);

  //  initShortSuccessions();

  createNodes();
  createArcs();

  // Set all status to authorized
  Tools::initVector2D(&dayShiftStatus_, nDays_, pScenario_->nShifts(), true);
  nPathsMin_ = 0;
}


void SubProblem::initSubproblemParam(int strategy) {
  maxRotationLength_ = pLiveNurse_->maxConsDaysWork();
  param_.rcsppMinNegativeLabels_ = static_cast<int>(round(
      param_.spMinColumnsRatioForIncrease_ * param_.nbMaxColumnsToAdd_));
  int nb_max_path = static_cast<int>(round(
      param_.spColumnsRatioForNumberPaths_ * param_.nbMaxColumnsToAdd_));
  param_.strategyLevel_ = strategy;
  switch (strategy) {
    // 0 -> [HeuristicMIP large search]
    //  short = all,
    //  min   = 0,
    //  max   = CD_max+3
    //
    case 0:param_.search_strategy_ = SP_BEST_FIRST;
      param_.rcsppMaxNegativeLabels_ = nb_max_path;
      param_.violateConsecutiveConstraints_ = true;
      param_.shortRotationsStrategy_ = 3;
      maxRotationLength_ += 3;
      break;

      // 1 -> [Exact legal only]
      //  short = first and last day,
      //  min   = CD_min,
      //  max   = CD_max
      //
    case 1:param_.search_strategy_ = SP_BREADTH_FIRST;
      param_.rcsppMaxNegativeLabels_ = -1;
      param_.violateConsecutiveConstraints_ = false;
      param_.shortRotationsStrategy_ = 2;
      maxRotationLength_ += 0;
      break;

      // 2 -> [Exact above legal]
      //  short = all,
      //  min   = 0,
      //  max   = CD_max+2
      //
    case 2:param_.search_strategy_ = SP_BREADTH_FIRST;
      param_.rcsppMaxNegativeLabels_ = -1;
      param_.violateConsecutiveConstraints_ = true;
      param_.shortRotationsStrategy_ = 3;
      maxRotationLength_ += 2;
      break;

      // 3 -> [Exact exhaustive search]
      //  short = all,
      //  min   = 0,
      //  max   = LARGE
      //
    case 3:param_.search_strategy_ = SP_BREADTH_FIRST;
      param_.rcsppMaxNegativeLabels_ = -1;
      param_.violateConsecutiveConstraints_ = true;
      param_.shortRotationsStrategy_ = 3;
      maxRotationLength_ =
          pLiveNurse_->pStateIni_->consDaysWorked_ + pLiveNurse_->nbDays_;
      break;
    default:  // UNKNOWN STRATEGY
      std::cout << "# Unknown strategy for the subproblem (" << strategy << ")"
                << std::endl;
      break;
  }
}

void SubProblem::updateParameters(bool useMoreTime) {
  // if backtracking -> restart at 0
  int diff = param_.strategyLevel_ - defaultStrategy_;
  if (useMoreTime || diff <= 1) initSubproblemParam(defaultStrategy_);
  else
    // otherwise, just try previous level if more than 2 levels
    initSubproblemParam(param_.strategyLevel_ - 1);
}

BoostRCSPPSolver *SubProblem::initRCSSPSolver() {
  return new BoostRCSPPSolver(&g_,
                              maxReducedCostBound_,
                              param_.verbose_,
                              param_.epsilon_,
                              param_.search_strategy_,
                              param_.rcsppMaxNegativeLabels_,
                              nullptr);
}

Penalties SubProblem::initPenalties() const {
  Penalties penalties;
  // Add each label
  // 0: CONS_DAYS
  penalties.addLabel(pLiveNurse()->minConsDaysWork(),
                     pLiveNurse()->maxConsDaysWork(),
                     pScenario_->weights().consDaysWork);

  // if a live nurse is defined (true when solving)
  if (pLiveNurse_) {
    // 1: DAYS
    penalties.addLabel(pLiveNurse_->totalShiftDurationResource_->getLb(),
                       pLiveNurse_->totalShiftDurationResource_->getUb(),
                       pLiveNurse_->totalShiftDurationResource_->getUbCost());
    // 2: WEEKEND
    penalties.addLabel(pLiveNurse_->totalWeekendResource_->getLb(),
                       pLiveNurse_->totalWeekendResource_->getUb(),
                       pLiveNurse_->totalWeekendResource_->getUbCost());
  }

  return penalties;
}

bool SubProblem::solve(bool initialSolve, bool relaxation) {
  bool ANS = false;
  vector<int> forbiddenArcs;
  while (!ANS) {
    // Maximum rotations length: update the bounds on the nodes if needed
    updatedMaxRotationLengthOnNodes(
        std::min(nDays_ + maxOngoingDaysWorked_,
            std::max(pLiveNurse()->maxConsDaysWork(), maxRotationLength_)));

    for (int a : forbiddenArcs) g().authorizeArc(a);
    forbiddenArcs.clear();

    if (!param_.violateConsecutiveConstraints_)
      forbiddenArcs = forbidViolationConsecutiveConstraints();

    ANS = SP::solve();

    if (param_.strategyLevel_ == maxSubproblemStrategyLevel_)
      break;
    // update strategy
    if (theSolutions_.size() < param_.rcsppMinNegativeLabels_)
      initSubproblemParam(param_.strategyLevel_+1);
  }
  return ANS;
}

// shortest path problem is to be solved
bool SubProblem::solveRCGraph(bool initialSolve, bool relaxation) {
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
  auto solver = std::unique_ptr<BoostRCSPPSolver>(initRCSSPSolver());

  std::vector<RCSolution>
      solutions = solver->solve(labels_, penalties, sinks, pScenario_);

  // 2 - Add back all forbidden edges
  g_.restoreForbiddenArcsToBoost(arcs_removed);

  // Extract the best reduced cost
  for (RCSolution &sol : solutions) {
    computeCost(nullptr, &sol);
    theSolutions_.push_back(sol);
    nPaths_++;
    nFound_++;
    if (sol.reducedCost() < bestReducedCost_)
      bestReducedCost_ = sol.reducedCost();
  }

  return !theSolutions_.empty();
}

// Initializes some cost vectors that depend on the nurse
void SubProblem::initStructuresForSolve() {
  // Start and End weekend costs
  Tools::initVector(&startWeekendCosts_, nDays_, .0);
  Tools::initVector(&endWeekendCosts_, nDays_, .0);

  if (pLiveNurse_->needCompleteWeekends()) {
    for (int k = 0; k < nDays_; k++) {
      if (pLiveNurse_->pContract_->isFirstWeekendDay(k))
        endWeekendCosts_[k] = pScenario_->weights().completeWeekend;
      else if (pLiveNurse_->pContract_->isLastWeekendDay(k))
        startWeekendCosts_[k] = pScenario_->weights().completeWeekend;
    }
  }

  // Preference costs.
  Tools::initVector2D(&preferencesCosts_, nDays_, pScenario_->nShifts(), .0);

  for (const auto &wish : pLiveNurse_->wishes())
    for (const PShift &pS : pScenario_->pShifts())
      preferencesCosts_[wish.first][pS->id] += wish.second.cost(pS);
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
                      pScenario_->nShiftTypes(),
                      nDays_,
                      0,
                      0,
                      -1);
  Tools::initVector3D(&principalToPrincipal_,
                      pScenario_->nShiftTypes(),
                      pScenario_->nShiftTypes(),
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
              c += historicalCost(arc_prop.pShifts.front()->id);
            // if rest shift, just add shift cost
            if (pg.shiftType() == 0)
              c += shiftCost(a, nullptr);
            else
              // otherwise, call startWorkCost method
              c += startWorkCost(a);
            g_.updateCost(a, c);
            // For an arc that starts on the first day, must update the
            // consumption based on the historical state
            if (k == minConsDays_ - 1)
              g_.updateConsumptions(a, startConsumption(k, arc_prop.pShifts));
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
  int shiftTypeIni = pLiveNurse_->pStateIni_->pShift_->type;
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
        cost += pScenario_->weights().consDaysOff;
    } else {
      // 2. The nurse was working
      // pay just penalty for min
      int diff = pLiveNurse_->minConsDaysWork() - nConsWorkIni;
      cost += std::max(.0, diff * pScenario_->weights().consDaysWork);

      int diff2 = pScenario_->minConsShiftsOfType(shiftTypeIni) - nConsShiftIni;
      cost += std::max(.0, diff2 * pScenario_->weights().consShifts);
    }
  } else {
    // otherwise, currently working
    // 1. The nurse was resting: pay more only if the rest is too short
    if (shiftTypeIni == 0) {
      int diffRest =
          pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      cost +=
          std::max(.0, diffRest * pScenario_->weights().consDaysOff);
    } else {
      // 2. The nurse was working
      // a. If the number of consecutive days worked has already exceeded the
      // max, subtract now the cost that will be added later
      int diffWork = nConsWorkIni - pLiveNurse()->maxConsDaysWork();
      cost -=
          std::max(.0, diffWork * pScenario_->weights().consDaysWork);

      // b. The nurse was working on a different shift: if too short,
      // add the corresponding cost
      int shiftType = pScenario_->shiftIDToShiftTypeID(currentShift);
      if (shiftTypeIni != shiftType) {
        int diff =
            pScenario_->minConsShiftsOfType(shiftTypeIni) - nConsShiftIni;
        cost += std::max(.0, diff * pScenario_->weights().consShifts);
      } else if (nConsShiftIni >=
          pScenario_->maxConsShiftsOfType(shiftTypeIni)) {
        // c. If working on the same shift type, need to update the
        // consecutive shift cost just if exceeding the max
        cost += pScenario_->weights().consShifts;
      }
    }
  }

  return cost;
}

std::vector<int> SubProblem::startConsumption(
    int day, const std::vector<PShift> &pShifts) const {
  if (pShifts.back()->isRest())
    return {0, 0, 0};

  int timeDuration = 0, size = 0;
  for (const PShift &pS : pShifts) {
    if (pS->isRest()) {
      timeDuration = 0;
      size = 0;
    } else {
      timeDuration += pS->duration;
      ++size;
    }
  }

  // set the consumption
  Weekend we;
  std::vector<int> c = {
      size,
      timeDuration,
      we.nWeekendsInInterval(day - size + 1, day)
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
  return canSuccStartHere(arc_prop.day, arc_prop.pShifts);
}

bool SubProblem::canSuccStartHere(
    int k, const std::vector<PShift> &shifts) const {
  // If the succession with the previous shift (day -1) is not allowed
  if (k == 0 && !shifts.front()->canSucceed(*pLiveNurse_->pStateIni_->pShift_))
    return false;
  // If some day-shift is forbidden
  for (const PShift &pS : shifts)
    if (!dayShiftStatus_[k++][pS->id])
      return false;
  return true;
}

// Forbids a day-shift couple
void SubProblem::forbidDayShift(int k, int s) {
  SP::forbidDayShift(k, s);

  int sh = pScenario_->shiftIDToShiftTypeID(s);
  principalGraphs_[sh].forbidDayShift(k, s);
}

// (re)Authorizes the day-shift couple
void SubProblem::authorizeDayShift(int k, int s) {
  SP::authorizeDayShift(k, s);

  int sh = pScenario_->shiftIDToShiftTypeID(s);
  principalGraphs_[sh].authorizeDayShift(k, s);
}

// Reset all authorizations to true
void SubProblem::resetAuthorizations() {
  SP::resetAuthorizations();

  // reset authorizations for all arcs and nodes
  g_.resetAuthorizations();
}

// forbid any arc that authorizes the violation of a consecutive constraint
vector<int> SubProblem::forbidViolationConsecutiveConstraints() {
  vector<int> forbiddenArcs;
  for (PrincipalGraph &pg : principalGraphs_) {
    vector<int> fArcs = pg.forbidViolationConsecutiveConstraints();
    forbiddenArcs = Tools::appendVectors(forbiddenArcs, fArcs);
  }
  return forbiddenArcs;
}

}  // namespace boostRCSPP
