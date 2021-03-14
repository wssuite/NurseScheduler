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

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "solvers/mp/RosterMP.h"

using std::stringstream;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;

#define DBG_AG

//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
SubProblem::SubProblem() :
    pScenario_(nullptr),
    nDays_(0),
    pContract_(nullptr),
    pLiveNurse_(nullptr),
    pCosts_(nullptr) {}

SubProblem::SubProblem(PScenario scenario,
                       int nDays,
                       PConstContract contract,
                       vector<State> *pInitState) :
    pScenario_(scenario), nDays_(nDays), pContract_(contract) {
  int max = pScenario_->maxDuration();
  maxTotalDuration_ = max * nDays;  // working everyday on the longest shift
  if (pInitState) init(*pInitState);
}

SubProblem::SubProblem(PScenario scenario,
                       int nDays,
                       PLiveNurse pNurse,
                       const SubproblemParam& param) :
    pScenario_(std::move(scenario)), nDays_(nDays),
    pContract_(pNurse->pContract_),
    pLiveNurse_(pNurse),
    param_(param) {
  int max = 0;
  for (auto s : pScenario_->pShifts())
    if (s->duration > max) max = s->duration;
  maxTotalDuration_ = max * nDays;  // working everyday on the longest shift
}

SubProblem::~SubProblem() {}

// Initialization function
void SubProblem::init(const vector<State> &initStates) {
  // Maximum number of consecutive days worked by a nurse ending at day -1
  maxOngoingDaysWorked_ = 0;
  for (const State &state : initStates)
    maxOngoingDaysWorked_ =
        std::max(state.consDaysWorked_, maxOngoingDaysWorked_);

  nFound_ = 0;
}

//--------------------------------------------
//
// Solve function
//
//--------------------------------------------

// Solve :
// Returns TRUE if negative reduced costs path were found;
// FALSE otherwise.
bool SubProblem::solve(PLiveNurse nurse,
                       const PDualCosts &costs,
                       const SubproblemParam &param,
                       const std::set<std::pair<int, int>> &forbiddenDayShifts,
                       double redCostBound) {
  bestReducedCost_ = 0;
  nFound_ = 0;
  param_ = param;
  maxReducedCostBound_ = redCostBound;  // Cost bound
  pLiveNurse_ = nurse;  // Store the new nurse
  pCosts_ = costs;  // Store the new cost

  // Maximum rotations length: update the bounds on the nodes if needed
  updatedMaxRotationLengthOnNodes(std::min(
      nDays_ + maxOngoingDaysWorked_, std::max(pContract_->maxConsDaysWork_,
                                               param_.maxRotationLength_)));

  resetAuthorizations();  // Reset authorizations
  resetSolutions();  // Delete all already existing solutions
  initStructuresForSolve();  // Initialize structures
  forbid(forbiddenDayShifts);  // Forbid arcs
  if (!param.violateConsecutiveConstraints_)
    forbidViolationConsecutiveConstraints();

  // Set to true if you want to display contract + preferences (for debug)
//  if (false) pLiveNurse_->printContractAndPreferences(pScenario_);

  timerSolve_.start();
  // Set of operation applied on the rcGraph to prepare the solving
  timerPresolve_.start();
  preprocess();
  timerPresolve_.stop();

  bool ANS = solveRCGraph();  // Solve the RCSPP on the rcGraph

  timerPostsolve_.start();
  postprocess();
  timerPostsolve_.stop();
  timerSolve_.stop();

  // and check that all solution respects forbidden shifts
#ifdef DBG
  for (RCSolution &sol : theSolutions_) {
    checkForbiddenDaysAndShifts(sol);
  }
#endif

  return ANS;
}

bool SubProblem::preprocess() {
  // update dual costs
  updateArcDualCosts();
  return true;
}

bool SubProblem::postprocess() {
  return true;
}

// Reset all solutions data (rotations, number of solutions, etc.)
void SubProblem::resetSolutions() {
  theSolutions_.clear();
  nPaths_ = 0;
}

//--------------------------------------------
//
// Functions for the pricing of the columns
//
//--------------------------------------------

// Initializes some cost vectors that depend on the nurse
void SubProblem::initStructuresForSolve() {
  // Start and End weekend costs
  Tools::initVector(&startWeekendCosts_, nDays_, .0);
  Tools::initVector(&endWeekendCosts_, nDays_, .0);

  if (pLiveNurse_->needCompleteWeekends()) {
    for (int k = 0; k < nDays_; k++) {
      if (Tools::isSaturday(k))
        endWeekendCosts_[k] = pScenario_->weights().WEIGHT_COMPLETE_WEEKEND;
      else if (Tools::isSunday(k))
        startWeekendCosts_[k] = pScenario_->weights().WEIGHT_COMPLETE_WEEKEND;
    }
  }



  // Preference costs.
  Tools::initVector2D(&preferencesCosts_, nDays_, pScenario_->nShifts(), .0);

  for (const auto &p : pLiveNurse_->wishesOff())
    for (const Wish &s : p.second)
      preferencesCosts_[p.first][s.shift] =
          pScenario_->weights().WEIGHT_PREFERENCES_OFF[s.level];
  for (const auto &p : pLiveNurse_->wishesOn())
    for (const Wish &s : p.second)
      preferencesCosts_[p.first][s.shift] =
          pScenario_->weights().WEIGHT_PREFERENCES_ON[s.level];
}

// Forbids the nodes that correspond to forbidden shifts
void SubProblem::forbid(const set<pair<int, int> > &forbiddenDayShifts) {
  for (const pair<int, int> &p : forbiddenDayShifts)
    forbidDayShift(p.first, p.second);
}

// Authorize the nodes that correspond to forbidden shifts
void SubProblem::authorize(const set<pair<int, int> > &forbiddenDayShifts) {
  for (const pair<int, int> &p : forbiddenDayShifts)
    authorizeDayShift(p.first, p.second);
}

// Generate random forbidden shifts
set<pair<int, int> > SubProblem::randomForbiddenShifts(int nbForbidden) {
  set<pair<int, int> > ans;
  for (int f = 0; f < nbForbidden; f++) {
    int k = Tools::randomInt(0, nDays_ - 1);
    int s = Tools::randomInt(1, pScenario_->nShifts() - 1);
    ans.insert(std::pair<int, int>(k, s));
  }
  return ans;
}

//--------------------------------------------
//
// PRINT FUNCTIONS
//
//--------------------------------------------

// Prints all rotations in the current list
void SubProblem::printAllSolutions() const {
  std::cout << "# HERE ARE ALL " << nPaths_
            << " ROTATIONS OF THE CURRENT SOLUTION LIST :" << std::endl;
  for (const RCSolution &sol : theSolutions_)
    std::cout << sol.toString();
  std::cout << "# " << std::endl;
}

// Print the list of currently forbidden day and shifts
void SubProblem::printForbiddenDayShift() const {
  std::cout << "# List of currently forbidden day-shift pairs :";
  bool anyForbidden = false;
  for (int k = 0; k < nDays_; k++) {
    bool alreadyStarted = false;
    for (int s = 1; s < pScenario_->nShifts(); s++) {
      if (isDayShiftForbidden(k, s)) {
        anyForbidden = true;
        if (!alreadyStarted) {
          std::cout << std::endl << "#      | Day " << k << " :";
          alreadyStarted = true;
        }
        std::cout << " " << pScenario_->shift(s).at(0);
      }
    }
  }
  if (!anyForbidden) std::cout << " NONE";
  std::cout << std::endl;
}

void SubProblem::checkForbiddenDaysAndShifts(const RCSolution &sol) const {
  int k = sol.firstDay();
  for (const PShift &pS : sol.pShifts())
    if (isDayShiftForbidden(k++, pS->id))
      Tools::throwError(
          "A RC solution uses the forbidden day %d and shift %s: %s",
          k - 1, pS->name.c_str(), sol.toString().c_str());
}
