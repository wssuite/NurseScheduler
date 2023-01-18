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
#include <utility>
#include <vector>

#include "solvers/mp/RosterMP.h"

using std::stringstream;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;


//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Dual Costs class
// update the dual values of every constraints based on the current solution
void DualCosts::updateDuals() {
  for (ConstraintMP* pC : pMaster_->columnConstraints())
    pC->updateDuals();
}

// update the dual values of every constraints randomly
void DualCosts::randomUpdateDuals(bool useInputData, int nPerturbations) {
  for (ConstraintMP* pC : pMaster_->columnConstraints())
    pC->randomUpdateDuals(useInputData, nPerturbations);
}

// return the dual cost of a stretch based on its consumption of
// every constraints
double DualCosts::getCost(
    int nurseNum,
    const Stretch &st,
    const PAbstractShift &prevS) const {
  double d = 0;
  for (ConstraintMP* pC : pMaster_->columnConstraints())
    d += pC->getDualCost(nurseNum, st, prevS);
  return d;
}

vector<double> DualCosts::getMaxDualValues() const {
  vector<double> maxDuals(pMaster_->nNurses());
  for (int n=0; n < maxDuals.size(); n++) {
    double d = 0, maxD = -XLARGE_SCORE;
    for (ConstraintMP* pC : pMaster_->columnConstraints()) {
      double v = pC->maxDualValue(n);
      d += v;
      if (maxD < v) maxD = v;
    }
    maxDuals[n] = d;
  }
  return maxDuals;
}

std::string DualCosts::toString() const {
  std::stringstream buff;
  for (ConstraintMP* pC : pMaster_->columnConstraints())
    buff << pC->toString();
  return buff.str();
}

std::string DualCosts::toString(int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  for (ConstraintMP* pC : pMaster_->columnConstraints())
    buff << pC->toString(nurseNum, st);
  return buff.str();
}

// Constructors and destructor
SubProblem::SubProblem() :
    pScenario_(nullptr),
    firstDayId_(0),
    nDays_(0),
    pLiveNurse_(nullptr),
    pCosts_(nullptr),
    rdm_(Tools::getANewRandomGenerator()),
    timerPresolve_("SP pre-solve"),
    timerSolve_("SP solve"),
    timerPostsolve_("SP post-solve") {}

SubProblem::SubProblem(PScenario scenario,
                       int firstDayId,
                       int nDays,
                       PLiveNurse pNurse,
                       const SubProblemParam& param) :
    pScenario_(std::move(scenario)),
    firstDayId_(firstDayId), nDays_(nDays),
    pLiveNurse_(std::move(pNurse)),
    param_(param),
    rdm_(Tools::getANewRandomGenerator()),
    timerPresolve_("SP pre-solve"),
    timerSolve_("SP solve"),
    timerPostsolve_("SP post-solve") {
  // working everyday on the longest shift
  maxTotalDuration_ = pScenario_->maxDuration() * nDays;
  init(*pScenario_->pInitialState());
}

SubProblem::~SubProblem() = default;

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
bool SubProblem::solve(const PDualCosts &costs,
                       const std::set<std::pair<int, int>> &forbiddenDayShifts,
                       double redCostBound,
                       bool relaxation) {
  // if not a new run, resolve now
  if (pCosts_ == costs) {
    if (param_.verbose_ >= 2)
      std::cout << "Resolving subproblem for nurse " << pLiveNurse_->name_
                << " (" << pLiveNurse_->num_ << ")" << std::endl;
    return solve(false, relaxation);
  }

  if (param_.verbose_ >= 2)
    std::cout << "Solving subproblem for nurse " << pLiveNurse_->name_
              << " (" << pLiveNurse_->num_ << ")" << std::endl;

  bestReducedCost_ = 0;
  nFound_ = 0;
  maxReducedCostBound_ = redCostBound;  // Cost bound
  pCosts_ = costs;  // Store the new cost

  resetAuthorizations();  // Reset authorizations
  resetSolutions();  // Delete all already existing solutions
  initStructuresForSolve();  // Initialize structures
  forbid(forbiddenDayShifts);  // Forbid arcs

  // Set to true if you want to display contract + preferences (for debug)
//  if (false) pLiveNurse_->printContractAndPreferences(pScenario_);

  return solve(true, relaxation);
}

bool SubProblem::solve(bool initialSolve, bool relaxation) {
  timerSolve_.start();
  if (initialSolve) {
    // Set of operation applied on the rcGraph to prepare the solving
    timerPresolve_.start();
    presolve();
    timerPresolve_.stop();
  }

  // Solve the RCSPP on the rcGraph
  bool ANS = solveRCGraph(initialSolve, relaxation);

  timerPostsolve_.start();
  postprocess();
  timerPostsolve_.stop();
  timerSolve_.stop();

  // and check that all solution respects forbidden shifts
#ifdef NS_DEBUG
  for (RCSolution &sol : theSolutions_) {
    checkForbiddenDaysAndShifts(sol);
  }
#endif

  return ANS;
}

bool SubProblem::presolve() {
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
        std::cout << " " << pScenario_->shiftName(s).at(0);
      }
    }
  }
  if (!anyForbidden) std::cout << " NONE";
  std::cout << std::endl;
}

void SubProblem::checkForbiddenDaysAndShifts(const RCSolution &sol) const {
  int k = sol.firstDayId();
  for (const PShift &pS : sol.pShifts())
    if (isDayShiftForbidden(k++, pS->id))
      Tools::throwError(
          "A RC solution uses the forbidden day %d and shift %s: %s",
          k - 1, pS->name.c_str(), sol.toString().c_str());
}

void SubProblem::computePreferencesCost(RCSolution *rcSol) const {
  /*
 * Compute preferencesCost
 */
  for (int k = rcSol->firstDayId(); k <= rcSol->lastDayId(); ++k) {
    double cost = pLiveNurse_->wishCostOfTheShift(k, rcSol->pShift(k));
    rcSol->addCost(cost, PREFERENCE_COST);
  }
}
