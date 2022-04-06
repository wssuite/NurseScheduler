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

#include "solvers/mp/RosterMP.h"

#include <map>
#include <memory>
#include <utility>

#include "solvers/mp/TreeManager.h"
#include "solvers/mp/sp/SubProblem.h"


//-----------------------------------------------------------------------------
//
//  S t r u c t   R o s t e r
//
//  A roster is a vector of shifts covering the whole horizon
//
//-----------------------------------------------------------------------------

// when branching on this pattern, this method add the corresponding forbidden
// shifts to the set.
// It will forbid any shifts on any days as the nurse already has a roster.
void RosterPattern::addForbiddenShifts(
    std::set<std::pair<int, int> > *forbidenShifts,
    int nbShifts,
    PDemand pDemand) const {
  // from the previous day to the day after the end of the rotation, forbid any
  // work shifts
  for (int day = firstDay(); day <= lastDay(); day++)
    for (int i = 1; i < nbShifts; ++i)
      forbidenShifts->insert(std::pair<int, int>(day, i));
}

void RosterPattern::checkReducedCost(const DualCosts &dualCosts,
                                     bool printBadPricing) const {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of the roster:
   ************************************************/
  double reducedCost = cost_ - dualCosts.getCost(nurseNum_, *this, nullptr);

  // Display: set to true if you want to display the details of the cost
  if (std::fabs(reducedCost_ - reducedCost) / (1 - reducedCost) > EPSILON) {
    // if do not print and not throwing an error
    if (!printBadPricing && reducedCost_ > reducedCost + EPSILON) return;

    std::cerr << "# " << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "Bad dual cost: " << reducedCost_ << " != " << reducedCost
              << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "#   | Base cost     : + " << cost_ << std::endl;
    std::cerr << costsToString();

    std::cerr << dualCosts.toString(nurseNum_, *this);
    std::cerr << toString();
    std::cerr << "# " << std::endl;

    // throw an error only when a significant misprice
    // Indeed, if the real reduced cost (reducedCost) is greater than the one
    // found by the pricing, some columns with a positive reduced cost could be
    // generated.
    // The reason why the other situation can arise is that some path in the
    // subproblem could under estimate the real cost. These paths won't be
    // found when the subproblems are solved at optimality, but could  be
    // present when using heuristics.
//    if (reducedCost_ < reducedCost + EPSILON)
      Tools::throwError("Invalid pricing of a roster.");
  }
}

//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r M P
//
// Build and solve the master problem of the column generation scheme with
// rosters
//
//-----------------------------------------------------------------------------

RosterMP::RosterMP(const PScenario& pScenario,
                   SolverType solver) :
    MasterProblem(pScenario, solver) {
  lagrangianBoundAvail_ = true;
}

RosterMP::~RosterMP() {
  delete assignmentConstraint_;
}

PPattern RosterMP::getPattern(MyVar *var) const {
  return std::make_shared<RosterPattern>(var, pScenario_);
}

// Main method to build the rostering problem for a given input
void RosterMP::build(const SolverParam &param) {
  /* Roster assignment constraints */
  assignmentConstraint_ = new RosterAssignmentConstraint(this);

  /* build the rest of the model */
  MasterProblem::build(param);

  /* Change the branching rule */
  auto *pRule = new RosterBranchingRule(this,
                                   dynamic_cast<RestTree *>(pTree()),
                                   "branching rule");
  pModel_->addBranchingRule(pRule);
}

// Provide an initial solution to the solver. If empty, add artificial columns
void RosterMP::initializeSolution(const std::vector<Roster> &solution) {
  // rosters are added for each nurse of the initial solution
  if (!solution.empty()) {
    const char *baseName("initialRoster");
    // build the roster of each nurse
    for (int i = 0; i < pScenario_->nNurses(); ++i) {
      RosterPattern pat(solution[i].pShifts(), i);
      computePatternCost(&pat);
      pModel_->addInitialColumn(createColumn(pat, baseName));
    }
  }
}

// Create a new rotation variable
// add the correct constraints and coefficients for the nurse i working
// on a rotation
// if s=-1, the nurse works on all shifts
MyVar *RosterMP::addColumn(int nurseNum, const RCSolution &solution) {
  // Build rotation from RCSolution
  RosterPattern pat(solution, nurseNum);
  pat.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  computePatternCost(&pat);
  DualCosts dualCosts(this);
  pat.checkReducedCost(dualCosts, pPricer()->isLastRunOptimal());
  checkIfPatternAlreadyPresent(pat);
#endif
  return createColumn(pat, "roster");
}

// get a reference to the restsPerDay_ for a Nurse
std::vector<MyVar *> RosterMP::getRestVarsPerDay(PLiveNurse pNurse,
                                                 int day) const {
  std::vector<MyVar *> restRosters;
  for (MyVar *var : pModel_->getActiveColumns()) {
    if (Pattern::nurseNum(var) == pNurse->num_
        && pScenario_->isRestShift(Pattern::shift(var, day)))
      restRosters.push_back(var);
  }
  return restRosters;
}

// compute the lagrangian bound
double RosterMP::computeLagrangianBound(double objVal) const {
  if (!pModel_->stab().stabCheckStoppingCriterion()) {
    std::cerr << "Cannot compute a lagrangian bound when stabilization "
                 "variables are present in the solution." << std::endl;
    return -LARGE_SCORE;
  }

  double sumRedCost = 0;
  for (double v : pPricer()->getLastMinReducedCosts())
    sumRedCost += v;
  return objVal + sumRedCost;
}

double RosterMP::getDaysCost() const {
  return getColumnsCost(TOTAL_WORK_COST);
}

double RosterMP::getWeekendCost() const {
  return getColumnsCost(TOTAL_WEEKEND_COST);
}
