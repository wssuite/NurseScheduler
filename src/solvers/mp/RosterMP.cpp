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

void RosterPattern::checkReducedCost(const PDualCosts &pCosts,
                                     bool printBadPricing) {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of the roster:
   ************************************************/
  double reducedCost = cost_ - pCosts->constant();
  for (int k = firstDay(); k <= lastDay(); ++k) {
    const PShift &pS = pShift(k);
    /* Working dual cost */
    if (pS->isWork())
      reducedCost -= pCosts->workedDayShiftCost(k, pS->id);
  }

  // Display: set to true if you want to display the details of the cost
  if (std::fabs(reducedCost_ - reducedCost) / (1 - reducedCost) > 1e-3) {
    // if do not print and not throwing an error
    if (!printBadPricing && reducedCost_ > reducedCost + 1e-3) return;

    std::cerr << "# " << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "Bad dual cost: " << reducedCost_ << " != " << reducedCost
              << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "#   | Base cost     : + " << cost_ << std::endl;
    std::cerr << costsToString();

    std::cerr << "#   | Constant: - " << pCosts->constant() << std::endl;
    for (int k = firstDay(); k <= lastDay(); ++k) {
      const PShift &pS = pShift(k);
      if (pS->isWork())
        std::cerr << "#   | Work day-shift " << k << ": - "
                  << pCosts->workedDayShiftCost(k, pS->id)
                  << std::endl;
    }
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
//    if (reducedCost_ < reducedCost + 1e-3)
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
    MasterProblem(pScenario, solver),
    assignmentCons_(pScenario->nNurses()) {
  lagrangianBoundAvail_ = true;
}

RosterMP::~RosterMP() = default;

PPattern RosterMP::getPattern(MyVar *var) const {
  return std::make_shared<RosterPattern>(var, pScenario_);
}

// Main method to build the rostering problem for a given input
void RosterMP::build(const SolverParam &param) {
  /* Roster assignment constraints */
  buildAssignmentCons(param);

  /* build the rest of the model */
  MasterProblem::build(param);

  /* Change the branching rule */
  pRule_ = new RosterBranchingRule(this,
                                   dynamic_cast<RestTree *>(pTree_),
                                   "branching rule");
  pModel_->addBranchingRule(pRule_);
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
      pModel_->addInitialColumn(addRoster(pat, baseName));
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
  PDualCosts costs = buildDualCosts(theLiveNurses_[nurseNum]);
  pat.checkReducedCost(costs, pPricer_->isLastRunOptimal());
  checkIfPatternAlreadyPresent(pat.getCompactPattern());
#endif
  return addRoster(pat, "roster", false);
}

MyVar *RosterMP::addRoster(const RosterPattern &roster,
                           const char *baseName,
                           bool coreVar) {
  // nurse index
  const int nurseNum = roster.nurseNum();

  // Column var, its name, and affected constraints with their coefficients
  MyVar *var;
  char name[255];
  std::vector<MyCons *> cons;
  std::vector<double> coeffs;

  /* Skills coverage constraints */
  addSkillsCoverageConsToCol(&cons, &coeffs, roster);

  /* add variable to model */
  snprintf(name, sizeof(name), "%s_N%d", baseName, nurseNum);
  if (coreVar) {
    pModel_->createPositiveVar(&var,
                               name,
                               roster.cost(),
                               roster.getCompactPattern());
    unsigned int i;
    for (i = 0; i < static_cast<int>(cons.size()); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  } else {
    /* Roster  assignment constraint s added only for real rosters
     * to be sure that the artificial variables can always be used to create
     * a feasible solution
     */
    addRosterConsToCol(&cons, &coeffs, nurseNum);

    pModel_->createIntColumn(&var,
                             name,
                             roster.cost(),
                             roster.getCompactPattern(),
                             roster.reducedCost(),
                             cons,
                             coeffs);
  }
  return var;
}

/* Build each set of constraints
 * Add also the coefficient of a column for each set
 */
void RosterMP::buildAssignmentCons(const SolverParam &param) {
  char name[255];
  // build the roster assignment constraint for each nurse
  for (int i = 0; i < pScenario_->nNurses(); i++) {
    snprintf(name, sizeof(name), "feasibilityAssignmentVar_N%d", i);
    MyVar *feasibilityVar;
    pModel_->createPositiveFeasibilityVar(&feasibilityVar, name);
    snprintf(name, sizeof(name), "assignmentCons_N%d", i);
    pModel_->createEQConsLinear(&assignmentCons_[i],
                                name,
                                1,
                                {feasibilityVar},
                                {1});
  }
}

int RosterMP::addRosterConsToCol(std::vector<MyCons *> *cons,
                                 std::vector<double> *coeffs,
                                 int i) {
  cons->push_back(assignmentCons_[i]);
  coeffs->push_back(1);
  return 1;
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
  if (!stabCheckStoppingCriterion()) {
    std::cerr << "Cannot compute a lagrangian bound when stabilization "
                 "variables are present in the solution." << std::endl;
    return -LARGE_SCORE;
  }

  double sumRedCost = 0;
  for (double v : pPricer_->getLastMinReducedCosts())
    sumRedCost += v;
  return objVal + sumRedCost;
}

double RosterMP::getDaysCost() const {
  return getColumnsCost(TOTAL_WORK_COST);
}

double RosterMP::getWeekendCost() const {
  return getColumnsCost(TOTAL_WEEKEND_COST);
}

/* retrieve the dual values */
double RosterMP::getConstantDualvalue(PLiveNurse pNurse) const {
  return pModel_->getDual(assignmentCons_[pNurse->num_]);
}
