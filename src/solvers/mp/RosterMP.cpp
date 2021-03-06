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

#include "RosterMP.h"

#include <map>
#include <memory>

#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
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

void RosterPattern::computeCost(const MasterProblem *pMaster,
                                const PLiveNurse &pNurse) {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
 * Compute all the costs of a roster:
 ************************************************/
  PScenario pScenario = pMaster->pScenario();

  /*
  * Compute resources costs
  */
  cost_ = 0;
  computeResourcesCosts(pMaster, *pNurse->pStateIni_);

  /*
   * Compute complete weekend
   */
  if (pNurse->needCompleteWeekends()) {
    int k = 0;
    bool rest = false;
    for (const PShift &pS : stretch_.pShifts()) {
      // on sunday, if complete weekend, it's either:
      // work on saturday (rest=false) and sunday
      // rest on saturday (rest=true) and sunday
      if (Tools::isSunday(k) && (rest ^ pS->isRest()))
        addCost(pScenario->weights().WEIGHT_COMPLETE_WEEKEND,
                COMPLETE_WEEKEND_COST);
      rest = pS->isRest();
      k++;
    }
  }

  /*
   * Compute preferencesCost
   */
  for (int k = stretch_.firstDay(); k <= stretch_.lastDay(); ++k) {
    int level = pNurse->wishesOffLevel(k, shift(k));
    if (level != -1)
      addCost(pScenario->weights().WEIGHT_PREFERENCES_OFF[level],
              PREFERENCE_COST);
    level = pNurse->wishesOnLevel(k, shift(k));
    if (level != -1)
      addCost(pScenario->weights().WEIGHT_PREFERENCES_ON[level],
              PREFERENCE_COST);
  }
}

void RosterPattern::checkReducedCost(const PDualCosts &pCosts,
                                     PScenario Scenario,
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
  if (abs(reducedCost_ - reducedCost) / (1 - reducedCost) > 1e-3) {
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
    std::cerr << toString(pCosts->nDays());
    std::cerr << "# " << std::endl;

    // throw an error only when a significant misprice
    // Indeed, if the real reduced cost (reducedCost) is greater than the one
    // found by the pricing, some columns with a positive reduced cost could be
    // generated.
    // The reason why the other situation can arise is that some path in the
    // subproblem could under estimate the real cost. These paths won't be
    // found when the subproblems are solved at optimality, but could  be
    // present when using heuristics.
    if (reducedCost_ < reducedCost + 1e-3)
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

RosterMP::RosterMP(PScenario pScenario,
                   PDemand pDemand,
                   PPreferences pPreferences,
                   std::vector<State> *pInitState,
                   SolverType solver) :
    MasterProblem(pScenario, pDemand, pPreferences, pInitState, solver),
    assignmentCons_(pScenario->nbNurses_) {
  lagrangianBoundAvail_ = true;
}

RosterMP::~RosterMP() {}

PPattern RosterMP::getPattern(MyVar *var) const {
  return std::make_shared<RosterPattern>(var->getPattern(), pScenario_);
}

// Main method to build the rostering problem for a given input
void RosterMP::build(const SolverParam &param) {
  /* Roster assignment constraints */
  buildAssignmentCons(param);

  /* build the rest of the model */
  MasterProblem::build(param);

  /* Change the branching rule */
  pRule_ = new RosterBranchingRule(this,
                                   static_cast<RestTree *>(pTree_),
                                   "branching rule");
  pModel_->addBranchingRule(pRule_);
}

// Provide an initial solution to the solver. If empty, add artificial columns
void RosterMP::initializeSolution(const std::vector<Roster> &solution) {
  // rosters are added for each nurse of the initial solution
  if (!solution.empty()) {
    const char *baseName("initialRoster");
    // build the roster of each nurse
    for (int i = 0; i < pScenario_->nbNurses_; ++i) {
      // load the roster of nurse i
      const Roster &roster = solution[i];
      // build the shift vector
      std::vector<int> shifts(nDays());
      for (int k = 0; k < nDays(); ++k)
        shifts[k] = roster.shift(k);
      RosterPattern pat(shifts, pScenario_, i);
      pat.computeCost(this, theLiveNurses_[i]);
      pModel_->addInitialColumn(addRoster(pat, baseName));
    }
  }
}

std::vector<PResource> RosterMP::createResources(
    const PLiveNurse &pN, std::map<int, CostType> *resourceCostType) const {
  std::vector<PResource> resources;
  PWeights pW = pScenario_->pWeights_;
  int k = 0;

  // initialize resource on the total number of working days
  resources.emplace_back(std::make_shared<SoftTotalShiftDurationResource>(
      pN->minTotalShifts(),
      pN->maxTotalShifts(),
      pW->WEIGHT_TOTAL_SHIFTS,
      pW->WEIGHT_TOTAL_SHIFTS,
      std::make_shared<AnyWorkShift>(),
      nDays()));
  if (resourceCostType) (*resourceCostType)[k] = DAYS_COST;
  resources.back()->setId(k++);

  // initialize resource on the total number of week-ends
  resources.emplace_back(std::make_shared<SoftTotalWeekendsResource>(
      pN->maxTotalWeekends(),
      pW->WEIGHT_TOTAL_WEEKENDS,
      nDays()));
  if (resourceCostType) (*resourceCostType)[k] = WEEKEND_COST;
  resources.back()->setId(k++);

  // initialize resource on the number of consecutive worked days
  resources.emplace_back(std::make_shared<SoftConsShiftResource>(
      pN->minConsDaysWork(),
      pN->maxConsDaysWork(),
      pW->WEIGHT_CONS_DAYS_WORK,
      pW->WEIGHT_CONS_DAYS_WORK,
      std::make_shared<AnyWorkShift>(),
      nDays()));
  if (resourceCostType) (*resourceCostType)[k] = CONS_WORK_COST;
  resources.back()->setId(k++);

  // initialize resources on the number of consecutive shifts of each type
  for (int st = 0; st < pScenario_->nbShiftsType_; st++) {
    shared_ptr<AbstractShift> absShift =
        std::make_shared<AnyOfTypeShift>(st, pScenario_->intToShiftType_[st]);
    if (absShift->isWork()) {
      resources.emplace_back(std::make_shared<SoftConsShiftResource>(
          pScenario_->minConsShiftsOf(st),
          pScenario_->maxConsShiftsOf(st),
          pW->WEIGHT_CONS_SHIFTS,
          pW->WEIGHT_CONS_SHIFTS,
          absShift,
          nDays()));
      if (resourceCostType) (*resourceCostType)[k] = CONS_SHIFTS_COST;
    } else if (absShift->isRest()) {
      resources.emplace_back(std::make_shared<SoftConsShiftResource>(
          pN->minConsDaysOff(),
          pN->maxConsDaysOff(),
          pW->WEIGHT_CONS_DAYS_WORK,
          pW->WEIGHT_CONS_DAYS_WORK,
          absShift,
          nDays()));
      if (resourceCostType) (*resourceCostType)[k] = CONS_REST_COST;
    }
    resources.back()->setId(k++);
  }

  return resources;
}

// Create a new rotation variable
// add the correct constraints and coefficients for the nurse i working
// on a rotation
// if s=-1, the nurse works on all shifts
MyVar *RosterMP::addColumn(int nurseNum, const RCSolution &solution) {
  // Build rotation from RCSolution
  RosterPattern pat(
      solution.shifts, pScenario_, nurseNum, DBL_MAX, solution.cost);
  pat.computeCost(this, theLiveNurses_[nurseNum]);
  pat.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  PDualCosts costs = buildDualCosts(theLiveNurses_[nurseNum]);
  pat.checkReducedCost(costs, pScenario_, pPricer_->isLastRunOptimal());
  std::vector<double> pattern = pat.getCompactPattern();
  checkIfPatternAlreadyPresent(pattern);
#endif
  return addRoster(pat, "roster", false);
}

MyVar *RosterMP::addRoster(const RosterPattern &roster,
                           const char *baseName,
                           bool coreVar) {
  // nurse index
  const int nurseNum = roster.nurseNum_;

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
                               roster.cost_,
                               roster.getCompactPattern());
    for (unsigned int i = 0; i < cons.size(); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  } else {
    /* Roster  assignment constraint s added only for real rosters
     * to be sure that the artificial variables can always be used to create
     * a feasible solution
     */
    addRosterConsToCol(&cons, &coeffs, nurseNum);

    pModel_->createIntColumn(&var,
                             name,
                             roster.cost_,
                             roster.getCompactPattern(),
                             roster.reducedCost_,
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
  for (int i = 0; i < pScenario_->nbNurses_; i++) {
    snprintf(name, sizeof(name), "feasibilityAssignmentVar_N%d", i);
    MyVar *feasibilityVar;
    pModel_->createPositiveFeasibilityVar(&feasibilityVar, name);
    snprintf(name, sizeof(name), "assignmentCons_N%d", i);
    pModel_->createEQConsLinear(&assignmentCons_[i],
                                name,
                                1,
                                {feasibilityVar},
                                {1});
    // STAB:Add stabilization variable
    if (param.isStabilization_) {
      snprintf(name, sizeof(name), "stabAssignment_%i", i);
      addStabVariables(param, name, assignmentCons_[i], true, true);
    }
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

// return the costs of all active columns associated to the type
double RosterMP::getColumnsCost(CostType costType) const {
  return getColumnsCost(costType, pModel_->getActiveColumns());
}

double RosterMP::getColumnsCost(CostType costType,
                                const std::vector<MyVar *> &vars) const {
  double cost = 0;
  for (MyVar *var : vars) {
    double value = pModel_->getVarValue(var);
    if (value > epsilon()) {
      RosterPattern ros(var->getPattern(), pScenario_);
      ros.computeCost(this, theLiveNurses_[ros.nurseNum_]);
      cost += ros.cost(costType) * value;
    }
  }
  return cost;
}

double RosterMP::getDaysCost() const {
  return getColumnsCost(DAYS_COST);
}

double RosterMP::getWeekendCost() const {
  return getColumnsCost(WEEKEND_COST);
}

/* retrieve the dual values */
double RosterMP::getConstantDualvalue(PLiveNurse pNurse) const {
  return pModel_->getDual(assignmentCons_[pNurse->num_]);
}
