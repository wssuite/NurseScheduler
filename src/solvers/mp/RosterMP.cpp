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

void RosterPattern::checkDefaultCost(const MasterProblem *pMaster,
                                     const PLiveNurse &pNurse) const {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the costs of a roster:
   ************************************************/
  PScenario pScenario = pMaster->pScenario();

  // initialize costs
  double cost = 0;

  /*
   * Compute initial cost
   */

  // initial state of the nurse
  int lastShiftType = pNurse->pStateIni_->shiftType_;
  int nbConsShifts = pNurse->pStateIni_->consShifts_;
  int nbConsDaysWorked = pNurse->pStateIni_->consDaysWorked_;
  int nbConsDaysOff = pNurse->pStateIni_->consDaysOff_;

  // 1. if the initial shift has already exceeded the max,
  // fix these values to the UB to not pay it twice
  if (lastShiftType > 0) {  // was working
    if (nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))
      nbConsShifts = pScenario->maxConsShiftsOf(lastShiftType);
    if (nbConsDaysWorked > pNurse->maxConsDaysWork())
      nbConsDaysWorked = pNurse->maxConsDaysWork();
  } else if (nbConsShifts > pNurse->maxConsDaysOff()) {
    nbConsDaysOff = pNurse->maxConsDaysOff();
  }

  // 2. compute the cost
  for (const PShift &pS : stretch_.pShifts()) {
    // a. same shift type -> increment the counters
    if (lastShiftType == pS->type) {
      if (pS->isWork()) {
        nbConsShifts++;
        nbConsDaysWorked++;
      } else {
        nbConsDaysOff++;
      }
      continue;
    }
    // b. different shift type -> add the corresponding cost
    // i) goes to rest
    if (pS->isRest()) {
      // compute cost
      cost += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
      cost += pNurse->consDaysCost(nbConsDaysWorked);
      // update counters
      nbConsShifts = 0;
      nbConsDaysWorked = 0;
      nbConsDaysOff = 1;
    } else if (lastShiftType == 0) {
      // ii) goes to work
      // compute cost
      cost += pNurse->consDaysOffCost(nbConsDaysOff);
      // update counters
      nbConsDaysOff = 0;
      nbConsDaysWorked = 1;
      nbConsShifts = 1;
    } else {
      // iii) continue to work on a different shift
      // compute cost
      cost += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
      // update counters
      nbConsShifts = 1;
      nbConsDaysWorked++;
    }
    // update
    lastShiftType = pS->type;
  }

  // pay the max for last day
  if (nbConsDaysOff > pNurse->maxConsDaysOff())
    cost += pNurse->consDaysOffCost(nbConsDaysOff);
  if (lastShiftType > 0
      && nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))
    cost += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
  if (nbConsDaysWorked > pNurse->maxConsDaysWork())
    cost += pNurse->consDaysCost(nbConsDaysWorked);

  /*
   * Compute preferencesCost
   */

  for (int k = firstDay(); k <= lastDay(); ++k) {
    int level = pNurse->wishesOffLevel(k, shift(k));
    if (level != -1)
      cost += pScenario->weights().WEIGHT_PREFERENCES_OFF[level];
    level = pNurse->wishesOnLevel(k, shift(k));
    if (level != -1)
      cost += pScenario->weights().WEIGHT_PREFERENCES_ON[level];
  }

  /*
   * Compute time duration and complete weekend
   */
  int nWeekends = 0;
  int k = 0;
  bool rest = false;
  for (const PShift &pS : stretch_.pShifts()) {
    if (pS->isWork()) {
      if (Tools::isSaturday(k)) {
        nWeekends++;
      } else if (rest && Tools::isSunday(k)) {
        nWeekends++;
        if (pNurse->needCompleteWeekends())
          cost += pScenario->weights().WEIGHT_COMPLETE_WEEKEND;
      }
      rest = false;
    } else {
      if (pNurse->needCompleteWeekends() && !rest && Tools::isSunday(k))
        cost += pScenario->weights().WEIGHT_COMPLETE_WEEKEND;
      rest = true;
    }
    k++;
  }

  if (pNurse->minTotalShifts() - duration() > 0)
    cost += pScenario->weights().WEIGHT_TOTAL_SHIFTS
        * (pNurse->minTotalShifts() - duration());
  if (duration() - pNurse->maxTotalShifts() > 0)
    cost += pScenario->weights().WEIGHT_TOTAL_SHIFTS
        * (duration() - pNurse->maxTotalShifts());
  cost += pNurse->totalWeekendCost(nWeekends);

  if (std::abs(cost - cost_) > 1e-3) {
    std::cerr << "# " << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: " << cost_ << " != " << cost
              << std::endl;
    std::cerr << "# " << std::endl;
    std::cerr << "#   | Base cost     : + " << cost_ << std::endl;
    std::cerr << costsToString();
    std::cerr << toString(pMaster->nDays());
    std::cerr << "# " << std::endl;
    Tools::throwError("defaultComputeCost and computeCost do not get "
                      "the same result for RosterPattern.");
  }
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

RosterMP::RosterMP(const PScenario& pScenario,
                   PDemand pDemand,
                   PPreferences pPreferences,
                   std::vector<State> *pInitState,
                   SolverType solver) :
    MasterProblem(pScenario, std::move(pDemand), std::move(pPreferences),
        pInitState, solver),
    assignmentCons_(pScenario->nNurses()) {
  lagrangianBoundAvail_ = true;
}

RosterMP::~RosterMP() = default;

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

std::map<PResource, CostType>
    RosterMP::defaultgeneratePResources(const PLiveNurse &pN) const {
  const Weights &w = pScenario_->weights();

  std::map<PResource, CostType> mResources = {
      // initialize resource on the total number of working days
      {std::make_shared<SoftTotalShiftDurationResource>(
          pN->minTotalShifts(),
          pN->maxTotalShifts(),
          w.WEIGHT_TOTAL_SHIFTS,
          w.WEIGHT_TOTAL_SHIFTS,
          std::make_shared<AnyWorkShift>(),
          nDays()), TOTAL_WORK_COST},
      // initialize resource on the total number of week-ends
      {std::make_shared<SoftTotalWeekendsResource>(
          pN->maxTotalWeekends(),
          w.WEIGHT_TOTAL_WEEKENDS,
          nDays()), TOTAL_WEEKEND_COST},
      // initialize resource on the number of consecutive worked days
      {std::make_shared<SoftConsShiftResource>(
          pN->minConsDaysWork(),
          pN->maxConsDaysWork(),
          w.WEIGHT_CONS_DAYS_WORK,
          w.WEIGHT_CONS_DAYS_WORK,
          std::make_shared<AnyWorkShift>(),
          nDays(),
          pN->pStateIni_->consDaysWorked_),
       CONS_WORK_COST}
  };

  // initialize resources on the number of consecutive shifts of each type
  for (int st = 0; st < pScenario_->nShiftTypes(); st++) {
    shared_ptr<AbstractShift> absShift =
        std::make_shared<AnyOfTypeShift>(st, pScenario_->shiftType(st));
    if (absShift->isWork()) {
      int consShiftsInitial =
          absShift->includes(*pN->pStateIni_->pShift_) ?
          pN->pStateIni_->consShifts_ : 0;
      mResources[std::make_shared<SoftConsShiftResource>(
          pScenario_->minConsShiftsOf(st),
          pScenario_->maxConsShiftsOf(st),
          w.WEIGHT_CONS_SHIFTS,
          w.WEIGHT_CONS_SHIFTS,
          absShift,
          nDays(),
          consShiftsInitial)] = CONS_SHIFTS_COST;
    } else if (absShift->isRest()) {
      mResources[std::make_shared<SoftConsShiftResource>(
          pN->minConsDaysOff(),
          pN->maxConsDaysOff(),
          w.WEIGHT_CONS_DAYS_WORK,
          w.WEIGHT_CONS_DAYS_WORK,
          absShift,
          nDays(),
          pN->pStateIni_->consDaysOff_)] = CONS_REST_COST;
    }
  }

  return mResources;
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
//  if (useDefaultResources())
//    pat.checkDefaultCost(this, theLiveNurses_[nurseNum]);
#ifdef DBG
  PDualCosts costs = buildDualCosts(theLiveNurses_[nurseNum]);
  pat.checkReducedCost(costs, pPricer_->isLastRunOptimal());
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
  return getColumnsCost(TOTAL_WORK_COST);
}

double RosterMP::getWeekendCost() const {
  return getColumnsCost(TOTAL_WEEKEND_COST);
}

/* retrieve the dual values */
double RosterMP::getConstantDualvalue(PLiveNurse pNurse) const {
  return pModel_->getDual(assignmentCons_[pNurse->num_]);
}
