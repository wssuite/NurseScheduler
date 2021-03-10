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

#include "RotationMP.h"

#include <memory>
#include <algorithm>
#include <utility>

#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"

using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

// when branching on this pattern, this method add the corresponding forbidden
// shifts to the set.
// It will forbid all the shifts that would be worked on a day that is already
// covered by this pattern.
// Moreover, there needs to be a resting day before and after each rotation,
// so the shifts can also be forbidden on these two days
// (if the rotation is not at an extremity of the horizon).
void RotationPattern::addForbiddenShifts(
    std::set<std::pair<int, int> > *forbiddenShifts,
    int nbShifts, PDemand pDemand) const {
  // from the previous day to the day after the end of the rotation,
  // forbid any work shifts
  for (int day = firstDay()-1; day <= lastDay()+1; day++) {
    if (day < pDemand->firstDay_) continue;
    if (day >= pDemand->firstDay_ + pDemand->nDays_) continue;
    for (int i = 1; i < nbShifts; ++i)
      forbiddenShifts->insert(pair<int, int>(day, i));
  }
}

void RotationPattern::computeCost(const MasterProblem *pMaster,
                                  const PLiveNurse &pNurse) {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  PScenario pScenario = pMaster->pScenario();

  /************************************************
   * Compute all the costs of a rotation:
   ************************************************/
  // if first day of the planning, check on the past, otherwise 0 (rest)
  int lastShiftType = (firstDay() == 0) ? pNurse->pStateIni_->shiftType_ : -1;
  // nbConsShift = number of consecutive shift
  // if first day of the planning, check on the past, otherwise 0
  int nbConsShifts = (firstDay() == 0) ? pNurse->pStateIni_->consShifts_ : 0;
  // consShiftCost = cost of be outside of the interval [min,max] of the
  // consecutive shifts
  consShiftsCost_ = 0;

  // nbConsWorked = number of consecutive worked days
  // if first day of the planning, check on the past , otherwise 0
  int nbConsDaysWorked =
      (firstDay() == 0) ? pNurse->pStateIni_->consDaysWorked_ : 0;
  // consWorkedCost = cost of be outside of the interval [min,max] of the
  // consecutive worked days
  consDaysWorkedCost_ = 0;

  // cost of not doing the whole weekend
  completeWeekendCost_ = 0;

  // preferencesCost = cost of not respecting preferences
  preferenceCost_ = 0;

  // initial resting cost
  initRestCost_ = 0;

  /*
   * Compute consShiftCost
   */

  // if the initial shift has already exceeded the max, substract now the cost
  // that will be readd later
  if ((firstDay() == 0) && (lastShiftType > 0) &&
      (nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))) {
    consShiftsCost_ -=
        (nbConsShifts - pScenario->maxConsShiftsOf(lastShiftType))
            * pScenario->weights().WEIGHT_CONS_SHIFTS;
  }

  for (int k = firstDay(); k <= lastDay(); ++k) {
    const PShift &pS = pShift(k);
    if (lastShiftType == pS->type) {
      nbConsShifts++;
      continue;
    }
    if (lastShiftType > 0) {
      int diff = max(pScenario->minConsShiftsOf(lastShiftType) - nbConsShifts,
                     nbConsShifts - pScenario->maxConsShiftsOf(lastShiftType));
      if (diff > 0) {
        consShiftsCost_ += diff * pScenario->weights().WEIGHT_CONS_SHIFTS;
      }
    }
    // initialize nbConsShifts and lastShift
    nbConsShifts = 1;
    lastShiftType = pS->type;
  }

  // compute consShiftsCost for the last shift
  if (lastShiftType > 0) {
    int diff = max((lastDay()+1 == pMaster->nDays()) ? 0 :
                   pScenario->minConsShiftsOf(lastShiftType) - nbConsShifts,
                   nbConsShifts - pScenario->maxConsShiftsOf(lastShiftType));
    if (diff > 0)
      consShiftsCost_ += diff * pScenario->weights().WEIGHT_CONS_SHIFTS;
  }


  /*
   * Compute consDaysWorkedCost
   */

  // if already worked too much
  double diffDays = nbConsDaysWorked - pNurse->pContract_->maxConsDaysWork_;
  consDaysWorkedCost_ =
      (diffDays > 0) ? -diffDays * pScenario->weights().WEIGHT_CONS_DAYS_WORK
                     : 0;

  nbConsDaysWorked += nDays();
  // check if nbConsDaysWorked < min, if finishes on last day, does not count
  if (nbConsDaysWorked < pNurse->minConsDaysWork()
      && lastDay() + 1 < pMaster->nDays())
    consDaysWorkedCost_ += (pNurse->minConsDaysWork() - nbConsDaysWorked)
        * pScenario->weights().WEIGHT_CONS_DAYS_WORK;
  else if (nbConsDaysWorked > pNurse->maxConsDaysWork())
    // check if nbConsDaysWorked > max
    consDaysWorkedCost_ += (nbConsDaysWorked - pNurse->maxConsDaysWork())
        * pScenario->weights().WEIGHT_CONS_DAYS_WORK;

  /*
   * Compute completeWeekendCost
   */
  if (pNurse->needCompleteWeekends()) {
    // if first day is a Sunday, the saturday is not worked
    if (Tools::isSunday(firstDay()))
      completeWeekendCost_ += pScenario->weights().WEIGHT_COMPLETE_WEEKEND;
    // if last day + 1 is a Sunday, the sunday is not worked
    if (Tools::isSunday(lastDay()+1))
      completeWeekendCost_ += pScenario->weights().WEIGHT_COMPLETE_WEEKEND;
  }

  /*
   * Compute preferencesCost
   */

  for (int k = firstDay(); k <= lastDay(); ++k) {
    int level = pNurse->wishesOffLevel(k, shift(k));
    if (level != -1)
      preferenceCost_ += pScenario->weights().WEIGHT_PREFERENCES_OFF[level];
    level = pNurse->wishesOnLevel(k, shift(k));
    if (level != -1)
      preferenceCost_ += pScenario->weights().WEIGHT_PREFERENCES_ON[level];
  }

  /*
   * Compute initial resting cost
   */

  if (firstDay() == 0 && pNurse->pStateIni_->shiftType_ == 0) {
    int diff = pNurse->minConsDaysOff() - pNurse->pStateIni_->consDaysOff_;
    initRestCost_ =
        (diff > 0) ? diff * pScenario->weights().WEIGHT_CONS_DAYS_OFF : 0;
  }

  /*
   * Compute the sum of the cost and stores it in cost_
   */
  if (false) {
    cout << "# Costs:" << endl;
    cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
    cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
    cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
    cout << "#       | Preferences       : " << preferenceCost_ << endl;
    cout << "#       | Initial rest      : " << initRestCost_ << endl;
    cout << "# " << endl;
  }

  cost_ = consShiftsCost_ + consDaysWorkedCost_ + completeWeekendCost_
      + preferenceCost_ + initRestCost_;
}

void RotationPattern::checkReducedCost(const PDualCosts &pCosts,
                                       bool printBadPricing) {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of a rotation:
   ************************************************/

  double reducedCost(cost_);

  /* Working dual cost */
  for (int k = firstDay(); k <= lastDay(); ++k)
    reducedCost -= pCosts->workedDayShiftCost(k, shift(k));
  /* Start working dual cost */
  reducedCost -= pCosts->startWorkCost(firstDay());
  /* Stop working dual cost */
  reducedCost -= pCosts->endWorkCost(lastDay());
  /* Working on weekend */
  if (Tools::isSunday(firstDay()))
    reducedCost -= pCosts->workedWeekendCost();
  for (int k = firstDay(); k <= lastDay(); ++k)
    if (Tools::isSaturday(k))
      reducedCost -= pCosts->workedWeekendCost();


  // Display: set to true if you want to display the details of the cost
  if (std::fabs(reducedCost_ - reducedCost) / (1 - reducedCost) > 1e-3) {
    // if do not print and not throwing an error
    if (!printBadPricing && reducedCost_ > reducedCost + 1e-3) return;

    cout << "# " << endl;
    cout << "# " << endl;
    cout << "# Bad dual cost: "
         << reducedCost_ << " != " << reducedCost << endl;
    cout << "# " << endl;
    cout << "#   | Base cost     : + " << cost_ << endl;

    cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
    cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
    cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
    cout << "#       | Preferences       : " << preferenceCost_ << endl;
    cout << "#       | Initial rest      : " << initRestCost_ << endl;

    for (int k = firstDay(); k <= lastDay(); ++k)
      cout << "#   | Work day-shift: - "
           << pCosts->workedDayShiftCost(k, shift(k)) << endl;
    cout << "#   | Start work    : - " << pCosts->startWorkCost(firstDay())
         << endl;
    cout << "#   | Finish Work   : - "
         << pCosts->endWorkCost(lastDay()) << endl;
    if (Tools::isSunday(firstDay()))
      cout << "#   | Weekends      : - " << pCosts->workedWeekendCost() << endl;
    for (int k = firstDay(); k <= lastDay(); ++k)
      if (Tools::isSaturday(k))
        cout << "#   | Weekends      : - "
             << pCosts->workedWeekendCost() << endl;
    std::cout << toString(pCosts->nDays());
    cout << "# " << endl;

    // throw an error only when a significant misprice
    // Indeed, if the real reduced cost (reducedCost) is greater than the one
    // found by the pricing, some columns with a positive reduced cost could be
    // generated.
    // The reason why the other situation can arise is that some path in the
    // subproblem could under estimate the real cost. These paths won't be
    // found when the subproblems are solved at optimality, but could  be
    // present when using heuristics.
    if (reducedCost_ < reducedCost + 1e-3)
      Tools::throwError("Invalid pricing of a rotation.");
  }
}


//-----------------------------------------------------------------------------
//
//  C l a s s   R o t a t i o n M P
//
// Build and solve the master problem of the column generation scheme with
// rotations
//
//-----------------------------------------------------------------------------

RotationMP::RotationMP(const PScenario& pScenario,
                       PDemand pDemand,
                       PPreferences pPreferences,
                       std::vector<State> *pInitState,
                       SolverType solver) :
    MasterProblem(pScenario, std::move(pDemand), std::move(pPreferences),
                  pInitState, solver),
    restsPerDay_(pScenario->nNurses()),
    restingVars_(pScenario->nNurses()),
    longRestingVars_(pScenario->nNurses()),
    minWorkedDaysVars_(pScenario->nNurses()),
    maxWorkedDaysVars_(pScenario->nNurses()),
    maxWorkedWeekendVars_(pScenario->nNurses()),
    minWorkedDaysAvgVars_(pScenario->nNurses()),
    maxWorkedDaysAvgVars_(pScenario->nNurses()),
    maxWorkedWeekendAvgVars_(pScenario_->nNurses()),
    minWorkedDaysContractAvgVars_(pScenario->nContracts()),
    maxWorkedDaysContractAvgVars_(pScenario->nContracts()),
    maxWorkedWeekendContractAvgVars_(pScenario_->nContracts()),

    restFlowCons_(pScenario->nNurses()),
    workFlowCons_(pScenario->nNurses()),
    minWorkedDaysCons_(pScenario->nNurses()),
    maxWorkedDaysCons_(pScenario->nNurses()),
    maxWorkedWeekendCons_(pScenario->nNurses()),
    minWorkedDaysAvgCons_(pScenario->nNurses()),
    maxWorkedDaysAvgCons_(pScenario->nNurses()),
    maxWorkedWeekendAvgCons_(pScenario_->nNurses()),
    minWorkedDaysContractAvgCons_(pScenario->nContracts()),
    maxWorkedDaysContractAvgCons_(pScenario->nContracts()),
    maxWorkedWeekendContractAvgCons_(pScenario_->nContracts()) {
  // initialize the vectors indicating whether the min/max total constraints
  // with averaged bounds are considered
  Tools::initVector(&isMinWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedWeekendAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMinWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedWeekendContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
}

RotationMP::~RotationMP() = default;

PPattern RotationMP::getPattern(MyVar *var) const {
  return std::make_shared<RotationPattern>(var->getPattern(), pScenario_);
}

// build the, possibly fractional, roster corresponding to the solution
// currently stored in the model
vector3D<double> RotationMP::fractionalRoster() const {
  vector3D<double> roster = MasterProblem::fractionalRoster();

  // add rest values
  for (auto &nurseRoster : roster)
    for (auto &dayRoster : nurseRoster) {
      // fill the rest shift
      double restValue = 1;
      for (int s = 1; s < nShifts(); ++s)
        restValue -= dayRoster[s];
      dayRoster[0] = restValue;
    }

  return roster;
}

// build the rostering problem
void RotationMP::build(const SolverParam &param) {
  /* Rotation constraints */
  buildRotationCons(param);

  /* Min/Max constraints */
  buildMinMaxCons(param);

  /* build the rest of the model */
  MasterProblem::build(param);
}

// Build the columns corresponding to the initial solution
void RotationMP::initializeSolution(const vector<Roster> &solution) {
  if (solution.empty()) return;

  // rotations are added for each nurse of the initial solution
  string baseName("initialRotation");
  // build the rotations of each nurse
  for (int i = 0; i < pScenario_->nNurses(); ++i) {
    // load the roster of nurse i
    Roster roster = solution[i];

    bool workedLastDay = false;
    int lastShift = 0;
    map<int, int> shifts;
    // build all the successive rotation of this nurse
    for (int k = 0; k < pDemand_->nDays_; ++k) {
      // shift=0 => rest
      int shift = roster.shift(k);
      // if work, insert the shift in the map
      if (shift > 0) {
        shifts[k] = shift;
        lastShift = shift;
        workedLastDay = true;
      } else if (shift < 0 && lastShift > 0) {
        shifts[k] = lastShift;
        workedLastDay = true;
      } else if (workedLastDay) {
        // if stop to work, build the rotation
        RotationPattern rotation(shifts, pScenario_, i);
        rotation.computeCost(this, theLiveNurses_[i]);
        pModel_->addInitialColumn(addRotation(rotation, baseName.c_str()));
        shifts.clear();
        lastShift = shift;
        workedLastDay = false;
      }
    }
    // if work on the last day, build the rotation
    if (workedLastDay) {
      RotationPattern rotation(shifts, pScenario_, i);
      rotation.computeCost(this, theLiveNurses_[i]);
      pModel_->addInitialColumn(addRotation(rotation, baseName.c_str()));
      shifts.clear();
    }
  }
}

PDualCosts RotationMP::buildDualCosts(PLiveNurse pNurse) const {
  return std::make_shared<RotationDualCosts>(
      getShiftsDualValues(pNurse),
      getStartWorkDualValues(pNurse),
      getEndWorkDualValues(pNurse),
      getWorkedWeekendDualValue(pNurse),
      getConstantDualvalue(pNurse));
}

vector2D<double> RotationMP::getShiftsDualValues(PLiveNurse pNurse) const {
  vector2D<double> dualValues = MasterProblem::getShiftsDualValues(pNurse);

  const int i = pNurse->num_;
  const int p = pNurse->pContract_->id_;

  /* Min/Max constraints */
  double minWorkedDays = pModel_->getDual(minWorkedDaysCons_[i], true);
  double maxWorkedDays = pModel_->getDual(maxWorkedDaysCons_[i], true);

  double minWorkedDaysAvg =
      isMinWorkedDaysAvgCons_[i] ?
      pModel_->getDual(minWorkedDaysAvgCons_[i], true) : 0.0;
  double maxWorkedDaysAvg =
      isMaxWorkedDaysAvgCons_[i] ?
      pModel_->getDual(maxWorkedDaysAvgCons_[i], true) : 0.0;

  double minWorkedDaysContractAvg =
      isMinWorkedDaysContractAvgCons_[p] ?
      pModel_->getDual(minWorkedDaysContractAvgCons_[p], true) : 0.0;
  double maxWorkedDaysContractAvg =
      isMaxWorkedDaysContractAvgCons_[p] ?
      pModel_->getDual(maxWorkedDaysContractAvgCons_[p], true) : 0.0;

  /* Min/Max constraints */
  double d = minWorkedDays + minWorkedDaysAvg + minWorkedDaysContractAvg;
  d += maxWorkedDays + maxWorkedDaysAvg + maxWorkedDaysContractAvg;

  for (int k = 0; k < nDays(); ++k) {
    vector<double> &dualValues2 = dualValues[k];
    for (int s = 1; s < pScenario_->nShifts(); ++s)
      // adjust the dual in function of the time duration of the shift
      dualValues2[s - 1] += d * pScenario_->duration(s);
  }

  return dualValues;
}

vector<double> RotationMP::getStartWorkDualValues(const PLiveNurse& pNurse)
const {
  int i = pNurse->num_;
  vector<double> dualValues(nDays());

  // get dual value associated to the source
  dualValues[0] = pModel_->getDual(restFlowCons_[i][0], true);
  // get dual values associated to the work flow constraints
  // don't take into account the last which is the sink
  for (int k = 1; k < nDays(); ++k)
    dualValues[k] = pModel_->getDual(workFlowCons_[i][k - 1], true);

  return dualValues;
}

vector<double> RotationMP::getEndWorkDualValues(const PLiveNurse& pNurse)
const {
  const int i = pNurse->num_;
  vector<double> dualValues(nDays());

  // get dual values associated to the work flow constraints
  // don't take into account the first which is the source
  // take into account the cost, if the last day worked is k
  for (int k = 0; k < nDays() - 1; ++k)
    dualValues[k] = -pModel_->getDual(restFlowCons_[i][k + 1], true);

  // get dual value associated to the sink
  dualValues[nDays() - 1] =
      pModel_->getDual(workFlowCons_[i][nDays() - 1], true);

  return dualValues;
}

double RotationMP::getWorkedWeekendDualValue(const PLiveNurse& pNurse) const {
  int id = pNurse->num_;
  double dualVal = pModel_->getDual(maxWorkedWeekendCons_[id], true);
  if (isMaxWorkedWeekendAvgCons_[id])
    dualVal += pModel_->getDual(maxWorkedWeekendAvgCons_[id], true);
  if (isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_])
    dualVal += pModel_->getDual(
        maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
  return dualVal;
}

PDualCosts RotationMP::buildRandomDualCosts(bool optimalDemandConsidered,
                                            int NDaysShifts) const {
  return std::make_shared<RotationDualCosts>(
      getRandomWorkedDualCosts(optimalDemandConsidered, NDaysShifts),
      Tools::randomDoubleVector(
          pDemand_->nDays_,
          0, 7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND),
      Tools::randomDoubleVector(
          pDemand_->nDays_,
          0,
          7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND),
      Tools::randomDouble(0, 2*pScenario_->weights().WEIGHT_TOTAL_WEEKENDS),
      Tools::randomDouble(-10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND,
                          10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND));
}

//------------------------------------------------------------------------------
// Build the variable of the rotation as well as all the affected constraints
// with their coefficients. if s=-1, the nurse i works on all shifts
//------------------------------------------------------------------------------
MyVar *RotationMP::addColumn(int nurseNum, const RCSolution &solution) {
  // Build rotation from RCSolution
  RotationPattern rotation
      (solution.firstDay, solution.shifts, pScenario_,
       nurseNum, DBL_MAX, solution.cost);
  rotation.computeCost(this, theLiveNurses_[nurseNum]);
  rotation.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  PDualCosts costs = buildDualCosts(theLiveNurses_[nurseNum]);
  rotation.checkReducedCost(costs, pPricer_->isLastRunOptimal());
  std::vector<double> pattern = rotation.getCompactPattern();
  checkIfPatternAlreadyPresent(pattern);
#endif
  return addRotation(rotation, "rotation", false);
}

MyVar *RotationMP::addRotation(const RotationPattern &rotation,
                               const char *baseName,
                               bool coreVar) {
  // nurse index
  const int nurseNum = rotation.nurseNum_;

  // Column var, its name, and affected constraints with their coefficients
  MyVar *var;
  char name[255];
  vector<MyCons *> cons;
  vector<double> coeffs;

  /* Min/Max constraints */
  int nbWeekends = Tools::nWeekendsInInterval(rotation.firstDay(),
                                              rotation.lastDay());
  addMinMaxConsToCol(&cons, &coeffs, nurseNum, rotation.duration(), nbWeekends);

  /* Skills coverage constraints */
  addSkillsCoverageConsToCol(&cons, &coeffs, rotation);

  snprintf(name, sizeof(name), "%s_N%d", baseName, nurseNum);
  if (coreVar) {
    pModel_->createPositiveVar(&var,
                               name,
                               rotation.cost_,
                               rotation.getCompactPattern());
    for (unsigned int i = 0; i < cons.size(); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  } else {
    /* Rotation constraints
      They are added only for real rotations to be sure that the artificial variables
      can always be used to create a feasible solution
     */
    addRotationConsToCol(&cons,
                         &coeffs,
                         nurseNum,
                         rotation.firstDay(),
                         true,
                         false);
    addRotationConsToCol(&cons,
                         &coeffs,
                         nurseNum,
                         rotation.lastDay(),
                         false,
                         true);

    pModel_->createIntColumn(&var,
                             name,
                             rotation.cost_,
                             rotation.getCompactPattern(),
                             rotation.reducedCost_,
                             cons,
                             coeffs);
  }
  return var;
}

/*
 * Rotation constraints
 */
void RotationMP::buildRotationCons(const SolverParam &param) {
  char name[255];
  // build the rotation network for each nurse
  for (int i = 0; i < pScenario_->nNurses(); i++) {
    int minConsDaysOff(theLiveNurses_[i]->minConsDaysOff()),
        maxConsDaysOff(theLiveNurses_[i]->maxConsDaysOff()),
        initConsDaysOff(theLiveNurses_[i]->pStateIni_->consDaysOff_);
    // =true if we have to compute a cost for resting days exceeding the
    // maximum allowed
    // =false otherwise
    bool const maxRest = (maxConsDaysOff < nDays() + initConsDaysOff);
    // number of long resting arcs as function of maxRest
    int const nbLongRestingArcs((maxRest) ? maxConsDaysOff : minConsDaysOff);
    // first day when a rest arc exists =
    // nbLongRestingArcs - number of consecutive worked days in the past
    int const firstRestArc(std::min(
        std::max(0, nbLongRestingArcs - initConsDaysOff),
        nDays() - 1));
    // first day when a restingVar exists: at minimun 1
    // if firstRestArc=0, the first resting arc is a longRestingVar
    int const indexStartRestArc = std::max(1, firstRestArc);
    // number of resting arcs
    int const nbRestingArcs(nDays() - indexStartRestArc);

    // initialize vectors
    vector<vector<MyVar *> > restsPerDay2(nDays());
    vector<MyVar *> restingVars2(nbRestingArcs);
    vector<vector<MyVar *> > longRestingVars2(nDays());
    vector<MyCons *> restFlowCons2(nDays());
    vector<MyCons *> workFlowCons2(nDays());
    vector<MyVar *> stabRestFlowPlus2(nDays());
    vector<MyVar *> stabRestFlowMinus2(nDays());
    vector<MyVar *> stabWorkFlowPlus2(nDays());
    vector<MyVar *> stabWorkFlowMinus2(nDays());

    /*****************************************
     * Creating arcs
     *****************************************/
    for (int k = 0; k < nDays(); ++k) {
      /*****************************************
       * first long resting arcs
       *****************************************/
      if (k == 0) {
        // number of min long resting arcs
        int nbMinRestArcs(std::max(0, minConsDaysOff - initConsDaysOff));
        // initialize cost
        double cost(nbMinRestArcs * pScenario_->weights().WEIGHT_CONS_DAYS_OFF);
        RotationPattern rot = computeInitStateRotation(theLiveNurses_[i]);

        // initialize vectors
        // Must have a minimum of one long resting arcs
        vector<MyVar *> longRestingVars3_0(indexStartRestArc);

        // create minRest arcs
        for (int l = 1; l <= nbMinRestArcs; ++l) {
          cost -= pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
          snprintf(name, sizeof(name), "longRestingVars_N%d_%d_%d", i, 0, l);
          pModel_->createPositiveVar(&longRestingVars3_0[l - 1],
                                     name,
                                     cost + rot.cost_,
                                     rot.getCompactPattern());
          initialStateVars_.push_back(longRestingVars3_0[l - 1]);
          // add this resting arc for each day of rest
          for (int k1 = 0; k1 < l; ++k1)
            restsPerDay2[k1].push_back(longRestingVars3_0[l - 1]);
        }

        // create maxRest arcs, if maxRest=true
        if (maxRest) {
          for (int l = 1 + nbMinRestArcs; l <= firstRestArc; ++l) {
            snprintf(name, sizeof(name), "longRestingVars_N%d_%d_%d", i, 0, l);
            pModel_->createPositiveVar(&longRestingVars3_0[l - 1],
                                       name,
                                       rot.cost_,
                                       rot.getCompactPattern());
            initialStateVars_.push_back(longRestingVars3_0[l - 1]);
            // add this resting arc for each day of rest
            for (int k1 = 0; k1 < l; ++k1)
              restsPerDay2[k1].push_back(longRestingVars3_0[l - 1]);
          }
        }

        // create the only resting arc (same as a short resting arcs)
        if (firstRestArc == 0) {
          snprintf(name, sizeof(name), "restingVars_N%d_%d_%d", i, 0, 1);
          pModel_->createPositiveVar(&longRestingVars3_0[0],
                                     name,
                                     (maxRest) ?
                                     pScenario_->weights().WEIGHT_CONS_DAYS_OFF
                                         + rot.cost_ : rot.cost_,
                                     rot.getCompactPattern());
          initialStateVars_.push_back(longRestingVars3_0[0]);
          // add this resting arc for the first day of rest
          restsPerDay2[0].push_back(longRestingVars3_0[0]);
        }
        // store vectors
        longRestingVars2[0] = longRestingVars3_0;
      } else {
        /*****************************************
         * long resting arcs without the first ones
         *****************************************/
        // n umber of long resting arcs
        // = min(nbLongRestingArcs, number of possible long resting arcs)
        int nbLongRestingArcs2(std::min(nbLongRestingArcs, nDays() - k));
        // initialize cost
        // if the arc finishes the last day, the cost is 0.
        // Indeed it will be computed on the next planning
        double cost =
            minConsDaysOff * pScenario_->weights().WEIGHT_CONS_DAYS_OFF;

        // initialize vectors
        vector<MyVar *> longRestingVars3(nbLongRestingArcs2);

        // create minRest arcs
        for (int l = 1; l <= minConsDaysOff; ++l) {
          bool doBreak = false;
          cost -= pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
          snprintf(name,
                   sizeof(name),
                   "longRestingVars_N%d_%d_%d",
                   i,
                   k,
                   k + l);
          // if arc ends before the last day: normal cost
          if (l < nDays() - k) {
            pModel_->createPositiveVar(&longRestingVars3[l - 1], name, cost);
          } else {
            // otherwise, arc finishes on last day
            // so: cost=0 and we break the loop
            pModel_->createPositiveVar(&longRestingVars3[l - 1], name, 0);
            doBreak = true;
          }
          // add this resting arc for each day of rest
          for (int k1 = k; k1 < k + l; ++k1)
            restsPerDay2[k1].push_back(longRestingVars3[l - 1]);
          if (doBreak)
            break;
        }
        // create maxRest arcs, if maxRest=true
        if (maxRest) {
          for (int l = 1 + minConsDaysOff; l <= maxConsDaysOff; ++l) {
            // if exceed last days, break
            if (l > nDays() - k)
              break;
            snprintf(name,
                     sizeof(name),
                     "longRestingVars_N%d_%d_%d",
                     i,
                     k,
                     k + l);
            pModel_->createPositiveVar(&longRestingVars3[l - 1], name, 0);
            // add this resting arc for each day of rest
            for (int k1 = k; k1 < k + l; ++k1)
              restsPerDay2[k1].push_back(longRestingVars3[l - 1]);
          }
        }
        // store vectors
        longRestingVars2[k] = longRestingVars3;
      }
      /*****************************************
       * short resting arcs
       *****************************************/
      if (k >= indexStartRestArc) {
        snprintf(name, sizeof(name), "restingVars_N%d_%d_%d", i, k, k + 1);
        pModel_->createPositiveVar(&restingVars2[k - indexStartRestArc],
                                   name,
                                   (maxRest)
                                   ? pScenario_->weights().WEIGHT_CONS_DAYS_OFF
                                   : 0);
        // add this resting arc for this day of rest
        restsPerDay2[k].push_back(restingVars2[k - indexStartRestArc]);
      }
    }

    /*****************************************
     * Resting nodes constraints
     *****************************************/
    for (int k = 0; k < nDays(); ++k) {
      vector<double> coeffs(longRestingVars2[k].size());
      for (unsigned int l = 0; l < longRestingVars2[k].size(); ++l)
        coeffs[l] = 1;
      snprintf(name, sizeof(name), "restingNodes_N%d_%d", i, k);
      // Create flow constraints. out flow = 1 if source node (k=0)
      pModel_->createEQConsLinear(&restFlowCons2[k], name, (k == 0) ? 1 : 0,
                                  longRestingVars2[k], coeffs);

      // STAB:Add stabilization variables
      if (param.isStabilization_) {
        snprintf(name, sizeof(name), "stabRestFlow2_%i", k);
        addStabVariables(param, name, restFlowCons2[k], true, true);
      }
    }

    /*****************************************
     * Working nodes constraints
     *****************************************/
    for (int k = 1; k <= nDays(); ++k) {
      // take the min between the number of long resting arcs and
      // the number of possible in arcs
      int nbLongRestingArcs2 = std::min(nbLongRestingArcs, k);

      vector<MyVar *> vars;
      vector<double> coeffs;
      // add long resting arcs
      for (int l = 0; l < nbLongRestingArcs2; ++l) {
        // if the long resting arc starts on the source node,
        // check if there exists such an arc
        if ((l > k - 1) || (l >= longRestingVars2[k - 1 - l].size()))
          break;
        vars.push_back(longRestingVars2[k - 1 - l][l]);
        // compute in-flow for the sink
        if (k == nDays()) coeffs.push_back(1);
        else  // compute out-flow
          coeffs.push_back(-1);
      }
      // add resting arcs
      if (k == indexStartRestArc) {
        // just 1 out, if first restingVar
        // compute out-flow
        vars.push_back(restingVars2[0]);
        coeffs.push_back(1);
      } else if (k == nDays()) {
        // just 1 in, if last resting arcs
        // compute in-flow for the sink
        vars.push_back(restingVars2[restingVars2.size() - 1]);
        coeffs.push_back(1);
      } else if (k > indexStartRestArc) {
        // 2 otherwise: 1 in and 1 out
        // compute out-flow
        vars.push_back(restingVars2[k - 1 - indexStartRestArc]);
        coeffs.push_back(-1);
        vars.push_back(restingVars2[k - indexStartRestArc]);
        coeffs.push_back(1);
      }
      snprintf(name, sizeof(name), "workingNodes_N%d_%d", i, k);

      // Create flow constraints. in flow = 1 if sink node
      // (k==pDemand_->nbDays_)
      pModel_->createEQConsLinear(&workFlowCons2[k - 1],
                                  name,
                                  (k == nDays()) ? 1 : 0,
                                  vars,
                                  coeffs);

      // STAB:Add stabilization variables
      if (param.isStabilization_) {
        snprintf(name, sizeof(name), "stabWorkFlow2_%i", k - 1);
        addStabVariables(param, name, workFlowCons2[k - 1], true, true);
      }
    }

    // store vectors
    restsPerDay_[i] = restsPerDay2;
    restingVars_[i] = restingVars2;
    longRestingVars_[i] = longRestingVars2;
    restFlowCons_[i] = restFlowCons2;
    workFlowCons_[i] = workFlowCons2;
  }
}

int RotationMP::addRotationConsToCol(vector<MyCons *> *cons,
                                     vector<double> *coeffs,
                                     int i,
                                     int k,
                                     bool firstDay,
                                     bool lastDay) {
  // check if the rotation starts on day k
  if (firstDay) {
    // compute out-flow
    coeffs->push_back(1.0);
    // add to source constraint
    if (k == 0)
      cons->push_back(restFlowCons_[i][0]);
    else
      // add to work node constraint
      cons->push_back(workFlowCons_[i][k - 1]);

    return 1;
  } else if (lastDay) {
    // check if the rotation finishes on day k
    // add to sink constraint
    // compute in-flow
    if (k == nDays() - 1) {
      coeffs->push_back(1.0);
      cons->push_back(workFlowCons_[i][nDays() - 1]);
    } else {
      // add to rest node constraint
      // compute out-flow
      coeffs->push_back(-1.0);
      cons->push_back(restFlowCons_[i][k + 1]);
    }

    return 1;
  }

  return 0;
}

/*
 * Min/Max constraints
 */
void RotationMP::buildMinMaxCons(const SolverParam &param) {
  char name[255];
  for (int i = 0; i < pScenario_->nNurses(); i++) {
    /* min worked days constraint */
    snprintf(name, sizeof(name), "minWorkedDaysVar_N%d", i);
    pModel_->createPositiveVar(&minWorkedDaysVars_[i],
                               name,
                               weightTotalShiftsMin_[i]);
    snprintf(name, sizeof(name), "minWorkedDaysCons_N%d", i);
    pModel_->createGEConsLinear(&minWorkedDaysCons_[i],
                                name,
                                minTotalShifts_[i],
                                {minWorkedDaysVars_[i]},
                                {1});
    // STAB: Add stabilization variable
    if (param.isStabilization_) {
      snprintf(name, sizeof(name), "stabMinWorkedDays_%i", i);
      addStabVariables(param, name, minWorkedDaysCons_[i], false, true);
    }

    /* max worked days constraint */
    snprintf(name, sizeof(name), "maxWorkedDaysVar_N%d", i);
    pModel_->createPositiveVar(&maxWorkedDaysVars_[i],
                               name,
                               weightTotalShiftsMax_[i]);
    snprintf(name, sizeof(name), "maxWorkedDaysCons_N%d", i);
    pModel_->createLEConsLinear(&maxWorkedDaysCons_[i],
                                name,
                                maxTotalShifts_[i],
                                {maxWorkedDaysVars_[i]},
                                {-1});
    // STAB: Add stabilization variable
    if (param.isStabilization_) {
      snprintf(name, sizeof(name), "stabMaxWorkedDays_%i", i);
      addStabVariables(param, name, maxWorkedDaysCons_[i], true, false);
    }

    // add constraints on the total number of shifts to satisfy bounds that
    // correspond to the global bounds averaged over the weeks
    //
    // STAB: not implemented there yet
    if (!minTotalShiftsAvg_.empty() && !maxTotalShiftsAvg_.empty()
        && !weightTotalShiftsAvg_.empty()) {
      // only add the constraint if is tighter than the already added constraint
      if (minTotalShiftsAvg_[i] > minTotalShifts_[i]) {
        snprintf(name, sizeof(name), "minWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&minWorkedDaysAvgVars_[i],
                                   name,
                                   weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "minWorkedDaysAvgCons_N%d", i);
        pModel_->createGEConsLinear(
            &minWorkedDaysAvgCons_[i],
            name,
            minTotalShiftsAvg_[i],
            {minWorkedDaysVars_[i], minWorkedDaysAvgVars_[i]},
            {1, 1});

        isMinWorkedDaysAvgCons_[i] = true;
      }

      if (maxTotalShiftsAvg_[i] < maxTotalShifts_[i]) {
        snprintf(name, sizeof(name), "maxWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&maxWorkedDaysAvgVars_[i],
                                   name,
                                   weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "maxWorkedDaysAvgCons_N%d", i);
        pModel_->createLEConsLinear(
            &maxWorkedDaysAvgCons_[i],
            name,
            maxTotalShiftsAvg_[i],
            {maxWorkedDaysVars_[i], maxWorkedDaysAvgVars_[i]},
            {-1, -1});

        isMaxWorkedDaysAvgCons_[i] = true;
      }
    }

    snprintf(name, sizeof(name), "maxWorkedWeekendVar_N%d", i);
    pModel_->createPositiveVar(&maxWorkedWeekendVars_[i],
                               name,
                               weightTotalWeekendsMax_[i]);

    snprintf(name, sizeof(name), "maxWorkedWeekendCons_N%d", i);
    pModel_->createLEConsLinear(&maxWorkedWeekendCons_[i],
                                name,
                                maxTotalWeekends_[i],
                                {maxWorkedWeekendVars_[i]},
                                {-1});

    // STAB:Add stabilization variables
    if (param.isStabilization_) {
      snprintf(name, sizeof(name), "stabMaxWorkedWeekend_%i", i);
      addStabVariables(param, name, maxWorkedWeekendCons_[i], true, false);
    }

    // STAB: not implemented there yet
    if (!maxTotalWeekendsAvg_.empty() && !weightTotalWeekendsAvg_.empty()
        && maxTotalWeekendsAvg_[i] < theLiveNurses_[i]->maxTotalWeekends()
            - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) {
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgVar_N%d", i);
      pModel_->createPositiveVar(&maxWorkedWeekendAvgVars_[i],
                                 name,
                                 weightTotalWeekendsAvg_[i]);
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgCons_N%d", i);
      pModel_->createLEConsLinear(
          &maxWorkedWeekendAvgCons_[i],
          name,
          maxTotalWeekendsAvg_[i]
              - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
          {maxWorkedWeekendVars_[i], maxWorkedWeekendAvgVars_[i]},
          {-1, -1});

      isMaxWorkedWeekendAvgCons_[i] = true;
    }
  }

  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    if (!minTotalShiftsContractAvg_.empty()
        && !maxTotalShiftsContractAvg_.empty()
        && !weightTotalShiftsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "minWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&minWorkedDaysContractAvgVars_[p],
                                 name,
                                 weightTotalShiftsContractAvg_[p]);
      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedDaysContractAvgVars_[p],
                                 name,
                                 weightTotalShiftsContractAvg_[p]);

      snprintf(name, sizeof(name), "minWorkedDaysContractAvgCons_P%d", p);
      pModel_->createGEConsLinear(&minWorkedDaysContractAvgCons_[p],
                                  name,
                                  minTotalShiftsContractAvg_[p],
                                  {minWorkedDaysContractAvgVars_[p]},
                                  {1});

      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgCons_P%d", p);
      pModel_->createLEConsLinear(&maxWorkedDaysContractAvgCons_[p],
                                  name,
                                  maxTotalShiftsContractAvg_[p],
                                  {maxWorkedDaysContractAvgVars_[p]},
                                  {-1});

      isMinWorkedDaysContractAvgCons_[p] = true;
      isMaxWorkedDaysContractAvgCons_[p] = true;
    }

    if (!maxTotalWeekendsContractAvg_.empty()
        && !weightTotalWeekendsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedWeekendContractAvgVars_[p],
                                 name,
                                 weightTotalWeekendsContractAvg_[p]);

      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgCons_C%d", p);
      pModel_->createLEConsLinear(&maxWorkedWeekendContractAvgCons_[p],
                                  name,
                                  maxTotalWeekendsContractAvg_[p],
                                  {maxWorkedWeekendContractAvgVars_[p]},
                                  {-1});

      isMaxWorkedWeekendContractAvgCons_[p] = true;
    }
  }
}

int RotationMP::addMinMaxConsToCol(vector<MyCons *> *cons,
                                   vector<double> *coeffs,
                                   int i,
                                   int nbDays,
                                   int nbWeekends) {
  int nbCons(0);
  int p = theLiveNurses_[i]->pContract_->id_;
  ++nbCons;
  cons->push_back(minWorkedDaysCons_[i]);
  coeffs->push_back(nbDays);
  ++nbCons;
  cons->push_back(maxWorkedDaysCons_[i]);
  coeffs->push_back(nbDays);
  if (isMinWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons->push_back(minWorkedDaysAvgCons_[i]);
    coeffs->push_back(nbDays);
  }
  if (isMaxWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons->push_back(maxWorkedDaysAvgCons_[i]);
    coeffs->push_back(nbDays);
  }
  if (isMinWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons->push_back(minWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(nbDays);
  }
  if (isMaxWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons->push_back(maxWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(nbDays);
  }

  if (nbWeekends) {
    ++nbCons;
    cons->push_back(maxWorkedWeekendCons_[i]);
    coeffs->push_back(nbWeekends);

    if (isMaxWorkedWeekendAvgCons_[i]) {
      ++nbCons;
      cons->push_back(maxWorkedWeekendAvgCons_[i]);
      coeffs->push_back(nbWeekends);
    }

    if (isMaxWorkedWeekendContractAvgCons_[p]) {
      ++nbCons;
      cons->push_back(maxWorkedWeekendContractAvgCons_[p]);
      coeffs->push_back(nbWeekends);
    }
  }

  return nbCons;
}

double RotationMP::getColumnsCost(CostType costType) const {
  double cost = 0;
  if (costType == CONS_REST_COST)
    return pModel_->getTotalCost(restingVars_)
        + pModel_->getTotalCost(longRestingVars_)
            // cost for empty rotation: rotation for initial state followed by
            // rest -> already included in longRestingVars_
        - getColumnsCost(costType, initialStateVars_)
            // just initial rest costs;
        + getColumnsCost(CONS_REST_COST, pModel_->getActiveColumns());

  cost = getColumnsCost(costType, pModel_->getActiveColumns());
  if (costType == ROTATION_COST)  // add rest costs + historical costs
    cost += pModel_->getTotalCost(restingVars_)
        + pModel_->getTotalCost(longRestingVars_);
  else  // add historical non resting costs
    cost += getColumnsCost(costType, initialStateVars_);

  return cost;
}

double RotationMP::getColumnsCost(CostType costType,
                                  const vector<MyVar *> &vars) const {
  double cost = 0;
  for (MyVar *var : vars) {
    double value = pModel_->getVarValue(var);
    if (value > epsilon()) {
      RotationPattern rot(var->getPattern(), pScenario_);
      rot.computeCost(this, theLiveNurses_[rot.nurseNum_]);
      switch (costType) {
        case CONS_SHIFTS_COST: cost += rot.consShiftsCost_ * value;
          break;
        case CONS_WORK_COST: cost += rot.consDaysWorkedCost_ * value;
          break;
        case COMPLETE_WEEKEND_COST: cost += rot.completeWeekendCost_ * value;
          break;
        case PREFERENCE_COST: cost += rot.preferenceCost_ * value;
          break;
        case CONS_REST_COST: cost += rot.initRestCost_ * value;
          break;
        default: cost += rot.cost_ * value;
          break;
      }
    }
  }
  return cost;
}

double RotationMP::getDaysCost() const {
  return pModel_->getTotalCost(minWorkedDaysVars_) +
      pModel_->getTotalCost(maxWorkedDaysVars_);
}

double RotationMP::getWeekendCost() const {
  return pModel_->getTotalCost(maxWorkedWeekendVars_);
}

RotationPattern RotationMP::computeInitStateRotation(const PLiveNurse& pNurse) {
  // initialize rotation
  RotationPattern rot = RotationPattern(map<int, int>(), nullptr, pNurse->num_);

  // compute cost for previous cons worked shifts and days
  int lastShiftType = pNurse->pStateIni_->shiftType_;
  if (lastShiftType > 0) {
    int nbConsWorkedDays = pNurse->pStateIni_->consDaysWorked_;
    int diff = pNurse->minConsDaysWork() - nbConsWorkedDays;
    rot.consDaysWorkedCost_ +=
        (diff > 0) ? diff * pScenario_->weights().WEIGHT_CONS_DAYS_WORK : 0;

    int nbConsShifts = pNurse->pStateIni_->consShifts_;
    int diff2 = pScenario_->minConsShiftsOf(lastShiftType) - nbConsShifts;
    rot.consShiftsCost_ +=
        (diff2 > 0) ? diff2 * pScenario_->weights().WEIGHT_CONS_SHIFTS : 0;
  }
  rot.cost_ = rot.consDaysWorkedCost_ + rot.consShiftsCost_;

  return rot;
}
