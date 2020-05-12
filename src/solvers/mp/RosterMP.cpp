//
// Created by antoine legrain on 2020-04-22.
//

#include "RosterMP.h"
#include "solvers/mp/TreeManager.h"


//-----------------------------------------------------------------------------
//
//  S t r u c t   R o s t e r
//
//  A roster is a vector of shifts covering the whole horizon
//
//-----------------------------------------------------------------------------

void RosterPattern::computeCost(PScenario pScenario, const std::vector<PLiveNurse>& liveNurses){
  //check if pNurse points to a nurse
  if(nurseId_ == -1)
    Tools::throwError("LiveNurse = NULL");

  PLiveNurse pNurse = liveNurses[nurseId_];

  /************************************************
   * Compute all the costs of a roster:
   ************************************************/

  // initialize costs
  consShiftsCost_ = 0;
  consDaysOffCost_ = 0;
  consDaysWorkedCost_ = 0;
  completeWeekendCost_ = 0;
  preferenceCost_ = 0;
  initRestCost_ =  0;
  initWorkCost_ = 0;
  minDaysCost_ = 0;
  maxDaysCost_ = 0;
  maxWeekendCost_ = 0;

  /*
   * Compute initial cost
   */

  // initial state of the nurse
  int lastShiftType = pNurse->pStateIni_->shiftType_;
  int nbConsShifts = pNurse->pStateIni_->consShifts_;
  int nbConsDaysWorked = pNurse->pStateIni_->consDaysWorked_;
  int nbConsDaysOff = pNurse->pStateIni_->consDaysOff_;

  // 1. if the initial shift has already exceeded the max, substract now the cost that will be readd later
  if(lastShiftType) { // was working
    if (nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))
      consShiftsCost_ -= (nbConsShifts-pScenario->maxConsShiftsOf(lastShiftType))*WEIGHT_CONS_SHIFTS;
  } else if(nbConsShifts > pNurse->maxConsDaysOff()) {
    consDaysOffCost_ -= (nbConsDaysOff-pNurse->maxConsDaysOff())*WEIGHT_CONS_DAYS_OFF;
  }

  // 2. compute the cost
  bool initialCost = true;
  for(int s: shifts_){
    // a. same shift type -> increment the counters
    int shiftType = pScenario->shiftIDToShiftTypeID_[s];
    if(lastShiftType == shiftType) {
      if (shiftType > 0) {
        nbConsShifts++;
        nbConsDaysWorked++;
      } else nbConsDaysOff++;
      initialCost = false;
      continue;
    }

    // b. different shift type -> add the corresponding cost
    // i) goes to rest
    if(shiftType == 0) {
      // compute cost
      consShiftsCost_ += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
      consDaysWorkedCost_ += pNurse->consDaysCost(nbConsDaysWorked);
      if(initialCost)
        initWorkCost_ = consShiftsCost_ + consDaysWorkedCost_;
      // update counters
      nbConsShifts = 0;
      nbConsDaysWorked = 0;
      nbConsDaysOff = 1;
    }
    // ii) goes to work
    else if (lastShiftType == 0){
      // compute cost
      consDaysOffCost_ += pNurse->consDaysOffCost(nbConsDaysOff);
      if(initialCost)
        initRestCost_ = consDaysOffCost_;
      // update counters
      nbConsDaysOff = 0;
      nbConsDaysWorked = 1;
      nbConsShifts = 1;
    }
    // iii) continue to work on a different shift
    else {
      // compute cost
      consShiftsCost_ += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
      // update counters
      nbConsShifts = 1;
      nbConsDaysWorked ++;
    }

    // update
    lastShiftType = shiftType;
    initialCost = false;
  }

  // pay the max for last day
  if(nbConsDaysOff > pNurse->maxConsDaysOff())
    consDaysOffCost_ += pNurse->consDaysOffCost(nbConsDaysOff);
  if(lastShiftType > 0 && nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))
    consShiftsCost_ += pScenario->consShiftTypeCost(lastShiftType, nbConsShifts);
  if(nbConsDaysWorked > pNurse->maxConsDaysWork())
    consDaysWorkedCost_ += pNurse->consDaysCost(nbConsDaysWorked);

  /*
   * Compute preferencesCost
   */

  for(int k=0; k<length_; ++k) {
    int  level = pNurse->wishesOffLevel(k, shifts_[k]);
    if (level != -1)
      preferenceCost_ += WEIGHT_PREFERENCES_OFF[level];
  }

  for(int k=firstDay_; k<firstDay_+length_; ++k) {
    int  level = pNurse->wishesOnLevel(k, shifts_[k]);
    if (level != -1)
      preferenceCost_ += WEIGHT_PREFERENCES_ON[level];
  }

  /*
   * Compute time duration and complete weekend
   */
  timeDuration_ = 0;
  nbWeekends_ = 0;
  int k = 0;
  bool rest = false;
  for (int s: shifts_) {
    timeDuration_ += pScenario->timeDurationToWork_[s];
    if(pScenario->isWorkShift(s)) {
      if (Tools::isSaturday(k))
        nbWeekends_ ++;
      else if (rest && Tools::isSunday(k)) {
        nbWeekends_ ++;
        if (pNurse->needCompleteWeekends())
          completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
      }
      rest = false;
    } else {
      if (pNurse->needCompleteWeekends() && !rest && Tools::isSunday(k))
        completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
      rest = true;
    }
    k++;
  }

  if(pNurse->minTotalShifts() - timeDuration_ > 0)
    minDaysCost_ = WEIGHT_TOTAL_SHIFTS * (pNurse->minTotalShifts() - timeDuration_);
  if(timeDuration_ - pNurse->maxTotalShifts() > 0)
    maxDaysCost_ = WEIGHT_TOTAL_SHIFTS * (timeDuration_ - pNurse->maxTotalShifts());
  maxWeekendCost_ = pNurse->totalWeekendCost(nbWeekends_);


  if(false){
    std::cout << "# Costs:" << std::endl;
    std::cout << costsToString();
    std::cout << "# " << std::endl;
  }

  /*
   * Compute the sum of the cost and stores it in cost_
   */
  cost_ = consShiftsCost_ + consDaysOffCost_ + consDaysWorkedCost_ + completeWeekendCost_
      + preferenceCost_ + minDaysCost_ + maxDaysCost_ + maxWeekendCost_;
}


void RosterPattern::checkDualCost(DualCosts& costs, PScenario pScenario){
  //check if pNurse points to a nurse
  if(nurseId_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of a rotation:
   ************************************************/

  double dualCost(cost_);

  bool rest = false;
  dualCost -= costs.constant();
  for(int k=0; k<length_; ++k) {
    int s = shifts_[k];// if working
    if(pScenario->isWorkShift(s)) {
      /* Working dual cost */
      dualCost -= costs.workedDayShiftCost(k, s);
      /* Start working dual cost */
      if(rest) dualCost -= costs.startWorkCost(k);
      /* Working on weekend */
      if(Tools::isSaturday(k))
        dualCost -= costs.workedWeekendCost();
      else if(rest && Tools::isSunday(k))
        dualCost -= costs.workedWeekendCost();
      rest = false;
    }
    // if on rest
    else {
      /* Stop working dual cost */
      if(k>0 && !rest) dualCost -= costs.endWorkCost(k-1);
      rest = true;
    }
  }
  /* Stop working dual cost */
  if(!rest) dualCost -= costs.endWorkCost(length_-1);

  // Display: set to true if you want to display the details of the cost

  if(abs(dualCost_ - dualCost) / (1 - dualCost)  > EPSILON ){
    std::cout << "# " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "Bad dual cost: " << dualCost_ << " != " << dualCost << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "#   | Base cost     : + " << cost_ << std::endl;
    std::cout << costsToString(false);


    std::cout << "#   | Constant: - " << costs.constant() << std::endl;
    for(int k=0; k<length_; ++k)
      if(pScenario->isWorkShift(shifts_[k]))
        std::cout << "#   | Work day-shift " << k << ": - " << costs.workedDayShiftCost(k, shifts_[k]) << std::endl;
    std::cout << toString(costs.nDays(), pScenario->shiftIDToShiftTypeID_);
    std::cout << "# " << std::endl;

  }
}

std::string RosterPattern::costsToString(bool initialCosts) const {
  std::stringstream rep;
  rep << "#       | Consecutive shifts  : " << consShiftsCost_ << std::endl;
  rep << "#       | Consecutive days off: " << consDaysOffCost_ << std::endl;
  rep << "#       | Consecutive days    : " << consDaysWorkedCost_ << std::endl;
  rep << "#       | Complete weekends   : " << completeWeekendCost_ << std::endl;
  rep << "#       | Preferences         : " << preferenceCost_ << std::endl;
  if(initialCosts) {
    rep << "#       | Initial rest        : " << initRestCost_ << std::endl;
    rep << "#       | Initial work        : " << initWorkCost_ << std::endl;
  }
  rep << "#       | Min worked days     : " << minDaysCost_ << std::endl;
  rep << "#       | Max worked days     : " << maxDaysCost_ << std::endl;
  rep << "#       | Max worked weekend  : " << maxWeekendCost_ << std::endl;
  return rep.str();
}

std::string RosterPattern::toString(int nbDays, std::vector<int> shiftIDToShiftTypeID) const {
  if(nbDays == -1) nbDays = firstDay_+length_;
  std::stringstream rep;
  rep << "#   | ROTATION: N=" << nurseId_ << "  cost=" << cost_ << "  dualCost=" << dualCost_ << std::endl;
  rep << "#               |";
  for(unsigned int i=0; i<shifts_.size(); i++){
    if(shifts_[i] < 1) rep << "\t|";
    else {
      int t = shifts_[i];
      if(t < shiftIDToShiftTypeID.size()) {
        rep << shiftIDToShiftTypeID[t] << ":";
        rep << shifts_[i] << "|";
      } else rep << " " << shifts_[i] << " |";
    }
  }
  rep << std::endl;
  return rep.str();
}


//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r M P
//
// Build and solve the master problem of the column generation scheme with rosters
//
//-----------------------------------------------------------------------------

RosterMP::RosterMP(PScenario pScenario, PDemand pDemand, PPreferences pPreferences, std::vector<State> *pInitState,
    MySolverType solver):
    MasterProblem(pScenario, pDemand, pPreferences, pInitState, solver),
    assignmentCons_(pScenario->nbNurses_) {
  lagrangianBoundAvail_ = true;
}

RosterMP::~RosterMP() {}

PPattern RosterMP::getPattern(const std::vector<double>& pattern) const {
  return std::make_shared<RosterPattern>(pattern);
}

// Main method to build the rostering problem for a given input
void RosterMP::build(const SolverParam& param) {
  /* Roster assignment constraints */
  buildAssignmentCons(param);

  /* build the rest of the model */
  MasterProblem::build(param);

  /* Change the branching rule */
  pRule_ = new RosterBranchingRule(this, (RestTree*)pTree_, "branching rule");
  pModel_->addBranchingRule(pRule_);
}

// Provide an initial solution to the solver. If empty, add artificial columns
void RosterMP::initializeSolution(const std::vector<Roster>& solution) {
  // rosters are added for each nurse of the initial solution
  if (solution.size() != 0) {
    const char* baseName("initialRoster");
    //build the roster of each nurse
    for (int i = 0; i < pScenario_->nbNurses_; ++i) {
      //load the roster of nurse i
      const Roster& roster = solution[i];
      // build the shift vector
      std::vector<int> shifts(getNbDays());
      for (int k = 0; k < getNbDays(); ++k)
        shifts[k] = roster.shift(k);
      RosterPattern pat(shifts, i);
      pat.computeCost(pScenario_, theLiveNurses_);
      pModel_->addInitialColumn(addRoster(pat, baseName));
    }
  }
  // Else, add roster to always be feasible
  else {
    const char* baseName("feasibilityRoster");
    std::vector<int> shifts;
    Tools::initVector(shifts, getNbDays(), -1); // -1 = work on all shifts
    for(int i=0; i<pScenario_->nbNurses_; ++i) {
      RosterPattern pat(shifts, i, LARGE_SCORE);
      addRoster(pat, baseName, true);
    }
  }
}

//Create a new rotation variable
//add the correct constraints and coefficients for the nurse i working on a rotation
//if s=-1, the nurse works on all shifts
MyVar* RosterMP::addColumn(int nurseId, const RCSolution& solution) {
  // Build rotation from RCSolution
  RosterPattern pat(solution.shifts, nurseId, DBL_MAX, solution.cost);
  pat.computeCost(pScenario_, theLiveNurses_);
  pat.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  DualCosts costs = buildDualCosts(theLiveNurses_[nurseId]);
  pat.checkDualCost(costs, pScenario_);
  std::vector<double> pattern = pat.getCompactPattern();
  checkIfPatternAlreadyPresent(pattern);
#endif
  return addRoster(pat, "roster", false);
}

MyVar* RosterMP::addRoster(const RosterPattern& roster, const char* baseName, bool coreVar) {
  //nurse index
  int nurseId = roster.nurseId_;

  //Column var, its name, and affected constraints with their coefficients
  MyVar* var;
  char name[255];
  std::vector<MyCons*> cons;
  std::vector<double> coeffs;

  /* Skills coverage constraints */
  addSkillsCoverageConsToCol(cons, coeffs, roster);

  /* add variable to model */
  sprintf(name, "%s_N%d_%ld",baseName , nurseId, roster.id_);
  if(coreVar){
    pModel_->createPositiveVar(&var, name, roster.cost_, roster.getCompactPattern());
    for(unsigned int i=0; i<cons.size(); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  }
  else {
    /* Roster  assignment constraint s added only for real rosters
     * to be sure that the artificial variables can always be used to create a feasible solution
     */
    addRosterConsToCol(cons, coeffs, nurseId);

    bool relaxed = false;
    bool rest = true;
    for (int k = 0; k < getNbDays(); ++k) {
      int shiftType = pScenario_->shiftIDToShiftTypeID_[roster.getShift(k)];
      if (rest && isRelaxDay(shiftType)) {
        pModel_->createPositiveColumn(&var, name, roster.cost_, roster.getCompactPattern(), roster.dualCost_, cons,
                                      coeffs);
        relaxed = true;
        Tools::throwError("Behaviour nos tested for relaxDay with roster.");
        break;
      }
      rest = (shiftType== 0);
    }
    if(!relaxed)
      pModel_->createIntColumn(&var, name, roster.cost_, roster.getCompactPattern(), roster.dualCost_, cons, coeffs);
  }
  return var;
}

/* Build each set of constraints - Add also the coefficient of a column for each set */
void RosterMP::buildAssignmentCons(const SolverParam &parameters) {
  char name[255];
  //build the roster assignment constraint for each nurse
  for(int i=0; i<pScenario_->nbNurses_; i++){
    sprintf(name, "feasibilityAssignmentVar_N%d", i);
    MyVar* feasibilityVar(nullptr);
    pModel_->createPositiveVar(&feasibilityVar, name, LARGE_SCORE);
    sprintf(name, "assignmentCons_N%d", i);
    pModel_->createEQConsLinear(&assignmentCons_[i], name, 1, {feasibilityVar}, {1});
  }
}

int RosterMP::addRosterConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i) {
  cons.push_back(assignmentCons_[i]);
  coeffs.push_back(1);
  return 1;
}

//get a reference to the restsPerDay_ for a Nurse
std::vector<MyVar*> RosterMP::getRestVarsPerDay(PLiveNurse pNurse, int day) const {
  std::vector<MyVar*> restRosters;
  for(MyVar* var: pModel_->getActiveColumns()) {
    RosterPattern pat(var->getPattern());
    if(pat.nurseId_ == pNurse->id_ && pScenario_->isRestShift(pat.getShift(day)))
      restRosters.push_back(var);
  }
  return restRosters;
}

// compute the lagrangian bound
double RosterMP::computeLagrangianBound(double objVal) const {
  double sumRedCost = 0;
  for(double v: pPricer_->getLastMinOptimalReducedCost())
    sumRedCost += v;
  return objVal + sumRedCost;
}

// return the costs of all active columns associated to the type
double RosterMP::getColumnsCost(CostType costType, bool justHistoricalCosts) const {
  return getColumnsCost(costType, pModel_->getActiveColumns(), justHistoricalCosts);
}

double RosterMP::getColumnsCost(CostType costType, const std::vector<MyVar*>& vars, bool justHistoricalCosts) const {
  if(justHistoricalCosts && costType != REST_COST && costType != TOTAL_COST)
    Tools::throwError("It's not possible to retrive the historical costs for a cost type "
                      "different than REST_COST or TOTAL_COST.");

  double cost = 0;
  for(MyVar* var: vars){
    double value = pModel_->getVarValue(var);
    if(value > EPSILON){
      RosterPattern ros(var->getPattern());
      ros.computeCost(pScenario_, theLiveNurses_);
      double c = 0;
      switch(costType){
        case CONS_SHIFTS_COST: c += ros.consShiftsCost_;
          break;
        case CONS_WORKED_DAYS_COST: c += ros.consDaysWorkedCost_;
          break;
        case COMPLETE_WEEKEND_COST: c += ros.completeWeekendCost_;
          break;
        case PREFERENCE_COST: c += ros.preferenceCost_;
          break;
        case MIN_DAYS_COST: c += ros.minDaysCost_;
          break;
        case MAX_DAYS_COST: c += ros.maxDaysCost_;
          break;
        case MAX_WEEKEND_COST: c += ros.maxWeekendCost_;
          break;
        case REST_COST: c += ros.initRestCost_;
          if(!justHistoricalCosts) c += ros.consDaysOffCost_;
          break;
        default: c += ros.initWorkCost_+ros.initRestCost_;
          if(!justHistoricalCosts) c += ros.cost_;
          break;
      }
      cost += c * value;
    }
  }
  return cost;
}

double RosterMP::getMinDaysCost() const {
  return getColumnsCost(MIN_DAYS_COST, false);
}

double RosterMP::getMaxDaysCost() const {
  return getColumnsCost(MAX_DAYS_COST, false);
}

double RosterMP::getMaxWeekendCost() const {
  return getColumnsCost(MAX_WEEKEND_COST, false);
}

/* retrieve the dual values */
double RosterMP::getConstantDualvalue(PLiveNurse pNurse) const {
  return pModel_->getDual(assignmentCons_[pNurse->id_]);
}
