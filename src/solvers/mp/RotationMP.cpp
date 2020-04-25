//
// Created by antoine legrain on 2020-04-15.
//

#include "RotationMP.h"
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

void RotationPattern::computeCost(Scenario* pScenario, Preferences* pPreferences, const vector<LiveNurse*>& liveNurses, int horizon){
  //check if pNurse points to a nurse
  if(nurseId_ == -1)
    Tools::throwError("LiveNurse = NULL");

  LiveNurse* pNurse = liveNurses[nurseId_];

  // compute time duration
  computeTimeDuration(pScenario);

  /************************************************
   * Compute all the costs of a rotation:
   ************************************************/
  //   double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_ ;

  //if first day of the planning, check on the past, otherwise 0 (rest)
  int lastShiftType = (firstDay_==0) ? pNurse->pStateIni_->shiftType_ : -1;
  //nbConsShift = number of consecutive shift
  //if first day of the planning, check on the past, otherwise 0
  int nbConsShifts = (firstDay_==0) ? pNurse->pStateIni_->consShifts_ : 0;
  //consShiftCost = cost of be outside of the interval [min,max] of the consecutives shifts
  consShiftsCost_ = 0;

  //nbConsWorked = number of consecutive worked days
  //if first day of the planning, check on the past , otherwise 0
  int nbConsDaysWorked = (firstDay_==0) ? pNurse->pStateIni_->consDaysWorked_ : 0;
  //consWorkedCost = cost of be outside of the interval [min,max] of the consecutives worked days
  consDaysWorkedCost_ = 0 ;

  //cost of not doing the whole weekend
  completeWeekendCost_ = 0;

  //preferencesCost = cost of not respecting preferences
  preferenceCost_ = 0;

  //initial resting cost
  initRestCost_ =  0;

  /*
   * Compute consShiftCost
   */

  // if the initial shift has already exceeded the max, substract now the cost that will be readd later
  if( (firstDay_==0) && (lastShiftType>0) &&
      (nbConsShifts > pScenario->maxConsShiftsOf(lastShiftType))){
    consShiftsCost_ -= (nbConsShifts-pScenario->maxConsShiftsOf(lastShiftType))*WEIGHT_CONS_SHIFTS;
  }

  for(int k=firstDay_; k<firstDay_+length_; ++k){
    int  shiftType = pScenario->shiftIDToShiftTypeID_[shifts_[k]];
    if(lastShiftType == shiftType){
      nbConsShifts ++;
      continue;
    }
    if(lastShiftType > 0){
      int diff = max(pScenario->minConsShiftsOf(lastShiftType) - nbConsShifts,
                     nbConsShifts-pScenario->maxConsShiftsOf(lastShiftType));
      if(diff>0) {
        consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
      }
    }
    //initialize nbConsShifts and lastShift
    nbConsShifts = 1;
    lastShiftType = shiftType;
  }

  //compute consShiftsCost for the last shift
  if(lastShiftType>0) {
    int diff = max((firstDay_ + length_ == horizon) ? 0 : pScenario->minConsShiftsOf(lastShiftType) - nbConsShifts,
                   nbConsShifts - pScenario->maxConsShiftsOf(lastShiftType));
    if (diff > 0)
      consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
  }


  /*
   * Compute consDaysWorkedCost
   */

  // if already worked too much
  double diffDays = nbConsDaysWorked - pNurse->pContract_->maxConsDaysWork_;
  consDaysWorkedCost_ = (diffDays > 0) ? - diffDays * WEIGHT_CONS_DAYS_WORK : 0 ;

  nbConsDaysWorked += length_;
  //check if nbConsDaysWorked < min, if finishes on last day, does not count
  if(nbConsDaysWorked < pNurse->minConsDaysWork() && firstDay_+length_ < horizon)
    consDaysWorkedCost_ += (pNurse->minConsDaysWork() - nbConsDaysWorked) * WEIGHT_CONS_DAYS_WORK;
    //check if nbConsDaysWorked > max
  else if(nbConsDaysWorked > pNurse->maxConsDaysWork())
    consDaysWorkedCost_ += (nbConsDaysWorked - pNurse->maxConsDaysWork()) * WEIGHT_CONS_DAYS_WORK;

  /*
   * Compute completeWeekendCost
   */
  if(pNurse->needCompleteWeekends()){
    //if first day is a Sunday, the saturday is not worked
    if(Tools::isSunday(firstDay_))
      completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
    //if last day + 1 is a Sunday, the sunday is not worked
    if(Tools::isSunday(firstDay_+length_))
      completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
  }

  /*
   * Compute preferencesCost
   */

  for(int k=firstDay_; k<firstDay_+length_; ++k) {
    int  level = pPreferences->wantsTheShiftOffLevel(nurseId_, k, shifts_[k]);
    if (level != -1)
      preferenceCost_ += WEIGHT_PREFERENCES_OFF[level];
  }

  for(int k=firstDay_; k<firstDay_+length_; ++k) {
    int  level = pPreferences->wantsTheShiftOnLevel(nurseId_, k, shifts_[k]);
    if (level != -1)
      preferenceCost_ += WEIGHT_PREFERENCES_ON[level];
  }

  /*
   * Compute initial resting cost
   */

  if(firstDay_==0 && pNurse->pStateIni_->shiftType_==0){
    int diff = pNurse->minConsDaysOff() - pNurse->pStateIni_->consDaysOff_;
    initRestCost_ = (diff > 0) ? diff*WEIGHT_CONS_DAYS_OFF : 0;
  }

  /*
   * Compute the sum of the cost and stores it in cost_
   */
  if(false){
    cout << "# Costs:" << endl;
    cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
    cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
    cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
    cout << "#       | Preferences       : " << preferenceCost_ << endl;
    cout << "#       | Initial rest      : " << initRestCost_ << endl;
    cout << "# " << endl;
  }

  cost_ = consShiftsCost_ + consDaysWorkedCost_ + completeWeekendCost_ + preferenceCost_ +  initRestCost_;
}


void RotationPattern::checkDualCost(DualCosts& costs){
  //check if pNurse points to a nurse
  if(nurseId_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of a rotation:
   ************************************************/

  double dualCost(cost_);

  /* Working dual cost */
  for(int k=firstDay_; k<firstDay_+length_; ++k)
    dualCost -= costs.workedDayShiftCost(k, shifts_[k]);
  /* Start working dual cost */
  dualCost -= costs.startWorkCost(firstDay_);
  /* Stop working dual cost */
  dualCost -= costs.endWorkCost(firstDay_+length_-1);
  /* Working on weekend */
  if(Tools::isSunday(firstDay_))
    dualCost -= costs.workedWeekendCost();
  for(int k=firstDay_; k<firstDay_+length_; ++k)
    if(Tools::isSaturday(k))
      dualCost -= costs.workedWeekendCost();


  // Display: set to true if you want to display the details of the cost

  if(abs(dualCost_ - dualCost) / (1 - dualCost) > EPSILON ){
    cout << "# " << endl;
    cout << "# " << endl;
    cout << "Bad dual cost: " << dualCost_ << " != " << dualCost << endl;
    cout << "# " << endl;
    cout << "#   | Base cost     : + " << cost_ << endl;

    cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
    cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
    cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
    cout << "#       | Preferences       : " << preferenceCost_ << endl;
    cout << "#       | Initial rest      : " << initRestCost_ << endl;

    for(int k=firstDay_; k<firstDay_+length_; ++k)
      cout << "#   | Work day-shift: - " << costs.workedDayShiftCost(k, shifts_[k]) << endl;
    cout << "#   | Start work    : - " << costs.startWorkCost(firstDay_) << endl;
    cout << "#   | Finish Work   : - " << costs.endWorkCost(firstDay_+length_-1) << endl;
    if(Tools::isSunday(firstDay_))
      cout << "#   | Weekends      : - " << costs.workedWeekendCost() << endl;
    for(int k=firstDay_; k<firstDay_+length_; ++k)
      if(Tools::isSaturday(k))
        cout << "#   | Weekends      : - " << costs.workedWeekendCost() << endl;
//		std::cout << "#   | ROTATION:" << "  cost=" << cost_ << "  dualCost=" << dualCost_ << "  firstDay=" << firstDay_ << "  length=" << length_ << std::endl;
    std::cout << toString(costs.nDays());
    cout << "# " << endl;

  }
}

std::string RotationPattern::toString(int nbDays, std::vector<int> shiftIDToShiftTypeID) const {
  if(nbDays == -1) nbDays = firstDay_+length_;
  std::stringstream rep;
  rep << "#   | ROTATION: N=" << nurseId_ << "  cost=" << cost_ << "  dualCost=" << dualCost_ << "  firstDay=" << firstDay_ << "  length=" << length_ << std::endl;
  rep << "#               |";
  std::vector<int> allTasks (nbDays);
  for(auto itTask = shifts_.begin(); itTask != shifts_.end(); ++itTask)
    allTasks[itTask->first] = itTask->second;
  for(unsigned int i=0; i<allTasks.size(); i++){
    if(allTasks[i] < 1) rep << "\t|";
    else {
      int t = allTasks[i];
      if(t < shiftIDToShiftTypeID.size()) {
        rep << shiftIDToShiftTypeID[t] << ":";
        rep << allTasks[i] << "|";
      } else rep << " " << allTasks[i] << " |";
    }
  }
  rep << std::endl;
  return rep.str();
}


//-----------------------------------------------------------------------------
//
//  C l a s s   R o t a t i o n M P
//
// Build and solve the master problem of the column generation scheme with rotations
//
//-----------------------------------------------------------------------------

RotationMP::RotationMP(Scenario* pScenario, Demand* pDemand, Preferences* pPreferences,
    std::vector<State> *pInitState, MySolverType solver) :
     MasterProblem(pScenario, pDemand, pPreferences, pInitState, solver),
     restsPerDay_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
     minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
     minWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedWeekendAvgVars_(pScenario_->nbNurses_),
     minWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedWeekendContractAvgVars_(pScenario_->nbContracts_),

     restFlowCons_(pScenario->nbNurses_),
     workFlowCons_(pScenario->nbNurses_),
     minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
     minWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedWeekendAvgCons_(pScenario_->nbNurses_),
     minWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedWeekendContractAvgCons_(pScenario_->nbContracts_),
     // STAB
     stabRestFlowPlus_(pScenario->nbNurses_), stabRestFlowMinus_(pScenario->nbNurses_),
     stabWorkFlowPlus_(pScenario->nbNurses_), stabWorkFlowMinus_(pScenario->nbNurses_),
     stabMinWorkedDaysPlus_(pScenario->nbNurses_), stabMaxWorkedDaysMinus_(pScenario->nbNurses_), stabMaxWorkedWeekendMinus_(pScenario->nbNurses_) {

  // initialize the vectors indicating whether the min/max total constraints
  // with averaged bounds are considered
  Tools::initVector(isMinWorkedDaysAvgCons_, pScenario_->nbNurses_, false);
  Tools::initVector(isMaxWorkedDaysAvgCons_, pScenario_->nbNurses_, false);
  Tools::initVector(isMaxWorkedWeekendAvgCons_, pScenario_->nbNurses_, false);
  Tools::initVector(isMinWorkedDaysContractAvgCons_, pScenario_->nbContracts_, false);
  Tools::initVector(isMaxWorkedDaysContractAvgCons_, pScenario_->nbContracts_, false);
  Tools::initVector(isMaxWorkedWeekendContractAvgCons_, pScenario_->nbContracts_, false);
}

RotationMP::~RotationMP() {}

PPattern RotationMP::getPattern(const std::vector<double>& pattern) const {
  return std::make_shared<RotationPattern>(pattern);
}

//build the rostering problem
void RotationMP::build(const SolverParam& param){
  /* Rotation constraints */
  buildRotationCons(param);

  /* Min/Max constraints */
  buildMinMaxCons(param);

  /* build the rest of the model */
  MasterProblem::build(param);
}

//initialize the rostering problem with one column to be feasible if there is no initial solution
//otherwise build the columns corresponding to the initial solution
void RotationMP::initializeSolution(const vector<Roster>& solution) {
  string baseName("initialRotation");
  //rotations are added for each nurse of the initial solution
  if (solution.size() != 0) {
    //build the rotations of each nurse
    for (int i = 0; i < pScenario_->nbNurses_; ++i) {
      //load the roster of nurse i
      Roster roster = solution[i];

      bool workedLastDay = false;
      int lastShift = 0;
      map<int, int> shifts;
      //build all the successive rotation of this nurse
      for (int k = 0; k < getNbDays(); ++k) {
        //shift=0 => rest
        int shift = roster.shift(k);
        //if work, insert the shift in the map
        if (shift > 0) {
          shifts.insert(pair<int, int>(k, shift));
          lastShift = shift;
          workedLastDay = true;
        } else if (shift < 0 && lastShift > 0) {
          shifts.insert(pair<int, int>(k, lastShift));
          workedLastDay = true;
        }
          //if stop to work, build the rotation
        else if (workedLastDay) {
          RotationPattern rotation(shifts, i);
          rotation.computeCost(pScenario_, pPreferences_, theLiveNurses_, getNbDays());
          pModel_->addActiveColumn(addRotation(rotation, baseName.c_str()));
          shifts.clear();
          lastShift = shift;
          workedLastDay = false;
        }
      }
      //if work on the last day, build the rotation
      if (workedLastDay) {
        RotationPattern rotation(shifts, i);
        rotation.computeCost(pScenario_, pPreferences_, theLiveNurses_, getNbDays());
        pModel_->addActiveColumn(addRotation(rotation, baseName.c_str()));
        shifts.clear();
      }
    }
  } else {
    /* We add initial rotations to be always feasible */
    std::string baseName("feasibilityRotation");
    //We add a column with 1 everywhere for each nurse to be always feasible
    //build a map of shift -1 everywhere
    map<int,int> shifts;
    for(int k=0; k<getNbDays(); ++k)
      shifts.insert(pair<int,int>( k , -1 ));

    for(int i=0; i<pScenario_->nbNurses_; ++i){
      // DBG: Compute the cost of artificial variables in accordance to the soft
      // constraints
      double artificialCost = WEIGHT_TOTAL_SHIFTS*pScenario_->nbShifts_*(pScenario_->nbDays()-pScenario_->maxTotalShiftsOf(i));
      artificialCost += WEIGHT_CONS_DAYS_WORK*pScenario_->nbShifts_*(pScenario_->nbDays()-pScenario_->maxConsDaysWorkOf(i));
      for (int s = 1; s < pScenario_->nbShifts_; s++) {
        artificialCost += WEIGHT_CONS_SHIFTS*(pScenario_->nbDays()-pScenario_->maxConsShiftsOfTypeOf(s));
      }
      RotationPattern rotation(shifts, i, LARGE_SCORE);// artificialCost);//
      addRotation(rotation, baseName.c_str(), true);
    }
  }
}

vector2D<double> RotationMP::getShiftsDualValues(LiveNurse*  pNurse) const {
  vector2D<double> dualValues = MasterProblem::getShiftsDualValues(pNurse);

  int i = pNurse->id_;
  int p = pNurse->pContract_->id_;

  /* Min/Max constraints */
  double minWorkedDays = pModel_->getDual(minWorkedDaysCons_[i], true);
  double maxWorkedDays = pModel_->getDual(maxWorkedDaysCons_[i], true);

  double minWorkedDaysAvg = isMinWorkedDaysAvgCons_[i] ? pModel_->getDual(minWorkedDaysAvgCons_[i], true):0.0;
  double maxWorkedDaysAvg = isMaxWorkedDaysAvgCons_[i] ? pModel_->getDual(maxWorkedDaysAvgCons_[i], true):0.0;

  double minWorkedDaysContractAvg = isMinWorkedDaysContractAvgCons_[p] ?
                                    pModel_->getDual(minWorkedDaysContractAvgCons_[p], true):0.0;
  double maxWorkedDaysContractAvg = isMaxWorkedDaysContractAvgCons_[p] ?
                                    pModel_->getDual(maxWorkedDaysContractAvgCons_[p], true):0.0;

  /* Min/Max constraints */
  double d = minWorkedDays + minWorkedDaysAvg + minWorkedDaysContractAvg;
  d += maxWorkedDays + maxWorkedDaysAvg + maxWorkedDaysContractAvg;


  for(int k=0; k<getNbDays(); ++k){
    vector<double>& dualValues2 = dualValues[k];
    for(int s=1; s<pScenario_->nbShifts_; ++s)
      // adjust the dual in function of the time duration of the shift
      dualValues2[s-1] += d * pScenario_->timeDurationToWork_[s];
  }

  return dualValues;
}

vector<double> RotationMP::getStartWorkDualValues(LiveNurse* pNurse) const {
  int i = pNurse->id_;
  vector<double> dualValues(getNbDays());

  //get dual value associated to the source
  dualValues[0] =  pModel_->getDual(restFlowCons_[i][0], true);
  //get dual values associated to the work flow constraints
  //don't take into account the last which is the sink
  for(int k=1; k<getNbDays(); ++k)
    dualValues[k] = pModel_->getDual(workFlowCons_[i][k-1], true);

  return dualValues;
}

vector<double> RotationMP::getEndWorkDualValues(LiveNurse* pNurse) const {
  int i = pNurse->id_;
  vector<double> dualValues(getNbDays());

  //get dual values associated to the work flow constraints
  //don't take into account the first which is the source
  //take into account the cost, if the last day worked is k
  for(int k=0; k<getNbDays()-1; ++k)
    dualValues[k] = -pModel_->getDual(restFlowCons_[i][k+1], true);

  //get dual value associated to the sink
  dualValues[getNbDays()-1] =
      pModel_->getDual(workFlowCons_[i][getNbDays()-1], true);

  return dualValues;
}

double RotationMP::getWorkedWeekendDualValue(LiveNurse* pNurse) const{
  int id = pNurse->id_;
  double dualVal = pModel_->getDual(maxWorkedWeekendCons_[id], true);
  if (isMaxWorkedWeekendAvgCons_[id]) {
    dualVal += pModel_->getDual(maxWorkedWeekendAvgCons_[id], true);
  }
  if (isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_]) {
    dualVal += pModel_->getDual(
        maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
  }

  return dualVal;
}

//------------------------------------------------------------------------------
// Build the variable of the rotation as well as all the affected constraints
// with their coefficients. if s=-1, the nurse i works on all shifts
//------------------------------------------------------------------------------
MyVar* RotationMP::addColumn(int nurseId, const RCSolution& solution) {
  // Build rotation from RCSolution
  RotationPattern rotation(solution.firstDay, solution.shifts, nurseId, DBL_MAX, solution.cost);
  rotation.computeCost(pScenario_, pPreferences_, theLiveNurses_, getNbDays());
  rotation.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  DualCosts costs = buildDualCosts(theLiveNurses_[nurseId]);
  rotation.checkDualCost(costs);
  std::vector<double> pattern = rotation.getCompactPattern();
  checkIfPatternAlreadyPresent(pattern);
#endif
  return addRotation(rotation, "rotation", false);
}

MyVar* RotationMP::addRotation(const RotationPattern& rotation, const char* baseName, bool coreVar){
  //nurse index
  int nurseId = rotation.nurseId_;

  //Column var, its name, and affected constraints with their coefficients
  MyVar* var;
  char name[255];
  vector<MyCons*> cons;
  vector<double> coeffs;

  /* Min/Max constraints */
  int nbWeekends = Tools::containsWeekend(rotation.firstDay_, rotation.firstDay_+rotation.length_-1);
  //addMinMaxConsToCol(cons, coeffs, nurseId, rotation.length_, nbWeekends);
  addMinMaxConsToCol(cons, coeffs, nurseId, rotation.timeDuration_, nbWeekends);   // pour prendre en compte les heures plutôt que les jours

  /* Skills coverage constraints */
  addSkillsCoverageConsToCol(cons, coeffs, rotation);

  sprintf(name, "%s_N%d_%ld",baseName , nurseId, rotation.id_);
  if(coreVar){
    // DBG
    // The artificial variables are taken out of the flow constraints to
    // allow them to be used even after branching
    // Otherwise, the problem might mistakenly appear infeasible after
    // branching
    // addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_, true, false);
    // addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_+rotation.length_-1, false, true);

    pModel_->createPositiveVar(&var, name, rotation.cost_, rotation.getCompactPattern());
    for(unsigned int i=0; i<cons.size(); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  }
  else {
    /* Rotation constraints
      They are added only for real rotations to be sure that the artificial variables
      can always be used to create a feasible solution
     */
    addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_, true, false);
    addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_+rotation.length_-1, false, true);

    if (this->isRelaxDay(rotation.firstDay_)) {
      pModel_->createPositiveColumn(&var, name, rotation.cost_, rotation.getCompactPattern(), rotation.dualCost_, cons, coeffs);
    }
    else {
      pModel_->createIntColumn(&var, name, rotation.cost_, rotation.getCompactPattern(), rotation.dualCost_, cons, coeffs);
    }
  }
  return var;
}

/*
 * Rotation constraints
 */
void RotationMP::buildRotationCons(const SolverParam& param){
  char name[255];
  //build the rotation network for each nurse
  for(int i=0; i<pScenario_->nbNurses_; i++){
    int minConsDaysOff(theLiveNurses_[i]->minConsDaysOff()),
        maxConsDaysOff(theLiveNurses_[i]->maxConsDaysOff()),
        initConsDaysOff(theLiveNurses_[i]->pStateIni_->consDaysOff_);
    //=true if we have to compute a cost for resting days exceeding the maximum allowed
    //=false otherwise
    bool const maxRest = (maxConsDaysOff < getNbDays() + initConsDaysOff);
    //number of long resting arcs as function of maxRest
    int const nbLongRestingArcs((maxRest) ? maxConsDaysOff : minConsDaysOff);
    //first day when a rest arc exists =
    //nbLongRestingArcs - number of consecutive worked days in the past
    int const firstRestArc( std::min( std::max( 0, nbLongRestingArcs - initConsDaysOff ), getNbDays()-1 ) );
    //first day when a restingVar exists: at minimun 1
    //if firstRestArc=0, the first resting arc is a longRestingVar
    int const indexStartRestArc = std::max(1, firstRestArc);
    //number of resting arcs
    int const nbRestingArcs( getNbDays()- indexStartRestArc );

    //initialize vectors
    vector< vector< MyVar* > > restsPerDay2(getNbDays());
    vector< MyVar* > restingVars2(nbRestingArcs);
    vector< vector<MyVar*> > longRestingVars2(getNbDays());
    vector<MyCons*> restFlowCons2(getNbDays());
    vector<MyCons*> workFlowCons2(getNbDays());
    vector<MyVar*> stabRestFlowPlus2(getNbDays());
    vector<MyVar*> stabRestFlowMinus2(getNbDays());
    vector<MyVar*> stabWorkFlowPlus2(getNbDays());
    vector<MyVar*> stabWorkFlowMinus2(getNbDays());

    /*****************************************
     * Creating arcs
     *****************************************/
    for(int k=0; k<getNbDays(); ++k){
      /*****************************************
       * first long resting arcs
       *****************************************/
      if(k==0){
        //number of min long resting arcs
        int nbMinRestArcs( std::max(0, minConsDaysOff - initConsDaysOff) );
        //initialize cost
        int cost (nbMinRestArcs * WEIGHT_CONS_DAYS_OFF);
        RotationPattern rot = computeInitStateRotation(theLiveNurses_[i]);

        //initialize vectors
        //Must have a minimum of one long resting arcs
        vector<MyVar*> longRestingVars3_0(indexStartRestArc);

        //create minRest arcs
        for(int l=1; l<=nbMinRestArcs; ++l){
          cost -= WEIGHT_CONS_DAYS_OFF;
          sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
          pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, cost+rot.cost_, rot.getCompactPattern());
          initialStateVars_.push_back(longRestingVars3_0[l-1]);
          //add this resting arc for each day of rest
          for(int k1=0; k1<l; ++k1)
            restsPerDay2[k1].push_back(longRestingVars3_0[l-1]);
        }

        //create maxRest arcs, if maxRest=true
        if(maxRest){
          for(int l=1+nbMinRestArcs; l<=firstRestArc; ++l){
            sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
            pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, rot.cost_, rot.getCompactPattern());
            initialStateVars_.push_back(longRestingVars3_0[l-1]);
            //add this resting arc for each day of rest
            for(int k1=0; k1<l; ++k1)
              restsPerDay2[k1].push_back(longRestingVars3_0[l-1]);
          }
        }

        //create the only resting arc (same as a short resting arcs)
        if(firstRestArc == 0){
          sprintf(name, "restingVars_N%d_%d_%d", i, 0, 1);
          pModel_->createPositiveVar(&longRestingVars3_0[0], name,
              (maxRest) ? WEIGHT_CONS_DAYS_OFF+rot.cost_ : rot.cost_, rot.getCompactPattern());
          initialStateVars_.push_back(longRestingVars3_0[0]);
          //add this resting arc for the first day of rest
          restsPerDay2[0].push_back(longRestingVars3_0[0]);
        }
        //store vectors
        longRestingVars2[0] = longRestingVars3_0;
      }
        /*****************************************
         * long resting arcs without the first ones
         *****************************************/
      else{
        //number of long resting arcs = min(nbLongRestingArcs, number of possible long resting arcs)
        int nbLongRestingArcs2( std::min(nbLongRestingArcs, getNbDays()-k) );
        //initialize cost
        //if the arc finishes the last day, the cost is 0. Indeed it will be computed on the next planning
        int cost = minConsDaysOff * WEIGHT_CONS_DAYS_OFF;

        //initialize vectors
        vector<MyVar*> longRestingVars3(nbLongRestingArcs2);

        //create minRest arcs
        for(int l=1; l<=minConsDaysOff; ++l){
          bool doBreak = false;
          cost -= WEIGHT_CONS_DAYS_OFF;
          sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
          //if arc ends before the last day: normal cost
          if(l < getNbDays()-k)
            pModel_->createPositiveVar(&longRestingVars3[l-1], name, cost);
            //otherwise, arc finishes on last day
            //so: cost=0 and we break the loop
          else{
            pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
            doBreak = true;
          }
          //add this resting arc for each day of rest
          for(int k1=k; k1<k+l; ++k1)
            restsPerDay2[k1].push_back(longRestingVars3[l-1]);
          if(doBreak)
            break;
        }
        //create maxRest arcs, if maxRest=true
        if(maxRest)
          for(int l=1+minConsDaysOff; l<=maxConsDaysOff; ++l){
            //if exceed last days, break
            if(l > getNbDays()-k)
              break;
            sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
            pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
            //add this resting arc for each day of rest
            for(int k1=k; k1<k+l; ++k1)
              restsPerDay2[k1].push_back(longRestingVars3[l-1]);
          }
        //store vectors
        longRestingVars2[k] = longRestingVars3;
      }
      /*****************************************
       * short resting arcs
       *****************************************/
      if(k>=indexStartRestArc){
        sprintf(name, "restingVars_N%d_%d_%d", i, k, k+1);
        pModel_->createPositiveVar(&restingVars2[k-indexStartRestArc], name, (maxRest) ? WEIGHT_CONS_DAYS_OFF : 0);
        //add this resting arc for this day of rest
        restsPerDay2[k].push_back(restingVars2[k-indexStartRestArc]);
      }
    }

    /*****************************************
     * Resting nodes constraints
     *****************************************/
    for(int k=0; k<getNbDays(); ++k){
      vector<double> coeffs(longRestingVars2[k].size());
      for(unsigned int l=0; l<longRestingVars2[k].size(); ++l)
        coeffs[l] = 1;
      sprintf(name, "restingNodes_N%d_%d", i, k);
      //Create flow constraints. out flow = 1 if source node (k=0)
      pModel_->createEQConsLinear(&restFlowCons2[k], name, (k==0) ? 1 : 0,
                                  longRestingVars2[k], coeffs);

      // STAB:Add stabilization variables
      //
      if (param.isStabilization_)	{
        sprintf(name,"stabRestFlowPlus2_%i",k);
        pModel_->createPositiveVar(&stabRestFlowPlus2[k],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
        sprintf(name,"stabRestFlowMinus2_%i",k);
        pModel_->createPositiveVar(&stabRestFlowMinus2[k],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
        pModel_->addCoefLinear(restFlowCons2[k],stabRestFlowPlus2[k],1.0);
        pModel_->addCoefLinear(restFlowCons2[k],stabRestFlowMinus2[k],-1.0);
      }

    }

    /*****************************************
     * Working nodes constraints
     *****************************************/
    for(int k=1; k<=getNbDays(); ++k){
      //take the min between the number of long resting arcs and the number of possible in arcs
      int nbLongRestingArcs2 = std::min(nbLongRestingArcs,k);

      vector<MyVar*> vars;
      vector<double> coeffs;
      //add long resting arcs
      for(int l=0; l<nbLongRestingArcs2; ++l){
        //if the long resting arc starts on the source node,
        //check if there exists such an arc
        if( (l > k-1) || (l >= (int) longRestingVars2[k-1-l].size()) )
          break;
        vars.push_back(longRestingVars2[k-1-l][l]);
        //compute in-flow for the sink
        if(k==getNbDays())
          coeffs.push_back(1);
          //compute out-flow
        else
          coeffs.push_back(-1);
      }
      //add resting arcs
      //just 1 out, if first restingVar
      //compute out-flow
      if(k==indexStartRestArc){
        vars.push_back(restingVars2[0]);
        coeffs.push_back(1);
      }
        //just 1 in, if last resting arcs
        //compute in-flow for the sink
      else if (k==getNbDays()){
        vars.push_back(restingVars2[restingVars2.size()-1]);
        coeffs.push_back(1);
      }
        //2 otherwise: 1 in and 1 out
        //compute out-flow
      else if(k>indexStartRestArc){
        vars.push_back(restingVars2[k-1-indexStartRestArc]);
        coeffs.push_back(-1);
        vars.push_back(restingVars2[k-indexStartRestArc]);
        coeffs.push_back(1);
      }
      sprintf(name, "workingNodes_N%d_%d", i, k);
      //Create flow constraints. in flow = 1 if sink node (k==getNbDays())
      pModel_->createEQConsLinear(&workFlowCons2[k-1], name, (k==getNbDays()) ? 1 : 0,
                                  vars, coeffs);

      // STAB:Add stabilization variables
      //
      if (param.isStabilization_) {
        sprintf(name,"stabWorkFlowPlus2_%i",k-1);
        pModel_->createPositiveVar(&stabWorkFlowPlus2[k-1],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
        sprintf(name,"stabWorkFlowMinus2_%i",k-1);
        pModel_->createPositiveVar(&stabWorkFlowMinus2[k-1],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
        pModel_->addCoefLinear(workFlowCons2[k-1],stabWorkFlowPlus2[k-1],1.0);
        pModel_->addCoefLinear(workFlowCons2[k-1],stabWorkFlowMinus2[k-1],-1.0);
      }
    }

    //store vectors
    restsPerDay_[i] = restsPerDay2;
    restingVars_[i] = restingVars2;
    longRestingVars_[i] = longRestingVars2;
    restFlowCons_[i] = restFlowCons2;
    workFlowCons_[i] = workFlowCons2;

    // STAB
    stabRestFlowPlus_[i] = stabRestFlowPlus2;
    stabRestFlowMinus_[i] = stabRestFlowMinus2;
    stabWorkFlowPlus_[i] = stabWorkFlowPlus2;
    stabWorkFlowMinus_[i] = stabWorkFlowMinus2;
  }
}

int RotationMP::addRotationConsToCol(vector<MyCons*>& cons, vector<double>& coeffs, int i, int k, bool firstDay, bool lastDay){
  //check if the rotation starts on day k
  if(firstDay){
    //compute out-flow
    coeffs.push_back(1.0);
    //add to source constraint
    if(k==0)
      cons.push_back(restFlowCons_[i][0]);
      //add to work node constraint
    else
      cons.push_back(workFlowCons_[i][k-1]);

    return 1;
  }

    //check if the rotation finishes on day k
  else if(lastDay){
    //add to sink constraint
    //compute in-flow
    if(k==getNbDays()-1){
      coeffs.push_back(1.0);
      cons.push_back(workFlowCons_[i][getNbDays()-1]);
    }
      //add to rest node constraint
      //compute out-flow
    else{
      coeffs.push_back(-1.0);
      cons.push_back(restFlowCons_[i][k+1]);
    }

    return 1;
  }

  return 0;
}

/*
 * Min/Max constraints
 */
void RotationMP::buildMinMaxCons(const SolverParam& param){
  char name[255];
  for(int i=0; i<pScenario_->nbNurses_; i++){
    sprintf(name, "minWorkedDaysVar_N%d", i);
    pModel_->createPositiveVar(&minWorkedDaysVars_[i], name, weightTotalShiftsMin_[i]);
    sprintf(name, "maxWorkedDaysVar_N%d", i);
    pModel_->createPositiveVar(&maxWorkedDaysVars_[i], name, weightTotalShiftsMax_[i]);

    sprintf(name, "minWorkedDaysCons_N%d", i);
    vector<MyVar*> vars1 = {minWorkedDaysVars_[i]};
    vector<double> coeffs1 = {1};
    pModel_->createGEConsLinear(&minWorkedDaysCons_[i], name, minTotalShifts_[i], vars1, coeffs1);

    // STAB:Add stabilization variable
    //
    if (param.isStabilization_) {
      sprintf(name,"stabMinWorkedDaysPlus_%i",i);
      pModel_->createPositiveVar(&stabMinWorkedDaysPlus_[i],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
      pModel_->addCoefLinear(minWorkedDaysCons_[i],stabMinWorkedDaysPlus_[i],1.0);
    }

    sprintf(name, "maxWorkedDaysCons_N%d", i);
    vector<MyVar*> vars2 = {maxWorkedDaysVars_[i]};
    vector<double> coeffs2 = {-1};
    pModel_->createLEConsLinear(&maxWorkedDaysCons_[i], name, maxTotalShifts_[i], vars2, coeffs2);

    // STAB:Add stabilization variable
    //
    if (param.isStabilization_) {
      sprintf(name,"stabMaxWorkedDaysMinus_%i",i);
      pModel_->createPositiveVar(&stabMaxWorkedDaysMinus_[i],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
      pModel_->addCoefLinear(maxWorkedDaysCons_[i],stabMaxWorkedDaysMinus_[i],-1.0);
    }

    // add constraints on the total number of shifts to satisfy bounds that
    // correspond to the global bounds averaged over the weeks
    //
    // STAB: not implemented there yet
    if (!minTotalShiftsAvg_.empty() && !maxTotalShiftsAvg_.empty() && !weightTotalShiftsAvg_.empty()) {

      // only add the constraint if is tighter than the already added constraint
      if (minTotalShiftsAvg_[i] > minTotalShifts_[i]) {
        sprintf(name, "minWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&minWorkedDaysAvgVars_[i], name, weightTotalShiftsAvg_[i]);

        sprintf(name, "minWorkedDaysAvgCons_N%d", i);
        vector<MyVar*> varsAvg1 = {minWorkedDaysVars_[i], minWorkedDaysAvgVars_[i]};
        vector<double> coeffsAvg1 = {1,1};
        pModel_->createGEConsLinear(&minWorkedDaysAvgCons_[i], name, minTotalShiftsAvg_[i], varsAvg1, coeffsAvg1);

        isMinWorkedDaysAvgCons_[i] = true;
      }

      if (maxTotalShiftsAvg_[i] < maxTotalShifts_[i]) {
        sprintf(name, "maxWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&maxWorkedDaysAvgVars_[i], name, weightTotalShiftsAvg_[i]);

        sprintf(name, "maxWorkedDaysAvgCons_N%d", i);
        vector<MyVar*> varsAvg2 = {maxWorkedDaysVars_[i],maxWorkedDaysAvgVars_[i]};
        vector<double> coeffsAvg2 = {-1,-1};
        pModel_->createLEConsLinear(&maxWorkedDaysAvgCons_[i], name, maxTotalShiftsAvg_[i], varsAvg2, coeffsAvg2);

        isMaxWorkedDaysAvgCons_[i] = true;
      }
    }

    sprintf(name, "maxWorkedWeekendVar_N%d", i);
    pModel_->createPositiveVar(&maxWorkedWeekendVars_[i], name, weightTotalWeekendsMax_[i]);

    sprintf(name, "maxWorkedWeekendCons_N%d", i);
    vector<MyVar*> vars3 = {maxWorkedWeekendVars_[i]};
    vector<double> coeffs3 = {-1};
    pModel_->createLEConsLinear(&maxWorkedWeekendCons_[i], name, maxTotalWeekends_[i],
                                vars3, coeffs3);

    // STAB:Add stabilization variable
    //
    if (param.isStabilization_) {
      sprintf(name,"stabMaxWorkedWeekendMinus_%i",i);
      pModel_->createPositiveVar(&stabMaxWorkedWeekendMinus_[i],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
      pModel_->addCoefLinear(maxWorkedWeekendCons_[i],stabMaxWorkedWeekendMinus_[i],-1.0);
    }

    // STAB: not implemented there yet
    if ( !maxTotalWeekendsAvg_.empty()  && !weightTotalWeekendsAvg_.empty()
         && maxTotalWeekendsAvg_[i] < theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) {

      sprintf(name, "maxWorkedWeekendAvgVar_N%d", i);
      pModel_->createPositiveVar(&maxWorkedWeekendAvgVars_[i], name, weightTotalWeekendsAvg_[i]);

      sprintf(name, "maxWorkedWeekendAvgCons_N%d", i);
      vector<MyVar*> varsAvg3 = {maxWorkedWeekendVars_[i],maxWorkedWeekendAvgVars_[i]};
      vector<double> coeffsAvg3 = {-1,-1};
      pModel_->createLEConsLinear(&maxWorkedWeekendAvgCons_[i], name, maxTotalWeekendsAvg_[i]- theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
                                  varsAvg3, coeffsAvg3);

      isMaxWorkedWeekendAvgCons_[i] = true;

    }

  }

  // WEEKEND CUTS
  // sprintf(name, "sumMaxWorkedWeekendCons");
  // std::vector<pVar> varsSum3;
  // std::vector<double> coeffsSum3;
  // for(int i=0; i<pScenario_->nbNurses_; i++) {
  // 	varsSum3.push_back(maxWorkedWeekendVars_[i]);
  // 	coeffsSum3.push_back(1);
  // }
  // pModel_->createGEConsLinear(&sumMaxWorkedWeekendCons_, name, 4, varsSum3, coeffsSum3);



  for(int p=0; p<pScenario_->nbContracts_; ++p){

    if(!minTotalShiftsContractAvg_.empty() && !maxTotalShiftsContractAvg_.empty()  && !weightTotalShiftsContractAvg_.empty()){
      sprintf(name, "minWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&minWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);
      sprintf(name, "maxWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);

      sprintf(name, "minWorkedDaysContractAvgCons_P%d", p);
      vector<MyVar*> vars1 = {minWorkedDaysContractAvgVars_[p]};
      vector<double> coeffs1 = {1};
      pModel_->createGEConsLinear(&minWorkedDaysContractAvgCons_[p], name, minTotalShiftsContractAvg_[p], vars1, coeffs1);

      sprintf(name, "maxWorkedDaysContractAvgCons_P%d", p);
      vector<MyVar*> vars2 = {maxWorkedDaysContractAvgVars_[p]};
      vector<double> coeffs2 = {-1};
      pModel_->createLEConsLinear(&maxWorkedDaysContractAvgCons_[p], name, maxTotalShiftsContractAvg_[p], vars2, coeffs2);

      isMinWorkedDaysContractAvgCons_[p] = true;
      isMaxWorkedDaysContractAvgCons_[p] = true;
    }

    if(!maxTotalWeekendsContractAvg_.empty()  && !weightTotalWeekendsContractAvg_.empty()){
      sprintf(name, "maxWorkedWeekendContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedWeekendContractAvgVars_[p], name, weightTotalWeekendsContractAvg_[p]);

      sprintf(name, "maxWorkedWeekendContractAvgCons_C%d", p);
      vector<MyVar*> varsAvg3 = {maxWorkedWeekendContractAvgVars_[p]};
      vector<double> coeffsAvg3 = {-1 };
      pModel_->createLEConsLinear(&maxWorkedWeekendContractAvgCons_[p], name, maxTotalWeekendsContractAvg_[p],
                                  varsAvg3, coeffsAvg3);

      isMaxWorkedWeekendContractAvgCons_[p] = true;
    }
  }
}

int RotationMP::addMinMaxConsToCol(vector<MyCons*>& cons, vector<double>& coeffs, int i, int nbDays, int nbWeekends){
  int nbCons(0);
  int p = theLiveNurses_[i]->pContract_->id_;
  ++nbCons;
  cons.push_back(minWorkedDaysCons_[i]);
  coeffs.push_back(nbDays);
  ++nbCons;
  cons.push_back(maxWorkedDaysCons_[i]);
  coeffs.push_back(nbDays);
  if (isMinWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons.push_back(minWorkedDaysAvgCons_[i]);
    coeffs.push_back(nbDays);
  }
  if (isMaxWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons.push_back(maxWorkedDaysAvgCons_[i]);
    coeffs.push_back(nbDays);
  }
  if (isMinWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons.push_back(minWorkedDaysContractAvgCons_[p]);
    coeffs.push_back(nbDays);
  }
  if (isMaxWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons.push_back(maxWorkedDaysContractAvgCons_[p]);
    coeffs.push_back(nbDays);
  }


  if(nbWeekends){
    ++nbCons;
    cons.push_back(maxWorkedWeekendCons_[i]);
    coeffs.push_back(nbWeekends);

    if (isMaxWorkedWeekendAvgCons_[i]) {
      ++nbCons;
      cons.push_back(maxWorkedWeekendAvgCons_[i]);
      coeffs.push_back(nbWeekends);
    }

    if (isMaxWorkedWeekendContractAvgCons_[p]) {
      ++nbCons;
      cons.push_back(maxWorkedWeekendContractAvgCons_[p]);
      coeffs.push_back(nbWeekends);
    }
  }



  return nbCons;
}

double RotationMP::getColumnsCost(CostType costType, bool justHistoricalCosts) const {
  double cost = 0;
  if(justHistoricalCosts) {
    if(costType == REST_COST || costType == TOTAL_COST)
      cost += getColumnsCost(REST_COST, pModel_->getActiveColumns()); // just initial rest costs
    if(costType != REST_COST)
      // cost for empty rotation: rotation for initial state followed by rest
      // -> already included in longRestingVars_
      cost += getColumnsCost(costType, initialStateVars_);
    return cost;
  }

  if(costType == REST_COST)
    return pModel_->getTotalCost(restingVars_)
            + pModel_->getTotalCost(longRestingVars_)
            // cost for empty rotation: rotation for initial state followed by rest
            // -> already included in longRestingVars_
            - getColumnsCost(costType, initialStateVars_)
            + getColumnsCost(REST_COST, pModel_->getActiveColumns()); // just initial rest costs;

  cost = getColumnsCost(costType, pModel_->getActiveColumns());
  if(costType == TOTAL_COST) // add rest costs + historical costs
    cost += pModel_->getTotalCost(restingVars_) + pModel_->getTotalCost(longRestingVars_);
  else // add historical non resting costs
    cost += getColumnsCost(costType, initialStateVars_);

  return cost;
}

double RotationMP::getColumnsCost(CostType costType, const vector<MyVar*>& vars) const {
  double cost = 0;
  for(MyVar* var: vars){
    double value = pModel_->getVarValue(var);
    if(value > EPSILON){
      RotationPattern rot(var->getPattern());
      rot.computeCost(pScenario_, pPreferences_, theLiveNurses_, getNbDays());
      switch(costType){
        case CONS_SHIFTS_COST: cost += rot.consShiftsCost_*value;
          break;
        case CONS_WORKED_DAYS_COST: cost += rot.consDaysWorkedCost_*value;
          break;
        case COMPLETE_WEEKEND_COST: cost += rot.completeWeekendCost_*value;
          break;
        case PREFERENCE_COST: cost += rot.preferenceCost_*value;
          break;
        case REST_COST: cost += rot.initRestCost_*value;
          break;
        default: cost += rot.cost_*value;
          //            if(!initStateRotation && rot.second.length_>0){
          //               rot.second.toString(getNbDays());
          //               pModel_->toString(rot.first);
          //            }
          break;
      }
    }
  }
  return cost;
}

double RotationMP::getMinDaysCost() const {
  return pModel_->getTotalCost(minWorkedDaysVars_);
}

double RotationMP::getMaxDaysCost() const {
  return pModel_->getTotalCost(maxWorkedDaysVars_);
}

double RotationMP::getMaxWeekendCost() const {
  return pModel_->getTotalCost(maxWorkedWeekendVars_);
}

RotationPattern RotationMP::computeInitStateRotation(LiveNurse* pNurse){
  //initialize rotation
  RotationPattern rot = RotationPattern(map<int,int>(), pNurse->id_);

  //compute cost for previous cons worked shifts and days
  int lastShiftType = pNurse->pStateIni_->shiftType_;
  if(lastShiftType>0){
    int nbConsWorkedDays = pNurse->pStateIni_->consDaysWorked_;
    int diff = pNurse->minConsDaysWork() - nbConsWorkedDays;
    rot.consDaysWorkedCost_ += (diff>0) ? diff*WEIGHT_CONS_DAYS_WORK : 0;

    int nbConsShifts = pNurse->pStateIni_->consShifts_;
    int diff2 = pScenario_->minConsShiftsOf(lastShiftType) - nbConsShifts;
    rot.consShiftsCost_ += (diff2>0) ? diff2*WEIGHT_CONS_SHIFTS : 0;
  }
  rot.cost_ = rot.consDaysWorkedCost_ + rot.consShiftsCost_;

  return rot;
}

// STAB
// Update all the upper bounds of the stabilization variables by multiplying
// them by an input factor
void RotationMP::stabUpdateBound(OsiSolverInterface* solver, double factor) {
  MasterProblem::stabUpdateBound(solver, factor);
  for(int i=0; i<pScenario_->nbNurses_; i++){
    // stabilization variables corresponding to the flow constraints
    for(int k=0; k<getNbDays(); ++k) {
      multiplyUbInSolver(stabRestFlowMinus_[i][k], solver, factor);
      multiplyUbInSolver(stabRestFlowPlus_[i][k], solver, factor);
      multiplyUbInSolver(stabWorkFlowMinus_[i][k], solver, factor);
      multiplyUbInSolver(stabWorkFlowPlus_[i][k], solver, factor);
    }

    // stabilization variables corresponding to the global constraints of
    // of the nurses
    multiplyUbInSolver(stabMinWorkedDaysPlus_[i], solver, factor);
    multiplyUbInSolver(stabMaxWorkedDaysMinus_[i], solver, factor);
    multiplyUbInSolver(stabMaxWorkedWeekendMinus_[i], solver, factor);
  }
}

// STAB
// Update all the costs of the stabilization variables to the values
// corresponding dual variables with a small margin in input
void RotationMP::stabUpdateCost(OsiSolverInterface* solver, double margin) {
  MasterProblem::stabUpdateCost(solver, margin);
  for(int i=0; i<pScenario_->nbNurses_; i++){
    // stabilization variables corresponding to the flow constraints
    for(int k=0; k<getNbDays(); ++k) {
      double restFlowDual = pModel_->getDual(restFlowCons_[i][k], true);
      double workFlowDual = pModel_->getDual(workFlowCons_[i][k], true);
      updateVarCostInSolver(stabRestFlowMinus_[i][k], solver, -restFlowDual+margin);
      updateVarCostInSolver(stabRestFlowPlus_[i][k], solver, restFlowDual+margin);
      updateVarCostInSolver(stabWorkFlowMinus_[i][k], solver, -workFlowDual+margin);
      updateVarCostInSolver(stabWorkFlowPlus_[i][k], solver, workFlowDual+margin);
    }

    // stabilization variables corresponding to the global constraints of
    // of the nurses
    double minWorkedDaysDual = pModel_->getDual(minWorkedDaysCons_[i], true);
    double maxWorkedDaysDual = pModel_->getDual(maxWorkedDaysCons_[i], true);
    double maxWorkedWeekendDual = pModel_->getDual(maxWorkedWeekendCons_[i], true);
    updateVarCostInSolver(stabMinWorkedDaysPlus_[i], solver, minWorkedDaysDual+margin);
    updateVarCostInSolver(stabMaxWorkedDaysMinus_[i], solver, -maxWorkedDaysDual+margin);
    updateVarCostInSolver(stabMaxWorkedWeekendMinus_[i], solver, -maxWorkedWeekendDual+margin);
  }
}

// STAB
// Check the stopping criterion of the relaxation solution specific to the
// the stabilization
// The point is that current solution can be infeasible if  stabilization
// variables are non zero
bool RotationMP::stabCheckStoppingCriterion() const {
  if (!pModel_->getParameters().isStabilization_)
    return true;

  for(int i=0; i<pScenario_->nbNurses_; i++)
    // stabilization variables corresponding to the flow constraints
    for(int k=0; k<getNbDays(); ++k)
      if (pModel_->getVarValue(stabRestFlowMinus_[i][k]) > EPSILON ||
          pModel_->getVarValue(stabRestFlowPlus_[i][k]) > EPSILON ||
          pModel_->getVarValue(stabWorkFlowMinus_[i][k]) > EPSILON ||
          pModel_->getVarValue(stabWorkFlowPlus_[i][k]) > EPSILON)
        return false;

  for(int i=0; i<pScenario_->nbNurses_; i++)
    // stabilization variables corresponding to the global constraints of
    // of the nurses
    if (pModel_->getVarValue(stabMinWorkedDaysPlus_[i]) > EPSILON ||
        pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]) >EPSILON ||
        pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]) > EPSILON)
      return false;

  return MasterProblem::stabCheckStoppingCriterion();
}

// STAB: compute the lagrangian bound
//
double RotationMP::computeLagrangianBound(double objVal,double sumRedCost) const {
  double bound = MasterProblem::computeLagrangianBound(objVal, sumRedCost);
  double stabSumCostValue = 0.0;
  if (pModel_->getParameters().isStabilization_) {
    for(int i=0; i<pScenario_->nbNurses_; i++){
      // stabilization variables corresponding to the flow constraints
      for(int k=0; k<getNbDays(); ++k) {
        stabSumCostValue += stabRestFlowPlus_[i][k]->getCost()*pModel_->getVarValue(stabRestFlowPlus_[i][k]);
        stabSumCostValue += stabRestFlowMinus_[i][k]->getCost()*pModel_->getVarValue(stabRestFlowMinus_[i][k]);
        stabSumCostValue += stabWorkFlowPlus_[i][k]->getCost()*pModel_->getVarValue(stabWorkFlowPlus_[i][k]);
        stabSumCostValue += stabWorkFlowMinus_[i][k]->getCost()*pModel_->getVarValue(stabWorkFlowMinus_[i][k]);
      }

      // stabilization variables corresponding to the global constraints of
      // of the nurses
      stabSumCostValue += stabMinWorkedDaysPlus_[i]->getCost()*pModel_->getVarValue(stabMinWorkedDaysPlus_[i]);
      stabSumCostValue += stabMaxWorkedDaysMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]);
      stabSumCostValue += stabMaxWorkedWeekendMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]);
    }
  }
  return bound-stabSumCostValue;
}

// STAB: reset the costs and bounds of the stabilization variables
//
void RotationMP::stabResetBoundAndCost(OsiSolverInterface* solver, const SolverParam& param) {
  MasterProblem::stabResetBoundAndCost(solver, param);

  for(int i=0; i<pScenario_->nbNurses_; i++){
    // stabilization variables corresponding to the flow constraints
    for(int k=0; k<getNbDays(); ++k) {
      updateVarCostInSolver(stabRestFlowMinus_[i][k], solver, -param.stabCostIni_+param.stabCostMargin_);
      updateVarCostInSolver(stabRestFlowPlus_[i][k], solver, param.stabCostIni_+param.stabCostMargin_);
      updateVarCostInSolver(stabWorkFlowMinus_[i][k], solver, -param.stabCostIni_+param.stabCostMargin_);
      updateVarCostInSolver(stabWorkFlowPlus_[i][k], solver, param.stabCostIni_+param.stabCostMargin_);
      updateVarUbInSolver(stabRestFlowMinus_[i][k], solver, param.stabBoundIni_);
      updateVarUbInSolver(stabRestFlowPlus_[i][k], solver, param.stabBoundIni_);
      updateVarUbInSolver(stabWorkFlowMinus_[i][k], solver, param.stabBoundIni_);
      updateVarUbInSolver(stabWorkFlowPlus_[i][k], solver, param.stabBoundIni_);
    }

    // stabilization variables corresponding to the global constraints of
    // of the nurses
    updateVarCostInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabCostIni_+param.stabCostMargin_);
    updateVarCostInSolver(stabMaxWorkedDaysMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
    updateVarCostInSolver(stabMaxWorkedWeekendMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
    updateVarUbInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabBoundIni_);
    updateVarUbInSolver(stabMaxWorkedDaysMinus_[i], solver, param.stabBoundIni_);
    updateVarUbInSolver(stabMaxWorkedWeekendMinus_[i], solver, param.stabBoundIni_);
  }
}