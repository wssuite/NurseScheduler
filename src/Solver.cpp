/*
* Solver.cpp
*
*  Created on: 22 d������������������c. 2014
*      Author: jeremy
*/


#include <math.h>
#include <cmath>
#include "Solver.h"


//-----------------------------------------------------------------------------
//
//  C l a s s   S t a t N u r s e C t
//
// The instances of this class gather the status of the constraints that relate
// to the nurses.
//
//-----------------------------------------------------------------------------

// Constructor and destructor
//
StatCtNurse::StatCtNurse() {}
StatCtNurse::~StatCtNurse() {}

// initialize the statuses
//
void StatCtNurse::init(int nbDays) {
  nbDays_ = nbDays;

  // initialize all the cost vectors
  Tools::initVector(&costConsShifts_, nbDays_);
  Tools::initVector(&costConsDays_, nbDays_);
  Tools::initVector(&costConsDaysOff_, nbDays_);
  Tools::initVector(&costPref_, nbDays_);
  Tools::initVector(&costWeekEnd_, nbDays_);

  // initialize the violation vector
  for (int day = 0; day < nbDays_; day++) violSuccShifts_.push_back(false);
  for (int day = 0; day < nbDays_; day++) violSkill_.push_back(false);
}


//-----------------------------------------------------------------------------
//
//  C l a s s   L i v e N u r s e
//
// A live nurse is a nurse whose characteristics can evolve depending on
// the demand and on the planning that is being built
// They are needed in the solvers to duplicate the static nurses and define new
// attribute that can be modified.
//
//-----------------------------------------------------------------------------



// Constructor
//
LiveNurse::LiveNurse(const Nurse& nurse, Scenario* pScenario, int nbDays, int firstDay,
State* pStateIni,	map<int,set<int> >* pWishesOff):
Nurse(nurse.id_, nurse.name_, nurse.nbSkills_, nurse.skills_, nurse.pContract_),
pScenario_(pScenario), nbDays_(nbDays), firstDay_(firstDay),
pStateIni_(pStateIni), pWishesOff_(pWishesOff), pPosition_(0),
minWorkDaysNoPenaltyConsDays_(-1), maxWorkDaysNoPenaltyConsDays_(-1),
minWorkDaysNoPenaltyTotalDays_(-1), maxWorkDaysNoPenaltyTotalDays_(-1),
minAvgWorkDaysNoPenaltyTotalDays_(-1), maxAvgWorkDaysNoPenaltyTotalDays_(-1) {

  roster_.init(nbDays, firstDay);
  statCt_.init(nbDays);

  // initialize the states at each day
  states_.push_back(*pStateIni);
  for (int day = 0; day < nbDays_; day++) {
    State nextState;
    nextState.addDayToState(states_[day], 0);
    states_.push_back(nextState);
  }
}

LiveNurse::~LiveNurse() { }

// returns true if the nurse wishes the day-shift off
//
bool LiveNurse::wishesOff(int day, int shift) const {
  map<int,set<int> >::iterator itM = pWishesOff_->find(day);
  // If the day is not in the wish-list, no possible violation
  if(itM == pWishesOff_->end())  {
    return false;
  }
  // no preference either in the wish-list for that day
  else if(itM->second.find(shift) == itM->second.end()) {
    return false;
  }
  else {
    return true;
  }
}

// returns true if the nurses reached the maximum number of consecutive worked
// days or is resting and did not reach the minimum number of resting days yet
// if consecutive number of shifts will only be reached by violating maximum
// number of worked days, go to rest only if consecutive working days penalty
// is the larger
//
bool LiveNurse::needRest(int day) {
  State state = states_[day];

  if (state.shift_>0 && state.consDaysWorked_ >= maxConsDaysWork()) {
    if (state.consShifts_ < pScenario_->minConsShifts_[state.shift_]) {
      return WEIGHT_CONS_SHIFTS > WEIGHT_CONS_DAYS_WORK ? false:true;
    }
    return true;
  }
  if (state.shift_==0 && state.consDaysOff_ <= minConsDaysOff()) {
    return true;
  }
  return false;
}

// returns true if the nurse needs to work one more day to reach the minimum
// number of consecutive working days or consecutive shifts
// if consecutive number of shifts will only be reached by violating maximum
// number of worked days, go to work only if consecutive shift penalty is
// the larger
bool LiveNurse::needWork(int day) {
  State state = states_[day];

  if (state.shift_>0) {
    if (state.consDaysWorked_ < minConsDaysWork()) {
      return true;
    }
    else if (state.consShifts_ < pScenario_->minConsShifts_[state.shift_]) {
      if (state.consDaysWorked_ >= maxConsDaysWork()) {
        return WEIGHT_CONS_SHIFTS > WEIGHT_CONS_DAYS_WORK ? true:false;
      }
      return true;
    }
  }

  return false;
}

// return true if the nurse is free to go to rest or work more without penalty
//==13776==    by 0x46A780: Solver::Solver(Scenario*, Demand*, Preferences*, std::vector<State, std::allocator<State> >*) (Solver.cpp:327)

bool LiveNurse::isFreeToChoose(int day) {
  return !needWork(day) && !needRest(day);
}

// check the satisfaction of the hard constraints and record the violations
// check the soft constraints and record the costs of the violations and the
// remaining margin for the satisfied ones.
//
void LiveNurse::checkConstraints(const Roster& roster,
  const vector<State>& states, StatCtNurse& stat) {
  // check the satisfaction of the hard constraints and record the violations
  //
  for (int day = 0; day < nbDays_; day++) {

    // Check that the nurse has the assigned skill
    //
    if (roster.shift(day)) {
      stat.violSkill_[day] = !hasSkill(roster.skill(day));
    }
    else {
      stat.violSkill_[day] = false;
    }

    // Check the forbidden successor constraint
    //
    int lastShift = states[day].shift_;   // last shift assigned to the nurse
    int thisShift = roster.shift(day);    // shift assigned on this day

    stat.violSuccShifts_[day] = pScenario_->isForbiddenSuccessor(thisShift,lastShift);
  }

  // check the soft constraints and record the costs of the violations and the
  // remaining margin for the satisfied ones.
  //
  for (int day = 1; day <= nbDays_; day++) {

    // shift assigned on the previous day
    int shift = states[day].shift_;
    int prevShift = states[day-1].shift_;

    // first look at consecutive working days or days off
    //
    int missingDays=0, extraDays=0;
    stat.costConsDays_[day-1] = 0;
    stat.costConsDaysOff_[day-1] = 0;

    // compute the violations of consecutive working days an
    if (shift) {
      if (prevShift == 0) {
        missingDays = minConsDaysOff()-states[day-1].consDaysOff_;
      }

      stat.costConsDaysOff_[day-1] += (missingDays>0) ? WEIGHT_CONS_DAYS_OFF*missingDays:0;
      stat.costConsDays_[day-1] += (states[day].consDaysWorked_>maxConsDaysWork()) ? WEIGHT_CONS_DAYS_WORK:0;
    }
    else {
      if (prevShift > 0) {
        missingDays =minConsDaysWork()-states[day-1].consDaysWorked_;
      }
      extraDays = states[day].consDaysOff_-maxConsDaysOff();

      stat.costConsDays_[day-1] += (missingDays>0) ? WEIGHT_CONS_DAYS_WORK*missingDays:0;
      stat.costConsDaysOff_[day-1] += (extraDays>0) ? WEIGHT_CONS_DAYS_OFF:0;
    }

    // check the consecutive same shifts
    //
    stat.costConsShifts_[day-1] = 0;
    int missingShifts = 0;

    // count the penalty for minimum consecutive shifts only for the previous day
    // when the new shift is different
    if (shift != prevShift && prevShift > 0)  {
      missingShifts = pScenario_->minConsShifts_[prevShift]-states[day-1].consShifts_;
      stat.costConsShifts_[day-1] += (missingShifts>0) ? WEIGHT_CONS_SHIFTS*missingShifts:0;
    }

    // count the penalty for maximum consecutive shifts when the shift is worked
    // the last day will then be counted
    if (shift > 0) {
      stat.costConsShifts_[day-1] +=
        (states[day].consShifts_>pScenario_->maxConsShifts_[shift]) ? WEIGHT_CONS_SHIFTS:0;
    }

    // check the preferences
    //
    map<int,set<int> >::iterator itM = pWishesOff_->find(day-1);
    // If the day is not in the wish-list, no possible violation
    if(itM == pWishesOff_->end())  {
      stat.costPref_[day-1] = 0;
    }
    // no preference either in the wish-list for that day
    else if(itM->second.find(shift) == itM->second.end()) {
      stat.costPref_[day-1] = 0;
    }
    else {
      stat.costPref_[day-1] = WEIGHT_PREFERENCES;
    }

    // check the complete week-end (only if the nurse requires them)
    // this cost is only assigned to the sundays
    //
    if ( Tools::isSunday(day-1) && needCompleteWeekends()) {
      if ( (shift > 0 && prevShift == 0) || ( shift == 0 && prevShift > 0 )) {
        stat.costWeekEnd_[day-1] = WEIGHT_COMPLETE_WEEKEND;
      }
    }

  } // end for day

  // get the costs due to total number of working days and week-ends
  //
  stat.costTotalDays_ = 0;
  stat.costTotalWeekEnds_ = 0;
  if (true) {//pScenario_->thisWeek() == pScenario_->nbWeeks_) {
    int missingDays=0, extraDays=0;
    missingDays = std::max(0, minTotalShifts() - states[nbDays_].totalDaysWorked_);
    extraDays = std::max(0, states[nbDays_].totalDaysWorked_-maxTotalShifts());
    stat.costTotalDays_ = WEIGHT_TOTAL_SHIFTS*(extraDays+missingDays);

    int extraWeekEnds = 0;
    extraWeekEnds = std::max(0, states[nbDays_].totalWeekendsWorked_-maxTotalWeekends());
    stat.costTotalWeekEnds_ = WEIGHT_TOTAL_WEEKENDS * extraWeekEnds;
  }
}

// Build States from the roster
//
void LiveNurse::buildStates(){
   for(int k=1; k<states_.size(); ++k)
      states_[k].addDayToState(states_[k-1], roster_.shift(k-1));
}

//-----------------------------------------------------------------------------
// Compute the maximum and minimum number of working days from the input
// current state and in the next nbDays without getting any penalty for
// consecutive working days/days-off
//-----------------------------------------------------------------------------

void LiveNurse::computeMinMaxDaysNoPenaltyConsDay(State* pCurrentState, int nbDays,
  int &minWorkDaysNoPenaltyConsDays, int &maxWorkDaysNoPenaltyConsDays) {
  minWorkDaysNoPenaltyConsDays = 0;
  maxWorkDaysNoPenaltyConsDays = 0;

  // first treat the history
  //
  // remaining number of days when building a maximum and minimum working days
  // periods
  int nbDaysMin = nbDays, nbDaysMax = nbDays;
  if (pCurrentState->shift_ > 0) {
    // if the last shift was not a rest
    // keep working until maximum consecutive number of shifts for maximum
    // working days or until minimum consecutive number of shifts for minimum
    // working days
    maxWorkDaysNoPenaltyConsDays = std::min(nbDays,std::max(0, maxConsDaysWork()-pCurrentState->consDaysWorked_));
    minWorkDaysNoPenaltyConsDays = std::min(nbDays,std::max(0, minConsDaysWork()-pCurrentState->consDaysWorked_));

    // perform minimum rest for maximum working days and maximum rest for
    // minimum working days
    nbDaysMax = nbDays - std::min(nbDays,maxWorkDaysNoPenaltyConsDays+minConsDaysOff());
    nbDaysMin = nbDays - std::min(nbDays,minWorkDaysNoPenaltyConsDays+maxConsDaysOff());
  }
  else if (pCurrentState->shift_ == 0) {
    // if the last shift was a rest
    // keep resting until minimum consecutive number of rests for maximum
    // working days or until maximum consecutive number of rests for minimum
    // working days
    nbDaysMax = nbDays - std::min(nbDays,std::max(0, minConsDaysOff()-pCurrentState->consDaysOff_));
    nbDaysMin = nbDays - std::min(nbDays,std::max(0, maxConsDaysOff()-pCurrentState->consDaysOff_));
  }

  // lengths of the stints maximizing and minimizing the percentage of working
  // days respectively
  //
  int lengthStintMax = maxConsDaysWork() + minConsDaysOff();
  int lengthStintMin = minConsDaysWork() + maxConsDaysOff();

  // perform as many of these stints as possible
  //
  maxWorkDaysNoPenaltyConsDays += (nbDaysMax/lengthStintMax)*maxConsDaysWork()
    + std::min(nbDaysMax%lengthStintMax, maxConsDaysWork());
  minWorkDaysNoPenaltyConsDays += (nbDaysMin/lengthStintMin)*minConsDaysWork()
    + std::min(nbDaysMin%lengthStintMax, minConsDaysWork());
}


//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r
//
//  Solves the offline problem
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

// Specific constructor
Solver::Solver(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState):
  pScenario_(pScenario),  pDemand_(pDemand),
  pPreferences_(pPreferences), pInitState_(pInitState),
  rdm_(Tools::getANewRandomGenerator()),
  totalCostUnderStaffing_(-1), maxTotalStaffNoPenalty_(-1),
  isPreprocessedSkills_(false), isPreprocessedNurses_(false), status_(UNSOLVED) {

  // initialize the preprocessed data of the skills
  for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
    maxStaffPerSkillNoPenalty_.push_back(-1.0);
    maxStaffPerSkillAvgWork_.push_back(-1.0);
    skillRarity_.push_back(1.0);
  }

  // copy the nurses in the live nurses vector
  for (int i = 0; i < pScenario_->nbNurses_; i++) {
    theLiveNurses_.push_back(
      new LiveNurse( (pScenario_->theNurses_[i]), pScenario_, pDemand_->nbDays_,
      pDemand_->firstDay_, &(*pInitState_)[i], &(pPreferences_->wishesOff_[i])  ) );
  }

  // initialize the minimum and maximum number of total working days
  for (int i = 0; i < pScenario_->nbNurses(); i++) {
     //defaault min and max
    minTotalShifts_.push_back(theLiveNurses_[i]->minTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_);
    maxTotalShifts_.push_back(theLiveNurses_[i]->maxTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_);
    maxTotalWeekends_.push_back(theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_);

    //compute global penalties
    //default penalties for the moment
    weightTotalShiftsMin_.push_back(WEIGHT_TOTAL_SHIFTS);
    weightTotalShiftsMax_.push_back(WEIGHT_TOTAL_SHIFTS);
    weightTotalWeekendsMax_.push_back(WEIGHT_TOTAL_WEEKENDS);
  }
}


// Destructor
Solver::~Solver(){
   for(LiveNurse* pNurse: theLiveNurses_)
      delete pNurse;
}

// Load a solution in the solver
//
void Solver::loadSolution(vector<Roster> &solution) {

  if (solution.size() != pScenario_->nbNurses_) {
    Tools::throwError("Solver::loadSolution: there is not one roster per nurse in the input solution!");
  }

  solution_ = solution;
  for (int n = 0; n < pScenario_->nbNurses_; n++) {
    LiveNurse* pNurse = theLiveNurses_[n];
    pNurse->roster(solution[n]);
    pNurse->buildStates();
  }
}


//------------------------------------------------------------------------
// Preprocess the nurses
// go through the nurses to collect data regarding the potential shift and
// skill coverage of the nurses
//-------------------------------------------------------------------------

void Solver::preprocessTheNurses() {
  // local variables for conciseness of the code
  //
  int nbDays = pDemand_->nbDays_;

  this->specifyNursePositions();
  if (!isPreprocessedSkills_) this->preprocessTheSkills();
  this->computeMinMaxDaysNoPenaltyTotalDays();
  this->computeMinMaxDaysNoPenaltyConsDays();


  maxTotalStaffNoPenalty_ = 0;
  maxTotalStaffAvgWork_ = 0;
  for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
    maxStaffPerSkillNoPenalty_[sk] = 0;
    maxStaffPerSkillAvgWork_[sk] = 0;
  }
  for (LiveNurse* pNurse: theLiveNurses_)	{

    // add the working maximum number of working days to the maximum staffing
    //
    maxTotalStaffNoPenalty_ +=
      std::min(pNurse->maxWorkDaysNoPenaltyConsDays_, pNurse->maxWorkDaysNoPenaltyTotalDays_);
    maxTotalStaffAvgWork_ += pNurse->maxAvgWorkDaysNoPenaltyTotalDays_;

    // the staffing per skill is very rough here since Nurses can have
    // multiple skills
    // we thus weight the covering power of each nurse for each skill according
    // to the rarities of the skills
    double totalRarity = 0;
    for (int sk:pNurse->skills_) {
      totalRarity += skillRarity_[sk];
    }
    double maxWorkNoPenalty = (double)std::min(pNurse->maxWorkDaysNoPenaltyConsDays_, pNurse->maxWorkDaysNoPenaltyTotalDays_);
    for (int sk: pNurse->skills_) {
      maxStaffPerSkillNoPenalty_[sk] += (skillRarity_[sk]/totalRarity)*maxWorkNoPenalty;
      maxStaffPerSkillAvgWork_[sk] +=  (skillRarity_[sk]/totalRarity)*pNurse->maxAvgWorkDaysNoPenaltyTotalDays_;
    }
  }


  // initialize to zero the satisfied demand
  //
  Tools::initVector3D(&satisfiedDemand_, nbDays,pScenario_->nbShifts_, pScenario_->nbSkills_);

  isPreprocessedNurses_ = true;
}

// Find the position of each nurse
//
void Solver::specifyNursePositions() {

  for (LiveNurse* pNurse: theLiveNurses_) {
    // the skills of the nurse need to be compared to the skills of each
    // existing position to determine the position of the nurse
    bool isPosition = true;
    for (int i = 0; i < pScenario_->nbPositions() ; i++)	{
      Position* pPosition = pScenario_->pPositions()[i];
      isPosition = true;
      if (pPosition->nbSkills_ == pNurse->nbSkills_) {
        for (int i = 0; i < pNurse->nbSkills_; i++) {
          if (pNurse->skills_[i] != pPosition->skills_[i])	{
            isPosition = false;
            break;
          }
        }
      }
      else isPosition = false;

      if (isPosition) {
        pNurse->pPosition_ = pPosition;
        break;
      }
    }
    if (!isPosition) {
      Tools::throwError("The nurse has no position!");
    }
  }
}

// Compute the maximum and minimum number of working days in the period of
// the demand without getting any penalty for the total number of working days
//
void Solver::computeMinMaxDaysNoPenaltyTotalDays() {

  // number of days that will be covered after the current demand
  int nbDaysFuture = 7*(pScenario_->nbWeeks()-pScenario_->thisWeek())-pDemand_->nbDays_;

  // For each contract, compute the maximum and minimum number of working days
  // that can be done after the current demand without ensuing penalties due to
  // the min/max numbers of working days and days off
  //
  map<string,int> minWorkDaysFutureNoPenaltyConsDays;
  map<string,int> maxWorkDaysFutureNoPenaltyConsDays;
  map<string, Contract*>::const_iterator itC;
  for (itC = pScenario_->contracts_.begin(); itC != pScenario_->contracts_.end(); itC++)	{

    // initialization
    string name = itC->first;
    Contract* pContract = itC->second;
    minWorkDaysFutureNoPenaltyConsDays[name] = 0;
    maxWorkDaysFutureNoPenaltyConsDays[name] = 0;

    // we compute the minimum number of working days by assuming that the
    // first next week can start with a rest and we compute the maximum number
    // of working days by assuming that the first next week can start with work
    //

    // Compute the length of the stint that minimimizes/maximizes the number of
    // working days and  that respects the min/max number of consecutive working
    // days and days off
    // (a stint is a work rotation followed by a break)
    int lengthStintMin = pContract->minConsDaysWork_ + pContract->maxConsDaysOff_;
    int lengthStintMax = pContract->maxConsDaysWork_ + pContract->minConsDaysOff_;

    // In each case, the nurse would perform the min/max stint until there is
    // not enough remaining days to complete a stint. It will then start with
    // rest for the min stints or with work for the max stint
    minWorkDaysFutureNoPenaltyConsDays[name] +=
      (nbDaysFuture/lengthStintMin)*pContract->minConsDaysWork_
      + std::max((nbDaysFuture%lengthStintMin)-pContract->maxConsDaysOff_,0);
    maxWorkDaysFutureNoPenaltyConsDays[name] +=
      (nbDaysFuture/lengthStintMax)*pContract->maxConsDaysWork_
      + std::min((nbDaysFuture%lengthStintMax),pContract->maxConsDaysWork_);
  }

  // For each nurse, deduce the interval [dayMin,dayMax] defined as:
  // if the nurse works d days, with d in [dayMin,dayMax], then there exists a
  // sequence of future stints  such that the nurse respects the bounds on
  // the total number of working days without penalty
  for (LiveNurse* pNurse:theLiveNurses_) {
    pNurse->minWorkDaysNoPenaltyTotalDays_ = std::max(0,
      pNurse->minTotalShifts()- pNurse->totalDaysWorked()
      - maxWorkDaysFutureNoPenaltyConsDays[pNurse->contractName()]);
    pNurse->maxWorkDaysNoPenaltyTotalDays_ = std::max(0,
      pNurse->maxTotalShifts()- pNurse->totalDaysWorked()
      - minWorkDaysFutureNoPenaltyConsDays[pNurse->contractName()]);
  }

  // Get the interval [avgMin,avgMax] defined as:
  // if the nurse works d days, with d in [avgMin,avgMax], then if it works
  // the same number of days every week, the nurse will respect the bounds on
  // the total number of working days without penalty
  for (LiveNurse* pNurse:theLiveNurses_) {
    double demandMinAvgPerWeek = pNurse->minTotalShifts()/(double) pScenario_->nbWeeks();
    double demandMaxAvgPerWeek = pNurse->maxTotalShifts()/(double) pScenario_->nbWeeks();
    double demandMinAvgUntilThisDemand = demandMinAvgPerWeek*(pScenario_->thisWeek()+pDemand_->nbDays_/7.0);
    double demandMaxAvgUntilThisDemand = demandMaxAvgPerWeek*(pScenario_->thisWeek()+pDemand_->nbDays_/7.0);

    pNurse->minAvgWorkDaysNoPenaltyTotalDays_ = demandMinAvgUntilThisDemand-pNurse->totalDaysWorked();
    pNurse->maxAvgWorkDaysNoPenaltyTotalDays_ = std::max(0.0,
      demandMaxAvgUntilThisDemand-pNurse->totalDaysWorked());
  }
}

// compute the maximum and minimum number of working days in the period of
// the demand without getting any penalty for the number of consecutive
// shifts
// RqJO: this neglects the constraint of complete week-ends and the
// preferences ; they should be added later
//
void Solver::computeMinMaxDaysNoPenaltyConsDays() {
  for (LiveNurse* pNurse: theLiveNurses_)	{
    pNurse->computeMinMaxDaysNoPenaltyConsDay(pNurse->pStateIni_, pDemand_->nbDays_,
      pNurse->minWorkDaysNoPenaltyConsDays_,pNurse->maxWorkDaysNoPenaltyConsDays_);
  }
}

//------------------------------------------------------------------------
// Preprocees the skills to get their rarity
// the value depends on the minimum demand for this skill, on the number
// of nurses that have the skill and on the number of skills per nurse
// that have the skill
//------------------------------------------------------------------------

void Solver::preprocessTheSkills() {

  // this vector will contain for each skill: the number of weighted nurses
  // that have the skill ; the weight of each nurse is its number of skills
  vector<double> nbNursesWeighted;

  for (int sk=0; sk < pScenario_->nbSkills_; sk++) {
    nbNursesWeighted.push_back(0.0);
    for (LiveNurse* pNurse : theLiveNurses_) {
      if (pNurse->hasSkill(sk)) {
        nbNursesWeighted[sk]+=pNurse->maxTotalShifts()/pow((double)pNurse->nbSkills_,2);
      }
    }
    // the skill rarity is the ratio of the the demand for the skill to the
    // weighted number of nurses that have the skill
    skillRarity_[sk] = (double)pDemand_->minPerSkill_[sk]/nbNursesWeighted[sk];
  }

  // update the rarities of the skills in the scenario
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    pScenario_->pPosition(p)->updateRarities(skillRarity_);
  }

  isPreprocessedSkills_=true;
}

//------------------------------------------------------------------------------
// Create the vector of sorted nurses
// The nurses are ordered according to their position and the nurses that have
// the same position are shuffled
//------------------------------------------------------------------------------

void Solver::sortShuffleTheNurses() {

  vector<Position*> positionsSorted;
  vector<vector<LiveNurse*>> nursePerPosition;

  // first, sort the position in the order in which the nurses should be treated
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    positionsSorted.push_back(pScenario_->pPosition(p));
  }
  std::sort(positionsSorted.begin(),positionsSorted.end(),comparePositions);

  // then organize the nurses depending on their position
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    vector<LiveNurse*> emptyVector;
    nursePerPosition.push_back(emptyVector);
  }
  for (int n=0; n < pScenario_->nbNurses_; n++) {
    LiveNurse* pNurse = theLiveNurses_[n];
    nursePerPosition[pNurse->pPosition()->id()].push_back(pNurse);
  }

  // shuffle the nurses that have the same position
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    std::random_shuffle(nursePerPosition[p].begin(),nursePerPosition[p].end());
  }

  // fill the sorted vector of live nurses
  theNursesSorted_.clear();
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    int id = positionsSorted[p]->id();
    for (int n=0; n < nursePerPosition[id].size(); n++) {
      theNursesSorted_.push_back(nursePerPosition[id][n]);
    }
  }

  // todo: think about sorting the nurses according to their contract, we might
  // want to use the full-time nurses first
}

//------------------------------------------------------------------------------
// Initialize the greedy by preprocessing all the input attributes and sorting
// the shifts, skills, nurses
//------------------------------------------------------------------------------

void Solver::preprocessData() {

  // Preprocess the attributes of the greedy solver
  // the result of the preprocessing will be very useful to sort the attributes
  // before greedily covering the demand
  //
  if (!pDemand_->isPreprocessed_) pDemand_->preprocessDemand();
  if (!isPreprocessedSkills_) this->preprocessTheSkills();
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  // sort the skills
  SkillSorter compareSkills(skillRarity_);
  std::stable_sort(skillsSorted_.begin(),skillsSorted_.end(),compareSkills);

  // sort the shifts (except the shift 0 which must always be rest)
  ShiftSorter compareShifts(pScenario_->nbForbiddenSuccessors_);
  std::stable_sort(shiftsSorted_.begin(),shiftsSorted_.end(),compareShifts);

  // sort the nurses
  sortShuffleTheNurses();
}


//-----------------------------------------------------------------------------
// Compute the weights o the violation of the min/max number of working days
// For now, the update depends only on the initial states and on the contract
// of the nurses, on the number of days on the demand, on the number of weeks
// already treated and on the number of weeks left
// The required data on the nurses is mostly computed in preprocessTheNurses
//-----------------------------------------------------------------------------

void Solver::setBoundsAndWeights(WeightStrategy strategy){
	switch (strategy){

	case MAX :
	case MEAN:
		computeWeightsTotalShiftsForPrimalDual(strategy);
		break;

	case RANDOMMEANMAX:
		if(rdm_()%2 == 0)
			computeWeightsTotalShiftsForPrimalDual(MEAN);
		else
			computeWeightsTotalShiftsForPrimalDual(MAX);
		break;


	case BOUNDRATIO:
		computeBoundsAccordingToDemandSize();
		break;

	case NO_STRAT:
		computeWeightsTotalShiftsForStochastic();
		break;

	default:
		Tools::throwError("Weight/bound strategy not defined");
		break;
	}
}


void Solver::computeWeightsTotalShiftsForStochastic() {

  // clear the vectors that are about to be filled
  minTotalShifts_.clear();
  weightTotalShiftsMin_.clear();
  maxTotalShifts_.clear();
  weightTotalShiftsMax_.clear();
  maxTotalWeekends_.clear();
  weightTotalWeekendsMax_.clear();

  minTotalShiftsAvg_.clear();
  maxTotalShiftsAvg_.clear();
  weightTotalShiftsAvg_.clear();
  maxTotalWeekendsAvg_.clear();
  weightTotalWeekendsAvg_.clear();

	// The nurses must be preprocessed to retrieve the information relative to the
	// past activity of the nurses and to their capacity to work more in the future
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  // The important value to infer the importance of respecting the strict constraints
  // on the total number of working days/week-ends is the remaining number of days/week-ends
  // after the demand currently treated
  int remainingDays = 7*pScenario_->nbWeeks()-7*(pScenario_->thisWeek())-pDemand_->nbDays_;
  double factorRemainingDays = (double) remainingDays/(double)(7*pScenario_->nbWeeks());


  //number of weekends before the end of the planning
  int remainingWeekends = pScenario_->nbWeeks() - pScenario_->thisWeek();
  //portion of these remaining weekends in this demand
  double factorRemainingWeekends =  (pDemand_->nbDays_/7)* 1.0 / remainingWeekends;
  //	double factorRemainingWeekends = (double)remainingWeekends/(double)pScenario_->nbWeeks();

  // Compute the non-penalized intervals and the associated penalties
  for (int n = 0; n < pScenario_->nbNurses(); n++) {
	  LiveNurse* pNurse =  theLiveNurses_[n];

	  // compute the interval that must be respected to have a chance of not paying
	  // penalties in the future
	  minTotalShifts_.push_back(pNurse->minWorkDaysNoPenaltyTotalDays_);
	  weightTotalShiftsMin_.push_back(WEIGHT_TOTAL_SHIFTS);
	  maxTotalShifts_.push_back(pNurse->maxWorkDaysNoPenaltyTotalDays_);
	  weightTotalShiftsMax_.push_back(WEIGHT_TOTAL_SHIFTS);

	  //penalize each worked weekend with primal-dual costs
	  maxTotalWeekends_.push_back(0);
	  weightTotalWeekendsMax_.push_back(WEIGHT_TOTAL_WEEKENDS * pNurse->pStateIni_->totalWeekendsWorked_ *
			  1.0 / pNurse->maxTotalWeekends());


	  if(weightTotalWeekendsMax_[n]>WEIGHT_TOTAL_WEEKENDS) {
		  weightTotalWeekendsMax_[n] = WEIGHT_TOTAL_WEEKENDS;
		  maxTotalWeekendsAvg_.push_back(pNurse->maxTotalWeekends());
		  weightTotalWeekendsAvg_.push_back(WEIGHT_TOTAL_WEEKENDS);
	  }
	  else{
		  //compute the proportion of weekends that can be worked in this demand without exceeding the max in the future
		  //round with a certain probability to the floor or the ceil
		  int remainingWeekendsToWork = pNurse->maxTotalWeekends() - pNurse->pStateIni_->totalWeekendsWorked_;
		  double numberOfAuthorizedWeekend = remainingWeekendsToWork * factorRemainingWeekends;
		  maxTotalWeekendsAvg_.push_back( Tools::roundWithProbability(numberOfAuthorizedWeekend) );
     // maxTotalWeekendsAvg_.push_back(numberOfAuthorizedWeekend);
		  weightTotalWeekendsAvg_.push_back(WEIGHT_TOTAL_WEEKENDS);
	  }


	  /* Essai Antoine */
	  // first compute the values relative to the average number of working days
	  // the interval is larger for the first weeks and the associated penalty is smaller
	  // minTotalShiftsAvg_.push_back((1.0-0.25*factorRemainingDays)*pNurse->minAvgWorkDaysNoPenaltyTotalDays_);
	  // maxTotalShiftsAvg_.push_back((1.0+0.25*factorRemainingDays)*pNurse->maxAvgWorkDaysNoPenaltyTotalDays_);
	  // weightTotalShiftsAvg_.push_back((1.0-factorRemainingDays)*(double)WEIGHT_TOTAL_SHIFTS);
    double factorMarginOnAvg = 0.25*factorRemainingDays*7.0/(double)pDemand_->nbDays_;
    factorMarginOnAvg = (pDemand_->nbDays_%14==0) ? 0.0:factorMarginOnAvg;
    minTotalShiftsAvg_.push_back((1.0-factorMarginOnAvg)*pNurse->minAvgWorkDaysNoPenaltyTotalDays_);
    //double maxAvg = (double)pDemand_->nbDays_/7.0*pNurse->maxTotalShifts()/(double)pScenario_->nbWeeks();
		maxTotalShiftsAvg_.push_back((1.0+factorMarginOnAvg)*pNurse->maxAvgWorkDaysNoPenaltyTotalDays_);
    weightTotalShiftsAvg_.push_back((double)WEIGHT_TOTAL_SHIFTS);



	  // Affichage des valeurs
	  // std::cout << std::endl;
	  // std::cout << "###############################################" << std::endl;
	  // std::cout << "# Nurse " << pNurse->name_ << std::endl;
	  // std::cout << "#   | Contract: " << pNurse->pContract_->name_ << std::endl;
	  // std::cout << "#   |           " << (*(pNurse->pContract_)) << std::endl;
	  // std::cout << "#   | History : " << pNurse->pStateIni_->totalDaysWorked_ << " days and " << pNurse->pStateIni_->totalWeekendsWorked_ << " weekends" << std::endl;
	  // std::cout << "# " << std::endl;
	  // std::cout << "# Costs / Bounds:" << std::endl;
	  // std::cout << "#   | Min total shifts: " << pNurse->minWorkDaysNoPenaltyTotalDays_ << std::endl;
	  // std::cout << "#   | Max total shifts: " << pNurse->maxWorkDaysNoPenaltyTotalDays_ << std::endl;
	  // std::cout << "#   | Weight total we : " << weightTotalWeekendsMax_.back() << std::endl;
	  // std::cout << "#   | Min tot shifts av " << minTotalShiftsAvg_.back() << std::endl;
	  // std::cout << "#   | Max tot shifts av " << maxTotalShiftsAvg_.back() << std::endl;
	  // std::cout << "#   | Weight tot sh avg " << weightTotalShiftsAvg_.back() << std::endl;
	  // std::cout << "###############################################" << std::endl;
	  // std::cout << std::endl;




	  /* Essai Sam (couts d'oportunite) */
	  //		maxTotalWeekends_[n] = 0;
	  //		double costOfWeekendForNurse;
	  //
	  //		int nWEAlreadyDoneByNurse = pNurse->pStateIni_->totalWeekendsWorked_;
	  //		int nWERemainingIncludingNow = pScenario_->nbWeeks_ - pScenario_->thisWeek();
	  //		int nWEMaxForNurse = pNurse->pContract_->maxTotalWeekends_;
	  //
	  //		if(nWEAlreadyDoneByNurse >= nWEMaxForNurse)
	  //			costOfWeekendForNurse = nWERemainingIncludingNow * WEIGHT_TOTAL_WEEKENDS;
	  //		else if(nWERemainingIncludingNow <= nWEMaxForNurse)
	  //			costOfWeekendForNurse = std::max( 0, nWEAlreadyDoneByNurse + nWERemainingIncludingNow - nWEMaxForNurse );
	  //		else
	  //			costOfWeekendForNurse = nWEAlreadyDoneByNurse;
	  //
	  //		weightTotalWeekendsMax_[n] = costOfWeekendForNurse ;

	}

//  getchar();
}

void Solver::computeWeightsTotalShiftsForPrimalDual(WeightStrategy strategy){
	if(strategy == NO_STRAT)
		Tools::throwError("Weight strategy not defined.");

   // clear the vectors that are about to be filled
   minTotalShiftsAvg_.clear();
   maxTotalShiftsAvg_.clear();
   weightTotalShiftsAvg_.clear();
   maxTotalWeekendsAvg_.clear();
   weightTotalWeekendsAvg_.clear();

   minTotalShiftsContractAvg_.clear();
   maxTotalShiftsContractAvg_.clear();
   weightTotalShiftsContractAvg_.clear();
   maxTotalWeekendsContractAvg_.clear();
   weightTotalWeekendsContractAvg_.clear();

   vector<double> maxPrimalDualCostForContractDays;
   vector<double> maxPrimalDualCostForContractWE;
   vector<double> meanPrimalDualCostForContractDays;
   vector<double> meanPrimalDualCostForContractWE;
   vector<int> nbNursesPerContract;

      // The nurses must be preprocessed to retrieve the information relative to the
      // past activity of the nurses and to their capacity to work more in the future
      if (!isPreprocessedNurses_) this->preprocessTheNurses();

      //initialize the non-penalized intervals and the associated penalties for each contract
      for(int p=0; p<pScenario_->nbContracts_; ++p){
         minTotalShiftsContractAvg_.push_back(0);
         maxTotalShiftsContractAvg_.push_back(0);
         weightTotalShiftsContractAvg_.push_back(WEIGHT_TOTAL_SHIFTS);
         maxTotalWeekendsContractAvg_.push_back(0);
         weightTotalWeekendsContractAvg_.push_back(WEIGHT_TOTAL_WEEKENDS);
         maxPrimalDualCostForContractDays.push_back(0);
         maxPrimalDualCostForContractWE.push_back(0);
         meanPrimalDualCostForContractDays.push_back(0);
         meanPrimalDualCostForContractWE.push_back(0);
         nbNursesPerContract.push_back(0);
      }

      // Compute the non-penalized intervals and the associated penalties for each contract
      for (int n = 0; n < pScenario_->nbNurses(); n++) {
         LiveNurse* pNurse =  theLiveNurses_[n];
         double ratio = WEIGHT_TOTAL_SHIFTS * pNurse->pStateIni_->totalDaysWorked_;

         // Deactivate minimum constraints
         minTotalShifts_[n] = 0;
         weightTotalShiftsMin_[n] = 0; //WEIGHT_TOTAL_SHIFTS - ratio * 1.0 / pNurse->minTotalShifts();

         // Activate maximum constraints (work) with primal-dual cost
         // Remark: these costs are always nonnegative here because, for the competition, nurses are always in shortage
         //         If this does not apply -> find another formulation
         //
         maxTotalShifts_[n] = 0;
         weightTotalShiftsMax_[n] = ratio * 1.0 / pNurse->maxTotalShifts(); 								// Primal-dual cost of max working days
//         weightTotalShiftsMax_[n] += - WEIGHT_TOTAL_SHIFTS + ratio * 1.0 / pNurse->minTotalShifts(); 		// Primal-dual cost of min working days
         if(weightTotalShiftsMax_[n]>WEIGHT_TOTAL_SHIFTS) weightTotalShiftsMax_[n] = WEIGHT_TOTAL_SHIFTS;	// Must not be higher that WEIGHT
         else if(weightTotalShiftsMax_[n]<0) weightTotalShiftsMax_[n] = 0;									// Must not be negative

         // Activate maximum constraints (WE) with primal-dual cost
         //
         maxTotalWeekends_[n] = 0;
         weightTotalWeekendsMax_[n] = WEIGHT_TOTAL_WEEKENDS * pNurse->pStateIni_->totalWeekendsWorked_ *
            1.0 / pNurse->maxTotalWeekends();
         if(weightTotalWeekendsMax_[n]>WEIGHT_TOTAL_WEEKENDS) weightTotalWeekendsMax_[n] = WEIGHT_TOTAL_WEEKENDS;



         // std::cout << weightTotalShiftsMin_[n] << " " << weightTotalShiftsMax_[n] << " " << weightTotalWeekendsMax_[n] << std::endl;

         // Update average days/weekends in contract
         //
         int p = pNurse->pContract_->id_;
         minTotalShiftsContractAvg_[p] += std::max(0, pNurse->minTotalShifts() - pNurse->pStateIni_->totalDaysWorked_);
         maxTotalShiftsContractAvg_[p] += std::max(0, pNurse->maxTotalShifts() - pNurse->pStateIni_->totalDaysWorked_);
         maxTotalWeekendsContractAvg_[p] += std::max(0, pNurse->maxTotalWeekends() - pNurse->pStateIni_->totalWeekendsWorked_);

         // Check for maximum primal-dual costs
         //
         maxPrimalDualCostForContractDays[p] = std::max( maxPrimalDualCostForContractDays[p] , weightTotalShiftsMax_[n]);
         maxPrimalDualCostForContractWE[p] = std::max( maxPrimalDualCostForContractWE[p] , weightTotalWeekendsMax_[n]);
         meanPrimalDualCostForContractDays[p] += weightTotalShiftsMax_[n];
         meanPrimalDualCostForContractWE[p] += weightTotalWeekendsMax_[n];
         nbNursesPerContract[p] ++;

      }

      //round the min/max values of the interval associated to the contract
      int remainingWeeks = pScenario_->nbWeeks() - pScenario_->thisWeek();
      int remainingDays = 7* remainingWeeks;
      int nDaysInHorizon = pDemand_->nbDays_;
      int nWeekendsInHorizon = (pDemand_->nbDays_+1) / 7;
      double ratioDays = nDaysInHorizon * (1.0/remainingDays);
      double ratioWeekends = nWeekendsInHorizon * (1.0/remainingWeeks);

      for(int p=0; p<pScenario_->nbContracts_; ++p){
         minTotalShiftsContractAvg_[p] = Tools::roundWithProbability( minTotalShiftsContractAvg_[p] * ratioDays);
         maxTotalShiftsContractAvg_[p] = Tools::roundWithProbability( maxTotalShiftsContractAvg_[p] * ratioDays);
         maxTotalWeekendsContractAvg_[p] = Tools::roundWithProbability( maxTotalWeekendsContractAvg_[p] * ratioWeekends);

         meanPrimalDualCostForContractDays[p] /= ((double)nbNursesPerContract[p]);
         meanPrimalDualCostForContractWE[p] /= ((double)nbNursesPerContract[p]);

         //to penalize some contracts vs the others
         const string contractName = pScenario_->intToContract_[p];
         Contract* pContract = (pScenario_->contracts_.at( contractName ));
         int lengthStintMin = pContract->minConsDaysWork_ + pContract->maxConsDaysOff_;
         //compute the number of days worked, if the nurse works the minimum days without penalties
         int minWorkDaysNoPenalty = pContract->minConsDaysWork_ * (7*pScenario_->nbWeeks_) / lengthStintMin;
         //compute the ratio between the min work days without penalties and the maximum number of shifts
         double ratioMinMax = minWorkDaysNoPenalty * 1.0 / pContract->maxTotalShifts_;
         double weightContract = 1 + (ratioMinMax)/10.0;
         weightTotalShiftsContractAvg_[p] *= weightContract;

         if(strategy == MAX){
        	 weightTotalShiftsContractAvg_[p] -= maxPrimalDualCostForContractDays[p];
        	 weightTotalWeekendsContractAvg_[p] -= maxPrimalDualCostForContractWE[p];
         }
         else if(strategy == MEAN){
        	 weightTotalShiftsContractAvg_[p] -= meanPrimalDualCostForContractDays[p];
        	 weightTotalWeekendsContractAvg_[p] -= meanPrimalDualCostForContractWE[p];
         }

         if(false){
        	 std::cout << "# " << std::endl;
        	 std::cout << "##################################################" << std::endl;
        	 const string contractName = pScenario_->intToContract_[p];
        	 std::cout << "# " << (*(pScenario_->contracts_.at( contractName ))) << std::endl;
           //std::cout << "# " << minWorkDaysNoPenalty << " " << ratioMinMax << ": " << weightContract << std::endl;
        	 std::cout << "# min/max : " << std::endl;
        	 std::cout << "#    | min total shifts: " << minTotalShiftsContractAvg_[p] << std::endl;
        	 std::cout << "#    | max total shifts: " << maxTotalShiftsContractAvg_[p] << std::endl;
        	 std::cout << "#    | max total we    : " << maxTotalWeekendsContractAvg_[p] << std::endl;
        	 std::cout << "# costs   : " << std::endl;
        	 std::cout << "#    | cost shift: " << weightTotalShiftsContractAvg_[p] << std::endl;
        	 std::cout << "#    | cost we   : " << weightTotalWeekendsContractAvg_[p] << std::endl;
        	 std::cout << "##################################################" << std::endl;
        	 std::cout << "# " << std::endl;
         }
      }
}


// Compute min/max bounds as the ratio : number of days in demand / total number of remaining days
//
void Solver::computeBoundsAccordingToDemandSize(){

	   // clear the vectors that are about to be filled
	   minTotalShiftsAvg_.clear();
	   maxTotalShiftsAvg_.clear();
	   weightTotalShiftsAvg_.clear();
	   maxTotalWeekendsAvg_.clear();
	   weightTotalWeekendsAvg_.clear();

	   minTotalShiftsContractAvg_.clear();
	   maxTotalShiftsContractAvg_.clear();
	   weightTotalShiftsContractAvg_.clear();
	   maxTotalWeekendsContractAvg_.clear();
	   weightTotalWeekendsContractAvg_.clear();


	   // The nurses must be preprocessed to retrieve the information relative to the
	   // past activity of the nurses and to their capacity to work more in the future
	   if (!isPreprocessedNurses_) this->preprocessTheNurses();

	   double ratio = pDemand_->nbDays_ / (7.0 * (pScenario_->nbWeeks_ - pScenario_->thisWeek()));
//	   std::cout << "# Ratio: " << ratio << std::endl;

	   // initialize the minimum and maximum number of total working days
	   for (int i = 0; i < pScenario_->nbNurses(); i++) {
		   //default min and max

		   //	     minTotalShifts_[i] = Tools::roundWithProbability((theLiveNurses_[i]->minTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_) * ratio);
		   //	     maxTotalShifts_[i] = Tools::roundWithProbability((theLiveNurses_[i]->maxTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_) * ratio);
		   //	     maxTotalWeekends_[i] = Tools::roundWithProbability((theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) * ratio);

		   minTotalShifts_[i] = std::floor((theLiveNurses_[i]->minTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_) * ratio);
		   maxTotalShifts_[i] = std::ceil((theLiveNurses_[i]->maxTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_) * ratio);
		   maxTotalWeekends_[i] = std::ceil((theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) * ratio);

//		   std::cout << "# " << theLiveNurses_[i]->pStateIni_->totalDaysWorked_ << " " << theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_ << "\t\t";
//		   std::cout << minTotalShifts_[i] << " " << maxTotalShifts_[i] << " " << maxTotalWeekends_[i] << std::endl;


	   }
}



//------------------------------------------------------------------------
// Compare two nurses based on their position
// the function is used to sort the nurses in ascending rank of their
// position
// if their positions have the same rank, then the smaller nurse is found
// by a lexicographic comparison of the rarity of the skills of the nurses
//------------------------------------------------------------------------

bool compareNurses(LiveNurse* n1, LiveNurse* n2) {
  return comparePositions(n1->pPosition(), n2->pPosition());
}

//------------------------------------------------------------------------
// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//------------------------------------------------------------------------

bool comparePositions(Position* p1, Position* p2) {
  if (p1->id() == p2->id()) {
    return false;
  }
  else if (p1->rank() == p2->rank()) {
    // the skillRarity vector is ALWAYS sorted in descending order, because the
    // updateRarities is the only setter for skillRarity and it sorts the vector
    for (int sk=0; sk<std::min(p1->nbSkills(),p2->nbSkills()); sk++) {
      if (p1->skillRarity(sk) != p2->skillRarity(sk)) {
        return p1->skillRarity(sk) > p2->skillRarity(sk);
      }
      return p1->nbSkills() > p2->nbSkills();
    }
  }
  else {
    return p1->rank() < p2->rank();
  }
  return true;
}


//------------------------------------------------------------------------
// Check the feasibility of the demand with these nurses
//------------------------------------------------------------------------

bool checkFeasibility() {
  return true;
}

// get the total cost of the current solution
// the solution is simply given by the roster of each nurse
double Solver::solutionCost(int nbDays) {
  double totalCost = 0.0;
  int nbNurses = pScenario_->nbNurses_;
  int nbShifts = pScenario_->nbShifts_, nbSkills = pScenario_->nbSkills_;

  // reset the satisfied demand to compute it from scratch
  Tools::initVector3D(&satisfiedDemand_,nbDays, nbShifts, nbSkills);

  // first add the individual cost of each nurse
  for (int n = 0; n < nbNurses; n++) {
    LiveNurse *pNurse = theLiveNurses_[n];
    pNurse->checkConstraints(pNurse->roster_, pNurse->states_, pNurse->statCt_);
    StatCtNurse stat = pNurse->statCt_;

    for (int day = 0; day < nbDays ; day++) {
      totalCost += stat.costConsDays_[day]+stat.costConsDaysOff_[day]+
        stat.costConsShifts_[day]+stat.costPref_[day]+stat.costWeekEnd_[day];

      if (pNurse->roster_.shift(day) > 0) {
        satisfiedDemand_[day][pNurse->roster_.shift(day)][pNurse->roster_.skill(day)]++;
      }
    }

    if ( pScenario_->thisWeek()+pScenario_->nbWeeksLoaded() == pScenario_->nbWeeks() ) {
      totalCost += stat.costTotalDays_+stat.costTotalWeekEnds_;
    }
  }
  std::cout << "Total cost due to individual soft constraints = " << totalCost << std::endl;

  // add the cost of non-optimal demand
  for(int day = 0; day < nbDays; day++) {
    for (int sh = 1; sh < nbShifts ; sh++) {
      for (int sk = 0; sk < nbSkills; sk++) {
        int missingStaff;
        missingStaff = std::max(0, pDemand_->optDemand_[day][sh][sk] - satisfiedDemand_[day][sh][sk]);
        totalCost += WEIGHT_OPTIMAL_DEMAND*missingStaff;
      }
    }
  }

  std::cout << "Total cost = " << totalCost << std::endl;


  return totalCost;
}

//------------------------------------------------
// Display functions
//------------------------------------------------

// Return the solution at day k
//
vector<Roster> Solver::getSolutionAtDay(int k){
	vector<Roster> ans;
	for(int i=0; i<solution_.size(); i++){
		Roster r = solution_[i];
		int nbDays = k+1;
		int firstDay = r.firstDay();
		vector<int> shifts, skills;
		for(int j=0; j<=k; j++){
			shifts.push_back(r.shift(j));
			skills.push_back(r.skill(j));
		}
		Roster s (nbDays, firstDay, shifts, skills);
		ans.push_back(s);
	}
	return ans;
}

// return the final states of the nurses
//
vector<State> Solver::getFinalStates() {
  vector<State> pFinalStates;
  int nbDays = pDemand_->nbDays_;
  for (LiveNurse* pNurse: theLiveNurses_) {
    pFinalStates.push_back(pNurse->state(nbDays));
  }

  return pFinalStates;
}

// return the states of the nurses at day k
//
vector<State> Solver::getStatesOfDay(int k) {
  vector<State> pStatesOfDayK;
  int nbDays = pDemand_->nbDays_;
  for (LiveNurse* pNurse: theLiveNurses_) {
	  pStatesOfDayK.push_back(pNurse->state(k+1));
  }
  return pStatesOfDayK;
}

// display the whole solution
//
string Solver::solutionToString() {
   return solutionToString(pDemand_->firstDay_, pDemand_->nbDays_, pScenario_->thisWeek());
}

// display the whole solution week by week for nbWeeks weeks.
//
vector<string> Solver::solutionToString(int nbWeeks) {
   vector<string> solutions;

   //build the solution for each week
   int firstDay = pDemand_->firstDay_;
   for(int w=0; w<nbWeeks; ++w){
      solutions.push_back(solutionToString(firstDay, 7, pScenario_->thisWeek()+w));
      firstDay += 7;
   }

   return solutions;
}

// display the solution between firstDay and firstDay+nbDays in the required format
//
string Solver::solutionToString(int firstDay, int nbDays, int firstWeek){
  std::stringstream rep;
  int nbNurses = pScenario_->nbNurses_;

  // write to stringstream that can then be printed in any output file
  // follow the template described by the competition
  rep << "SOLUTION" << std::endl;
  rep << firstWeek << " " << pScenario_->name_ << std::endl;
  rep << std::endl;

  // compute the total number of assignments that are not rests
  // if no shift is assigned to a nurse on given, it still counts
  int nbAssignments = 0;
  for (int n = 0; n < nbNurses; n ++) {
    for (int day = firstDay; day < firstDay+nbDays; day++){
      if (theLiveNurses_[n]->roster_.shift(day) > 0) nbAssignments++;
    }
  }

  rep << "ASSIGNMENTS = " << nbAssignments << std::endl;
  for (int n = 0; n < nbNurses; n ++) {
    for (int day = firstDay; day < firstDay+nbDays; day++){
      int shift = theLiveNurses_[n]->roster_.shift(day);
      if (shift != 0) {
        rep << theLiveNurses_[n]->name_ << " "<< Tools::intToDay(day) << " ";
        if (shift > 0) {
          int skill = theLiveNurses_[n]->roster_.skill(day);
          rep << pScenario_->intToShift_[shift] << " ";
          rep << pScenario_->intToSkill_[skill] << std::endl;
        }
        else { // shift < 0 so no shif is assigned
          rep << "Unassigned TBD"<< std::endl;
          LiveNurse* pNurse = theLiveNurses_[n];
          std::cout  << "This shift " << pNurse->states_[day+1].shift_ << std::endl;
        }
      }
    }
  }

  return rep.str();
}

// display the whole solution in a more readable format and append advanced
// information on the solution quality
//
string Solver::solutionToLogString() {
  std::stringstream rep;
  int nbNurses = pScenario_->nbNurses_, nbShifts = pScenario_->nbShifts_;
  int nbSkills = pScenario_->nbSkills_;
  int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nbDays_;

  rep << "Complete shift schedule" << std::endl << std::endl;
  rep << "\t\t\t";
  for (int day = firstDay; day < firstDay+nbDays; day++) {
    rep << "| " << Tools::intToDay(day).at(0) << " ";
  }
  rep << "|" << std::endl;
  rep << "-------------------------------------"<< std::endl;

  for (int n = 0; n < nbNurses; n ++) {
    LiveNurse* pNurse = theLiveNurses_[n];
    rep << pNurse->name_ << "\t";
    for (int day = firstDay; day < firstDay+nbDays; day++){
      int shift = pNurse->roster_.shift(day);
      if (shift != 0) {
        rep << "| " <<  pScenario_->intToShift_[shift].at(0) << " ";
      }
      else {
        rep << "| - ";
      }

      if(Tools::isSunday(day)) rep << "| ";
    }
    rep << "|" << std::endl;
  }
  rep << std::endl;

  // compute the total cost and in the mean time update the structures of each
  // live nurse that contains all the required information on soft and hard
  // constraints satisfaction
  //
  double totalCost = solutionCost();

  // store temporarily the data that is about to be written
  //
  int violMinCover = 0, violReqSkill = 0, violForbiddenSucc = 0;
  double costOptCover = 0, costTotalDays = 0, costTotalWeekEnds = 0;
  double costConsDays = 0, costConsDaysOff=0, costConsShifts = 0, costPref = 0, costWeekEnds = 0;

  // constraints related to the demand
  for (int day = firstDay; day < firstDay+nbDays; day++) {
    for (int sh = 1; sh < nbShifts; sh++) {
      for (int sk = 0; sk < nbSkills; sk++) {
        violMinCover += std::max(0,pDemand_->minDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]);
        costOptCover += WEIGHT_OPTIMAL_DEMAND
          * std::max(0,pDemand_->optDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]);
//        if(pDemand_->minDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]>0)
//           std::cout << day << " " << sh  << " " << sk << " " << pDemand_->minDemand_[day][sh][sk] << " " << satisfiedDemand_[day][sh][sk] << std::endl;
//        else if(pDemand_->optDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]>0)
//                   std::cout << day << " " << sh  << " " << sk << " " << satisfiedDemand_[day][sh][sk] << " " << pDemand_->optDemand_[day][sh][sk] << std::endl;
//        else if(pDemand_->optDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]<0)
//                           std::cout << "*" << day << " " << sh  << " " << sk << " " << satisfiedDemand_[day][sh][sk] << " " << pDemand_->optDemand_[day][sh][sk] << std::endl;
      }
    }
  }


  for (int n = 0; n < nbNurses; n ++) {
    LiveNurse* pNurse = theLiveNurses_[n];

    if ( pScenario_->thisWeek()+pScenario_->nbWeeksLoaded() == pScenario_->nbWeeks() ) {
      costTotalDays += pNurse->statCt_.costTotalDays_;
      costTotalWeekEnds += pNurse->statCt_.costTotalWeekEnds_;
    }

    for (int day = firstDay; day < firstDay+nbDays; day++){
      // record the violations
      int skill = pNurse->roster_.skill(day);
      int shift = pNurse->roster_.shift(day);
      int prevShift = pNurse->states_[day].shift_;
      violReqSkill += shift == 0 ? 0 : (pNurse->hasSkill(skill)? 0:1);
      violForbiddenSucc += pScenario_->isForbiddenSuccessor(shift, prevShift)? 1: 0;

      // the other costs per soft constraint can be read from the stat structure
      costConsDays += pNurse->statCt_.costConsDays_[day];
      costConsDaysOff += pNurse->statCt_.costConsDaysOff_[day];
      costConsShifts += pNurse->statCt_.costConsShifts_[day];
      costPref += pNurse->statCt_.costPref_[day];
      costWeekEnds += pNurse->statCt_.costWeekEnd_[day];

    }
  }

  // write the status of hard and soft constraints
  //
  rep << "Hard constraints violations\n";
  rep << "---------------------------\n";
  rep << "Minimal coverage constraints: " << violMinCover << std::endl;
  rep << "Required skill constraints: " << violReqSkill << std::endl;
  rep << "Illegal shift type succession constraints: " << violForbiddenSucc << std::endl;
  rep << "Single assignment per day: 0" << std::endl;

  rep << "\nCost per constraint type\n";
  rep << "------------------------\n";
  rep << "Total assignment constraints: " << costTotalDays << std::endl;
  //rep << "Consecutive constraints: " << costConsDays+costConsShifts << std::endl;
  rep << "Consecutive working days constraints: " << costConsDays << std::endl;
  rep << "Consecutive days off constraints: " << costConsDaysOff << std::endl;
  rep << "Consecutive shifts constraints: " << costConsShifts << std::endl;
  rep << "Preferences: " << costPref << std::endl;
  rep << "Max working weekend: " << costTotalWeekEnds << std::endl;
  rep << "Complete weekends: " << costWeekEnds << std::endl;
  rep << "Optimal coverage constraints: " << costOptCover << std::endl;

  rep << "\n---------------------------\n";
  rep << "\nTotal cost: " << totalCost << std::endl;

  return rep.str();
}
