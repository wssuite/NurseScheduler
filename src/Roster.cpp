#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"
#include "Roster.h"


// Constructor: initialize planning from nothing
//
Roster::Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map<int,std::set<int>>* pWishesOff, const State& initialState):
    nbDays_(nbDays), firstDay_(firstDay),
    pScenario_(pScenario), pNurse_(pNurse), pWishesOff_(pWishesOff) {

  // initialize the roster with a rest week
  for (int i = 0; i < nbDays; i++) {
    tasks_.push_back({0,0});
  }

  // initialize the states at each day
  states_.push_back(initialState);
  for (int day = 0; day < nbDays_; day++) {
    State nextState;
    nextState.addDayToState(states_[day], 0);
  }

  // initialize all the cost vectors
  Tools::initVector(&costConsShifts_, nbDays_);
  Tools::initVector(&costConsDays_, nbDays_);
  Tools::initVector(&costPreferences_, nbDays_);
  Tools::initVector(&costCompleteWeekEnd_, nbDays_);

  // initialize the violation vector
  for (int day = 0; day < nbDays_; day++) violationSuccShifts_.push_back(false);
}

// Constructor: initialize planning from an input set of tasks for the nurse
//
Roster::Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
std::map<int,std::set<int>>* pWishesOff, const State& initialState,
vector<task> inputTasks):
nbDays_(nbDays), firstDay_(firstDay), pScenario_(pScenario),
pNurse_(pNurse), pWishesOff_(pWishesOff), tasks_(inputTasks) {

  // initialize the states at each day
  states_.push_back(initialState);
  for (int day = 0; day < nbDays_; day++) {
    State nextState;
    nextState.addDayToState(states_[day], tasks_[day].shift);
    states_.push_back(nextState);
  }

  // initialize all the cost vectors
  Tools::initVector(&costConsShifts_, nbDays_);
  Tools::initVector(&costConsDays_, nbDays_);
  Tools::initVector(&costPreferences_, nbDays_);
  Tools::initVector(&costCompleteWeekEnd_, nbDays_);

  // initialize the violation vector
  for (int day = 0; day < nbDays_; day++) violationSuccShifts_.push_back(false);

  // check the soft and hard constraints
  checkSoftConstraints();
  checkHardConstraints();

}

// Destructor
Roster::~Roster(){
}

// check the satisfaction of the hard constraints and record the violations
//
void Roster::checkHardConstraints() {

  for (int day = 0; day < nbDays_; day++) {

    // last shift assigned to the nurse
    int lastShift = states_[day].shift_;
    // shift assigned on this day
    int thisShift = tasks_[day].shift;

    for (int nbShift = 0; nbShift < pScenario_->nbForbiddenSuccessors_[lastShift]; nbShift++) {
       if (thisShift == pScenario_->forbiddenSuccessors_[lastShift][nbShift])  {
         violationSuccShifts_[day] = true;
         break;
       }
    }
  }

}

// check the soft constraints and record the costs of the violations and the
// remaining margin for the satisfied ones.
//
void Roster::checkSoftConstraints() {
  for (int day = 1; day <= nbDays_; day++) {

    // shift assigned on this day
    int shift = states_[day].shift_;

    // first look at consecutive working days or days off
    //
    if (switchOff_[day])  {
      int missingDays=0, extraDays=0;

      // compute the violations of consecutive working days
      if (shift) {
        missingDays = pNurse_->minConsDaysWork_-states_[day].consDaysWorked_;
        extraDays = states_[day].consDaysWorked_-pNurse_->maxConsDaysWork_;

        costConsDays_[day] = (missingDays>0) ? WEIGHT_CONS_DAYS_WORK*missingDays:0;
        costConsDays_[day] += (extraDays>0) ? WEIGHT_CONS_DAYS_WORK*extraDays:0;
      }
      // compute the violations of consecutive days off
      else {
        missingDays = pNurse_->minConsDaysOff_-states_[day].consDaysOff_;
        extraDays = states_[day].consDaysOff_-pNurse_->maxConsDaysOff_;

        costConsDays_[day] = (missingDays>0) ? WEIGHT_CONS_DAYS_OFF*missingDays:0;
        costConsDays_[day] += (extraDays>0) ? WEIGHT_CONS_DAYS_OFF*extraDays:0;
      }
    }
    else  {
      costConsDays_[day] = 0;
    }

    // check the consecutive same shifts
    //
    if (switchShift_[day])  {
      int missingShifts = 0, extraShifts = 0;

      // it only makes sense if the nurse is working that day
      if (shift) {
        missingShifts = pScenario_->minConsShifts_[shift]-states_[day].consShifts_;
        extraShifts =  states_[day].consShifts_-pScenario_->maxConsShifts_[shift];

        costConsShifts_[day] = (missingShifts>0) ? WEIGHT_CONS_SHIFTS*missingShifts:0;
        costConsShifts_[day] += (extraShifts>0) ? WEIGHT_CONS_SHIFTS*extraShifts:0;
      }
    }
    else {
      costConsShifts_[day] = 0;
    }

    // check the preferences
    //
    map<int,set<int> >::iterator itM = pWishesOff_->find(day);
    // If the day is not in the wish-list, no possible violation
    if(itM == pWishesOff_->end())  {
      costPreferences_[day] = 0;
    }
    // no preference either in the wish-list for that day
    else if(itM->second.find(shift) == itM->second.end()) {
      costPreferences_[day] = 0;
    }
    else {
      costPreferences_[day] = WEIGHT_PREFERENCES;
    }

    // check the complete week-end (only if the nurse requires them)
    // this cost is only assigned to the sundays
    //
    if ( (day+this->firstDay_)%7 == 6 && pNurse_->isCompleteWeekEnds_) {
      if (states_[day].consDaysWorked_==1 || states_[day].consDaysOff_==1) {
        costCompleteWeekEnd_[day] = WEIGHT_COMPLETE_WEEKEND;
      }
    }

  } // end for day

}
