/*
* Solver.cpp
*
*  Created on: 22 déc. 2014
*      Author: jeremy
*/


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
pStateIni_(pStateIni), pWishesOff_(pWishesOff) {

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
    stat.violSkill_[day] = hasSkill(roster.skill(day));

    // Check the forbidden successor constraint
    //
    int lastShift = states[day].shift_;   // last shift assigned to the nurse
    int thisShift = roster.shift(day);    // shift assigned on this day

    stat.violSuccShifts_[day] = pScenario_->isForbiddenSuccessor(thisShift,lastShift);
  }

  // initialize the
  // check the soft constraints and record the costs of the violations and the
  // remaining margin for the satisfied ones.
  //
  for (int day = 1; day <= nbDays_; day++) {

    // shift assigned on this day
    int shift = states[day].shift_;

    // first look at consecutive working days or days off
    //
    if (roster.switchOff(day))  {
      int missingDays=0, extraDays=0;

      // compute the violations of consecutive working days
      if (shift) {
        missingDays = minConsDaysWork()-states[day].consDaysWorked_;
        extraDays = states[day].consDaysWorked_-maxConsDaysWork();

        stat.costConsDays_[day] = (missingDays>0) ? WEIGHT_CONS_DAYS_WORK*missingDays:0;
        stat.costConsDays_[day] += (extraDays>0) ? WEIGHT_CONS_DAYS_WORK*extraDays:0;
      }
      // compute the violations of consecutive days off
      else {
        missingDays =minConsDaysOff()-states[day].consDaysOff_;
        extraDays = states[day].consDaysOff_-maxConsDaysOff();

        stat.costConsDays_[day] = (missingDays>0) ? WEIGHT_CONS_DAYS_OFF*missingDays:0;
        stat.costConsDays_[day] += (extraDays>0) ? WEIGHT_CONS_DAYS_OFF*extraDays:0;
      }
    }
    else  {
      stat.costConsDays_[day] = 0;
    }

    // check the consecutive same shifts
    //
    if (roster.switchShift(day))  {
      int missingShifts = 0, extraShifts = 0;

      // it only makes sense if the nurse is working that day
      if (shift) {
        missingShifts = pScenario_->minConsShifts_[shift]-states[day].consShifts_;
        extraShifts =  states[day].consShifts_-pScenario_->maxConsShifts_[shift];

        stat.costConsShifts_[day] = (missingShifts>0) ? WEIGHT_CONS_SHIFTS*missingShifts:0;
        stat.costConsShifts_[day] += (extraShifts>0) ? WEIGHT_CONS_SHIFTS*extraShifts:0;
      }
    }
    else {
      stat.costConsShifts_[day] = 0;
    }

    // check the preferences
    //
    map<int,set<int> >::iterator itM = pWishesOff_->find(day);
    // If the day is not in the wish-list, no possible violation
    if(itM == pWishesOff_->end())  {
      stat.costPref_[day] = 0;
    }
    // no preference either in the wish-list for that day
    else if(itM->second.find(shift) == itM->second.end()) {
      stat.costPref_[day] = 0;
    }
    else {
      stat.costPref_[day] = WEIGHT_PREFERENCES;
    }

    // check the complete week-end (only if the nurse requires them)
    // this cost is only assigned to the sundays
    //
    if ( (day+this->firstDay_)%7 == 6 && needCompleteWeekends()) {
      if (states[day].consDaysWorked_==1 || states[day].consDaysOff_==1) {
        stat.costWeekEnd_[day] = WEIGHT_COMPLETE_WEEKEND;
      }
    }

  } // end for day

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
  pPreferences_(pPreferences), pInitState_(pInitState) {

    // initialize the preprocessed data of the skills
    for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
      maxStaffPerSkill_.push_back(0);
      maxStaffPerSkillNoPenalty_.push_back(0);
      skillRarity_.push_back(0);
    }

    // copy the nurses in the live nurses vector
    for (int i = 0; i < pScenario_->nbNurses_; i++) {
      theLiveNurses_.push_back(
        new LiveNurse( (pScenario_->theNurses_[i]), pScenario_, pDemand_->nbDays_,
        pDemand_->firstDay_, &(*pInitState_)[i], &(pPreferences_->wishesOff_[i])  ) );
    }


  }

// Destructor
Solver::~Solver(){}


//------------------------------------------------
// Preprocess the data
//------------------------------------------------

// go through the nurses to collect data regarding the potential shift and
// skill coverage of the nurses
//
void Solver::preprocessTheNurses() {
  // local variables for conciseness of the code
  //
  int nbNurses = pScenario_->nbNurses_, nbSkills = pScenario_->nbSkills_,
    nbDays = pDemand_->nbDays_;

  maxTotalStaff_ = 0;
  for (int n = 0; n < nbNurses; n++)	{
    LiveNurse* pNurse = theLiveNurses_[n];
    pNurse->maxWorkDays_ = 0;
    pNurse->minWorkDays_ = 0;

    // compute the maximum and minimum number of working days in the period of
    // the demand without getting any penalty for the number of consecutive
    // shifts
    // RqJO: this neglects the constraint of complete week-ends and the
    // preferences ; they should be added later
    //

    // first treat the history
    //
    // remaining number of days when building a maximum and minimum working days
    // periods
    int nbDaysMax=nbDays, nbDaysMin = nbDays;
    State* pState = &(pInitState_->at(n));
    if (pState->shift_ != 0) {
      // if the last shift was not a rest
      // keep working until maximum consecutive number of shifts for maximum
      // working days or until minimum consecutive number of shifts for minimum
      // working days
      pNurse->maxWorkDays_ = std::max(0, pNurse->maxConsDaysWork()-pState->consDaysWorked_);
      pNurse->minWorkDays_ = std::max(0, pNurse->minConsDaysWork()-pState->consDaysWorked_);

      // perform minimum rest for maximum working days and maximum rest for
      // minimum working days
      nbDaysMax = nbDays - pNurse->maxWorkDays_ - pNurse->minConsDaysOff();
      nbDaysMin = nbDays - pNurse->minWorkDays_ - pNurse->maxConsDaysOff();
    }
    else {
      // if the last shift was a rest
      // keep resting until minimum consecutive number of rests for maximum
      // working days or until maximum consecutive number of rests for minimum
      // working days
      nbDaysMax = nbDays - std::max(0, pNurse->minConsDaysOff()-pState->consDaysOff_);
      nbDaysMin = nbDays - std::max(0, pNurse->maxConsDaysWork()-pState->consDaysOff_);
    }

    // lengths of the stints maximizing and minimizing the percentage of working
    // days respectively
    //
    int lengthWorkStint = pNurse->maxConsDaysWork() + pNurse->minConsDaysOff();
    int lengthRestStint = pNurse->minConsDaysWork() + pNurse->maxConsDaysOff();

    // perform as many of these stints as possible
    //
    pNurse->maxWorkDays_ += (nbDaysMax/lengthWorkStint)*pNurse->maxConsDaysWork()+
      + std::min(nbDaysMax%lengthWorkStint, pNurse->maxConsDaysWork());
    pNurse->minWorkDays_ = (nbDaysMin/lengthRestStint)*pNurse->minConsDaysWork()+
      + std::min(nbDaysMin%lengthWorkStint, pNurse->minConsDaysWork());

    // add the working maximum number of working days to the maximum staffing
    //
    maxTotalStaffNoPenalty_ += pNurse->maxWorkDays_;
    maxTotalStaff_ += nbDays;
    for (int i = 0; i < pNurse->nbSkills_; i++) {
      // RqJO: the staffing per skill is very rough here since Nurses can have
      // multiple skills. A better data structure should be found.
      int sk = pNurse->skills_[i];
      maxStaffPerSkill_[sk] += nbDays;
      maxStaffPerSkillNoPenalty_[sk] += pNurse->maxWorkDays_;
    }
  }

  for (int sk = 0; sk < nbSkills; sk++) {
    std::cout << "Skill " << sk << " : " << maxStaffPerSkillNoPenalty_[sk] << std::endl;
  }

  // Compute the
}
