/*
* Solver.cpp
*
*  Created on: 22 déc. 2014
*      Author: jeremy
*/


#include "Solver.h"


// Specific constructor
Solver::Solver(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState):
  pScenario_(pScenario),  pDemand_(pDemand),
  pPreferences_(pPreferences), pInitState_(pInitState) {

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

  maxTotalStaff_ == nbNurses * nbDays;
  for (int n = 0; n < nbNurses; n++)	{
    maxWorkDays_.push_back(0);
    minWorkDays_.push_back(0);
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
      maxWorkDays_[n] = std::max(0, pScenario_->maxConsDaysWorkOf(n)-pState->consDaysWorked_);
      minWorkDays_[n] = std::max(0, pScenario_->minConsDaysWorkOf(n)-pState->consDaysWorked_);

      // perform minimum rest for maximum working days and maximum rest for
      // minimum working days
      nbDaysMax = nbDays - maxWorkDays_[n] - pScenario_->minConsDaysOffOf(n);
      nbDaysMin = nbDays - minWorkDays_[n] - pScenario_->maxConsDaysOffOf(n);
    }
    else {
      // if the last shift was a rest
      // keep resting until minimum consecutive number of rests for maximum
      // working days or until maximum consecutive number of rests for minimum
      // working days
      nbDaysMax = nbDays - std::max(0, pScenario_->minConsDaysOffOf(n)-pState->consDaysOff_);
      nbDaysMin = nbDays - std::max(0, pScenario_->maxConsDaysWorkOf(n)-pState->consDaysOff_);
    }

    // lengths of the stints maximizing and minimizing the percentage of working
    // days respectively
    //
    int lengthWorkStint = pScenario_->maxConsDaysWorkOf(n)+pScenario_->minConsDaysOffOf(n);
    int lengthRestStint = pScenario_->minConsDaysWorkOf(n)+pScenario_->maxConsDaysOffOf(n);

    // perform as many of these stints as possible
    //
    maxWorkDays_[n] += (nbDaysMax/lengthWorkStint)*pScenario_->maxConsDaysWorkOf(n)+
      + std::min(nbDaysMax%lengthWorkStint, pScenario_->maxConsDaysWorkOf(n));
    minWorkDays_[n] = (nbDaysMin/lengthRestStint)*pScenario_->minConsDaysWorkOf(n)+
      + std::min(nbDaysMin%lengthWorkStint, pScenario_->minConsDaysWorkOf(n));

    // add the working maximum number of working days to the maximum staffing
    //
    maxTotalStaffNoPenalty_ += maxWorkDays_[n];
    maxTotalStaff_ += nbDays;
    for (int i = 0; i < pScenario_->theNurses_[n].nbSkills_; i++) {
      int sk = pScenario_->theNurses_[n].skills_[i];
      maxStaffPerSkill_[sk] += nbDays;
      maxStaffPerSkillNoPenalty_[sk] += maxWorkDays_[n];
      }
  }

  // Compute the
}
