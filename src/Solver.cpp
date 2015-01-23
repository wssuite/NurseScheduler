/*
* Solver.cpp
*
*  Created on: 22 déc. 2014
*      Author: jeremy
*/


#include "Solver.h"


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
LiveNurse::LiveNurse(const Nurse& nurse):
Nurse(nurse.id_, nurse.name_, nurse.nbSkills_, nurse.skills_, nurse.pContract_){
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
      theLiveNurses_.push_back(new LiveNurse( (pScenario_->theNurses_[i]) ) );
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
