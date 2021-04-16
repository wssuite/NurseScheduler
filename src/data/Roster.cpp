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

#include "data/Roster.h"
#include "tools/Tools.h"

using std::vector;
using std::map;
using std::string;
using std::pair;

// Constructor: initialize planning from nothing
//
Roster::Roster(int nbDays, int firstDay, const PShift &pSDefault) {
  // initialize the roster with a rest week
  init(nbDays, firstDay, pSDefault);
}

// Constructor: initialize planning from an input set of shifts for the nurse
//
Roster::Roster(int nbDays, int firstDay, const vector<PShift> &shifts) :
    nbDays_(nbDays), firstDay_(firstDay), pShifts_(shifts) {
  // initialize the states at each day
  // states_.push_back(initialState);
  // for (int day = 0; day < nbDays_; day++) {
  //   State nextState;
  //   int newShiftType = pScenario->shiftIDToShiftTypeID_[shifts_[day]];
  //   nextState.addDayToState(states_[day], newShiftType, shifts_[day],
  //   pScenario_->timeDurationToWork_[shifts_[day]);
  //   states_.push_back(nextState);
  // }
}

// Constructor: initialize planning from an input set of shifts and skills
//
Roster::Roster(int nbDays,
               int firstDay,
               const std::vector<PShift> &shifts,
               const std::vector<int> &skills) :
    nbDays_(nbDays), firstDay_(firstDay), pShifts_(shifts) {
  // set the skill assignment
  for (int day = 0; day < nbDays_; day++) skills_.push_back(skills[day]);
}

// Destructor
Roster::~Roster() {}

// initialize the roster
//
void Roster::init(
    int nbDays, int firstDay, const PShift &pSDefault, int skillDefault) {
  nbDays_ = nbDays;
  firstDay_ = firstDay;
  pShifts_.resize(nbDays, pSDefault);
  skills_.resize(nbDays, skillDefault);
}

// re-inialize the roster
//
void Roster::reset(const PShift &pSDefault) {
  skills_.clear();
  pShifts_.clear();
  init(nbDays_, firstDay_, pSDefault);
}

// get a vector of consecutive states that will result from applying the
// the roster from a given initial state
//
vector<State> Roster::getStates(const State &stateIni, PScenario pScenario) {
  vector<State> states;

  // initialize the states at each day
  states.push_back(stateIni);
  for (int day = 0; day < nbDays_; day++) {
    State nextState;
    nextState.addDayToState(states[day], pShifts_[day]);
    states.push_back(nextState);
  }

  return states;
}

// assign a task at on a given day
//
void Roster::assignTask(int day, const PShift &pS, int skill) {
  pShifts_[day] = pS;
  skills_[day] = skill;
}

// add a roster at the end of the roster
//
void Roster::push_back(const Roster &roster) {
  if (!this->nbDays_) {
    this->copy(roster);
  } else if (!Tools::isSunday(firstDay_ + nbDays_ - 1)) {
    Tools::throwError("Roster::push_back: "
                      "The last day of the current roster is not a sunday!");
  } else {
    nbDays_ += roster.nbDays();
    for (int day = 0; day < roster.nbDays(); day++) {
      pShifts_.push_back(roster.pShift(day));
      skills_.push_back(roster.skill(day));
    }
  }

  if (!Tools::isSunday(firstDay_ + nbDays_ - 1)) {
    Tools::throwError("Roster::push_back: "
                      "The last day of the updated roster is not a sunday!");
  }
}

// copy the input roster
//
void Roster::copy(const Roster &roster) {
  firstDay_ = roster.firstDay_;
  nbDays_ = roster.nbDays_;
  pShifts_ = roster.pShifts_;
  skills_ = roster.skills_;
}
