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
Roster::Roster(int nbDays, int firstDay) {
  // initialize the roster with a rest week
  init(nbDays, firstDay);
}

// Constructor: initialize planning from an input set of shifts for the nurse
//
Roster::Roster(int nbDays, int firstDay, const vector<int> &shifts) :
    nbDays_(nbDays), firstDay_(firstDay), shifts_(shifts) {
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
               const std::vector<int> &shifts,
               const std::vector<int> &skills) :
    nbDays_(nbDays), firstDay_(firstDay), shifts_(shifts) {
  // set the skill assignment
  for (int day = 0; day < nbDays_; day++) skills_.push_back(skills[day]);
}

// Destructor
Roster::~Roster() {}

// initialize the roster
//
void Roster::init(int nbDays, int firstDay, int skillDefault) {
  nbDays_ = nbDays;
  firstDay_ = firstDay;

  // initialize the vectors of skills and shifts
  for (int day = 0; day < nbDays_; day++) skills_.push_back(skillDefault);
  for (int day = 0; day < nbDays_; day++) shifts_.push_back(0);
}

// re-inialize the roster
//
void Roster::reset() {
  skills_.clear();
  shifts_.clear();
  init(nbDays_, firstDay_);
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
    int shiftType = pScenario->shiftIDToShiftTypeID_[shifts_[day]];
    nextState.addDayToState(states[day],
                            shiftType,
                            shifts_[day],
                            pScenario_->timeDurationToWork_[shifts_[day]]);
    states.push_back(nextState);
  }

  return states;
}

// assign a task at on a given day
//
void Roster::assignTask(int day, int shift, int skill) {
  shifts_[day] = shift;
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
      shifts_.push_back(roster.shift(day));
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
  firstDay_ = roster.firstDay();
  nbDays_ = roster.nbDays();
  if (!shifts_.empty()) {
    shifts_.clear();
    skills_.clear();
  }
  for (int day = 0; day < nbDays_; day++) {
    shifts_.push_back(roster.shift(day));
    skills_.push_back(roster.skill(day));
  }
}
