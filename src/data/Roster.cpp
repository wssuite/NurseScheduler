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

#include <utility>
#include "tools/Tools.h"

using std::vector;
using std::map;
using std::string;
using std::pair;

// Constructor: initialize planning from an input set of shifts and skills
//
Roster::Roster(int firstDay,
               std::vector<PShift> shifts,
               const std::vector<int> &skills) :
    Stretch(firstDay, std::move(shifts)) {
  // set the skill assignment
  skills_ = skills;
  skills_.resize(nDays());
}

// Destructor
Roster::~Roster() = default;

// initialize the roster
void Roster::init(
    int firstDay, int nDays, const PShift &pSDefault, int skillDefault) {
  Stretch::init(firstDay, nDays, pSDefault);
  skills_.clear();
  skills_.resize(nDays, skillDefault);
}

// re-initialize the roster
void Roster::reset(const PShift &pSDefault, int skillDefault) {
  init(firstDayId_, nDays(), pSDefault, skillDefault);
}

// assign a task at on a given day
//
void Roster::assignTask(int day, const PShift &pS, int skill) {
  assignShift(day, pS);
  skills_[day] = skill;
}

// add a roster at the end of the roster
//
void Roster::pushBack(const Roster &roster) {
  if (!nDays()) {
    this->copy(roster);
  } else {
    Stretch::pushBack(roster);
    skills_.insert(skills_.end(), roster.skills_.begin(), roster.skills_.end());
  }
}

// copy the input roster
//
void Roster::copy(const Roster &roster) {
  Stretch::copy(roster);
  skills_ = roster.skills_;
}

// get a vector of consecutive states that will result from applying the
// the roster from a given initial state
//
vector<State> Roster::getStates(
    const State &stateIni, const PScenario &pScenario) {
  vector<State> states;

  // initialize the states at each day
  states.push_back(stateIni);
  for (int day = 0; day < nDays(); day++) {
    State nextState;
    nextState.addDayToState(states[day], pShifts_[day]);
    states.push_back(nextState);
  }

  return states;
}
