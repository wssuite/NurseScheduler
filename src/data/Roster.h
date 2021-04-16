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

#ifndef SRC_DATA_ROSTER_H_
#define SRC_DATA_ROSTER_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "tools/Tools.h"
#include "data/Scenario.h"
#include "data/Nurse.h"

//-----------------------------------------------------------------------------
// struct task:
// a task is a shift performed on a given day with a given skill
//
//-----------------------------------------------------------------------------
typedef struct { int shift; int skill; } task;

//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r
//
//  Schedule of a single nurse
//
//-----------------------------------------------------------------------------

class Roster {
 public:
  // Default constructor
  //
  Roster() : nbDays_(0), firstDay_(0) {}

  // Constructor form no particular planning
  //
  Roster(int nbDays, int firstDay, const PShift &pSDefault);

  // Constructor: initialize planning from an input set of shifts for the nurse
  //
  Roster(int nbDays, int firstDay, const std::vector<PShift> &shifts);

  // Constructor: initialize planning from an input set of shifts and skills
  //
  Roster(int nbDays,
         int firstDay,
         const std::vector<PShift> &shifts,
         const std::vector<int> &skills);

  // Destructor
  ~Roster();

 private:
  // number of days in the roster and index of the first day
  //
  int nbDays_, firstDay_;

  // pointer to the scenario, the nurse under consideration and her wishes in
  // terms of days off
  // (the key of the map is the day and the value is the set of wishes)
  //
  PScenario pScenario_;

  // vector containing for each day the shift assigned to the nurse
  // the vector contains exactly one element per day
  // the shift 0 corresponds to a rest
  //
  std::vector<PShift> pShifts_;

  // vector containing for each day the shift assigned to the nurse
  // the vector contains exactly one element per day
  // if the nurse is resting, the skill has no importance
  //
  std::vector<int> skills_;

 public:
  // Basic getters
  int firstDay() const { return firstDay_; }
  int nbDays() const { return nbDays_; }
  const PShift & pShift(int day) const { return pShifts_[day]; }
  const vector<PShift> & pShifts() const { return pShifts_; }
  int skill(int day) const { return skills_[day]; }
  const vector<int> &skills() const { return skills_; }

  // initialize the roster
  //
  void init(
      int nbDays, int firstDay, const PShift &pSDefault, int shiftDefault = 0);

  // re-inialize the roster
  //
  void reset(const PShift &pSDefault);

  // get a vector of consecutive states that will result from applying the
  // the roster from a given initial state
  //
  std::vector<State> getStates(const State &pStateIni, PScenario pScenario);

  // assign a task at on a given day
  //
  void assignTask(int day, const PShift &pS, int skill = 0);

  // add a roster at the end of the roster
  //
  void push_back(const Roster &roster);

  // copy the input roster
  //
  void copy(const Roster &roster);
};

#endif  // SRC_DATA_ROSTER_H_
