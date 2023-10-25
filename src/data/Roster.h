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
//
//  C l a s s   R o s t e r
//
//  Schedule of a single nurse : it's a stretch with some skills
//
//-----------------------------------------------------------------------------
class Roster : public Stretch {
 public:
  // Default constructor
  Roster() : Stretch(), cost_(INFEAS_COST) {}

  // Constructor: initialize planning from an input set of shifts and skills
  Roster(int firstDay,
         std::vector<PShift> shifts,
         const std::vector<int> &skills);

  // Destructor
  ~Roster();

 private:
  // vector containing for each day the shift assigned to the nurse
  // the vector contains exactly one element per day
  // if the nurse is resting, the skill has no importance
  std::vector<int> skills_;
  double cost_;
  std::map<int, double> costPerType_;

 public:
  // Basic getters
  int skill(int day) const { return skills_[day]; }
  const vector<int> &skills() const { return skills_; }

  // re-initialize the roster
  void init(int firstDay, int nDays,
            const PShift &pSDefault, int skillDefault = 0);

  // clear the shifts and skills and fill it with these default values
  void reset(const PShift &pSDefault, int skillDefault = 0);

  // assign a task at on a given day
  void assignTask(int day, const PShift &pS, int skill = 0);

  // add a roster at the end of the roster
  void pushBack(const Roster &roster);

  // copy the input roster
  void copy(const Roster &roster);

  double cost() const;
  void cost(double cost);

  const std::map<int, double> &costs() const;
  void costs(const std::map<int, double> &costPerType);
  void cost(int type, double cost);

  // get a vector of consecutive states that will result from applying the
  // the roster from a given initial state
  std::vector<State> states(
      const State &pStateIni, const PScenario &pScenario);
};

#endif  // SRC_DATA_ROSTER_H_
