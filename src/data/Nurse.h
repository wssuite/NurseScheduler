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

#ifndef SRC_DATA_NURSE_H_
#define SRC_DATA_NURSE_H_

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "tools/Tools.h"
#include "Scenario.h"

class Scenario;

//-----------------------------------------------------------------------------
//
//  C l a s s   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------
class Contract {
 public:
  // Id of the contract and index in the vector intToContract_
  //
  const int id_;

  // Name of the contract
  //
  const std::string name_;

  // Minimum and maximum total number of shifts over the time period
  //
  const int minTotalShifts_, maxTotalShifts_;

  // Minimum and maximum number of consecutive days worked
  //
  const int minConsDaysWork_, maxConsDaysWork_;

  // Minimum and maximum number of consecutive days off
  //
  const int minConsDaysOff_, maxConsDaysOff_;

  // Maximum number of weekends worked, and complete weekend constraint
  //
  const int maxTotalWeekends_;
  const int needCompleteWeekends_;

  // weights
  const PWeights weights_;

  // Constructor and Destructor
  //
  Contract(int id, std::string name, int minTotalShifts, int maxTotalShifts,
           int minConsDaysWork, int maxConsDaysWork,
           int minConsDaysOff, int maxConsDaysOff,
           int maxTotalWeekends, int needCompleteWeekends, PWeights weights) :
      id_(id),
      name_(name),
      minTotalShifts_(minTotalShifts),
      maxTotalShifts_(maxTotalShifts),
      minConsDaysWork_(minConsDaysWork),
      maxConsDaysWork_(maxConsDaysWork),
      minConsDaysOff_(minConsDaysOff),
      maxConsDaysOff_(maxConsDaysOff),
      maxTotalWeekends_(maxTotalWeekends),
      needCompleteWeekends_(needCompleteWeekends),
      weights_(weights) {}

  // Cost function for consecutive identical shifts
  //
  double consDaysCost(int n) const;
  double consDaysOffCost(int n) const;
  double totalShiftCost(int n) const;
  double totalWeekendCost(int n) const;

  // Display methods: toString + override operator<< (easier)
  //
  std::string toString();
  friend std::ostream &operator<<(std::ostream &outs, Contract obj) {
    return outs << obj.toString();
  }
};

//-----------------------------------------------------------------------------
//
//  C l a s s   P o s i t i o n
//
//  A position (job) is a set of skills that a nurse may possess
//
//-----------------------------------------------------------------------------

class Position {
 public:
  // Constructor and Destructor
  //
  Position(int index, std::vector<int> skills);

  ~Position() {}

 public:
  // Index of the position
  //
  const int id_;

  // Vector of skills for this position.
  // For simplicity, the skill indices are sorted.
  //
  const std::vector<int> skills_;

 private:
  // Positions that are below and above this one in the hierarchy
  // this is deduced from the dominance criterion implemented in compare()
  //
  std::vector<PPosition> positionsBelow_;
  std::vector<PPosition> positionsAbove_;
  int nBelow_, nAbove_;

  // Rarity of the skills that appear in this position
  //
  std::vector<double> skillRarity_;

  // Rank of the position with regard to the dominance criterion in compare()
  // rank i contains all the positions that are dominated only by positions
  // with a rank smaller than i (the smaller rank is 0)
  //
  int rank_;

 public:
  // basic getters
  //
  int id() const { return id_; }
  int nSkills() const { return skills_.size(); }
  int skill(int sk) const { return skills_[sk]; }
  const std::vector<int> &skills() const { return skills_; }
  int nBelow() const { return nBelow_; }
  int nAbove() const { return nAbove_; }
  PPosition positionsBelow(int i) const { return positionsBelow_[i]; }
  PPosition positionsAbove(int i) const { return positionsAbove_[i]; }
  double skillRarity(int sk) const { return skillRarity_[sk]; }
  int rank() const { return rank_; }

  // basic setters
  //
  void rank(int i) { rank_ = i; }

  // Display method: toString
  //
  std::string toString() const;

  // Compare this position with the input position
  // The dominance criterion is that a position p1 with skills sk1 dominates p2
  // with skills sk2 if and only if (sk1 contains sk2) and sk1 has more skills
  // than sk2
  // The function returns 1 if this position dominates, -1 if it is dominated
  // and 0 if there is no dominance
  //
  int compare(const Position &p);

  // returns true if the position shares at least one skill
  // with the input position
  //
  bool shareSkill(const Position &p) const;

  // set positions above and below
  //
  void addBelow(PPosition pPosition);
  void addAbove(PPosition pPosition);

  // reset the list of positions below and above
  //
  void resetAbove();
  void resetBelow();

  // update the rarity of the skills
  // the input is the vector of the rarity of all the skills
  // the vector is sorted without record of the corresponding skill because it
  // is used only to compare two positions with the same rank
  //
  void updateRarities(std::vector<double> allRarities);
};



//-----------------------------------------------------------------------------
//
//  C l a s s   N u r s e
//
//  Class that contains all the attributes describing the characteristics and
//  the planning of each nurse
//
//-----------------------------------------------------------------------------

class Nurse {
 public:
  // Constructor and destructor
  Nurse(int id,
        std::string name,
        int nShifts,
        std::vector<int> skills,
        std::vector<int> availableShifts,
        PConstContract contract);

  Nurse(int id, const Nurse &nurse);

  ~Nurse();

  // the constant attributes of the nurses are public
 public:
  //-----------------------------------------------------------------------------
  // Constant characteristics of the nurses (no set method)
  //-----------------------------------------------------------------------------
  // Id of the nurse within a scenario
  // (=index in the vector<Nurse> theNurse of the Scenario)
  //
  const int num_;

  // name of the nurse
  //
  const std::string name_;

  // number of skills and vector of the skills indices
  // for simplicity, the vector of skills is sorted
  //
  const std::vector<int> skills_;

  // available shifts
  const std::vector<int> availableShifts_;

  // Her contract type
  //
  PConstContract pContract_;

 protected:
  //-----------------------------------------------------------------------------
  // Other constant characteristics of the nurses that could not be set in the
  // constructor
  // (only getters for these fields)
  //-----------------------------------------------------------------------------
  std::vector<bool> hasSkill_, isAvailableShifts_;

 public:
  // Basic getters
  //
  int nSkills() const { return skills_.size(); }
  int minTotalShifts() const { return pContract_->minTotalShifts_; }
  int maxTotalShifts() const { return pContract_->maxTotalShifts_; }
  int minConsDaysWork() const { return pContract_->minConsDaysWork_; }
  int maxConsDaysWork() const { return pContract_->maxConsDaysWork_; }
  int minConsDaysOff() const { return pContract_->minConsDaysOff_; }
  int maxConsDaysOff() const { return pContract_->maxConsDaysOff_; }
  int maxTotalWeekends() const { return pContract_->maxTotalWeekends_; }
  int needCompleteWeekends() const { return pContract_->needCompleteWeekends_; }

  double consDaysCost(int n) const { return pContract_->consDaysCost(n); }
  double consDaysOffCost(int n) const { return pContract_->consDaysOffCost(n); }
  double totalShiftCost(int n) const { return pContract_->totalShiftCost(n); }
  double totalWeekendCost(int n) const {
    return pContract_->totalWeekendCost(n);
  }

  // Avanced getters
  //
  bool hasSkill(int skill) const;
  std::string contractName() { return pContract_->name_; }
  bool isShiftAvailable(int s) const { return isAvailableShifts_[s]; }

  // Display methods: toString
  //
  std::string toString() const;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   P r e f e r e n c e s
//
//  Describes the preferences of a nurse for a certain period of time
//  They are given as a vector (entry = nurseNum).
//  Each element is a map<int,set<int>> whose keys are the days,
//  and values are the sets of wished shift(s) OFF on that day.
//
//-----------------------------------------------------------------------------

struct Wish {
  int shift;
  PREF_LEVEL level;
};

class Preferences {
 public:
  // Constructor and destructor
  Preferences();
  ~Preferences();

  // Constructor with initialization to a given number of nurses
  Preferences(int nbNurses, int nbDays, int nbShifts);

  // Initialization with a map corresponding to the input nurses
  // and no wished Shift-Off.
  Preferences(const std::vector<PNurse> &pNurses, int nbDays, int nbShifts);

  static int wishLevel(const std::map<int, std::vector<Wish> > &wishes,
                       int day,
                       int shift);

 protected:
  // Number of nurses
  //
  int nbNurses_;

  // Number of days considered in that case
  //
  int nbDays_;

  // Total number of possible shifts
  //
  int nbShifts_;

  // For each nurse, maps the day to the set of shifts
  // that he/she wants to have off
  //
  std::map<int, std::map<int, std::vector<Wish> > > wishesOff_;
  std::map<int, std::map<int, std::vector<Wish> > > wishesOn_;

 public:
  // For a given day, and a given shift, adds it to the wish-list for OFF-SHIFT
  void addShiftOff(int nurse, int day, int shift, PREF_LEVEL level);
  void addShiftOn(int nurse, int day, int shift, PREF_LEVEL level);

  // Adds the whole day to the wish-list
  void addDayOff(int nurse, int day, PREF_LEVEL level);
  void addDayOn(int nurse, int day, PREF_LEVEL level);

  const std::map<int, std::vector<Wish>> &nurseWishesOff(int id) const {
    return wishesOff_.at(id);
  }
  const std::map<int, std::vector<Wish>> &nurseWishesOn(int id) const {
    return wishesOn_.at(id);
  }

  // True if the nurses wants that shift off
  bool wantsTheShiftOff(int nurse, int day, int shift) const;
  bool wantsTheShiftOn(int nurse, int day, int shift) const;

  // Returns level if the nurse wants that shift off : -1 otherwise
  int wantsTheShiftOffLevel(int nurseNum, int day, int shift) const;
  int wantsTheShiftOnLevel(int nurseNum, int day, int shift) const;

  // True if the nurses wants the whole day off
  bool wantsTheDayOff(int nurse, int day) const;
  bool wantsTheDayOn(int nurse, int day) const;

  // Total number of shifts off that the nurse wants
  int howManyShiftsOff(int nurse) const;
  int howManyShiftsOn(int nurse) const;

  // Number of whole days off that the nurse wants
  int howManyDaysOff(int nurse, int dayMin, int dayMax) const;
  int howManyDaysOn(int nurse, int dayMin, int dayMax) const;

  // add another week preferences at the end of the current one
  //
  void pushBack(PPreferences pPref);

  // Keep the preferences relative to the days in [begin,end)
  PPreferences keep(int begin, int end);

  // Remove the preferences relative to the nbDays first days
  PPreferences removeNFirstDays(int nbDays);

  // Display methods: toString + override operator<< (easier)
  //
  std::string toString(PScenario pScenario = nullptr) const;
  friend std::ostream &operator<<(std::ostream &outs,
                                  const Preferences &obj) {
    return outs << obj.toString();
  }
};

#endif  // SRC_DATA_NURSE_H_
