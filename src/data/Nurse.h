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
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "tools/Tools.h"
#include "data/Scenario.h"

class Scenario;

//-----------------------------------------------------------------------------
//
//  C l a s s   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------
class BaseResource {
 public:
  virtual BaseResource* clone() const = 0;
};
typedef std::shared_ptr<BaseResource> PBaseResource;

class Contract {
 public:
  // Id of the contract and index in the vector intToContract_
  //
  const int id_;

  // Name of the contract
  const std::string name_;

  // Minimum and maximum total number of shifts over the time period
  const int minTotalShifts_ = 0, maxTotalShifts_ = 99;
  const double costMinTotalShifts_ = 20, costMaxTotalShifts_ = 20;

  // Minimum and maximum number of consecutive days worked
  const int minConsDaysWork_ = 0, maxConsDaysWork_ = 99;
  const double costMinConsDaysWork_ = 30, costMaxConsDaysWork_ = 30;

  // Minimum and maximum number of consecutive days off
  const int minConsDaysOff_ = 0, maxConsDaysOff_ = 99;
  const double costMinConsDaysOff_ = 30, costMaxConsDaysOff_ = 30;

  // weekend constraints
  const Weekend weekend_;
  const int maxTotalWeekends_ = 99;
  const double costMaxTotalWeekends_ = 30;
  const bool needCompleteWeekends_ = false;
  const double costCompleteWeekends_ = 30;
  const bool identShiftTypesDuringWeekend_ = false;
  const bool noNightShiftBeforeFreeWeekend_ = false;

  // true if the resting periods must contain at least two days after a
  // series of night shifts
  const bool twoFreeDaysAfterNightShifts_ = false;

  // true if nurses can be assigned to a skill they do not possess
  const bool alternativeSkill_ = false;
  const double costAlternativeSkill_ = 0;

  // true if nurses can be assigned to a shift they cannot perform
  const bool alternativeShift_ = false;
  const double costAlternativeShift_ = 0;

  // vector of forbidden patterns
  vector<int> forbiddenPatternIds_;

  // vector of resources
  vector<PBaseResource> pBaseResources_;

  // Constructor and Destructor
  //
  Contract(int id,
           std::string name,
           int minTotalShifts,
           int maxTotalShifts,
           int minConsDaysWork,
           int maxConsDaysWork,
           int minConsDaysOff,
           int maxConsDaysOff,
           int maxTotalWeekends,
           int needCompleteWeekends,
           vector<PBaseResource> pBaseResources = {},
           DayOfWeek firstWeekendDayId = SATURDAY,
           DayOfWeek lastWeekendDayId = SUNDAY,
           bool identShiftTypesDuringWeekend = false,
           bool noNightShiftBeforeFreeWeekend = false,
           bool twoFreeDaysAfterNightShifts = false,
           bool alternativeSkill = false,
           double costAlternativeSkill = 0,
           bool alternativeShift = false,
           double costAlternativeShift = 0) :
      id_(id),
      name_(std::move(name)),
      minTotalShifts_(minTotalShifts),
      maxTotalShifts_(maxTotalShifts),
      minConsDaysWork_(minConsDaysWork),
      maxConsDaysWork_(maxConsDaysWork),
      minConsDaysOff_(minConsDaysOff),
      maxConsDaysOff_(maxConsDaysOff),
      weekend_(firstWeekendDayId, lastWeekendDayId),
      maxTotalWeekends_(maxTotalWeekends),
      needCompleteWeekends_(needCompleteWeekends),
      identShiftTypesDuringWeekend_(identShiftTypesDuringWeekend),
      noNightShiftBeforeFreeWeekend_(noNightShiftBeforeFreeWeekend),
      twoFreeDaysAfterNightShifts_(twoFreeDaysAfterNightShifts),
      alternativeSkill_(alternativeSkill),
      costAlternativeSkill_(costAlternativeSkill),
      alternativeShift_(alternativeShift),
      costAlternativeShift_(costAlternativeShift),
      pBaseResources_(std::move(pBaseResources)) {}

  Contract(int id,
           std::string name,
           vector<PBaseResource> pBaseResources = {},
           DayOfWeek firstWeekendDayId = SATURDAY,
           DayOfWeek lastWeekendDayId = SUNDAY,
           vector<int> patternIds = {},
           bool alternativeSkill = false,
           double costAlternativeSkill = 0,
           bool alternativeShift = false,
           double costAlternativeShift = 0) :
      id_(id),
      name_(std::move(name)),
      weekend_(firstWeekendDayId, lastWeekendDayId),
      alternativeSkill_(alternativeSkill),
      costAlternativeSkill_(costAlternativeSkill),
      alternativeShift_(alternativeShift),
      costAlternativeShift_(costAlternativeShift),
      forbiddenPatternIds_(std::move(patternIds)),
      pBaseResources_(std::move(pBaseResources)) {}

  // Information about the weekend
  bool isWeekend(const AbstractDay& d) const {
    return weekend_.isWeekend(d);
  }

  bool isWeekend(int dayId) const {
    return weekend_.isWeekend(dayId);
  }

  bool isFirstWeekendDay(int dayId) const {
    return weekend_.isFirstWeekendDay(dayId);
  }

  bool isLastWeekendDay(int dayId) const {
    return weekend_.isLastWeekendDay(dayId);
  }

  // Cost function that recover total penalty for a given number of
  // consecutive or total assignments
  double consDaysCost(int ndays) const;
  double consDaysOffCost(int ndays) const;
  double totalShiftCost(int ndays) const;
  double totalWeekendCost(int ndays) const;

  // Add a resource to the vector of PResources
  void addResource(const PBaseResource& pR) { pBaseResources_.push_back(pR); }

  // Display methods: toStringINRC2 + override operator<< (easier)
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
class Nurse;

class Position {
 public:
  // Constructor and Destructor
  //
  Position(int index,
           std::vector<int> skills,
           std::vector<int> alternativeSkills = {},
           double alternativeSkillsCost = 0);
  Position(int index,
           std::vector<int> skills,
           bool addOtherSkillsAsAlternative,
           double alternativeSkillsCost,
           int nSkills);
  Position(int index, const Nurse &nurse);

  ~Position() = default;

 public:
  // Index of the position
  //
  const int id_;

  // Vector of skills for this position.
  // For simplicity, the skill indices are sorted.
  //
  const std::vector<int> skills_;
  const std::vector<int> alternativeSkills_;
  const double alternativeSkillsCost_;

 private:
  // union of the sklls and alternative skills
  std::vector<int> allSkills_;

  // Positions that are below and above this one in the hierarchy
  // this is deduced from the dominance criterion implemented in compare()
  //
  std::vector<PPosition> positionsBelow_;
  std::vector<PPosition> positionsAbove_;

  // Rarity of the skills that appear in this position
  //
  std::vector<double> skillRarity_;

  // Rank of the position with regard to the dominance criterion in compare()
  // rank i contains all the positions that are dominated only by positions
  // with a rank smaller than i (the smaller rank is 0)
  //
  int rank_;

  // initalize rarity and check if skills are sorted
  void init();

 public:
  // basic getters
  //
  int id() const { return id_; }
  int nSkills() const { return skills_.size(); }
  int skill(int sk) const { return skills_[sk]; }
  const std::vector<int> &skills() const { return skills_; }
  int nAltSkills() const { return alternativeSkills_.size(); }
  const std::vector<int> &allSkills() const { return allSkills_; }
  int nBelow() const { return positionsBelow_.size(); }
  int nAbove() const { return positionsAbove_.size(); }
  PPosition positionsBelow(int i) const { return positionsBelow_[i]; }
  PPosition positionsAbove(int i) const { return positionsAbove_[i]; }
  double skillRarity(int sk) const { return skillRarity_[sk]; }
  int rank() const { return rank_; }

  // basic setters
  //
  void rank(int i) { rank_ = i; }

  // Display method: toStringINRC2
  //
  std::string toString() const;

  bool dominate(const Position &p) const;

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

  // return true  if this position corresponds to the one of the nurse
  //
  bool isNursePosition(const Nurse &nurse);

  // set positions above and below
  //
  void addBelow(const PPosition& pPosition);
  void addAbove(const PPosition& pPosition);

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
        int nSkills,
        std::vector<int> skills,
        std::vector<int> availableShifts,
        PContract contract);

  Nurse(int id,
        std::string name,
        int nShifts,
        int nSkills,
        std::vector<int> skills,
        std::vector<int> alternativeSkills,
        std::vector<int> availableShifts,
        PContract contract);

  Nurse(int id,
        std::string name,
        int nShifts,
        int nSkills,
        std::vector<int> skills,
        bool addOtherSkillsAsAlternative,
        std::vector<int> availableShifts,
        PContract contract);

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
  const std::vector<int> alternativeSkills_;

  // available shifts and alternative ones
  const std::vector<int> availableShifts_, alternativeShifts_;

  // Their contract type
  PContract pContract_;

 protected:
  //-----------------------------------------------------------------------------
  // Other constant characteristics of the nurses that could not be set in the
  // constructor
  // (only getters for these fields)
  //-----------------------------------------------------------------------------
  std::vector<bool> hasSkill_, isAvailOrAltShift_;

  // Resources of the nurse
  //
  vector<PBaseResource> pBaseResources_;

  // initialize vectors and checked they are sorted
  void init();


 public:
  // Basic getters
  //
  int nSkills() const { return skills_.size(); }
  int nAltSkills() const { return alternativeSkills_.size(); }
  int minTotalShifts() const { return pContract_->minTotalShifts_; }
  int maxTotalShifts() const { return pContract_->maxTotalShifts_; }
  int minConsDaysWork() const { return pContract_->minConsDaysWork_; }
  int maxConsDaysWork() const { return pContract_->maxConsDaysWork_; }
  int minConsDaysOff() const { return pContract_->minConsDaysOff_; }
  int maxConsDaysOff() const { return pContract_->maxConsDaysOff_; }
  int maxTotalWeekends() const { return pContract_->maxTotalWeekends_; }
  bool needCompleteWeekends() const {
    return  pContract_->needCompleteWeekends_; }
  bool alternativeSkillCategory() const {
    return pContract_->alternativeSkill_; }
  bool identShiftTypesDuringWeekend() const {
    return pContract_->identShiftTypesDuringWeekend_;  }
  bool noNightShiftBeforeFreeWeekEnd_() const {
    return pContract_->noNightShiftBeforeFreeWeekend_; }
  bool twoFreeDaysAfterNightShifts() const {
    return pContract_->twoFreeDaysAfterNightShifts_; }

  double consDaysCost(int n) const { return pContract_->consDaysCost(n); }
  double consDaysOffCost(int n) const { return pContract_->consDaysOffCost(n); }
  double totalShiftCost(int n) const { return pContract_->totalShiftCost(n); }
  double totalWeekendCost(int n) const {
    return pContract_->totalWeekendCost(n);
  }
  PContract pContract() const { return pContract_; }

  // Avanced getters
  //
  bool hasSkill(int skill) const;
  int contractId() const { return pContract_->id_; }
  std::string contractName() const { return pContract_->name_; }
  bool isShiftAvailOrAlt(int s) const { return isAvailOrAltShift_[s]; }
  bool isShiftNotAvailNorAlt(int s) const { return !isShiftAvailOrAlt(s); }

  // Add a resource to the vector
  void addBaseResource(const PBaseResource& pR);

  // Display methods: toStringINRC2
  //
  std::string toString() const;
};

// Compare nurses in ascending order of their ids
//
bool compareNursesById(const PNurse& n1, const PNurse& n2);

//-----------------------------------------------------------------------------
//
//  C l a s s   P r e f e r e n c e s
//
//  Describes the preferences of a nurse for a certain period of time
//  They are given as a vector (entry = nurseNum).
//  Each element is a map<int,set<int>> whose keys are the days,
//  and values are the sets of wished shift(s) OFF or ON on that day.
//
//-----------------------------------------------------------------------------

class Wish {
  std::vector<PAbstractShift> pAShiftsOff_, pAShiftsOn_;
  std::vector<double> costsOff_, costsOn_;

 public:
  Wish() = default;
  Wish(PAbstractShift pAS, double cost, bool off) {
    if (off) {
      pAShiftsOff_.push_back(std::move(pAS));
      costsOff_.push_back(cost);
    } else {
      pAShiftsOn_.push_back(std::move(pAS));
      costsOn_.push_back(cost);
    }
  }

  bool off() const { return !pAShiftsOff_.empty(); }

  bool dayOff() const {
    for (const auto &pAS : pAShiftsOff_)
      if (pAS->isAnyWork()) return true;
    return false;
  }

  bool on() const { return !pAShiftsOn_.empty(); }

  bool dayOn() const {
    for (const auto &pAS : pAShiftsOn_)
      if (pAS->isAnyWork()) return true;
    return false;
  }

  void add(const Wish &wish) {
    pAShiftsOff_.insert(pAShiftsOff_.end(),
                        wish.pAShiftsOff_.begin(), wish.pAShiftsOff_.end());
    costsOff_.insert(costsOff_.end(),
                     wish.costsOff_.begin(), wish.costsOff_.end());
    pAShiftsOn_.insert(pAShiftsOn_.end(),
                       wish.pAShiftsOn_.begin(), wish.pAShiftsOn_.end());
    costsOn_.insert(costsOn_.end(),
                    wish.costsOn_.begin(), wish.costsOn_.end());
  }

  // Penalize if the wish is violated
  // return the associated cost (default 0 if nothing is penalized)
  double cost(const PShift &pShift) const {
    double cost = 0;
    auto itC = costsOff_.begin();
    // OFF: add the costs of the violated wishes
    for (const auto &pAS : pAShiftsOff_) {
      if (pAS->includes(*pShift)) cost += *itC;
      ++itC;
    }
    // ON: add the costs of the violated wishes
    itC = costsOn_.begin();
    for (const auto &pAS : pAShiftsOn_) {
      if (!pAS->includes(*pShift)) cost += *itC;
      ++itC;
    }
    return cost;
  }

  std::string toString() const {
    std::stringstream rep;
    if (!pAShiftsOff_.empty()) {
      rep << "Off:";
      auto itC = costsOff_.begin();
      for (const auto &pAS : pAShiftsOff_)
         rep << (itC == costsOff_.begin() ? "" : ",") << pAS->name
             << "(" << *itC++ << ")";
    }
    if (!pAShiftsOn_.empty()) {
      rep << "On:";
      auto itC = costsOn_.begin();
      for (const auto &pAS : pAShiftsOn_)
        rep << (itC == costsOn_.begin() ? "" : ",") << pAS->name
            << "(" << *itC++ << ")";
    }
    return rep.str();
  }
};

class Preferences {
 public:
  // Constructor and destructor
  Preferences() = default;
  ~Preferences() = default;

  // Constructor with initialization to a given number of nurses
  Preferences(int nbNurses, int nbDays, int nbShifts);

  // Initialization with a map corresponding to the input nurses
  // and no wished Shift-Off.
  Preferences(const std::vector<PNurse> &pNurses, int nbDays, int nbShifts);

  // Returns the cost of the nurse wish for the shift
  static double wishCostOfTheShift(const std::map<int, Wish> &wishes,
                                   int day,
                                   const PShift &pS);

 protected:
  // Number of days considered in that case
  int nbDays_;

  // Total number of possible shifts
  int nbShifts_;

  // For each nurse, maps the day to the abstract shift that
  // he/she wants to have off/on
  std::map<int, std::map<int, Wish>> wishes_;

 public:
  // For a given day, and a given shift, adds it to the wish-list
  const Wish& addShift(int nurse, int day, const Wish &wish);

  const Wish& addShiftOff(
      int nurse, int day, const PAbstractShift &pShift, double cost);
  const Wish& addShiftOn(
      int nurse, int day, const PAbstractShift &pShift, double cost);

  // Adds the whole day to the wish-list
  const Wish& addDayOff(int nurse, int day, double cost);
  const Wish& addDayOn(int nurse, int day, double cost);

  const std::map<int, Wish> &nurseWishes(int id) const {
    return wishes_.at(id);
  }

  // Returns true if the shift respect the nurse wish
  bool wishTheShift(int nurse, int day, const PShift &pShift) const;

  // Returns the cost of the nurse wish for the shift
  double wishCostOfTheShift(
      int nurseNum, int day, const PShift &pShift) const;

  // Total number of shifts off that the nurse wants
  int howManyShiftsOff(int nurse) const;
  int howManyShiftsOn(int nurse) const;

  // Number of whole days off that the nurse wants
  int howManyDaysOff(int nurse, int dayMin, int dayMax) const;
  int howManyDaysOn(int nurse, int dayMin, int dayMax) const;

  // add another week preferences at the end of the current one
  void pushBack(const PPreferences& pPref);

  // Keep the preferences relative to the days in [begin,end)
  PPreferences keep(int begin, int end) const;

  // Keep the preferences relative to the nurses
  PPreferences keep(const std::vector<PNurse> &pNurses) const;

  // Remove the preferences relative to the nbDays first days
  PPreferences removeNFirstDays(int nbDays);

  // Display methods: toStringINRC2 + override operator<< (easier)
  //
  std::string toString() const;
  friend std::ostream &operator<<(std::ostream &outs,
                                  const Preferences &obj) {
    return outs << obj.toString();
  }

  static string toString(const std::map<int, Wish> &wishes);
};

#endif  // SRC_DATA_NURSE_H_
