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

#include "data/Nurse.h"

#include <iostream>
#include <algorithm>
#include <memory>
#include <utility>
#include <functional>


using std::vector;
using std::map;
using std::string;
using std::pair;


// Compare nurses in ascending order of their ids
//
bool compareNursesById(const PNurse& n1, const PNurse& n2) {
  return (n1->num_ < n2->num_);
}

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------

// Cost function for consecutive identical shifts
//
double Contract::consDaysCost(int ndays) const {
  if (minConsDaysWork_ - ndays > 0)
    return (costMinConsDaysWork_ * (minConsDaysWork_ - ndays));
  if (ndays - maxConsDaysWork_ > 0)
    return (costMaxConsDaysWork_ * (ndays - maxConsDaysWork_));
  return 0;
}

double Contract::consDaysOffCost(int ndays) const {
  if (minConsDaysOff_ - ndays > 0)
    return (costMinConsDaysWork_ * (minConsDaysOff_ - ndays));
  if (ndays - maxConsDaysOff_ > 0)
    return (costMaxConsDaysWork_ * (ndays - maxConsDaysOff_));
  return 0;
}

// Print method
//
std::string Contract::toString() {
  std::stringstream rep;
  rep << name_ << "  -  ";
  rep << "Tot:" << minTotalShifts_ << "<" << maxTotalShifts_ << "  |  ";
  rep << "Work:" << minConsDaysWork_ << "<" << maxConsDaysWork_ << "  |  ";
  rep << "Rest:" << minConsDaysOff_ << "<" << maxConsDaysOff_ << "  |  ";
  rep << "WE:" << maxTotalWeekends_ << "  |  ";
  if (!needCompleteWeekends_) rep << "NOT";
  rep << "complete";
  return rep.str();
}


//-----------------------------------------------------------------------------
//
//  C l a s s   P o s i t i o n
//
//  A position (job) is a set of skills that a nurse may possess
//
//-----------------------------------------------------------------------------

SkillsSet::SkillsSet(
        int nSkills, std::vector<int> skills, const Contract &contract):
skills_(std::move(skills)), hasSkill_(nSkills, false) {
  int i = 0;
  for (int altSk : contract.alternativeSkills_) {
    // if not found as a skill, add alternative skill
    if (std::find(skills_.begin(), skills_.end(), altSk) == skills_.end()) {
      double c = contract.costAlternativeSkills_.at(i);
      // add only if soft
      if (isSoftCost(c)) {
        alternativeSkills_.push_back(altSk);
        alternativeSkillsCosts_.push_back(
                contract.costAlternativeSkills_.at(i));
      }
    }
    i++;
  }

  // Verify that the vector of skills is sorted
  for (auto it = skills_.begin(); it+1 < skills_.end(); ++it)
    if (*it >= *(it + 1))
      Tools::throwError("The skills in a nurse are not sorted "
                        "or some skill is repeated!");
  // Verify that the vector of skills is sorted
  for (auto it = alternativeSkills_.begin();
       it+1 < alternativeSkills_.end(); ++it)
    if (*it >= *(it + 1))
      Tools::throwError("The alternative skills in a nurse are not sorted "
                        "or some skill is repeated!");

  allSkills_ = Tools::appendVectors(skills_, alternativeSkills_);
  std::sort(allSkills_.begin(), allSkills_.end());

  // build hasSkill_
  for (int sk : allSkills_)
    hasSkill_[sk] = true;
}

bool SkillsSet::operator==(const SkillsSet &skillsSet) const {
  if (nSkills() == skillsSet.nSkills() &&
      nAltSkills() == skillsSet.nAltSkills()) {
    // check if same skills
    if (Tools::includes(skills_, skillsSet.skills_) &&
        Tools::includes(alternativeSkills_, skillsSet.alternativeSkills_)) {
      // check if same cost
      for (int i = 0; i < nAltSkills(); i++)
        if (abs(alternativeSkillsCosts_[i] -
                skillsSet.skillCost(alternativeSkills_[i])) > 1e-3)
          return false;
      return true;
    }
  }
  return false;
}

// Constructor and Destructor
Position::Position(int index,
                   int nSkills,
                   std::vector<int> skills,
                   const Contract &contract) :
        SkillsSet(nSkills, std::move(skills), contract),
        id_(index), rank_(0) {
  init();
}

Position::Position(int index, SkillsSet skillsSet) :
        SkillsSet(std::move(skillsSet)),
        id_(index), rank_(0) {
  init();
}

// initialize rarity and check if skills are sorted
void Position::init() {
  for (int sk = 0; sk < nSkills(); sk++)
    skillRarity_.push_back(1.0);
}

// Set positions above and below
//
void Position::addBelow(const PPosition& pPosition) {
  positionsBelow_.push_back(pPosition);
}
void Position::addAbove(const PPosition& pPosition) {
  positionsAbove_.push_back(pPosition);
}

// Print method
//
std::string Position::toString() const {
  std::stringstream rep;
  rep << id_ << ": ";
  if (id_ < 10) rep << " ";
  if (id_ < 100) rep << " ";

  // print the rank and the list of skills
  rep << ": rank = " << rank_;
  rep << " ; skills = [ ";
  for (int sk : allSkills()) rep << sk << " ";
  rep << "]\t";

  // print the rank
  // print the list of positions above and below this one
  if (!positionsBelow_.empty()) {
    rep << std::endl;
    rep << "#\t\t\t\t\tdominates positions:      ";
    for (const auto& position : positionsBelow_) {
      rep << "\t" << position->id_;
    }
  }
  if (!positionsAbove_.empty()) {
    rep << std::endl;
    rep << "#\t\t\t\t\tis dominated by positions:";
    for (const auto& position : positionsAbove_) {
      rep << "\t" << position->id_;
    }
  }
  return rep.str();
}

// update the rarity of the skills
// the input is the vector of the rarity of all the skills
// the vector is sorted without record of the corresponding skill because it
// is used only to compare two positions with the same rank
//
void Position::updateRarities(std::vector<double> allRarities) {
  for (int sk = 0; sk < this->nSkills(); sk++) {
    skillRarity_[sk] = allRarities[skill(sk)];
  }
  std::sort(skillRarity_.begin(), skillRarity_.end(), std::greater<>());
}

// Check if skills of p are included in the position skills and
// if all the skills of p are also in all the position skills
bool Position::dominate(const Position &p) const {
  if (!Tools::includes(skills(), p.skills()))
    return false;
  // otherwise look at alternative skills
  // check if all skills of p are included in allSkills
  return Tools::includes(allSkills(), p.altSkills());
}

// Compare this position with the input position
// The dominance criterion is that a position p1 with skills sk1 dominates p2
// with skills sk2 if and only if (sk1 contains sk2) and sk1 has more skills
// than sk2
// The function returns 1 if this position dominates, -1 if it is dominated
// and 0 if there is no dominance
//
int Position::compare(const Position &p) {
  // no possible dominance if both positions have as many skills
  if (p.nAllSkills() == this->nAllSkills()) {
    return 0;
  } else if (p.nAllSkills() > this->nAllSkills()) {
    // only p can dominate if it has more skills
    // the comparison of the two skill lists is based on the fact that they are
    // both sorted
    if (p.dominate(*this))
      return -1;
    return 0;
  } else if (this->nAllSkills() > p.nAllSkills()) {
    // only this position can dominate if it has more skills
    if (dominate(p))
      return 1;
    return 0;
  }
  return 0;
}

// returns true if the position shares at least
// one skill with the input position
bool Position::shareSkill(const Position &p) const {
  for (int s : allSkills())
    if (find(p.allSkills().begin(), p.allSkills().end(), s) !=
        p.allSkills().end())
      return true;
  return false;
}

// return true if this position corresponds to the one of the nurse
bool Position::isNursePosition(const Nurse &nurse) {
  return *this == nurse;
}

// reset the list of positions below and above
//
void Position::resetAbove() {
  positionsAbove_.clear();
}
void Position::resetBelow() {
  positionsBelow_.clear();
}

//------------------------------------------------------------------------
// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//------------------------------------------------------------------------
bool comparePositions(const Position &p1, const Position &p2) {
  if (p1.id() == p2.id()) {
    return false;
  } else if (p1.rank() == p2.rank()) {
    // the skillRarity vector is ALWAYS sorted in descending order, because the
    // updateRarities is the only setter for skillRarity and it sorts the vector
    for (int sk = 0; sk < std::min(p1.nAllSkills(), p2.nAllSkills()); sk++)
      if (p1.skillRarity(sk) != p2.skillRarity(sk))
        return p1.skillRarity(sk) > p2.skillRarity(sk);
    return p1.nSkills() > p2.nSkills();
  } else {
    return p1.rank() < p2.rank();
  }
  return true;
}


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   P r e f e r e n c e
//
//-----------------------------------------------------------------------------

// Initialization with a map of size nNurses with no wished Shift-Off.
Preferences::Preferences(int nbNurses, int nbDays, int nbShifts):
    nbDays_(nbDays), nbShifts_(nbShifts) {
  for (int n=0; n < nbNurses; ++n)
    wishes_[n];
}

// Initialization with a map corresponding to the input nurses
// and no wished Shift-Off.
Preferences::Preferences(
    const vector<PNurse> &nurses, int nbDays, int nbShifts):
    nbDays_(nbDays), nbShifts_(nbShifts) {
  for (const PNurse &pN : nurses)
    wishes_[pN->num_];
}

// Returns the cost of the nurse wish for the shift
double Preferences::wishCostOfTheShift(const std::map<int, Wish> &wishes,
                                       int day,
                                       const PShift &pS) {
  // If the day is not in the wish-list, return 0
  auto itM = wishes.find(day);
  if (itM == wishes.end())
    return 0;
  // return the cost of the shift for the wish
  return itM->second.cost(pS);
}


// add a wish for the nurse on a given day
const Wish& Preferences::addShift(int nurseNum,
                                  int day,
                                  const Wish &wish) {
  wishes_.at(nurseNum)[day].add(wish);
  return wishes_.at(nurseNum)[day];
}

// Add a wished day-shift off for a nurse
const Wish& Preferences::addShiftOff(int nurseNum,
                                     int day,
                                     const PAbstractShift &pAShift,
                                     double cost) {
  return addShift(nurseNum, day, Wish(pAShift, true, cost));
}

// Returns the cost of the nurse wish for the shift
double Preferences::wishCostOfTheShift(
    int nurseNum, int day, const PShift &pShift) const {
  return wishCostOfTheShift(wishes_.at(nurseNum), day, pShift);
}

// Total number of shifts off that the nurse wants
int Preferences::howManyShiftsOff(int nurseNum) const {
  int nOff = 0;
  for (const auto &p : wishes_.at(nurseNum))
    nOff += p.second.off();
  return nOff;
}

// Number of whole days off that the nurse wants
int Preferences::howManyDaysOff(int nurseNum, int dayMin, int dayMax) const {
  int nbDayOff = 0;
  // look at every wishes of the nurse
  for (const auto &p : wishes_.at(nurseNum))
    nbDayOff += (p.second.dayOff() && (p.first >= dayMin)
                 && (p.first <= dayMax));
  return nbDayOff;
}

//////////////////////////////////////////////////////////////////

// Add a wished day-shift on for a nurse
const Wish& Preferences::addShiftOn(int nurseNum,
                                    int day,
                                    const PAbstractShift &pAShift,
                                    double cost) {
  return addShift(nurseNum, day, Wish(pAShift, false, cost));
}

// Total number of shifts on that the nurse wants
int Preferences::howManyShiftsOn(int nurseNum) const {
  int nOn = 0;
  for (const auto &p : wishes_.at(nurseNum))
    nOn += p.second.on();
  return nOn;
}

// Number of whole days on that the nurse wants
int Preferences::howManyDaysOn(int nurseNum, int dayMin, int dayMax) const {
  int nbDayOn = 0;
  // look at every wishes of the nurse
  for (const auto &p : wishes_.at(nurseNum))
    nbDayOn += (p.second.dayOn() && (p.first >= dayMin) && (p.first <= dayMax));
  return nbDayOn;
}

// add another week preferences at the end of the current one
void Preferences::pushBack(const PPreferences& pPref) {
  // check if same scenario
  if ((nbShifts_ != pPref->nbShifts_)) {
    string error = "Preferences are not compatible";
    Tools::throwError(error.c_str());
  }

  // update the whishes
  for (const auto &p : wishes_)
    for (const auto &wish : pPref->wishes_.at(p.first))
      addShift(p.first, nbDays_ + wish.first, wish.second);

  // update the number of days
  nbDays_ += pPref->nbDays_;
}

// create a new preferences that contains the current preferences and
// another week preferences at the end
PPreferences Preferences::append(const PPreferences& pPref) const {
  PPreferences pNewPref = std::make_shared<Preferences>(*this);
  pNewPref->pushBack(pPref);
  return pNewPref;
}

// K the preferences relative to the nbDays first days
PPreferences Preferences::keep(int begin, int end) const {
  PPreferences pPref = std::make_shared<Preferences>();
  for (const auto &p : wishes_)
    for (const auto &wish : wishes_.at(p.first))
      if (wish.first >= begin && wish.first < end)
        pPref->addShift(p.first, wish.first, wish.second);

  return pPref;
}

// Keep the preferences relative to the nurses
PPreferences Preferences::keep(const std::vector<PNurse> &pNurses) const {
  PPreferences pPref = std::make_shared<Preferences>();
  for (const PNurse& pNurse : pNurses)
    pPref->wishes_[pNurse->num_] = wishes_.at(pNurse->num_);
  return pPref;
}

// Remove the preferences relative to the nbDays first days
PPreferences Preferences::removeNFirstDays(int nbDays) {
  return keep(nbDays, nbDays_);
}

// Display method: toString()
//
string Preferences::toString() const {
  std::stringstream rep;
  rep << "# Preferences:" << std::endl;
  for (const auto &p : wishes_) {
    rep << "#       " << p.first << ": " << toString(p.second);
  }
  rep << std::endl;

  return rep.str();
}

string Preferences::toString(const std::map<int, Wish> &wishes) {
  std::stringstream rep;
  for (const auto &p : wishes)
    rep << p.first << "  ->  " << p.second.toString() << "   ";
  rep << std::endl;
  return rep.str();
}



//-----------------------------------------------------------------------------
//
//  C l a s s   N u r s e
//
//  Class that contains all the attributes describing the characteristics and
//  the planning of each nurse
//
//-----------------------------------------------------------------------------


// Constructor and destructor
// Note : need both with const Contract and (non-const) Contract
// because non-const is used in our code, and const is needed so that
// we can override the operator= and have vector<Nurse>.
// We need to override it because vector members should have some properties
// (assignable a.o., which implies non-const)
//
Nurse::Nurse(int id,
             string name,
             int nSkills,
             vector<int> skills,
             std::vector<int> availableShifts,
             PContract pContract) :
    SkillsSet(nSkills, std::move(skills), *pContract),
    num_(id),
    name_(std::move(name)),
    availableShifts_(std::move(availableShifts)),
    pContract_(std::move(pContract)) {
  init();
}

Nurse::Nurse(int id, const Nurse &nurse) :
    SkillsSet(nurse),
    num_(id),
    name_(nurse.name_),
    availableShifts_(nurse.availableShifts_),
    pContract_(nurse.pContract_),
    pBaseResources_(nurse.pBaseResources_),
    isShiftAvail_(nurse.isShiftAvail_) {}

Nurse::~Nurse() = default;

void Nurse::init() {
  auto it = ++availableShifts_.begin();
  for (; it < availableShifts_.end(); ++it)
    if (*it < *(it - 1))
      Tools::throwError("The available shifts in a nurse are not sorted "
                        "or at least a shift is repeated!");

  isShiftAvail_.resize(availableShifts_.back()+1);
  for (int i : availableShifts_)
    isShiftAvail_[i] = true;

  pBaseResources_ = pContract_->pBaseResources_;
}

// Print method
//
string Nurse::toString() const {
  std::stringstream rep;
  rep << num_ << ": ";
  if (num_ < 10) rep << " ";
  if (num_ < 100) rep << " ";
  rep << name_ << "\t" << nSkills() << " [ ";
  for (int sk : skills()) rep << sk << " ";
  rep << "]\t";
  if (nSkills() == 1) rep << "\t";
  rep << name_ << "\t" << nAltSkills() << " [ ";
  for (int sk : altSkills()) rep << sk << " ";
  rep << "]\t";
  if (nAltSkills() == 1) rep << "\t";
  rep << pContract_->name_;
  return rep.str();
}

void Nurse::addBaseResource(const PBaseResource& pR) {
  pBaseResources_.push_back(pR);
}
