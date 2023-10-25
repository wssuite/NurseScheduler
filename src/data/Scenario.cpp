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

#include "Scenario.h"

#include <algorithm>
#include <utility>
#include <sstream>

#include "data/Nurse.h"
#include "tools/Tools.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//-----------------------------------------------------------------------------

// return a state for day d with None shift
State State::noneState(const ShiftsFactory &sF, int d) {
  return State(sF.pNoneShift()->pIncludedShifts().front(), d);
}

// Destructor
State::~State() = default;

// Function that appends a new day worked on a given shiftType to an input state
// to update this state
// RqJO: I slghtly modified the method to take into account the possibility to
// add in the state that no task has been assigned on this day
//
void State::addDayToState(const State &prevState, const PShift &pS) {
  //  int  timeWorked = 1;               // days worked
  //  int  timeWorked = timeDurationToWork_[newShift];      // hours worked

  // Total shifts worked if it is a worked day
  totalTimeWorked_ = prevState.totalTimeWorked_ + pS->duration;
  totalDaysWorked_ = prevState.totalDaysWorked_ + pS->isWork();

  // index of the previous shift
  int prevShiftType = -1;
  if (prevState.pShift_)
    prevShiftType = prevState.pShift_->type;

  // Treat the case in which no shift is assigned to the nurse on this day
  if (pS->type < 0) {
    totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
    consShifts_ = prevState.consShifts_;
    consDaysWorked_ = prevState.consDaysWorked_;
    consDaysOff_ = prevState.consDaysOff_;
    pShift_ = pS;

    // increment the day index
    dayId_ = prevState.dayId_ + 1;
    return;
  }
  if (prevShiftType >= 0) {
    // Total weekends worked:
    // +1 IF : new day is a Sunday and the nurse works
    // on prevState.pAShift_ or newShift
    if (Day::isFirstDayOfWeek(dayId_)) {
      if (pS->isWork() || prevState.pShift_->isWork()) {
        totalWeekendsWorked_ = prevState.totalWeekendsWorked_ + 1;
        consWeekendWorked_ = prevState.consWeekendWorked_ + 1;
        consWeekendOff_ = 0;
      } else {
        totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
        consWeekendOff_ = prevState.consWeekendOff_ + 1;
        consWeekendWorked_ = 0;
      }
    } else {
      totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
      consWeekendWorked_ = prevState.consWeekendWorked_;
      consWeekendOff_ = prevState.consWeekendOff_;
    }

    // Consecutives : +1 iff it is the same as the previous one
    consShifts_ = (pS->type == prevShiftType) ? prevState.consShifts_ + 1 : 1;

    // Consecutive Days Worked :
    // +1 if the new one is worked (!=0), 0 if it is a rest (==0)
    consDaysWorked_ = pS->isWork() ? (prevState.consDaysWorked_ + 1) : 0;

    // Consecutive Days off :
    // +1 if the new one is off (==0), 0 if it is worked (!=0)
    consDaysOff_ = pS->isWork() ? 0 : (prevState.consDaysOff_ + 1);
  } else {  // the previous shift was not assigned but this one is
    if (pS->isWork()) {
      totalTimeWorked_ = prevState.totalTimeWorked_ + pS->duration
          + (prevState.consDaysWorked_ > 0);
      totalDaysWorked_ = prevState.totalDaysWorked_ + 1
          + (prevState.consDaysWorked_ > 0);
      if (Day::isFirstDayOfWeek(dayId_)) {
        totalWeekendsWorked_ = prevState.totalWeekendsWorked_ + 1;
        consWeekendWorked_ = prevState.consWeekendWorked_ + 1;
        consWeekendOff_ = 0;
      } else {
        totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
        consWeekendWorked_ = prevState.consWeekendWorked_;
        consWeekendOff_ = prevState.consWeekendOff_;
      }
      consDaysWorked_ =
          (prevState.consDaysWorked_ > 0) ? (prevState.consDaysWorked_ + 1
              - prevShiftType) : 1;
      consShifts_ = 1;
      consDaysOff_ = 0;
    } else {
      if (Day::isFirstDayOfWeek(dayId_)) {
        consWeekendOff_ = prevState.consWeekendOff_ + 1;
        consWeekendWorked_ = 0;
      } else {
        consWeekendWorked_ = prevState.consWeekendWorked_;
        consWeekendOff_ = prevState.consWeekendOff_;
      }
      totalTimeWorked_ = prevState.totalTimeWorked_;
      totalDaysWorked_ = prevState.totalDaysWorked_;
      totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
      consDaysWorked_ = 0;
      consShifts_ = 0;
      consDaysOff_ = (prevState.consDaysOff_ > 0) ? (prevState.consDaysOff_ + 1
          - prevShiftType) : 1;
    }
  }

  // update shift
  pShift_ = pS;
  // increment the day index
  dayId_ = prevState.dayId_ + 1;  // increment the day index
}

// set each total to 0. Just keep the consecutive counters
void State::resetTotal() {
  totalTimeWorked_ = 0;
  totalDaysWorked_ = 0;
  totalWeekendsWorked_ = 0;
}

// Display method: toString
string State::toString(const ShiftsFactory& shiftsFactory) const {
  if (pShift_->type >= 0)
    return toString(shiftsFactory.pAnyTypeShift(pShift_->type)->name);
  return toString(pShift_->name);
}

string State::toString(const string &shiftType) const {
  std::stringstream rep;
  rep << totalTimeWorked_ << " " << totalDaysWorked_ << " "
      << totalWeekendsWorked_ << " "
      << (shiftType.empty() ? std::to_string(pShift_->type) : shiftType) << " "
      << consShifts_ << " " << consDaysWorked_ << " "
      << consDaysOff_ << std::endl;
  return rep.str();
}



//-----------------------------------------------------------------------------
//
//  C l a s s   S c e n a r i o
//
//  Class that contains all the attributes describing the scenario
//
//-----------------------------------------------------------------------------

vector2D<int> forbiddenShiftTypeSuccessors(
    const vector<PShift> &pShifts,
    const vector2D<int> &shiftTypeIDToShiftID,
    const vector<int> &shiftIDToShiftTypeID) {
  vector2D<int> forbidSuccs(shiftTypeIDToShiftID.size());
  int nShiftType = shiftTypeIDToShiftID.size();
  for (int st1 = 0; st1 < nShiftType; st1++)
    for (int st2 = 0; st2 < nShiftType; st2++) {
      if (st1 == st2) continue;
      int nShiftToFound2 = shiftTypeIDToShiftID[st2].size();
      bool allFound = true;
      for (int s1 : shiftTypeIDToShiftID[st1]) {
        const PShift pS1 = pShifts.at(s1);
        int n = 0;
        for (const PShift &pS2 : pShifts)
          if (!pS2->canSucceed(*pS1) && pS2->type == st2)
            n++;
        if (n != nShiftToFound2) {
          allFound = false;
          break;
        }
      }
      if (allFound) forbidSuccs[st1].push_back(st2);
    }
  return forbidSuccs;
}

// Constructor and destructor
//
Scenario::Scenario(std::string name,
                   int nbWeeks,
                   vector<std::string> intToSkill,
                   std::map<std::string, int> skillToInt,
                   vector<PShift> pShifts,
                   vector<std::string> intToShiftType,
                   vector<int> minConsShiftType,
                   vector<int> maxConsShiftType,
                   vector<PContract> contracts,
                   std::vector<PNurse> theNurses,
                   PWeights weights,
                   bool isINRC2) :
    name_(std::move(name)),
    nWeeks_(nbWeeks),
    intToSkill_(std::move(intToSkill)),
    skillToInt_(std::move(skillToInt)),
    pShifts_(std::move(pShifts)),
    shiftsFactory_(pShifts_, intToShiftType),
    theContracts(std::move(contracts)),
    pWeights_(std::move(weights)),
    pNurses_(std::move(theNurses)),
    minConsShiftType_(std::move(minConsShiftType)),
    maxConsShiftType_(std::move(maxConsShiftType)),
    nShiftOffRequests_(0),
    nShiftOnRequests_(0) {
  // To make sure that it is modified later when reading the history data file
  thisWeek_ = -1;
  init();
}

Scenario::Scenario(std::string name,
                   int nbWeeks,
                   vector<std::string> intToSkill,
                   std::map<std::string, int> skillToInt,
                   vector<PShift> pShifts,
                   vector<std::string> intToShiftType,
                   vector<PContract> contracts,
                   std::vector<PNurse> theNurses,
                   string header) :
        name_(std::move(name)),
        nWeeks_(nbWeeks),
        intToSkill_(std::move(intToSkill)),
        skillToInt_(std::move(skillToInt)),
        pShifts_(std::move(pShifts)),
        shiftsFactory_(pShifts_, intToShiftType),
        theContracts(std::move(contracts)),
        pNurses_(std::move(theNurses)),
        nShiftOffRequests_(0),
        nShiftOnRequests_(0),
        thisWeek_(0) {
  init();
}

// Hybrid copy constructor : this is only called when constructing
// a new scenario that copies most parameters
// from the input scenario but for only a subgroup of nurses
//
Scenario::Scenario(const PScenario &pScenario,
                   const vector<PNurse> &theNurses,
                   vector<PDemand> pDemands,
                   PPreferences pPreferences) :
    Scenario(*pScenario) {
  // change nurses and reset positions
  pNurses_ = theNurses;
  pPositions_.clear();
  nursesPerPosition_.clear();
  // Preprocess the vector of nurses
  // This creates the positions
  this->preprocessTheNurses();

  // The nurses are already preprocessed at this stage
  // Load the input week demand and preferences
  this->linkWithDemand(std::move(pDemands));
  this->linkWithPreferences(std::move(pPreferences));
}

Scenario::~Scenario() = default;

void Scenario::init() {
  // build shiftToInt_, shiftTypeToInt_
  for (int i = 0; i < nShifts(); i++)
    shiftToInt_[pShift(i)->name] = i;
  for (int i = 0; i < shiftsFactory_.pAnyTypeShifts().size(); i++)
    shiftTypeToInt_[pAnyTypeShift(i)->name] = i;

  // check the contracts index are correctly set
  for (int n=0; n < nContracts(); n++)
    if (theContracts.at(n)->id_ != n)
      Tools::throwError("Contract %s has the index %d different than "
                        "its position %d",
                        theContracts.at(n)->name_.c_str(),
                        theContracts.at(n)->id_, n);

  // match nurse name to their index
  for (int n = 0; n < nNurses(); n++)
    nurseNameToInt_[pNurses_[n]->name_] = n;

  // Preprocess the vector of nurses
  // This creates the positions
  this->preprocessTheNurses();
}

void Scenario::setPreferences(PPreferences pPreferences) {
  pPreferences_ = std::move(pPreferences);
}

// return true if the shift shNext is a forbidden successor of sh
//
//  bool Scenario::isForbiddenSuccessor(int shNext, int shLast) {
//    if (shLast <= 0) return false;
//
//    for (int i = 0; i < nbForbiddenSuccessors_[shLast]; i++) {
//      if (shNext == forbiddenSuccessors_[shLast][i]) {
//        return true;
//      }
//    }
//    return false;
//  }

// return true if the shift shNext is
// a forbidden successor of shift shLast (via types)
bool Scenario::isForbiddenSuccessorShift_Shift(int shNext, int shLast) const {
  if (shLast <= 0) return false;

  const PShift pSLast = pShift(shLast), pSNext = pShift(shNext);
  pSNext->canSucceed(*pSLast);
  return !pSNext->canSucceed(*pSLast);
}

// return true if the shiftType shTypeNext is
// a forbidden successor of shiftType shTypeLast
bool Scenario::isForbiddenSuccessorShiftType_ShiftType(int shTypeNext,
                                                       int shTypeLast) const {
  if (shTypeLast <= 0) return false;

  const vector<int> &lastSucc = pShiftsOfType(shTypeLast).front()->successors;
  int sNext = pShiftsOfType(shTypeNext).front()->id;
  return std::find(lastSucc.begin(), lastSucc.end(), sNext) == lastSucc.end();
}

// return the min/max consecutive shifts of the same type as the argument

int Scenario::minConsShifts(int whichShift) const {
  if (minConsShiftType_.empty())
    Tools::throwError("minConsShiftType is not defined in Scenario.");
  return minConsShiftsOfType(pShift(whichShift)->type);
}

int Scenario::maxConsShifts(int whichShift) const {
  if (maxConsShiftType_.empty())
    Tools::throwError("maxConsShiftType is not defined in Scenario.");
  return maxConsShiftsOfType(pShift(whichShift)->type);
}

int Scenario::minConsShiftsOfType(int whichShiftType) const {
  if (whichShiftType == 0)
    Tools::throwError(
        "Behavior not defined for a rest shift. Please use minConsDaysOffOf.");
  if (minConsShiftType_.empty())
    Tools::throwError("minConsShiftType is not defined in Scenario.");
  return minConsShiftType_[whichShiftType];
}

int Scenario::maxConsShiftsOfType(int whichShiftType) const {
  if (whichShiftType == 0)
    Tools::throwError(
        "Behavior not defined for a rest shift. Please use maxConsDaysOffOf.");
  if (maxConsShiftType_.empty())
    Tools::throwError("maxConsShiftType is not defined in Scenario.");
  return maxConsShiftType_[whichShiftType];
}

// update the scenario to treat a new week
//
void Scenario::updateNewWeek(vector<PDemand> pDemands,
                             PPreferences pPreferences,
                             const vector<State> &initialStates) {
  // set the demand, preferences and initial states
  this->linkWithDemand(std::move(pDemands));
  this->linkWithPreferences(std::move(pPreferences));
  this->setInitialState(initialStates);

  // update the index of the week
  thisWeek_++;
}

void Scenario::linkWithPreferences(PPreferences pPreferences) {
  pPreferences_ = std::move(pPreferences);
  nShiftOffRequests_ = 0;
  for (const PNurse &nurse : pNurses_)
    nShiftOffRequests_ += pPreferences_->howManyShiftsOff(nurse->num_);
  nShiftOnRequests_ = 0;
  for (const PNurse &nurse : pNurses_)
    nShiftOnRequests_ += pPreferences_->howManyShiftsOn(nurse->num_);
}

void Scenario::pushBack(
        const vector<PDemand> &pDemands, PPreferences pPreferences) {
  for (int d = 0; d < pDemands_.size(); d++)
    pDemands_[d] = pDemands_[d]->append(pDemands[d]);
  pPreferences_ = pPreferences_->append(pPreferences);
  for (const PNurse &nurse : pNurses_)
    nShiftOffRequests_ += pPreferences->howManyShiftsOff(nurse->num_);
  for (const PNurse &nurse : pNurses_)
    nShiftOnRequests_ += pPreferences->howManyShiftsOn(nurse->num_);
}

//------------------------------------------------
// Display functions
//------------------------------------------------
string Scenario::toString() const {
  std::stringstream rep;
  rep << "######################################"
         "######################################"
      << std::endl;
  rep << "##############################"
         "    Scenario    "
         "##############################"
      << std::endl;
  rep << "######################################"
         "######################################"
      << std::endl;
  rep << "# " << std::endl;
  rep << "# NAME             \t= " << name_ << std::endl;
  rep << "# NUMBER_OF_WEEKS  \t= " << nWeeks_ << std::endl;
  rep << "# " << std::endl;
  rep << "# SKILLS           \t= " << nSkills() << std::endl;
  for (int i = 0; i < nSkills(); i++) {
    rep << "#                  \t= " << i << ":" << intToSkill_[i] << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# SHIFTS           \t= " << nShifts() - 1 << std::endl;
  for (const auto &pS : this->pShifts()) {
    if (pS->isRest()) continue;
    rep << "#                  \t= ";
    rep << pS->id << ":" << pS->name;
    if (!minConsShiftType_.empty())
      rep << " \t(" << minConsShifts(pS->id) << ","
          << maxConsShifts(pS->id) << ")";
    rep << std::endl;
  }

  rep << "# " << std::endl;
  rep << "# CONTRACTS        " << std::endl;
  for (const auto &contract : theContracts) {
    rep << "#\t\t\t" << *(contract) << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# NURSES           \t= " << nNurses() << std::endl;
  for (const auto &pN : pNurses_)
    rep << "#\t\t\t" << pN->toString() << std::endl;
  rep << "# " << std::endl;
  rep << "# POSITIONS        \t= " << nPositions() << std::endl;
  for (const auto &pP : pPositions_)
    rep << "#\t\t\t" << pP->toString() << std::endl;
  if (!pDemands_.empty()) {
    // write the demand using the member method toString
    // do not write the preprocessed information at this stage
    for (const auto &pD : pDemands_) {
      rep << pD->toString(false) << std::endl;
      rep << "# " << std::endl;
    }

    // write the preferences
    rep << "# " << std::endl;
    rep << pPreferences_->toString();
  }
  if (thisWeek_ > -1) {
    rep << "# " << std::endl;
    rep << "# INITIAL STATE    \t= WEEK Nb " << thisWeek_ << std::endl;
    for (int n = 0; n < nNurses(); n++) {
      rep << "#\t\t\t" << pNurses_[n]->name_ << " ";
      rep << initialState_[n].toString(shiftsFactory_);
    }
  }
  rep << "########################################"
         "####################################"
      << std::endl;
  return rep.str();
}


//------------------------------------------------
// Preprocess functions
//------------------------------------------------

// presolve the nurses to get the types
//
void Scenario::preprocessTheNurses() {
  if (!pPositions_.empty())
    Tools::throwError("The nurse preprocessing is run for the second time!");

  // Go through the nurses, and create their positions when it has not already
  // been done
  for (const PNurse &nurse : pNurses_) {
    // go through every existing position to see if the position of this nurse
    // has already been created
    bool positionExists = false;
    for (const PPosition &pos : pPositions_) {
      positionExists = pos->isNursePosition(*nurse);
      if (positionExists) break;
    }

    // create the position if it doesn't exist
    if (!positionExists)
      pPositions_.push_back(std::make_shared<Position>(nPositions(), *nurse));
  }

  // build the list of position dominance
  for (int i = 0; i < nPositions(); i++) {
    for (int j = i + 1; j < nPositions(); j++) {
      if (pPositions_[i]->compare(*pPositions_[j]) == 1) {
        pPositions_[i]->addBelow(pPositions_[j]);
        pPositions_[j]->addAbove(pPositions_[i]);
      }
      if (pPositions_[i]->compare(*pPositions_[j]) == -1) {
        pPositions_[i]->addAbove(pPositions_[j]);
        pPositions_[j]->addBelow(pPositions_[i]);
      }
    }
  }

  // compute the rank of each position
  vector<bool> isRanked;
  bool isAllRanked = false;
  for (const auto &pP : pPositions_) {
    isRanked.push_back(false);
    pP->rank(0);
    pP->resetAbove();
    pP->resetBelow();
  }
  while (!isAllRanked) {
    for (int i = 0; i < nPositions(); i++) {
      if (!isRanked[i]) {
        int rankIni = pPositions_[i]->rank();

        // go through the positions above the considered position to increment
        // its rank
        if (pPositions_[i]->nAbove()) {
          for (int j = 0; j < pPositions_[i]->nAbove(); j++) {
            int currentRank = pPositions_[i]->rank();
            int newRank = pPositions_[j]->rank() + 1;
            pPositions_[i]->rank(std::max(currentRank, newRank));
          }
        }

        // the position is treated when the rank is not modified
        // by the loop above
        if (pPositions_[i]->rank() == rankIni) isRanked[i] = true;
      }
    }

    // check if all the positions are ranked
    isAllRanked = true;
    for (int i = 0; i < nPositions(); i++) {
      if (!isRanked[i]) {
        isAllRanked = false;
        break;
      }
    }
  }

  // Organize the nurses by position
  for (int i = 0; i < nPositions(); i++) {
    vector<PNurse> nursesInThisPosition;
    nursesPerPosition_.push_back(nursesInThisPosition);
  }

  for (const PNurse &pNurse : pNurses_) {
    // the skills of the nurse need to be compared to the skills of each
    // existing position to determine the position of the nurse
    bool isPosition = false;
    for (int i = 0; i < nPositions(); i++) {
      PPosition pPosition = pPositions_[i];
      isPosition = (*pNurse == *pPosition);
      if (*pNurse == *pPosition) {
        nursesPerPosition_[pPosition->id()].push_back(pNurse);
        break;
      }
    }
    if (!isPosition)
      Tools::throwError("The nurse has no position!");
  }
}

// compute the connected components of the positions rcspp
// (one edge between two positions indicate that they share a skill)
//
void Scenario::computeConnectedPositions() {
  vector<PPosition> pRemainingPositions(pPositions_);

  // First build the connected components of positions
  while (!pRemainingPositions.empty()) {
    PPosition pPos = pRemainingPositions.back();
    pRemainingPositions.pop_back();
    vector<PPosition> connectedPositions;
    connectedPositions.push_back(pPos);

    vector<PPosition> nextPositions;
    nextPositions.push_back(pPos);

    while (!nextPositions.empty()) {
      // treat the next position in the waiting list
      pPos = nextPositions.back();
      nextPositions.pop_back();

      // add all the positions with a common skill in the connected component
      int nbRemaining = pRemainingPositions.size();
      int ind = 0;
      while (ind < nbRemaining) {
        if (pPos->shareSkill(*pRemainingPositions[ind])) {
          // update the connected component, the waiting list and the list of
          // positions which are not in any connected component
          connectedPositions.push_back(pRemainingPositions[ind]);
          nextPositions.push_back(pRemainingPositions[ind]);
          pRemainingPositions.erase(pRemainingPositions.begin() + ind);
          nbRemaining--;
        } else {
          ind++;
        }
      }
    }
    componentsOfConnectedPositions_.push_back(connectedPositions);
  }

  // Get the nurses that belong to each component
  for (const auto &component : componentsOfConnectedPositions_) {
    vector<PNurse> pNursesInThisComponent;

    for (const PPosition &p : component)
      for (const auto &pN : nursesPerPosition_[p->id()])
        pNursesInThisComponent.push_back(pN);

    std::stable_sort(pNursesInThisComponent.begin(),
                     pNursesInThisComponent.end(),
                     compareNursesById);
    nursesPerConnectedComponentOfPositions_.push_back(pNursesInThisComponent);
  }
}

