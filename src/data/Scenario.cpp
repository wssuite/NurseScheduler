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
#include <string>

#include "solvers/Solver.h"
#include "data/Nurse.h"
#include "tools/MyTools.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//-----------------------------------------------------------------------------

// Destructor
State::~State() {}

//// Updates the state if a new day is worked on shiftType of newShiftType
// void State::addNewDay(int newShiftType) {
//
//  // Total shifts worked if it is a worked day
//  totalTimeWorked_ += (newShiftType ? 1 : 0);
//
//  // Total weekends worked :
//  // +1 IF : new day is a Sunday and the nurse works on pShift_ or newShift
//  if (Tools::isSunday(dayId_ - 1) and (newShiftType or shiftType_))
//    totalWeekendsWorked_++;
//
//  // Consecutives : +1 iff it is the same as the previous one
//  consShifts_ = (shiftType_ == newShiftType) ? (consShifts_ + 1) : 1;
//
//  // Consecutive Days Worked :
//  // +1 if the new one is worked (!=0), 0 if it is a rest (==0)
//  consDaysWorked_ = shiftType_ ? (consDaysWorked_ + 1) : 0;
//
//  // Current shiftType worked : updated with the new one
//  shiftType_ = newShiftType;
//
//  // increment the day index
//  dayId_++;
//}


// Function that appends a new day worked on a given shift to an input state
//

// void State::addDayToState(const State& prevState, int newShift,
// const PScenario pScenario)   {
//   int newShiftType = pScenario->shiftIDToShiftTypeID_[newShift];

//   addDayToState(prevState, newShiftType);
// }

// Function that appends a new day worked on a given shiftType to an input state
// to update this state
// RqJO: I slghtly modified the method to take into account the possibility to
// add in the state that no task has been assigned on this day
//
void State::addDayToState(const State &prevState,
                          int newShiftType,
                          int newShift,
                          int timeWorked) {
  //  int  timeWorked = 1;               // days worked
  //  int  timeWorked = timeDurationToWork_[newShift];      // hours worked

  // Total shifts worked if it is a worked day
  totalTimeWorked_ =
      prevState.totalTimeWorked_ + (newShiftType > 0 ? timeWorked : 0);

  // index of the previous shift
  int prevShiftType = prevState.shiftType_;

  // Treat the case in which no shift is assigned to the nurse on this day
  if (newShiftType < 0) {
    totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
    consShifts_ = prevState.consShifts_;
    consDaysWorked_ = prevState.consDaysWorked_;
    consDaysOff_ = prevState.consDaysOff_;

    shiftType_ = prevShiftType < 0 ? prevShiftType - 1 : -1;
    shift_ = 0;
  } else if (prevShiftType >= 0) {
    // Total weekends worked:
    // +1 IF : new day is a Sunday and the nurse works
    // on prevState.pShift_ or newShift
    if (Tools::isSunday(dayId_ - 1) && (newShiftType || prevState.shiftType_))
      totalWeekendsWorked_ = prevState.totalWeekendsWorked_ + 1;
    else
      totalWeekendsWorked_ = prevState.totalWeekendsWorked_;

    // Consecutives : +1 iff it is the same as the previous one
    consShifts_ = (newShiftType && newShiftType == prevState.shiftType_) ?
                  prevState.consShifts_ + 1 : (newShiftType ? 1 : 0);

    // Consecutive Days Worked :
    // +1 if the new one is worked (!=0), 0 if it is a rest (==0)
    consDaysWorked_ = newShiftType ? (prevState.consDaysWorked_ + 1) : 0;

    // Consecutive Days off :
    // +1 if the new one is off (==0), 0 if it is worked (!=0)
    consDaysOff_ = newShiftType ? 0 : (prevState.consDaysOff_ + 1);

    shiftType_ = newShiftType;
    shift_ = newShift;
  } else {  // the previous shift was not assigned but this one is
    if (newShiftType > 0) {
      totalTimeWorked_ = prevState.totalTimeWorked_ + timeWorked
          + (prevState.consDaysWorked_ > 0 ? (-prevShiftType)
                                           : 0);  // SERGEB ??????
      totalWeekendsWorked_ =
          Tools::isSunday(dayId_) ? prevState.totalWeekendsWorked_ + 1
                                  : prevState.totalWeekendsWorked_;
      consDaysWorked_ =
          (prevState.consDaysWorked_ > 0) ? (prevState.consDaysWorked_ + 1
              - prevShiftType) : 1;
      consShifts_ = 1;
      consDaysOff_ = 0;
      shiftType_ = newShiftType;
      shift_ = newShift;
    } else {
      totalTimeWorked_ = prevState.totalTimeWorked_;
      totalWeekendsWorked_ = prevState.totalWeekendsWorked_;
      consDaysWorked_ = 0;
      consShifts_ = 0;
      consDaysOff_ = (prevState.consDaysOff_ > 0) ? (prevState.consDaysOff_ + 1
          - prevShiftType) : 1;
      shiftType_ = newShiftType;
      shift_ = newShift;
    }
  }

  // increment the day index
  dayId_ = prevState.dayId_ + 1;
}

// Display method: toString
//
string State::toString() {
  std::stringstream rep;
  rep << totalTimeWorked_ << " " << totalWeekendsWorked_ << " " << shiftType_
      << " ";
  if (shiftType_) rep << consShifts_ << " " << consDaysWorked_;
  else
    rep << "0 0";
  if (shiftType_) rep << " 0";
  else
    rep << " " << consShifts_;
  rep << std::endl;
  return rep.str();
}



//-----------------------------------------------------------------------------
//
//  C l a s s   S c e n a r i o
//
//  Class that contains all the attributes describing the scenario
//
//-----------------------------------------------------------------------------

// Constructor and destructor
//
Scenario::Scenario(string name,
                   int nbWeeks,
                   int nbSkills,
                   vector<string> intToSkill,
                   map<string, int> skillToInt,
                   int nbShifts,
                   vector<string> intToShift,
                   map<string, int> shiftToInt,
                   vector<int> hoursToWork,
                   vector<int> shiftIDToShiftTypeID,
                   int nbShiftsType,
                   vector<string> intToShiftType,
                   map<string, int> shiftTypeToInt,
                   vector2D<int> shiftTypeIDToShiftID,
                   vector<int> minConsShiftType,
                   vector<int> maxConsShiftType,
                   vector<int> nbForbiddenSuccessors,
                   vector2D<int> forbiddenSuccessors,
                   vector<PShift> pShifts,
                   int nbContracts,
                   vector<string> intToContract,
                   map<string, PConstContract> contracts,
                   int nbNurses,
                   vector<PNurse> theNurses,
                   map<string, int> nurseNameToInt,
                   PWeights weights) :
    name_(std::move(name)),
    nbWeeks_(nbWeeks),
    nbSkills_(nbSkills),
    intToSkill_(std::move(intToSkill)),
    skillToInt_(std::move(skillToInt)),
    nbShifts_(nbShifts),
    intToShift_(std::move(intToShift)),
    shiftToInt_(std::move(shiftToInt)),
    timeDurationToWork_(std::move(hoursToWork)),
    shiftIDToShiftTypeID_(std::move(shiftIDToShiftTypeID)),
    nbShiftsType_(nbShiftsType),
    intToShiftType_(std::move(intToShiftType)),
    shiftTypeToInt_(std::move(shiftTypeToInt)),
    shiftTypeIDToShiftID_(std::move(shiftTypeIDToShiftID)),
    pShifts_(std::move(pShifts)),
    nbContracts_(nbContracts),
    intToContract_(std::move(intToContract)),
    pContracts_(std::move(contracts)),
    nbNurses_(nbNurses),
    theNurses_(std::move(theNurses)),
    nurseNameToInt_(std::move(nurseNameToInt)),
    pWeights_(std::move(weights)),
    minConsShiftType_(std::move(minConsShiftType)),
    maxConsShiftType_(std::move(maxConsShiftType)),
    nbForbiddenSuccessors_(std::move(nbForbiddenSuccessors)),
    forbiddenSuccessors_(std::move(forbiddenSuccessors)),
    pWeekDemand_(0),
    nbShiftOffRequests_(0),
    nbShiftOnRequests_(0),
    nbPositions_(0) {
  // To make sure that it is modified later when reading the history data file
  //
  thisWeek_ = -1;
  nbWeeksLoaded_ = 1;

  // Preprocess the vector of nurses
  // This creates the positions
  //
  this->preprocessTheNurses();
}

// Hybrid copy constructor : this is only called when constructing
// a new scenario that copies most parameters
// from the input scenario but for only a subgroup of nurses
//
Scenario::Scenario(PScenario pScenario,
                   const vector<PNurse> &theNurses,
                   PDemand pDemand,
                   PPreferences pWeekPreferences) :
    name_(pScenario->name_),
    nbWeeks_(pScenario->nbWeeks_),
    nbSkills_(pScenario->nbSkills_),
    intToSkill_(pScenario->intToSkill_),
    skillToInt_(pScenario->skillToInt_),
    nbShifts_(pScenario->nbShifts_),
    intToShift_(pScenario->intToShift_),
    shiftToInt_(pScenario->shiftToInt_),
    timeDurationToWork_(pScenario->timeDurationToWork_),
    shiftIDToShiftTypeID_(pScenario->shiftIDToShiftTypeID_),
    nbShiftsType_(pScenario->nbShiftsType_),
    intToShiftType_(pScenario->intToShiftType_),
    shiftTypeToInt_(pScenario->shiftTypeToInt_),
    shiftTypeIDToShiftID_(pScenario->shiftTypeIDToShiftID_),
    pShifts_(pScenario->pShifts_),
    nbContracts_(pScenario->nbContracts_),
    intToContract_(pScenario->intToContract_),
    pContracts_(pScenario->pContracts_),
    nbNurses_(theNurses.size()),
    theNurses_(theNurses),
    nurseNameToInt_(pScenario->nurseNameToInt_),
    pWeights_(pScenario->pWeights_),
    minConsShiftType_(pScenario->minConsShiftType_),
    maxConsShiftType_(pScenario->maxConsShiftType_),
    nbForbiddenSuccessors_(pScenario->nbForbiddenSuccessors_),
    forbiddenSuccessors_(pScenario->forbiddenSuccessors_),
    pWeekDemand_(nullptr),
    nbShiftOffRequests_(0),
    nbShiftOnRequests_(0),
    pWeekPreferences_(nullptr),
    thisWeek_(pScenario->thisWeek()),
    nbWeeksLoaded_(pScenario->nbWeeksLoaded()),
    nbPositions_(0),
    nursesPerPosition_(0) {
  // Preprocess the vector of nurses
  // This creates the positions
  //
  this->preprocessTheNurses();

  // The nurses are already preprocessed at this stage
  // Load the input week demand and preferences
  //
  this->linkWithDemand(pDemand);
  this->linkWithPreferences(pWeekPreferences);
}

Scenario::~Scenario() = default;

void Scenario::setWeekPreferences(PPreferences weekPreferences) {
  pWeekPreferences_ = weekPreferences;
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
bool Scenario::isForbiddenSuccessorShift_Shift(int shNext, int shLast) {
  if (shLast <= 0) return false;

  int shTypeNext = shiftIDToShiftTypeID_[shNext];
  int shTypeLast = shiftIDToShiftTypeID_[shLast];

  return isForbiddenSuccessorShiftType_ShiftType(shTypeNext, shTypeLast);
}

// return true if the shift shNext is
// a forbidden successor of shiftType shTypeLast
bool Scenario::isForbiddenSuccessorShift_ShiftType(int shNext, int shTypeLast) {
  if (shTypeLast <= 0) return false;

  int shTypeNext = shiftIDToShiftTypeID_[shNext];

  return isForbiddenSuccessorShiftType_ShiftType(shTypeNext, shTypeLast);
}

// return true if the shiftType shTypeNext is
// a forbidden successor of shiftType shTypeLast
bool Scenario::isForbiddenSuccessorShiftType_Shift(int shTypeNext, int shLast) {
  if (shLast <= 0) return false;

  int shTypeLast = shiftIDToShiftTypeID_[shLast];

  return isForbiddenSuccessorShiftType_ShiftType(shTypeNext, shTypeLast);
}

// return true if the shiftType shTypeNext is
// a forbidden successor of shiftType shTypeLast
bool Scenario::isForbiddenSuccessorShiftType_ShiftType(int shTypeNext,
                                                       int shTypeLast) {
  if (shTypeLast <= 0) return false;

  for (int i = 0; i < nbForbiddenSuccessors_[shTypeLast]; i++) {
    if (shTypeNext == forbiddenSuccessors_[shTypeLast][i]) {
      return true;
    }
  }
  return false;
}

const int Scenario::minTotalShiftsOf(int whichNurse) const {
  return theNurses_[whichNurse]->minTotalShifts();
}
int Scenario::maxTotalShiftsOf(int whichNurse) const {
  return theNurses_[whichNurse]->maxTotalShifts();
}
int Scenario::minConsDaysWorkOf(int whichNurse) const {
  return theNurses_[whichNurse]->minConsDaysWork();
}
int Scenario::maxConsDaysWorkOf(int whichNurse) const {
  return theNurses_[whichNurse]->maxConsDaysWork();
}
int Scenario::minConsDaysOffOf(int whichNurse) const {
  return theNurses_[whichNurse]->maxConsDaysOff();
}
int Scenario::maxConsDaysOffOf(int whichNurse) const {
  return theNurses_[whichNurse]->maxConsDaysOff();
}
int Scenario::maxTotalWeekendsOf(int whichNurse) const {
  return theNurses_[whichNurse]->maxTotalWeekends();
}
bool Scenario::isCompleteWeekendsOf(int whichNurse) const {
  return theNurses_[whichNurse]->needCompleteWeekends();
}

// return the min/max consecutive shifts of the same type as the argument

int Scenario::minConsShiftsOfTypeOf(int whichShift) {
  int shiftType = shiftIDToShiftTypeID_[whichShift];
  return minConsShiftsOf(shiftType);
}

int Scenario::maxConsShiftsOfTypeOf(int whichShift) {
  int shiftType = shiftIDToShiftTypeID_[whichShift];
  return maxConsShiftsOf(shiftType);
}

int Scenario::minConsShiftsOf(int whichShiftType) {
  if (whichShiftType == 0)
    Tools::throwError(
        "Behavior not defined for a rest shift. Please use minConsDaysOffOf.");
  return minConsShiftType_[whichShiftType];
}

int Scenario::maxConsShiftsOf(int whichShiftType) {
  if (whichShiftType == 0)
    Tools::throwError(
        "Behavior not defined for a rest shift. Please use maxConsDaysOffOf.");
  return maxConsShiftType_[whichShiftType];
}

// update the scenario to treat a new week
//
void Scenario::updateNewWeek(PDemand pDemand,
                             PPreferences pPreferences,
                             const vector<State> &initialStates) {
  // set the demand, preferences and initial states
  this->linkWithDemand(pDemand);
  this->linkWithPreferences(pPreferences);
  this->setInitialState(initialStates);

  // update the index of the week
  thisWeek_++;
}

void Scenario::linkWithPreferences(PPreferences pPreferences) {
  pWeekPreferences_ = pPreferences;
  nbShiftOffRequests_ = 0;
  for (PNurse nurse : theNurses_)
    nbShiftOffRequests_ += pWeekPreferences_->howManyShiftsOff(nurse->num_);
  nbShiftOnRequests_ = 0;
  for (PNurse nurse : theNurses_)
    nbShiftOnRequests_ += pWeekPreferences_->howManyShiftsOn(nurse->num_);
}

//------------------------------------------------
// Display functions
//------------------------------------------------

// Display methods: toString + override operator<< (easier)
//
string Scenario::toString() {
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
  rep << "# NUMBER_OF_WEEKS  \t= " << nbWeeks_ << std::endl;
  rep << "# " << std::endl;
  rep << "# SKILLS           \t= " << nbSkills_ << std::endl;
  for (int i = 0; i < nbSkills_; i++) {
    rep << "#                  \t= " << i << ":" << intToSkill_[i] << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# SHIFTS           \t= " << nbShifts_ - 1 << std::endl;
  for (int i = 0; i < nbShifts_; i++) {
    if (isRestShift(i)) continue;
    rep << "#                  \t= ";
    rep << i << ":" << intToShift_[i];
    if (i > 0)
      rep << " \t(" << minConsShiftsOfTypeOf(i) << ","
          << maxConsShiftsOfTypeOf(i) << ")";
    rep << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# FORBIDDEN        " << std::endl;
  for (int i = 0; i < nbShifts_; i++) {
    rep << "#\t\t\t" << intToShift_[i] << "\t-> ";
    for (int j = 0; j < nbForbiddenSuccessors_[i]; j++) {
      rep << intToShift_[forbiddenSuccessors_[i][j]] << " ";
    }
    rep << std::endl;
  }
  rep << "# CONTRACTS        " << std::endl;
  for (const auto &contract : pContracts_) {
    rep << "#\t\t\t" << *(contract.second) << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# NURSES           \t= " << nbNurses_ << std::endl;
  for (int i = 0; i < nbNurses_; i++) {
    rep << "#\t\t\t" << theNurses_[i]->toString() << std::endl;
  }
  rep << "# " << std::endl;
  rep << "# POSITIONS        \t= " << nbPositions_ << std::endl;
  for (int i = 0; i < nbPositions_; i++) {
    rep << "#\t\t\t" << pPositions_[i]->toString() << std::endl;
  }
  if (weekName_ != "") {
    // write the demand using the member method toString
    // do not write the preprocessed information at this stage
    //
    rep << pWeekDemand_->toString(false) << std::endl;

    // write the preferences
    //
    rep << "# " << std::endl;
    rep << "# WISHED SHIFTS OFF" << std::endl;
    for (PNurse nurse : theNurses_) {
      // Display only if the nurse has preferences
      const map<int, vector<Wish> >
          &prefNurse = pWeekPreferences_->nurseWishesOff(nurse->num_);
      if (!prefNurse.empty()) {
        rep << "#\t\t\t" << nurse->num_ << "\t" << nurse->name_ << "\t";
        for (const auto &itWishlist : prefNurse) {
          rep << Tools::intToDay(itWishlist.first) << ": ";
          bool first = true;
          for (const auto &itShift : itWishlist.second) {
            if (first) first = false;
            else
              rep << ",";
            rep << intToShift_[itShift.shift];
          }
          rep << "    ";
        }
        rep << std::endl;
      }
    }
    for (PNurse nurse : theNurses_) {
      // Display only if the nurse has preferences
      const map<int, vector<Wish> >
          &prefNurse = pWeekPreferences_->nurseWishesOn(nurse->num_);
      if (!prefNurse.empty()) {
        rep << "#\t\t\t" << nurse->num_ << "\t" << nurse->name_ << "\t";
        for (const auto &itWishlist : prefNurse) {
          rep << Tools::intToDay(itWishlist.first) << ": ";
          bool first = true;
          for (const auto &itShift : itWishlist.second) {
            if (first) first = false;
            else
              rep << ",";
            rep << intToShift_[itShift.shift];
          }
          rep << "    ";
        }
        rep << std::endl;
      }
    }
  }
  if (thisWeek_ > -1) {
    rep << "# " << std::endl;
    rep << "# INITIAL STATE    \t= WEEK Nb " << thisWeek_ << std::endl;
    for (int n = 0; n < nbNurses_; n++) {
      rep << "#\t\t\t" << theNurses_[n]->name_ << " ";
      State s = initialState_[n];
      rep << s.totalTimeWorked_ << " " << s.totalWeekendsWorked_ << " "
          << intToShiftType_[s.shiftType_] << " ";
      if (s.shiftType_) {
        rep << s.consShifts_ << " " << s.consDaysWorked_;
      } else {
        rep << "0 0";
      }
      if (s.shiftType_) rep << " 0";
      else
        rep << " " << s.consShifts_;
      rep << std::endl;
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

// preprocess the nurses to get the types
//
void Scenario::preprocessTheNurses() {
  if (nbPositions_) {
    Tools::throwError("The nurse preprocessing is run for the second time!");
  }

  // Go through the nurses, and create their positions when it has not already
  // been done
  //
  for (PNurse nurse : theNurses_) {
    bool positionExists = nbPositions_ != 0;
    int nbSkills = nurse->nbSkills_;
    vector<int> skills = nurse->skills_;

    // go through every existing position to see if the position of this nurse
    // has already been created
    for (PPosition pos : pPositions_) {
      positionExists = true;
      if (pos->nbSkills_ == nbSkills) {
        for (int i = 0; i < nbSkills; i++) {
          if (skills[i] != pos->skills_[i]) {
            positionExists = false;
            break;
          }
        }
      } else {
        positionExists = false;
      }
      if (positionExists) break;
    }

    // create the position if if doesn't exist
    if (!positionExists) {
      pPositions_.emplace_back(std::make_shared<Position>(nbPositions_,
                                                          nbSkills,
                                                          skills));
      nbPositions_++;
    }
  }

  // build the list of position dominance
  for (int i = 0; i < nbPositions_; i++) {
    for (int j = i + 1; j < nbPositions_; j++) {
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
  for (int i = 0; i < nbPositions_; i++) {
    isRanked.push_back(false);
    pPositions_[i]->rank(0);
    pPositions_[i]->resetAbove();
    pPositions_[i]->resetBelow();
  }
  while (!isAllRanked) {
    for (int i = 0; i < nbPositions_; i++) {
      if (!isRanked[i]) {
        int rankIni = pPositions_[i]->rank();

        // go through the positions above the considered position to increment
        // its rank
        if (pPositions_[i]->nbAbove()) {
          for (int j = 0; j < pPositions_[i]->nbAbove(); j++) {
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
    for (int i = 0; i < nbPositions_; i++) {
      if (!isRanked[i]) {
        isAllRanked = false;
        break;
      }
    }
  }

  // Organize the nurses by position
  for (int i = 0; i < nbPositions(); i++) {
    vector<PNurse> nursesInThisPosition;
    nursesPerPosition_.push_back(nursesInThisPosition);
  }
  for (PNurse nurse : theNurses_) {
    // the skills of the nurse need to be compared to the skills of each
    // existing position to determine the position of the nurse
    bool isPosition = true;
    for (int i = 0; i < nbPositions(); i++) {
      PPosition pPosition = pPositions_[i];
      isPosition = true;
      if (pPosition->nbSkills_ == nurse->nbSkills_) {
        for (int i = 0; i < nurse->nbSkills_; i++) {
          if (nurse->skills_[i] != pPosition->skills_[i]) {
            isPosition = false;
            break;
          }
        }
      } else {
        isPosition = false;
      }

      if (isPosition) {
        nursesPerPosition_[pPosition->id()].push_back(nurse);
        break;
      }
    }
    if (!isPosition) {
      Tools::throwError("The nurse has no position!");
    }
  }
}

// Compare nurses in ascending order of their ideas
//
bool compareNursesById(PNurse n1, PNurse n2) {
  return (n1->num_ < n2->num_);
}

bool comparePairSecond(pair<int, int> p1, pair<int, int> p2) {
  return p1.second < p2.second;
}

// compute the connected components of the positions rcspp
// (one edge between two positions indicate that they share a skill)
//
void Scenario::computeConnectedPositions() {
  vector<PPosition> pRemainingPositions(pPositions_);

  PPosition pPos = pRemainingPositions.back();
  pRemainingPositions.pop_back();


  // First build the connected components of positions
  while (!pRemainingPositions.empty()) {
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

    // load the next remaining position
    if (!pRemainingPositions.empty()) {
      pPos = pRemainingPositions.back();
      pRemainingPositions.pop_back();
    }
  }

  // Get the nurses that belong to each component
  for (unsigned int c = 0; c < componentsOfConnectedPositions_.size(); c++) {
    vector<PNurse> pNursesInThisComponent;

    for (PPosition p : componentsOfConnectedPositions_[c])
      for (unsigned int i = 0; i < nursesPerPosition_[p->id()].size(); i++)
        // nursesInThisComponent.push_back(nursesPerPosition_[p->id()][i]);
        pNursesInThisComponent.push_back(nursesPerPosition_[p->id()][i]);

    std::stable_sort(pNursesInThisComponent.begin(),
                     pNursesInThisComponent.end(),
                     compareNursesById);
    nursesPerConnectedComponentOfPositions_.push_back(pNursesInThisComponent);
  }
}
