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

#ifndef SRC_DATA_SCENARIO_H_
#define SRC_DATA_SCENARIO_H_

#include <memory>
#include <map>
#include <vector>
#include <string>

#include "tools/Tools.h"
#include "data/Shift.h"
#include "Demand.h"

// the penalties for violating the soft constraints on the nurses' schedules
// are in the problem definition
// they are set as static constant values in case they need to be shared with
// other classes (e.g. solvers)
//
enum PREF_LEVEL { WEAK = 0, MODERATE = 1, STRONG = 2, COMPULSORY = 3 };
const std::map<PREF_LEVEL, std::string>
    levelsToString = {{WEAK, "WEAK"}, {MODERATE, "MODERATE"},
                      {STRONG, "STRONG"}, {COMPULSORY, "COMPULSORY"}};
struct Weights {
  Weights() {}
  Weights(double weightOptimalDemand,
          double weightConsShifts,
          const double weightConsDaysWork,
          const double weightConsDaysOff,
          const std::vector<double> &weightPreferencesOff,
          const std::vector<double> &weightPreferencesOn,
          const double weightCompleteWeekend,
          const double weightTotalShifts,
          const double weightTotalWeekends) :
      WEIGHT_OPTIMAL_DEMAND(weightOptimalDemand),
      WEIGHT_CONS_SHIFTS(weightConsShifts),
      WEIGHT_CONS_DAYS_WORK(weightConsDaysWork),
      WEIGHT_CONS_DAYS_OFF(weightConsDaysOff),
      WEIGHT_PREFERENCES_OFF(weightPreferencesOff),
      WEIGHT_PREFERENCES_ON(weightPreferencesOn),
      WEIGHT_COMPLETE_WEEKEND(weightCompleteWeekend),
      WEIGHT_TOTAL_SHIFTS(weightTotalShifts),
      WEIGHT_TOTAL_WEEKENDS(weightTotalWeekends) {}

  const double WEIGHT_OPTIMAL_DEMAND = 30;
  const double WEIGHT_CONS_SHIFTS = 15;
  const double WEIGHT_CONS_DAYS_WORK = 30;
  const double WEIGHT_CONS_DAYS_OFF = 30;
  const std::vector<double> WEIGHT_PREFERENCES_OFF = {10, 20, 50, 1000};
  const std::vector<double> WEIGHT_PREFERENCES_ON = {-10, -20, -50, 1000};
  const double WEIGHT_COMPLETE_WEEKEND = 30;
  const double WEIGHT_TOTAL_SHIFTS = 20;
  const double WEIGHT_TOTAL_WEEKENDS = 30;
};

typedef std::shared_ptr<Weights> PWeights;
class Scenario;
typedef std::shared_ptr<Scenario> PScenario;
class Nurse;
typedef std::shared_ptr<Nurse> PNurse;
class Contract;
typedef std::shared_ptr<const Contract> PConstContract;
class Position;
typedef std::shared_ptr<Position> PPosition;
class Preferences;
typedef std::shared_ptr<Preferences> PPreferences;

//-----------------------------------------------------------------------------
//
//  C l a s s   S t a t e
//
//  Describes the current (or initial) state of a nurse at D-day
//
//-----------------------------------------------------------------------------
// TODO(JO): the state contains the number of week-ends that have been worked
//  before, but if the day is a weekend day, we need to know if the week-end
//  is already counted in this total
class State {
 public:
  // Constructor and Destructor
  State() : dayId_(0),
            totalTimeWorked_(0),
            totalWeekendsWorked_(0),
            consDaysWorked_(0),
            consShifts_(0),
            consDaysOff_(0),
            consWeekendWorked_(0),
            consWeekendOff_(0),
            shiftType_(0),
            pShift_(nullptr) {}
  ~State();

  // Constructor with attributes
  State(int dayId,
        int totalTimeWorked,
        int totalWeekendsWorked,
        int consDaysWorked,
        int consShifts,
        int consDaysOff,
        int shiftType,
        int shift) :
      dayId_(dayId),
      totalTimeWorked_(totalTimeWorked),
      totalWeekendsWorked_(totalWeekendsWorked),
      consDaysWorked_(consDaysWorked),
      consShifts_(consShifts),
      consDaysOff_(consDaysOff),
      consWeekendWorked_(0),
      consWeekendOff_(0),
      shiftType_(shiftType),
      shift_(shift),
      pShift_(std::make_shared<Shift>(shift, shiftType)) {}

  // Function that appends a new day worked on a given shiftType
  // to the previous ones
  //
  // void addNewDay(int newShiftType);

  // Function that appends a new day worked on a given shiftType
  // to an input state to update this state
  //
  void addDayToState(const State &prevState,
                     int newShiftType,
                     int newShift,
                     int timeWorked);


  // Function that appends a new day worked on a given shift
  // to an input state
  //
  // void addDayToState(const State& prevState, int newShift,
  // const PScenario pScenario);


  // Display methods: toString + override operator<< (easier)
  //
  std::string toString();
  friend std::ostream &operator<<(std::ostream &outs, State obj) {
    return outs << obj.toString();
  }

 public:
  // Index of the day in the planning horizon
  // WARNING : THE FIRST DAY IS ALWAYS SUPPOSED TO BE A MONDAY !!!!!!!!!!!!!
  //           If it may not be the case, the code should be modified,
  //           namely when counting the weekends worked
  //
  int dayId_;

  // Total nummber of days and weekends worked
  //
  int totalTimeWorked_, totalWeekendsWorked_;

  // number of consecutive days worked ending at D,
  // and of consecutive days worked on the same shiftType
  // ending at D (including RESTSHIFT = 0) and shiftType worked on D-Day.
  //
  int consDaysWorked_, consShifts_, consDaysOff_,
      consWeekendWorked_, consWeekendOff_;

  // Type of shift worked on D-Day. It can be a rest shift (=0).
  // A negative value -d means that the nurse has not been assigned a task for
  // the last d days
  //
  int shiftType_;

  int shift_;

  PShift pShift_;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   S c e n a r i o
//
//  Class that contains all the attributes describing the scenario
//
//-----------------------------------------------------------------------------

class Scenario {
 public:
  // Constructor and destructor
  //
  Scenario(std::string name,
           int nbWeeks,
           int nbSkills,
           std::vector<std::string> intToSkill,
           std::map<std::string, int> skillToInt,
           int nbShifts,
           std::vector<std::string> intToShift,
           std::map<std::string, int> shiftToInt,
           std::vector<int> hoursToWork,
           std::vector<int> shiftIDToShiftTypeID,
           int nbShiftsType,
           std::vector<std::string> intToShiftType,
           std::map<std::string, int> shiftTypeToInt,
           std::vector<std::vector<int> > shiftTypeIDToShiftID,
           std::vector<int> minConsShiftsType,
           std::vector<int> maxConsShiftsType,
           std::vector<int> nbForbiddenSuccessors,
           vector2D<int> forbiddenSuccessors,
           vector<PShift> pShifts,
           int nbContracts,
           std::vector<std::string> intToContract,
           std::map<std::string, PConstContract> contracts,
           int nbNurses,
           std::vector<PNurse> theNurses,
           std::map<std::string, int> nurseNameToInt,
           PWeights weights);

  // Hybrid copy constructor : this is only called when constructing
  // a new scenario that copies most parameters from the input scenario
  // but for only a subgroup of nurses.
  //
  Scenario(PScenario pScenario,
           const std::vector<PNurse> &theNurses,
           PDemand pDemand,
           PPreferences pWeekPreferences);

  ~Scenario();

  // constant attributes are public
  // TODO(JO): the above member variables are all public, but we still
  //  declared public getters, why do that?
 public:
  // name of the scenario
  //
  const std::string name_;

  // total number of weeks and current week being planned
  //
  const int nbWeeks_;

  // number of skills, a std::map and a std::vector matching the name
  // of each skill to an index and reversely
  //
  const int nbSkills_;
  const std::vector<std::string> intToSkill_;
  const std::map<std::string, int> skillToInt_;

  // number of shifts, a std::map and a std::vector matching the name
  // of each shift to an index and reversely

  const int nbShifts_;
  const std::vector<std::string> intToShift_;
  const std::map<std::string, int> shiftToInt_;
  const std::vector<int> timeDurationToWork_, shiftIDToShiftTypeID_;

  bool isRestShift(int shift) const {
    return shiftIDToShiftTypeID_[shift] == 0;
  }
  bool isWorkShift(int shift) const { return !isRestShift(shift); }
  bool isAnyShift(int shift) const { return shift == -1; }
  bool isRestShiftType(int st) const {
    return st == 0;
  }
  bool isWorkShiftType(int st) const { return !isRestShiftType(st); }

  // number of typeshifts, a std::map and a std::vector matching the name of
  // each type shift to an index and,
  // reversely minimum and maximum number consecutive assignments
  // for each shift, and penalty for violations of these bounds.
  //
  const int nbShiftsType_;
  const std::vector<std::string> intToShiftType_;
  const std::map<std::string, int> shiftTypeToInt_;
  const vector2D<int> shiftTypeIDToShiftID_;
  const vector<PShift> pShifts_;

  // std::vector of possible contract types
  //
  const int nbContracts_;
  const std::vector<std::string> intToContract_;
  const std::map<std::string, PConstContract> pContracts_;

  // number of nurses, and std::vector of all the nurses
  //
  const int nbNurses_;
  const std::vector<PNurse> theNurses_;
  std::map<std::string, int> nurseNameToInt_;

  // weights of the cost functions
  const PWeights pWeights_;

 private:
  const std::vector<int> minConsShiftType_, maxConsShiftType_;

  // for each shift, the number of forbidden successors and a table containing
  // the indices of these forbidden successors
  //
  const std::vector<int> nbForbiddenSuccessors_;
  const vector2D<int> forbiddenSuccessors_;

  //------------------------------------------------
  // From the Week data file
  //------------------------------------------------
  // Name of the week
  std::string weekName_;
  // Current week demand for each DAY, SHIFT, and SKILL
  //
  PDemand pWeekDemand_ = nullptr;

  // Shift off requests : Preferences for each nurse :
  // which (day,shift) do they want off ?
  //
  int nbShiftOffRequests_;
  int nbShiftOnRequests_;
  PPreferences pWeekPreferences_ = nullptr;
  //------------------------------------------------


  //------------------------------------------------
  // From the History data file
  //------------------------------------------------
  // Initial historical state of the nurses
  //
  std::vector<State> initialState_;
  // range of the weeks that are being scheduled
  //
  int thisWeek_;
  int nbWeeksLoaded_;
  //------------------------------------------------


  //------------------------------------------------
  // From the custom file
  //------------------------------------------------
  //------------------------------------------------

  //------------------------------------------------
  // From the preprocessing of the nurses
  //------------------------------------------------
  // std::vector of existing positions
  //
  int nbPositions_;
  std::vector<PPosition> pPositions_;
  vector2D<PNurse> nursesPerPosition_;
  vector2D<PPosition> componentsOfConnectedPositions_;
  vector2D<PNurse> nursesPerConnectedComponentOfPositions_;


  //------------------------------------------------

 public:
  //------------------------------------------------
  // Getters and setters
  //------------------------------------------------

  // getters for the private class attributes
  //
  int nbWeeks() { return nbWeeks_; }
  int thisWeek() { return thisWeek_; }
  int nbWeeksLoaded() { return nbWeeksLoaded_; }
  std::string weekName() { return weekName_; }
  PDemand pWeekDemand() { return pWeekDemand_; }
  int nbShifts() { return nbShifts_; }
  int nbShiftOffRequests() { return nbShiftOffRequests_; }
  int nbShiftOnRequests() { return nbShiftOnRequests_; }
  PPreferences pWeekPreferences() { return pWeekPreferences_; }
  std::vector<State> *pInitialState() { return &initialState_; }
  int nbSkills() { return nbSkills_; }
  int nbPositions() { return nbPositions_; }
  const std::vector<PPosition> &pPositions() const { return pPositions_; }
  PPosition pPosition(int p) const { return pPositions_[p]; }
  int nbNurses() { return nbNurses_; }
  int nbOfConnectedComponentsOfPositions() {
    return componentsOfConnectedPositions_.size();
  }
  const std::vector<PPosition> &componentOfConnectedPositions(int c) const {
    return componentsOfConnectedPositions_[c];
  }
  const std::vector<PNurse> &
  nursesInConnectedComponentOfPositions(int c) const {
    return nursesPerConnectedComponentOfPositions_[c];
  }
  const Weights &weights() const { return *pWeights_; }

  int nbForbiddenSuccessorsShift(int shift) {
    int shiftType = shiftIDToShiftTypeID_[shift];
    return nbForbiddenSuccessors_[shiftType];
  }

  int nbForbiddenSuccessorsShiftType(int shiftType) {
    return nbForbiddenSuccessors_[shiftType];
  }

  const std::vector<int> &nbForbiddenSuccessors() {
    return nbForbiddenSuccessors_;
  }

  // getter for the maximum number of consecutive worked days
  // before the planning horizon
  //
  int maxConDaysWorkedInHistory() {
    int ANS = 0;
    for (auto p : initialState_) {
      if ((p.consDaysWorked_ > ANS) && (p.shiftType_ > 0))
        ANS = p.consDaysWorked_;
    }
    return ANS;
  }

  // getters for the attributes of the nurses
  //
  const int minTotalShiftsOf(int whichNurse) const;
  int maxTotalShiftsOf(int whichNurse) const;
  int minConsDaysWorkOf(int whichNurse) const;
  int maxConsDaysWorkOf(int whichNurse) const;
  int minConsDaysOffOf(int whichNurse) const;
  int maxConsDaysOffOf(int whichNurse) const;
  int maxTotalWeekendsOf(int whichNurse) const;
  bool isCompleteWeekendsOf(int whichNurse) const;

  // getters for consecutive type of shifts

  int minConsShiftsOfTypeOf(int whichShift);
  int maxConsShiftsOfTypeOf(int whichShift);

  int minConsShiftsOf(int whichShiftType);
  int maxConsShiftsOf(int whichShiftType);

  // Cost function for consecutive identical shifts
  //
  double consShiftCost(int sh, int n) {
    if (minConsShiftsOfTypeOf(sh) - n > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (minConsShiftsOfTypeOf(sh) - n));
    if (n - maxConsShiftsOfTypeOf(sh) > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (n - maxConsShiftsOfTypeOf(sh)));
    return 0;
  }

  double consShiftTypeCost(int sh, int n) {
    if (minConsShiftsOf(sh) - n > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (minConsShiftsOf(sh) - n));
    if (n - maxConsShiftsOf(sh) > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (n - maxConsShiftsOf(sh)));
    return 0;
  }

  // getters for the attribute of the demand
  //
  int firstDay() { return pWeekDemand_->firstDay_; }
  int nbDays() { return pWeekDemand_->nDays_; }

  // Setters to class attributes

  // when reading the week file (Demand and preferences)
  //
  void setWeekName(std::string weekName) { weekName_ = weekName; }
  void setWeekDemand(PDemand pDemand) {
    pWeekDemand_ = pDemand;
  }
  void setTNbShiftOffRequests(int nbShiftOffRequests) {
    nbShiftOffRequests_ = nbShiftOffRequests;
  }
  void setTNbShiftOnRequests(int nbShiftOnRequests) {
    nbShiftOnRequests_ = nbShiftOnRequests;
  }
  void setWeekPreferences(PPreferences weekPreferences);

  // when reading the history file
  //
  void setThisWeek(int thisWeek) { thisWeek_ = thisWeek; }
  void addAWeek() { ++nbWeeksLoaded_; }
  void setInitialState(const std::vector<State> &initialState) {
    initialState_ = initialState;
  }

  // return true if the shift shNext is a forbidden successor of shLast
  //
  bool isForbiddenSuccessorShift_Shift(int shNext, int shLast);
  bool isForbiddenSuccessorShift_ShiftType(int shNext, int shTypeLast);
  bool isForbiddenSuccessorShiftType_Shift(int shTypeNext, int shLast);
  bool isForbiddenSuccessorShiftType_ShiftType(int shTypeNext, int shTypeLast);

  // update the scenario to treat a new week
  //
  void updateNewWeek(PDemand pDemand,
                     PPreferences pPreferences,
                     const std::vector<State> &initialStates);

  // Link the scenario with the Demand and the Preferences
  //
  void linkWithDemand(PDemand pDemand) {
    weekName_ = pDemand->name_;
    pWeekDemand_ = pDemand;
  }

  void linkWithPreferences(PPreferences pPreferences);

  //------------------------------------------------
  // Display functions
  //------------------------------------------------

  // display the whole scenario
  //
  std::string toString();

  //------------------------------------------------
  // Preprocess functions
  //------------------------------------------------

  // preprocess the nurses to get the types
  //
  void preprocessTheNurses();

  // compute the connected components of the positions rcspp
  // (one edge between two positions indicate that they share a skill)
  //
  void computeConnectedPositions();
};

#endif  // SRC_DATA_SCENARIO_H_
