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

static const int REST_SHIFT_ID = 0;

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
            pShift_(nullptr) {}
  ~State();

  // Constructor with attributes
  State(int dayId,
        int totalTimeWorked,
        int totalWeekendsWorked,
        int consDaysWorked,
        int consShifts,
        int consDaysOff,
        const PShift &pShift) :
      dayId_(dayId),
      totalTimeWorked_(totalTimeWorked),
      totalWeekendsWorked_(totalWeekendsWorked),
      consDaysWorked_(consDaysWorked),
      consShifts_(consShifts),
      consDaysOff_(consDaysOff),
      consWeekendWorked_(Tools::nWeekendsInInterval(dayId - consDaysWorked + 1,
                                                    dayId)),
      consWeekendOff_(Tools::nWeekendsInInterval(dayId - consDaysOff + 1,
                                                 dayId)),
      pShift_(pShift) {}

  // Function that appends a new day worked on a given shiftType
  // to the previous ones
  //
  // void addNewDay(int newShiftType);

  // Function that appends a new day worked on a given shiftType
  // to an input state to update this state
  //
  void addDayToState(const State &prevState, const PShift &pS);


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

 private:
  // name of the scenario
  //
  const std::string name_;

  // total number of weeks and current week being planned
  //
  const int nWeeks_;

  // number of skills, a std::map and a std::vector matching the name
  // of each skill to an index and reversely
  //
  const int nSkills_;
  const std::vector<std::string> intToSkill_;
  const std::map<std::string, int> skillToInt_;

  // number of shifts, a std::map and a std::vector matching the name
  // of each shift to an index and reversely

  const int nShifts_;
  const std::vector<std::string> intToShift_;
  const std::map<std::string, int> shiftToInt_;
  const std::vector<int> timeDurationToWork_, shiftIDToShiftTypeID_;

  // number of typeshifts, a std::map and a std::vector matching the name of
  // each type shift to an index and,
  // reversely minimum and maximum number consecutive assignments
  // for each shift, and penalty for violations of these bounds.
  //
  const int nShiftTypes_;
  const std::vector<std::string> intToShiftType_;
  const std::map<std::string, int> shiftTypeToInt_;
  const vector2D<int> shiftTypeIDToShiftID_;
  const vector<PShift> pShifts_;

  // std::vector of possible contract types
  //
  const int nContracts_;
  const std::vector<std::string> intToContract_;
  const std::map<std::string, PConstContract> pContracts_;

  // weights of the cost functions
  const PWeights pWeights_;

  // number of nurses, and std::vector of all the nurses
  //
  int nNurses_;
  std::vector<PNurse> pNurses_;
  std::map<std::string, int> nurseNameToInt_;

  const std::vector<int> minConsShiftType_, maxConsShiftType_;

  // for each shift, the number of forbidden successors and a table containing
  // the indices of these forbidden successors
  //
  const std::vector<int> nForbiddenSuccessors_;
  const vector2D<int> forbiddenSuccessors_, forbiddenShiftTypeSuccessors_;

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
  int nShiftOffRequests_;
  int nShiftOnRequests_;
  PPreferences pWeekPreferences_ = nullptr;
  //------------------------------------------------


  //------------------------------------------------
  // From the History data file
  //------------------------------------------------
  // Initial historical state of the nurses
  //
  std::vector<State> initialState_;
  // True when the rosters created must be cyclic
  // When activated, the initial states are all default values
  bool cyclic_ = false;

  // range of the weeks that are being scheduled
  //
  int thisWeek_;
  int nWeeksLoaded_;
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
  int nPositions_;
  std::vector<PPosition> pPositions_;
  vector2D<PNurse> nursesPerPosition_;
  vector2D<PPosition> componentsOfConnectedPositions_;
  vector2D<PNurse> nursesPerConnectedComponentOfPositions_;


  //------------------------------------------------

 public:
  //------------------------------------------------
  // Getters and setters
  //------------------------------------------------

  bool isRestShift(int shift) const {
    return shiftIDToShiftTypeID_[shift] == 0;
  }
  bool isWorkShift(int shift) const { return !isRestShift(shift); }
  bool isAnyShift(int shift) const { return shift == -1; }
  bool isRestShiftType(int st) const {
    return st == 0;
  }
  bool isWorkShiftType(int st) const { return !isRestShiftType(st); }
  int duration(int s) const { return timeDurationToWork_[s]; }
  int maxDuration() const {
    return *max_element(timeDurationToWork_.begin(), timeDurationToWork_.end());
  }

  const std::string &shift(int i) const { return intToShift_[i]; }
  int shift(const std::string &s) const { return shiftToInt_.at(s); }
  const PShift &pRestShift() const { return pShifts_[REST_SHIFT_ID]; }
  const PShift &pAnyWorkShift() const {
    for (const auto &pS : pShifts_)
      if (pS->isWork())
        return pS;
    Tools::throwError("There is no work shift defined.");
    return pShifts_.back();
  }
  const PShift &pAnyWorkShift(const std::vector<int> &shifts) const {
    for (int s : shifts)
      if (pShift(s)->isWork())
        return pShift(s);
    Tools::throwError(
        "There is no work shift defined within the given vector.");
    return pShift(shifts.back());
  }
  const PShift &pShift(int s) const { return pShifts_[s]; }
  const PShift &pShift(const std::string &s) const {
    return pShifts_[shift(s)];
  }
  const vector<PShift> &pShifts() const { return pShifts_; }

  const std::string &skill(int i) const { return intToSkill_[i]; }
  int skill(const std::string &s) const { return skillToInt_.at(s); }

  const std::string &shiftType(int i) const { return intToShiftType_[i]; }
  int shiftType(const std::string &s) const { return shiftTypeToInt_.at(s); }

  const std::string &contract(int c) const { return intToContract_[c]; }
  const std::map<string, PConstContract> &pContracts() const {
    return pContracts_;
  }
  const PConstContract &pContract(const string &c) const {
    return pContracts_.at(c);
  }
  const PConstContract &pContract(int c) const {
    return pContracts_.at(contract(c));
  }

  int shiftIDToShiftTypeID(int s) const { return shiftIDToShiftTypeID_[s]; }

  const std::vector<int> &shiftTypeIDToShiftID(int st) const {
    return shiftTypeIDToShiftID_[st];
  }

  // getters for the private class attributes
  //
  const std::string &name() const { return name_; }
  int nWeeks() const { return nWeeks_; }
  int thisWeek() const { return thisWeek_; }
  int nWeeksLoaded() const { return nWeeksLoaded_; }
  const std::string &weekName() const { return weekName_; }
  PDemand pWeekDemand() { return pWeekDemand_; }
  int nShifts() const { return nShifts_; }
  int nShiftTypes() const { return nShiftTypes_; }
  int nShiftOffRequests() const { return nShiftOffRequests_; }
  int nShiftOnRequests() const { return nShiftOnRequests_; }
  PPreferences pWeekPreferences() { return pWeekPreferences_; }
  std::vector<State> *pInitialState() { return &initialState_; }
  int nSkills() const { return nSkills_; }
  const std::map<string, int> &skillsToInt() const { return skillToInt_; }
  int nContracts() const { return nContracts_; }
  int nPositions() const { return nPositions_; }
  const std::vector<PPosition> &pPositions() const { return pPositions_; }
  PPosition pPosition(int p) const { return pPositions_[p]; }
  int nNurses() { return nNurses_; }
  const PNurse &pNurse(int nurseNum) { return pNurses_[nurseNum]; }
  const vector<PNurse> &pNurses() { return pNurses_; }
  int nurse(const std::string &name) const {
    return nurseNameToInt_.at(name);
  }
  int nbOfConnectedComponentsOfPositions() const {
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

  const std::vector<int> &nForbiddenSuccessors() const {
    return nForbiddenSuccessors_;
  }

  // getter for the maximum number of consecutive worked days
  // before the planning horizon
  //
  int maxConDaysWorkedInHistory() const {
    int ANS = 0;
    for (auto p : initialState_) {
      if ((p.consDaysWorked_ > ANS) && (p.pShift_->isWork()))
        ANS = p.consDaysWorked_;
    }
    return ANS;
  }

  void enableCyclic() {
    initialState_.resize(nNurses_);
    cyclic_ = true;
  }

  bool isCyclic() const {
    return cyclic_;
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

  int minConsShifts(int whichShift) const;
  int maxConsShifts(int whichShift) const;

  int minConsShiftsOfType(int whichShiftType) const;
  int maxConsShiftsOfType(int whichShiftType) const;

  // Cost function for consecutive identical shifts
  //
  double consShiftCost(int sh, int n) const {
    if (minConsShifts(sh) - n > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (minConsShifts(sh) - n));
    if (n - maxConsShifts(sh) > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (n - maxConsShifts(sh)));
    return 0;
  }

  double consShiftTypeCost(int sh, int n) const {
    if (minConsShiftsOfType(sh) - n > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (minConsShiftsOfType(sh) - n));
    if (n - maxConsShiftsOfType(sh) > 0)
      return (pWeights_->WEIGHT_CONS_SHIFTS * (n - maxConsShiftsOfType(sh)));
    return 0;
  }

  // getters for the attribute of the demand
  //
  int firstDay() const { return pWeekDemand_->firstDay_; }
  int horizon() const { return pWeekDemand_->nDays_; }

  // Setters to class attributes

  // when reading the week file (Demand and preferences)
  //
  void setWeekName(std::string weekName) { weekName_ = weekName; }
  void setWeekDemand(PDemand pDemand) {
    pWeekDemand_ = pDemand;
  }
  void setTNbShiftOffRequests(int nbShiftOffRequests) {
    nShiftOffRequests_ = nbShiftOffRequests;
  }
  void setTNbShiftOnRequests(int nbShiftOnRequests) {
    nShiftOnRequests_ = nbShiftOnRequests;
  }
  void setWeekPreferences(PPreferences weekPreferences);

  // when reading the history file
  //
  void setThisWeek(int thisWeek) { thisWeek_ = thisWeek; }
  void addAWeek() { ++nWeeksLoaded_; }
  void setInitialState(const std::vector<State> &initialState) {
    initialState_ = initialState;
  }

  // return true if the shift shNext is a forbidden successor of shLast
  //
  bool isForbiddenSuccessorShift_Shift(int shNext, int shLast) const;
  bool isForbiddenSuccessorShiftType_ShiftType(
      int shTypeNext, int shTypeLast) const;

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
  std::string toString() const;

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
