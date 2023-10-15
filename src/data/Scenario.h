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
#include <utility>
#include <vector>
#include <string>

#include "tools/Tools.h"
#include "data/Shift.h"
#include "data/Demand.h"

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
  Weights() = default;
  Weights(double underCoverage,
          double weightAlternativeSkills,
          double weightConsShifts,
          const double weightConsDaysWork,
          const double weightConsDaysOff,
          std::vector<double> weightPreferences,
          const double weightCompleteWeekend,
          const double weightTotalShifts,
          const double weightTotalWeekends,
          double overCoverage = -1) :
      underCoverage(underCoverage),
      overCoverage(overCoverage),
      alternativeSkills(weightAlternativeSkills),
      consShifts(weightConsShifts),
      consDaysWork(weightConsDaysWork),
      consDaysOff(weightConsDaysOff),
      preferences(std::move(weightPreferences)),
      completeWeekend(weightCompleteWeekend),
      totalShifts(weightTotalShifts),
      totalWeekends(weightTotalWeekends) {}

  const double underCoverage = 30;
  const double overCoverage = -1;
  const double alternativeSkills = 20;
  const double consShifts = 15;
  const double consDaysWork = 30;
  const double consDaysOff = 30;
  const std::vector<double> preferences = {10, 20, 50, 1000};
  const double completeWeekend = 30;
  const double totalShifts = 20;
  const double totalWeekends = 30;
};

typedef std::shared_ptr<Weights> PWeights;
class Scenario;
typedef std::shared_ptr<Scenario> PScenario;
class Nurse;
typedef std::shared_ptr<Nurse> PNurse;
class Contract;
typedef std::shared_ptr<Contract> PContract;
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
// TODO(JO): the state contains the number of weekends that have been worked
//  before, but if the day is a weekend day, we need to know if the weekend
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
        PShift pShift) :
      dayId_(dayId),
      totalTimeWorked_(totalTimeWorked),
      totalWeekendsWorked_(totalWeekendsWorked),
      consDaysWorked_(consDaysWorked),
      consShifts_(consShifts),
      consDaysOff_(consDaysOff),
      pShift_(std::move(pShift)) {
    Weekend we(SATURDAY, SUNDAY);
    consWeekendWorked_ =
        we.nWeekendsInInterval(dayId - consDaysWorked + 1, dayId);
    consWeekendOff_ = we.nWeekendsInInterval(dayId - consDaysOff + 1, dayId);
  }

  // Function that appends a new day worked on a given shiftType
  // to an input state to update this state
  //
  void addDayToState(const State &prevState, const PShift &pS);

  // set each total to 0. Just keep the consecutive counters
  void resetTotal();


  // Function that appends a new day worked on a given shift
  // to an input state
  //
  // void addDayToState(const State& prevState, int newShift,
  // const PScenario pScenario);


  // Display methods: toStringINRC2 + override operator<< (easier)
  //
  std::string toString() const;
  friend std::ostream &operator<<(std::ostream &outs, const State &obj) {
    return outs << obj.toString();
  }

 public:
  // Index of the day in the planning horizon
  // WARNING : THE FIRST DAY IS ALWAYS SUPPOSED TO BE A MONDAY !!!!!!!!!!!!!
  //           If it may not be the case, the code should be modified,
  //           namely when counting the weekends worked
  //
  int dayId_;

  // Total number of days and weekends worked
  int totalTimeWorked_, totalWeekendsWorked_;

  // number of consecutive days worked
  int consDaysWorked_, consShifts_, consDaysOff_,
      consWeekendWorked_, consWeekendOff_;

  // Type of shift worked on D-Day. It can be a rest shift (=0).
  // A negative value -d means that the nurse has not been assigned a task for
  // the last d days
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
           vector<std::string> intToSkill,
           std::map<std::string, int> skillToInt,
           vector<PShift> pShifts,
           vector<std::string> intToShiftType,
           vector<int> minConsShiftsType,
           vector<int> maxConsShiftsType,
           int nbContracts,
           vector<PContract> contracts,
           int nbNurses,
           std::vector<PNurse> theNurses,
           std::map<std::string, int> nurseNameToInt,
           PWeights weights,
           string header,
           bool isINRC,
           bool isINRC2);

  // Hybrid copy constructor : this is only called when constructing
  // a new scenario that copies most parameters from the input scenario
  // but for only a subgroup of nurses.
  //
  Scenario(const PScenario &pScenario,
           const std::vector<PNurse> &theNurses,
           PDemand pDemand,
           PPreferences pWeekPreferences);

  ~Scenario();

  // first day of the horizon and
  const Day firstDay_ = Day(0);

  // type of instance
  const bool isINRC2_ = true;
  const bool isINRC_ = false;

 private:
  // name of the scenario
  //
  const std::string name_;
  // starting date of the horizon
  std::tm startDate_ = std::tm();

  // total number of weeks
  const int nWeeks_;

  // number of skills, a std::map and a std::vector matching the name
  // of each skill to an index and reversely
  const int nSkills_;
  const std::vector<std::string> intToSkill_;
  const std::map<std::string, int> skillToInt_;

  // number of shifts, a std::map and a std::vector matching the name
  // of each shift to an index and reversely
  const int nShifts_;
  std::map<std::string, int> shiftToInt_;
  const vector<PShift> pShifts_;

  // Shifts factory creating the shifts hierarchy
  ShiftsFactory shiftsFactory_;

  // number of typeshifts, a std::map and a std::vector matching the name of
  // each type shift to an index and,
  // reversely minimum and maximum number consecutive assignments
  // for each shift, and penalty for violations of these bounds.
  int nShiftTypes_;
  std::map<std::string, int> shiftTypeToInt_;

  // std::vector of possible contract types
  //
  const int nContracts_;
  const std::vector<PContract> theContracts;

  // weights of the cost functions
  const PWeights pWeights_;

  // number of nurses, and std::vector of all the nurses
  int nNurses_;
  std::vector<PNurse> pNurses_;
  std::map<std::string, int> nurseNameToInt_;

  const std::vector<int> minConsShiftType_, maxConsShiftType_;

  //------------------------------------------------
  // From the Week data file
  //------------------------------------------------
  // Name of the week
  std::string weekName_;
  // Current week demand for each DAY, SHIFT, and SKILL
  //
  PDemand pDemand_ = nullptr;

  // Shift off requests : Preferences for each nurse :
  // which (day,shift) do they want off ?
  //
  int nShiftOffRequests_;
  int nShiftOnRequests_;
  PPreferences pPreferences_ = nullptr;
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

  // True when using shift types for nodes in the graph sub problems
  bool shiftTypeGraph_ = false;

  // range of the weeks that are being scheduled
  //
  int thisWeek_;

  //------------------------------------------------


  //------------------------------------------------
  // From the custom file
  //------------------------------------------------
  //------------------------------------------------
  const std::string header_;

 private:
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

  const std::tm &startDate() { return startDate_; }
  void setStartDate(std::tm date) { startDate_ = date; }
  bool isRestShift(int shift) const {
    return pShifts_.at(shift)->isRest();
  }
  bool isWorkShift(int shift) const { return !isRestShift(shift); }
  int duration(int s) const { return pShift(s)->duration; }
  double maxDuration() const {
    auto it = max_element(pShifts_.begin(), pShifts_.end(),
                          [](const PShift &pS1, const PShift &pS2) {
                            return pS1->duration < pS2->duration;
                          });
    return (*it)->duration;
  }

  const std::string &shiftName(int i) const { return pShift(i)->name; }
  int shift(const std::string &s) const { return shiftToInt_.at(s); }

  const ShiftsFactory &shiftsFactory() const { return shiftsFactory_; }

  const PShift &pShift(int s) const { return pShifts_.at(s); }
  const PShift &pShift(const std::string &s) const {
    return pShifts_[shift(s)];
  }
  const vector<PShift> &pShifts() const { return pShifts_; }

  const PShift &pRestShift() const {
    return shiftsFactory_.pAnyRestShift()->pIncludedShifts().front();
  }

  const std::string &skillName(int i) const { return intToSkill_[i]; }
  int skillId(const std::string &s) const { return skillToInt_.at(s); }

  const std::string &shiftType(int i) const {
    return shiftsFactory_.pAnyTypeShift(i)->name;
  }
  int shiftType(const std::string &s) const { return shiftTypeToInt_.at(s); }

  const std::vector<PContract> &pContracts() const {
    return theContracts;
  }
  const PContract &pContract(int c) const {
    return theContracts[c];
  }

  int shiftIDToShiftTypeID(int s) const { return pShift(s)->type; }

  const PAbstractShift &pAnyTypeShift(int st) const {
    return shiftsFactory_.pAnyTypeShift(st);
  }

  const std::vector<PShift> &pShiftsOfType(int st) const {
    return pAnyTypeShift(st)->pIncludedShifts();
  }

  const std::vector<PShift> &pWorkShifts() const {
    return shiftsFactory_.pAnyWorkShift()->pIncludedShifts();
  }

  // getters for the private class attributes
  //
  const std::string &name() const { return name_; }
  int nWeeks() const { return nWeeks_; }
  int thisWeek() const { return thisWeek_; }
  const std::string &weekName() const { return weekName_; }
  PDemand pDemand() { return pDemand_; }
  int nShifts() const { return nShifts_; }
  int nShiftTypes() const { return nShiftTypes_; }
  int nShiftOffRequests() const { return nShiftOffRequests_; }
  int nShiftOnRequests() const { return nShiftOnRequests_; }
  PPreferences pWeekPreferences() { return pPreferences_; }
  std::vector<State> *pInitialState() { return &initialState_; }
  int nSkills() const { return nSkills_; }
  int nContracts() const { return nContracts_; }
  int nPositions() const { return nPositions_; }
  const std::vector<PPosition> &pPositions() const { return pPositions_; }
  PPosition pPosition(int p) const { return pPositions_[p]; }
  int nNurses() const { return nNurses_; }
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

  // getter for the maximum number of consecutive worked days
  // before the planning horizon
  //
  int maxConDaysWorkedInHistory() const {
    int ANS = 0;
    for (const auto &p : initialState_) {
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

  void enableShiftTypeGraph() {
    shiftTypeGraph_ = true;
  }

  bool shiftTypeGraph() const {
    return shiftTypeGraph_;
  }

  // getters for consecutive type of shifts
  int minConsShifts(int whichShift) const;
  int maxConsShifts(int whichShift) const;
  int minConsShiftsOfType(int whichShiftType) const;
  int maxConsShiftsOfType(int whichShiftType) const;

  // getter for custom file attributes

  string getHeader() const {
    return header_;
  }

  // Cost function for consecutive identical shifts
  //
  double consShiftCost(int sh, int n) const {
    if (minConsShifts(sh) - n > 0)
      return (pWeights_->consShifts * (minConsShifts(sh) - n));
    if (n - maxConsShifts(sh) > 0)
      return (pWeights_->consShifts * (n - maxConsShifts(sh)));
    return 0;
  }

  double consShiftTypeCost(int sh, int n) const {
    if (minConsShiftsOfType(sh) - n > 0)
      return (pWeights_->consShifts * (minConsShiftsOfType(sh) - n));
    if (n - maxConsShiftsOfType(sh) > 0)
      return (pWeights_->consShifts * (n - maxConsShiftsOfType(sh)));
    return 0;
  }

  // getters for the attribute of the demand
  //
  int firstDayIdOfDemand() const { return pDemand_->firstDayId_; }
  int nDays() const { return pDemand_->nDays_; }

  // Setters to class attributes

  // when reading the week file (Demand and preferences)
  //
  void setWeekName(std::string weekName) { weekName_ = std::move(weekName); }
  void setDemand(PDemand pDemand) {
    pDemand_ = std::move(pDemand);
  }
  void setTNbShiftOffRequests(int nbShiftOffRequests) {
    nShiftOffRequests_ = nbShiftOffRequests;
  }
  void setTNbShiftOnRequests(int nbShiftOnRequests) {
    nShiftOnRequests_ = nbShiftOnRequests;
  }
  void setPreferences(PPreferences pPreferences);

  // when reading the history file
  //
  void setThisWeek(int thisWeek) { thisWeek_ = thisWeek; }
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
    pDemand_ = pDemand;
  }

  void linkWithPreferences(PPreferences pPreferences);

  void pushBack(PDemand pDemand, PPreferences pPreferences);

  //------------------------------------------------
  // Display functions
  //------------------------------------------------

  // display the whole scenario
  //
  std::string toStringINRC2() const;

  //------------------------------------------------
  // Preprocess functions
  //------------------------------------------------

  // presolve the nurses to get the types
  //
  void preprocessTheNurses();

  // compute the connected components of the positions rcspp
  // (one edge between two positions indicate that they share a skill)
  //
  void computeConnectedPositions();
};

#endif  // SRC_DATA_SCENARIO_H_
