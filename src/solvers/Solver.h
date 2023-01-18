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

#ifndef SRC_SOLVERS_SOLVER_H_
#define SRC_SOLVERS_SOLVER_H_

#include <algorithm>
#include <memory>
#include <map>
#include <utility>
#include <vector>
#include <string>

#include "tools/Tools.h"
#include "tools/InputPaths.h"
#include "data/Nurse.h"
#include "data/Roster.h"
#include "data/Scenario.h"
#include "Parameters.h"
#include "solvers/DynamicWeights.h"


//-----------------------------------------------------------------------------
//
//  C l a s s   S t a t N u r s e C t
//
// The instances of this class gather the status of the constraints that relate
// to the nurses.
//
//-----------------------------------------------------------------------------

class StatCtNurse {
 public:
  // Constructor and destructor
  StatCtNurse() = default;
  ~StatCtNurse();

  // number of days of the associated demand
  int nbDays_{};

  // costs for the violation of soft constraints
  std::vector<double> costConsDays_;
  std::vector<double> costConsDaysOff_;
  std::vector<double> costConsShifts_;
  std::vector<double> costPref_;
  std::vector<double> costWeekEnd_;
  double costTotalDays_, costMissingDays_, costExceedingDays_;
  double costTotalWeekEnds_;

  // vector of booleans equal to true if the corresponding hard contraint is
  // violated on each day
  std::vector<bool> violSuccShifts_;  // forbidden successive shifts
  std::vector<bool> violSkill_;  // missing required skill

 public:
  // initialize the statuses
  void init(int nbDays);
};


//-----------------------------------------------------------------------------
//
// C l a s s  S k i l l S o r t e r
//
// This class is a function object used only to compare two skills
// the function is used to sort the skills in descending order of rarity
// (we want to treat the rarest skill first)
//
//-----------------------------------------------------------------------------

class SkillSorter {
 public:
  // take the field to sort by in the constructor
  explicit SkillSorter(std::vector<double> skillRarity)
      : skillRarity_(std::move(skillRarity)) {}
  bool operator()(const int sk1, const int sk2) {
    return skillRarity_[sk1] > skillRarity_[sk2];
  }
 private:
  std::vector<double> skillRarity_;
};

//-----------------------------------------------------------------------------
//
// C l a s s  S h i f t S o r t e r
//
// This class is a function object used only to compare two shifts
// the function is used to sort the skills in ascending order of the number of
// forbidden successors
// (we want to treat these shifts first)
//
//-----------------------------------------------------------------------------

class ShiftSorter {
 public:
  // take the field to sort by in the constructor
  explicit ShiftSorter(std::vector<PShift> pShifts, bool reverse = false)
      : reverse_(reverse),
        pShifts_(std::move(pShifts)) {}
  bool operator()(const int s1, const int s2) {
    if (reverse_)
      return pShifts_[s1]->successors.size() < pShifts_[s2]->successors.size();
    return pShifts_[s1]->successors.size() > pShifts_[s2]->successors.size();
  }
 private:
  bool reverse_;
  std::vector<PShift> pShifts_;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   L i v e N u r s e
//
// A live nurse is a nurse whose characteristics can evolve depending on
// the demand and on the planning that is being built
// They are needed in the solvers to duplicate the static nurses and define new
// attribute that can be modified.
//
// The attributes are left public, because they are meant to be modified at will
// by the solver, and because the live nurses are protected in the solver
// with no get or set method
//
//-----------------------------------------------------------------------------
class RCSolution;
class Resource;
typedef shared_ptr<Resource> PResource;
class SoftTotalShiftDurationResource;
class SoftTotalWeekendsResource;

class LiveNurse : public Nurse {
 public:
  // Constructor and destructor
//  LiveNurse(const Nurse &nurse, PScenario pScenario, int nbDays, int firstDay,
//            State *pStateIni, PPreferences pPreferences);
  LiveNurse(const Nurse &nurse, PScenario pScenario, int nbDays, int firstDay,
            State *pStateIni, PPreferences pPreferences, int num);
  ~LiveNurse();

  //----------------------------------------------------------------------------
  // Pointers to background data
  //----------------------------------------------------------------------------

  // Scenario under consideration
  PScenario pScenario_;

  //----------------------------------------------------------------------------
  // Data of the the particular period the live nurse is going to work
  //----------------------------------------------------------------------------
  int nbDays_, firstDay_;
  // WARNING: used to identify the original nurse when having sub scenario
  // You should normally use the field num_ of the nurse class
  const int nurseNum_;

  // Initial state
  State *pStateIni_;

  // Wishes of days off
  PPreferences pPreferences_;

  //----------------------------------------------------------------------------
  // Planning data
  //----------------------------------------------------------------------------

  // the current roster assigned to the nurse and the associated status of the
  // nurse constraints
  Roster roster_;
  StatCtNurse statCt_;

  // a vector of rosters with no penalty and a maximum number of worked days
  std::vector<Roster> maxFreeRosters_;

  // vector containing for each day the state of the nurse
  // the size is the number of days of the roster plus one, since the initial
  // and the final states are of importance
  std::vector<State> states_;

  // position of the nurse: this field is deduced from the list of skills
  PPosition pPosition_;

  // vector of resources
  vector<PResource> pResources_;
  // Use these two resources to update their values (LB, UB, costs)
  // for the dynamic scheduler
  std::shared_ptr<SoftTotalShiftDurationResource> totalShiftDurationResource_;
  std::shared_ptr<SoftTotalWeekendsResource> totalWeekendResource_;

  vector<PResource>& pResources() { return pResources_; }

  //----------------------------------------------------------------------------
  // Informative data
  //----------------------------------------------------------------------------

  // maximum and minimum number of working days for each nurse in the period of
  // the demand without getting any penalty for consecutive shifts
  // RqJO: this neglects the constraint of complete weekends and the
  // preferences ; they should be added later
  int minWorkDaysNoPenaltyConsDays_, maxWorkDaysNoPenaltyConsDays_;

  // maximum and minimum number of working days for each nurse in the period of
  // the demand without being sure to get penalty due to the total number of
  // working days
  int minWorkDaysNoPenaltyTotalDays_, maxWorkDaysNoPenaltyTotalDays_;

  // minimum and maximum average number of days that can be worked per week
  // without getting penalty to the total number of working days
  double minAvgWorkDaysNoPenaltyTotalDays_, maxAvgWorkDaysNoPenaltyTotalDays_;

  // basic getters
  PPosition pPosition() const { return pPosition_; }
  State state(int day) { return states_[day]; }

  // advanced getters
  int totalTimeWorked() const { return pStateIni_->totalTimeWorked_; }
  int totalWeekendsWorked() const { return pStateIni_->totalWeekendsWorked_; }

  // basic setters
  void roster(const Roster &inputRoster) { roster_ = inputRoster; }

  //----------------------------------------------------------------------------
  // Methods that relate to the future capacity of a nurse
  //----------------------------------------------------------------------------

  // Compute the maximum and minimum number of working days from the input
  // current state until the input lastDay without getting any penalty for
  // consecutive working days/days-off
  std::pair<int, int> computeMinMaxDaysNoPenaltyConsDay(
      State *pCurrentState, int lastDay);


  //----------------------------------------------------------------------------
  // Methods that relate to the rosters of a nurse
  //----------------------------------------------------------------------------
  const std::map<int, Wish> &wishes() const;

  // returns the cost of the nurse wish for the shift on the day
  double wishCostOfTheShift(int day, const PShift &pShift) const;

  // returns true if the nurses reached the maximum number of consecutive worked
  // days or is resting and did not reach the minimum number of resting days yet
  // if consecutive number of shifts will only be reached by violating maximum
  // number of worked days, go to rest only if consecutive working days penalty
  // is the larger
  bool needRest(int day);

  // returns true if the nurse needs to work one more day to reach the minimum
  // number of consecutive working days or consecutive shifts
  // if consecutive number of shifts will only be reached by violating maximum
  // number of worked days, go to work only if consecutive shift penalty is
  // the larger
  bool needWork(int day);

  // return true if the nurse is free to go to rest or work more without penalty
  bool isFreeToChoose(int day);

  // check the satisfaction of the hard constraints and record the violations
  // for the input roster and resulting states until nDays
  void checkConstraints(const Roster &roster,
                        const std::vector<State> &states,
                        int nbDays,
                        bool payExcessImmediately,
                        StatCtNurse *stat);

  // Build States from the roster
  void buildStates();

  // Print the contract type + preferences
  void printContractAndPreferences(const PScenario& pScenario) const;
};
typedef std::shared_ptr<LiveNurse> PLiveNurse;

// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//
bool comparePositions(const PPosition& p1, const PPosition& p2);

// Compare two nurses based on their position
// the function is used to sort the nurses in ascending rank of their
// position
// if their positions have the same rank, then the smaller nurse is found
// by a lexicographic comparison of the rarity of the skills of the nurses
//
bool compareNurses(const PLiveNurse& n1, const PLiveNurse& n2);

//-----------------------------------------------------------------------------
//
// C l a s s   S o l v e r
//
// Solves the offline problem
// From a given problem (number of weeks, nurses, etc.), can compute a
// solution.

//-----------------------------------------------------------------------------

class Solver {
 public:
  // Generic constructor and destructor
//  Solver() {}
  virtual ~Solver();

  // Specific constructor
  explicit Solver(const PScenario& pScenario);

  static Solver * newSolver(
      PScenario pScenario, Algorithm algo, SPType spType, SolverType sType);

  // Main method to solve the rostering problem for a given input and an
  // initial solution
  virtual double solve(const std::vector<Roster> &solution = {}) {
    return DBL_MAX;
  }

  // Main method to solve the rostering problem for a given input and an
  // initial solution and parameters
  virtual double solve(const SolverParam &parameters,
                       const std::vector<Roster> &solution = {}) {
    param_ = parameters;
    return solve(solution);
  }

  // Resolve the problem with another demand and keep the same preferences
  virtual double resolve(PDemand pDemand,
                         const SolverParam &parameters,
                         const std::vector<Roster> &solution = {}) {
    pDemand_ = std::move(pDemand);
    return solve(parameters, solution);
  }

  // if a solution, always integer
  // method is virtual if storing non integer solutions
  virtual bool isSolutionInteger() const {
    return !solution_.empty();
  }

  // add a solver that could be used to compute a solution
  void attachHeuristic(Solver *pSolver) {
    pHeuristics_.push_back(pSolver);
  }

  const std::vector<Solver*>& pHeuristics() const { return pHeuristics_; }

  // Should be protected (and not private) because Solver will have subclasses
 protected:
  //-----------------------------------------------------------------------------
  // Inputs of the solver: they are all recorded as pointers
  //-----------------------------------------------------------------------------

  // Recall the "const" attributes as pointers : Scenario informations
  PScenario pScenario_;

  // Minimum and optimum demand for each day, shift and skill
  PDemand pDemand_;

  // Preferences of the nurses (that vector must be of same length and in the
  // same order as the nurses)
  PPreferences pPreferences_;

  // pointer to the state of each nurse at the beginning of the time horizon
  std::vector<State> *pInitState_;

  // Timer started at the creation of the solver and stopped at destruction
  Tools::Timer timerTotal_;

  // current parameters of the solver (change with each solve)
  SolverParam param_;

  // store the job running the solver if it's the case
  Tools::Job job_;

  //-----------------------------------------------------------------------------
  // Manipulated data
  //-----------------------------------------------------------------------------

  // vector of LiveNurses. Initially a copy of the scenario nurses, they may
  // then be preprocessed and get enw attributes
  const std::vector<PLiveNurse> theLiveNurses_;

  // Dynamic weights. It is used for the INRC2 to manage the dynamic aspects
  DynamicWeights dynamicWeights_;

  // list of solvers that will be used as heuristics
  std::vector<Solver*> pHeuristics_;

  //-----------------------------------------------------------------------------
  // Outputs of the solver
  //-----------------------------------------------------------------------------

  // Status of the solver
  Status status_ = UNSOLVED;

  // a solution is a vector of rosters, one for each nurse
  // it is recorded in a vector (roster i in the vector corresponds to nurse i)
  double solutionCost_{};
  std::vector<Roster> solution_;

  // Objective value of the current solution
  // Warning: this value may not be updated every time it should be
  double objValue_ = XLARGE_SCORE;

  // staffing in the solution : a 3D vector that contains the number of nurses
  //  for each triple (day,shift,skill)
  vector3D<int> satisfiedDemand_;

  // total cost under-staffing cost and under staffing cost for each triple
  // (day,shift,skill)
  int totalCostUnderStaffing_;
  vector3D<int> costUnderStaffing_;

  // vectors of nurses, skills and shifts that shall be sorted before running
  // the greedy algorithms
  std::vector<PLiveNurse> theNursesSorted_;

 public:
  // check the type of instance we are solving
  bool isINRC2() { return pScenario_->isINRC2_; }
  bool isINRC() { return pScenario_->isINRC_; }

  double epsilon() const { return param_.epsilon_; }

  const SolverParam & parameters() const { return param_; }

  virtual void setParameters(const SolverParam &param) {
    param_ = param;
  }

  const DynamicWeights& getDynamicWeights() const {
    return dynamicWeights_;
  }

  void setDynamicWeightsStrategy(WeightStrategy strategy) {
    dynamicWeights_.setStrategy(strategy);
  }

  void attachJob(Tools::Job job) {
    job_ = job;
  }

  Tools::Job getJob() const { return job_; }

  //------------------------------------------------
  // Solution with rolling horizon process
  //------------------------------------------------
 protected:
  // list of days for which the integrity constraints would be relaxed
  // (The solver will not branch on these days)
  bool isPartialRelaxDays_ = false;
  std::vector<bool> isRelaxDay_;

  // list of available days and shifts. By default, all
  bool isPartialAvailable_ = false;
  vector3D<bool> nursesAvailabilities_;

  // list of nurses whose roster is fixed in current resolution
  bool isPartialFixNurses_ = false;
  std::vector<bool> isFixNurse_;

 public:
  bool isPartialRelaxed() const { return isPartialRelaxDays_; }
  bool isRelaxDay(int day) const {
    return isPartialRelaxDays_ && isRelaxDay_[day];
  }

  bool isPartialAvailable() const { return isPartialAvailable_; }
  bool isNurseAvailableOnDayShift(int nurseNum, int day, int shift) const {
    return !isPartialAvailable_ ||
           nursesAvailabilities_[nurseNum][day][shift];
  }

  bool isPartialFixNurses() const { return isPartialFixNurses_; }
  bool isFixNurse(int n) const {
    return isPartialFixNurses_ && isFixNurse_[n];
  }

  // relax/unrelax the integrality constraints of the variables corresponding
  // to input days
  virtual void relaxDays(const std::vector<bool> &isRelax) {}
  virtual void unrelaxDays(const std::vector<bool> &isUnrelax) {}

  // set availability for the days in fixDays based on
  // the current solution
  virtual void fixAvailabilityBasedOnSolution(
      const std::vector<bool> &fixDays,
      const std::vector<Roster> &solution = {}) {
    Tools::throwError(
        "Solver::fixAvailabilityBasedOnSolution() not implemented");
  }

  // set availability for the days before fixBefore (<=) based on
  // the current solution
  virtual void fixAvailabilityBasedOnSolution(
      int fixBefore, const std::vector<Roster> &solution = {}) {
    std::vector<bool> fixDays(nDays(), false);
    for (int k = 0; k <= fixBefore; ++k) fixDays[k] = true;
    fixAvailabilityBasedOnSolution(fixDays, solution);
  }

  // set the availability for each nurse
  virtual void nursesAvailabilities(
      const vector3D<bool> &nursesAvailabilities) {
    if (nursesAvailabilities.size() != nNurses())
      Tools::throwError("The input vector does not have the right size "
                        "for the method nursesAvailabilities.");
    isPartialAvailable_ = true;
    nursesAvailabilities_ = nursesAvailabilities;
  }

  virtual void resetNursesAvailabilities() {
    isPartialAvailable_ = false;
    nursesAvailabilities_.clear();
  }

  // fix/unfix all the variables corresponding to the input vector of nurses
  virtual void fixNurses(const std::vector<bool> &isFixNurse) {}
  virtual void unfixNurses(const std::vector<bool> &isUnfixNurse) {}

  // Solve the problem with a method that allows for a warm start
  virtual double rollingSolve(
      const SolverParam &parameters,
      int firstDay,
      const std::vector<Roster> &solution) {
    return 0.0;
  }

  // Special solve function for LNS
  // It is a priori the same as a regular, but it might be modified if needed
  virtual double LNSSolve(const SolverParam &parameters) {
    return 0.0;
  }

  // Solve the problem using a decomposition of the set nurses by connected
  // components of the rcspp of positions
  virtual double solveByConnectedPositions() { return 0.0; }

  //------------------------------------------------
  // Preprocess the data
  //------------------------------------------------

  // total potential staffing with and without penalty
  int maxTotalStaffNoPenalty_;
  int maxTotalStaffAvgWork_{};

  // potential staffing for each skill, with and without penalt
  std::vector<double> maxStaffPerSkillNoPenalty_;
  std::vector<double> maxStaffPerSkillAvgWork_;

  // rarity of the skills
  // it may depend on how many nurses have a skill and what the demand for this
  // skill is
  std::vector<double> skillRarity_;

  // indicators related to the preprocessing
  bool isPreprocessedSkills_;
  bool isPreprocessedNurses_;

  // Load a solution in the solver and build the states of the live nurses
  void loadSolution(const std::vector<Roster> &solution);

  // copy the solver solution to this solver
  void copySolution(Solver *pSolver);

  //------------------------------------------------
  // Preprocess functions
  //------------------------------------------------

  // go through the nurses to collect data regarding the potential shift and
  // skill coverage of the nurses
  void preprocessTheNurses();

  // Find the position of each nurse
  void specifyNursePositions();

  // compute the maximum and minimum number of working days in the period of
  // the demand without getting any penalty for the total number of working days
  void computeMinMaxDaysNoPenaltyTotalDays();

  // compute the maximum and minimum number of working days in the period of
  // the demand without getting any penalty for the number of consecutive
  // shifts
  // RqJO: this neglects the constraint of complete weekends and the
  // preferences ; they should be added later
  void computeMinMaxDaysNoPenaltyConsDays();

  // Compute the weights o the violation of the min/max number of working days
  // For now, the update depends only on the initial states and on the contract
  // of the nurses, on the number of days on the demand, on the number of weeks
  // already treated and on the number of weeks left
  // The required data on the nurses is mostly computed in preprocessTheNurses
  void boundsAndWeights(WeightStrategy strategy);

  // set the bounds based on the initial state (resources already consumed)
  void computeBoundsAccordingToInitialState();

  // Compute min/max bounds as the ratio and based on
  // the initial state (resources already consumed):
  // number of days in demand / total number of remaining days
  void computeBoundsAccordingToDemandSize();

  // Compute min/max bounds as the ratio and based on
  // the initial state (resources already consumed)
  void computeBoundsAccordingToRatio(double ratio);

  // preprocees the skills to get their rarity
  // the value depends on the demand for this skill, on the number of nurses
  // that have the skill and on the number of skills per nurse that have the
  // skill
  void preprocessTheSkills();

  // compute the rarity indicator for each skill
  void skillsRarity();

  // Create the vector of sorted nurses
  // The nurses are ordered according to their position and the nurses that have
  // the same position are shuffled
  void sortShuffleTheNurses();

  // Initialize the greedy by preprocessing all the input attributes and sorting
  // the shifts, skills, nurses
  void preprocessData();

  //------------------------------------------------
  // Postprocess functions
  //------------------------------------------------

  // check the feasibility of the demand with these nurses
  virtual bool checkFeasibility();

  // build the, possibly fractional, roster corresponding to the solution
  // currently stored in the model
  virtual vector3D<double> fractionalRoster() const { return {}; }

  // count the fraction of current solution that is integer
  double computeFractionOfIntegerInCurrentSolution() const;

// // compute the fractional penalty due to weekend that should be paid if in a
// // model with complete plannings
//  double computeFractionalWeekendPenalty();

  // get the total cost of the current solution
  // the solution is simply given by the roster of each nurse
  double computeSolutionCost(int nbDays, bool payExcessImmediately = true);

  double computeSolutionCost(bool payExcessImmediately = true) {
    return computeSolutionCost(pDemand_->nDays_, payExcessImmediately);
  }

  // get aggregate information on the solution and write them in a string
  std::string solutionStatisticsToString(int nbDays);

  //------------------------------------------------
  // Display functions
  //------------------------------------------------

  // return the status of the solution
  Status status(bool removeTIME_LIMIT = false) const {
    if (removeTIME_LIMIT && status_ == TIME_LIMIT)
      return isSolutionInteger() ? FEASIBLE : INFEASIBLE;
    return status_;
  }
  void status(Status status) { status_ = status; }

  // return/set solution_
  const std::vector<Roster> &solution() const { return solution_; }
  void addRosterToSolution(const Roster &roster) {
    solution_.push_back(roster);
  }
  void solution(const std::vector<Roster> &solution) {
    solution_ = solution;
  }

  double objValue() const { return objValue_; }

  // get the timer
  const Tools::Timer& timerTotal() const { return timerTotal_; }

  void resetTimerTotal() {
    timerTotal_.reset();
  }

  // return the solution, but only for the k first days
  std::vector<Roster> solutionAtDay(int k);

  // convert the internal solution of a solver into a interpretable one
  // return false if fails, true otherwise
  virtual bool storeSolution() { return false; }
  virtual std::string costsConstraintsToString() const { return ""; }

  // return the final states of the nurses
  std::vector<State> finalStates();

  // Returns the states(k+1) since the states start at 0
  // (hence, the state at the end of day k is state(k+1)
  std::vector<State> statesOfDay(int k);

  // display the whole solution in the required format
  std::string solutionToString();
  void solutionToTxt(string outdir);
  void solutionToXmlINRC(string outdir = "");
  std::string solutionToSolINRC();

  // display the whole solution week by week for nbWeeks weeks in the required
  // format
  std::vector<std::string> solutionToString(int nbWeeks);

  // display the solution between firstDay and firstDay+nbDays in the required
  // format
  std::string solutionToString(int firstDay, int nbDays, int firstWeek);

  // display the solution in a more readable format and append advanced
  // information on the solution quality
  std::string writeResourceCosts();
  std::string solutionToLogString();
  // if resetInitialState, do not count previous totals in the initial state
  string writeResourceCostsPerNurse(bool resetInitialState = false);
  string writeResourceCostsINRC2();

  PScenario pScenario() const { return pScenario_; }

  const std::vector<PLiveNurse> &pLiveNurses() const {
    return theLiveNurses_;
  }
  const std::vector<PLiveNurse> &sortedLiveNurses() const {
    return theNursesSorted_;
  }

  std::vector<State> *pInitialStates() const { return pInitState_; }

 public:
  // Returns the number of days over which the solver solves the problem
  int firstDayId() const { return pDemand_->firstDayId_; }
  int nDays() const { return pDemand_->nDays_; }

  // Returns the number of nurses
  int nNurses() const { return theLiveNurses_.size(); }

  // Returns the number of shifts
  int nShifts() const { return pDemand_->nShifts_; }

  // Returns the demand
  PDemand pDemand() const { return pDemand_; }

  virtual double LB() const { return -XLARGE_SCORE; }

  // Extend the rosters in the solution with the days covered by the input
  // solution
  void extendSolution(std::vector<Roster> solutionExtension);

  // Print the current best solution
  virtual std::string currentSolToString() const { return ""; }

  // When a solution of multiple consecutive weeks is available,
  // display the complete solution in the log and write the solution of
  // the weeks separately
  bool displaySolutionMultipleWeeks(bool addNumNursesInFileName = false,
                                    std::string outDir = "");

  // set the function to override
  // the default resource generation function defaultgeneratePResources
  // in the master problem
  void setGeneratePResourcesFunction(
      const std::function<std::vector<PResource>(const PLiveNurse &)> &f) {
    generatePResourcesFunc_ = f;
    useDefaultResources_ = !f;
  }

  bool useDefaultResources() const {
    return useDefaultResources_;
  }

  std::function<std::vector<PResource>(const PLiveNurse &)>
      generatePResourcesFunc() const {
    return generatePResourcesFunc_;
  }

 protected:
  bool useDefaultResources_ = true;

 private:
  // Functions to generate the resources for a given nurse for the subproblem
  // if this function is defined, it will override any default function
  std::function<std::vector<PResource>(const PLiveNurse &)>
      generatePResourcesFunc_ = nullptr;
};

#endif  // SRC_SOLVERS_SOLVER_H_
