/*
 * Solver.h
 *
 *  Created on: 22 déc. 2014
 *      Author: jeremy
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "MyTools.h"
#include "Nurse.h"
#include "Roster.h"
#include "Scenario.h"
#include "SolverInput.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S t a t N u r s e C t
//
// The instances of this class gather the status of the constraints that relate
// to the nurses.
//
//-----------------------------------------------------------------------------

class StatCtNurse{

public:
   // Constructor and destructor
   //
   StatCtNurse();
   ~StatCtNurse();

   // number of days of the associated demand
   //
   int nbDays_;

   // costs for the violation of soft constraints
   //
   vector<int> costConsDays_;
   vector<int> costConsDaysOff_;
   vector<int> costConsShifts_;
   vector<int> costPref_;
   vector<int> costWeekEnd_;
   int costTotalDays_;
   int costTotalWeekEnds_;

   // vector of booleans equal to true if the corresponding hard contraint is
   // violated on each day
   //
   vector<bool> violSuccShifts_; // forbidden successive shifts
   vector<bool> violSkill_; // missing required skill

public:
   // initialize the statuses
   //
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
   SkillSorter (const vector<double> skillRarity) : skillRarity_(skillRarity) {}
   bool operator() (const int sk1, const int sk2) {
      return skillRarity_[sk1] > skillRarity_[sk2];
   }
private:
   vector<double> skillRarity_;
};

//-----------------------------------------------------------------------------
//
// C l a s s  S h i f t S o r t e r
//
// This class is a function object used only to compare two shifts
// the function is used to sort the skills in ascending ordrer of the number of
// forbidden successors
// (we want to treat these shifts first)
//
//-----------------------------------------------------------------------------

class ShiftSorter {
public:
   // take the field to sort by in the constructor
   ShiftSorter (const vector<int> nbForbiddenSuccessors) : nbForbiddenSuccessors_(nbForbiddenSuccessors), reverse_(false) {}
   ShiftSorter (const vector<int> nbForbiddenSuccessors, bool reverse) : nbForbiddenSuccessors_(nbForbiddenSuccessors), reverse_(reverse) {}
   bool operator() (const int sh1, const int sh2) {
      return (reverse_?-1:1)*nbForbiddenSuccessors_[sh1] < (reverse_?-1:1)*nbForbiddenSuccessors_[sh2];
   }
private:
   bool reverse_;
   vector<int> nbForbiddenSuccessors_;
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
class LiveNurse : public Nurse {

public:

   // Constructor and destructor
   //
   LiveNurse(const Nurse& nurse, Scenario* pScenario, int nbDays, int firstDay,
      State* pStateIni,	map<int,set<int> >* pWishesOff);
   ~LiveNurse();

public:

   //----------------------------------------------------------------------------
   // Pointers to background data
   //----------------------------------------------------------------------------

   // Scenario under consideration
   Scenario* pScenario_;

   //----------------------------------------------------------------------------
   // Data of the the particular period the live nurse is going to work
   //----------------------------------------------------------------------------
   int nbDays_, firstDay_;

   // Initial state
   State* pStateIni_;

   // Wishes of days off
   map<int,set<int> >* pWishesOff_;

   //----------------------------------------------------------------------------
   // Planning data
   //----------------------------------------------------------------------------

   // the current roster assigned to the nurse and the associated status of the
   // nurse constraints
   //
   Roster roster_;
   StatCtNurse statCt_;

   // a vector of rosters with no penalty and a maximum number of worked days
   //
   vector<Roster> maxFreeRosters_;

   // vector containing for each day the state of the nurse
   // the size is the number of days of the roster plus one, since the initial
   // and the final states are of importance
   //
   vector<State> states_;

   // position of the nurse: this field is deduced from the list of skills
   //
   Position* pPosition_;

   //----------------------------------------------------------------------------
   // Informative data
   //----------------------------------------------------------------------------

   // maximum and minimum number of working days for each nurse in the period of
   // the demand without getting any penalty for consecutive shifts
   // RqJO: this neglects the constraint of complete week-ends and the
   // preferences ; they should be added later
   //
   int minWorkDaysNoPenaltyConsDays_, maxWorkDaysNoPenaltyConsDays_;

   // maximum and minimum number of working days for each nurse in the period of
   // the demand without being sure to get penalty due to the total number of
   // working days
   //
   int minWorkDaysNoPenaltyTotalDays_, maxWorkDaysNoPenaltyTotalDays_;

   // minimum and maximum average number of days that can be worked per week
   // without getting penalty to the total number of working days
   //
   double minAvgWorkDaysNoPenaltyTotalDays_, maxAvgWorkDaysNoPenaltyTotalDays_;

public:
  // basic getters
  //
  Position* pPosition() const {return pPosition_;}
  State state(int day) {return states_[day];}

  // advanced getters
  //
  int totalDaysWorked() {return pStateIni_->totalDaysWorked_;}
  int totalWeekendsWorked() {return pStateIni_->totalWeekendsWorked_;}

  // basic setters
  //
  void roster(Roster &inputRoster) {roster_ = inputRoster;}

  //----------------------------------------------------------------------------
  // Methods that relate to the future capacity of a nurse
  //----------------------------------------------------------------------------

  // Compute the maximum and minimum number of working days from the input
  // current state until the input lastDay without getting any penalty for
  // consecutive working days/days-off
  //
  void computeMinMaxDaysNoPenaltyConsDay(State* pCurrentState, int lastDay,
    int &minWorkDaysNoPenaltyConsDays, int &maxWorkDaysNoPenaltyConsDays);


   //----------------------------------------------------------------------------
   // Methods that relate to the rosters of a nurse
   //----------------------------------------------------------------------------

   // assign a task at on a given day and update the states of the nurse
   //
   void assignTask(task t, int day);

   // returns true if the nurse wishes the day-shift off
   //
   bool wishesOff(int day, int shift) const;

   // returns true if the nurses reached the maximum number of consecutive worked
   // days or is resting and did not reach the minimum number of resting days yet
   // if consecutive number of shifts will only be reached by violating maximum
   // number of worked days, go to rest only if consecutive working days penalty
   // is the the larger
   //
   bool needRest(int day);

   // returns true if the nurse needs to work one more day to reach the minimum
   // number of consecutive working days or consecutive shifts
   // if consecutive number of shifts will only be reached by violating maximum
   // number of worked days, go to work only if consecutive shift penalty is
   // the larger
   bool needWork(int day);

   // return true if the nurse is free to go to rest or work more without penalty
   //
   bool isFreeToChoose(int day);

   // check the satisfaction of the hard constraints and record the violations
   // for the input roster and resulting states
   //
   void checkConstraints(const Roster& roster, const vector<State>& states, StatCtNurse& stat);

   // Build States from the roster
   //
   void buildStates();


};


// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//
bool comparePositions(Position* p1, Position* p2);

// Compare two nurses based on their position
// the function is used to sort the nurses in ascending rank of their
// position
// if their positions have the same rank, then the smaller nurse is found
// by a lexicographic comparison of the rarity of the skills of the nurses
//
bool compareNurses(LiveNurse* n1, LiveNurse* n2);


//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r
//
//  Solves the offline problem
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

enum Algorithm{GREEDY, GENCOL, STOCHASTIC_GREEDY, STOCHASTIC_GENCOL, NONE};
enum Status{UNSOLVED,FEASIBLE,INFEASIBLE,OPTIMAL};

class Solver{

public:

   // Generic constructor and destructor
   Solver() {}
   virtual ~Solver();

   // Specific constructor
   Solver(Scenario* pScenario, Demand* pDemand,
      Preferences* pPreferences, vector<State>* pInitState);

   Solver(Scenario* pScenario, Demand* pDemand,
      Preferences* pPreferences, vector<State>* pInitState,
      vector<double> minTotalShifts, vector<double> maxTotalShifts,
      vector<double> minTotalShiftsAvg, vector<double> maxTotalShiftsAvg, vector<double> weightTotalShiftsAvg,
      vector<double> maxTotalWeekendsAvg, vector<double> weightTotalWeekendsAvg);

   // Main method to solve the rostering problem for a given input and an initial solution
   virtual double solve(vector<Roster> solution = {}) { return DBL_MAX;}

   // Main method to evaluate an initial state for a given input and an initial solution
   //same as solve if not redefine
   virtual double evaluate(vector<Roster> solution = {}) {
      return solve(solution);
   }

   // Should be protected (and not private) because Solver will have subclasses
protected:

   //-----------------------------------------------------------------------------
   // Inputs of the solver: they are all recorded as pointers
   //-----------------------------------------------------------------------------

   // Recall the "const" attributes as pointers : Scenario informations
   //
   Scenario* pScenario_;

   // Minimum and optimum demand for each day, shift and skill
   //
   Demand* pDemand_;

   // Preferences of the nurses (that vector must be of same length and in the
   // same order as the nurses)
   //
   Preferences* pPreferences_;

   // pointer to the state of each nurse at the beginning of the time horizon
   //
   vector<State>* pInitState_;

   //-----------------------------------------------------------------------------
   // Manipulated data
   //-----------------------------------------------------------------------------

   // vector of LiveNurses. Initially a copy of the scenario nurses, they may
   // then be preprocessed and get enw attributes
   //
   vector<LiveNurse*> theLiveNurses_;

   // Preprocessed minimum and maximum number of working days on all the weeks
   //
   vector<double> minTotalShifts_;
   vector<double> maxTotalShifts_;

   // Interval inside of which there is no penalty for the total number of
   // working days (for each nurse)
   // This interval is computed from the max/min number of working days averaged
   // over the number of remaining weeks
   vector<double> minTotalShiftsAvg_;
   vector<double> maxTotalShiftsAvg_;

   // Penalties for values outside of [minTotalShiftsAvg_,maxTotalShiftsAvg_]
   vector<double> weightTotalShiftsAvg_;

   // Number of worked week-ends below which there is no penalty for the
   // total number of working week-ends
   // This interval is computed from the max number of working week-ends averaged
   // over the number of remaining weeks
   vector<double> maxTotalWeekendsAvg_;

   // Penalties for the number of working weekends on the current period
   // (for each nurse)
   vector<double> weightTotalWeekendsAvg_;

   //Penalties
   vector<double> weightTotalShiftsMin_, weightTotalShiftsMax_, weightTotalWeekendsMax_;

   //-----------------------------------------------------------------------------
   // Outputs of the solver
   //-----------------------------------------------------------------------------

   // a solution is a vector of rosters, one for each nurse
   // it is recorded in a vector (roster i in the vector corresponds to nurse i)
   //
   vector<Roster> solution_;

   // staffing in the solution : a 3D vector that contains the number of nurses
   //  for each triple (day,shift,skill)
   //
   vector3D satisfiedDemand_;

   // total cost under-staffing cost and under staffing cost for each triple
   // (day,shift,skill)
   //
   int totalCostUnderStaffing_;
   vector3D costUnderStaffing_;

   // vectors of nurses, skills and shifts that shall be sorted before running
   // the greedy algorithms
   //
   vector<LiveNurse*> theNursesSorted_;
   vector<int> shiftsSorted_;
   vector<int> skillsSorted_;

public:

   //------------------------------------------------
   // Preprocess the data
   //------------------------------------------------

   // total potential staffing with and without penalty
   //
   int maxTotalStaffNoPenalty_;
   int maxTotalStaffAvgWork_;

   // potential staffing for each skill, with and without penalt
   vector<double> maxStaffPerSkillNoPenalty_;
   vector<double> maxStaffPerSkillAvgWork_;

   // rarity of the skills
   // it may depend on how many nurses have a skill and what the demand for this
   // skill is
   vector<double> skillRarity_;

   // indicators related to the preprocessing
   //
   bool isPreprocessedSkills_;
   bool isPreprocessedNurses_;

   // Status of the solver
   //
   Status status_;


public:

  // Load a solution in the solver and build the states of the live nurses
  //
  void loadSolution(vector<Roster> &solution);

   //------------------------------------------------
   // Preprocess functions
   //------------------------------------------------

   // go through the nurses to collect data regarding the potential shift and
   // skill coverage of the nurses
   //
   void preprocessTheNurses();

   // Find the position of each nurse
   //
   void specifyNursePositions();

   // compute the maximum and minimum number of working days in the period of
   // the demand without getting any penalty for the total number of working days
   //
   void computeMinMaxDaysNoPenaltyTotalDays();

   // compute the maximum and minimum number of working days in the period of
   // the demand without getting any penalty for the number of consecutive
   // shifts
   // RqJO: this neglects the constraint of complete week-ends and the
   // preferences ; they should be added later
   //
   void computeMinMaxDaysNoPenaltyConsDays();

   // Compute the weights o the violation of the min/max number of working days
   // For now, the update depends only on the initial states and on the contract
   // of the nurses, on the number of days on the demand, on the number of weeks
   // already treated and on the number of weeks left
   // The required data on the nurses is mostly computed in preprocessTheNurses
   //
   void computeWeightsTotalShiftsForStochastic();

   void computeWeightsTotalShiftsForPrimalDual();

   // preprocees the skills to get their rarity
   // the value depends on the demand for this skill, on the number of nurses
   // that have the skill and on the number of skills per nurse that have the
   // skill
   //
   void preprocessTheSkills();

   // compute the rarity indicator for each skill
   //
   void getSkillsRarity();

   // Create the vector of sorted nurses
   // The nurses are ordered according to their position and the nurses that have
   // the same position are shuffled
   //
   void sortShuffleTheNurses();

   // Initialize the greedy by preprocessing all the input attributes and sorting
   // the shifts, skills, nurses
   //
   void preprocessData();

   //------------------------------------------------
   // Postprocess functions
   //------------------------------------------------

   // check the feasibility of the demand with these nurses
   //
   bool checkFeasibility();

   // get the total cost of the current solution
   // the solution is simply given by the roster of each nurse
   double solutionCost();

   //------------------------------------------------
   // Display functions
   //------------------------------------------------

   // return the status of the solution
   //
   Status getStatus() {return status_;}

   // return solution_
   //
   vector<Roster> getSolution() { return solution_; }

   // return the final states of the nurses
   //
   vector<State> getFinalStates();

   // display the whole solution in the required format
   //
   string solutionToString();

   // display the whole solution week by week for nbWeeks weeks in the required format
   //
   vector<string> solutionToString(int nbWeeks);

   // display the solution between firstDay and firstDay+nbDays in the required format
   //
   string solutionToString(int firstDay, int nbDays,  int firstWeek);

   // display the solution in a more readable format and append advanced
   // information on the solution quality
   //
   string solutionToLogString();

};


#endif /* SOLVER_H_ */
