/*
 * Solver.h
 *
 *  Created on: 22 dec. 2014
 *      Author: jeremy
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "tools/MyTools.h"
#include "tools/InputPaths.h"
#include "data/Nurse.h"
#include "data/Roster.h"
#include "data/Scenario.h"
#include "solvers/SolverInput.h"

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
	std::vector<int> costConsDays_;
	std::vector<int> costConsDaysOff_;
	std::vector<int> costConsShifts_;
	std::vector<int> costPref_;
	std::vector<int> costWeekEnd_;
	int costTotalDays_;
	int costTotalWeekEnds_;

	// margins with respect to the global constraints
	//
	int deltaTotalDays_;
	int deltaWeekEnds_;

	// vector of booleans equal to true if the corresponding hard contraint is
	// violated on each day
	//
	std::vector<bool> violSuccShifts_; // forbidden successive shifts
	std::vector<bool> violSkill_; // missing required skill

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
	SkillSorter (const std::vector<double> skillRarity) : skillRarity_(skillRarity) {}
	bool operator() (const int sk1, const int sk2) {
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
// the function is used to sort the skills in ascending ordrer of the number of
// forbidden successors
// (we want to treat these shifts first)
//
//-----------------------------------------------------------------------------

class ShiftSorter {
public:
	// take the field to sort by in the constructor
	ShiftSorter (const std::vector<int> nbForbiddenSuccessors) : reverse_(false), nbForbiddenSuccessors_(nbForbiddenSuccessors) {}
	ShiftSorter (const std::vector<int> nbForbiddenSuccessors, bool reverse) : reverse_(reverse), nbForbiddenSuccessors_(nbForbiddenSuccessors) {}
	bool operator() (const int sh1, const int sh2) {
		return (reverse_?-1:1)*nbForbiddenSuccessors_[sh1] < (reverse_?-1:1)*nbForbiddenSuccessors_[sh2];
	}
private:
	bool reverse_;
	std::vector<int> nbForbiddenSuccessors_;
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
class LiveNurse;
typedef std::shared_ptr<LiveNurse> PLiveNurse;
class LiveNurse : public Nurse {

public:

	// Constructor and destructor
	//
	LiveNurse(const Nurse& nurse, PScenario pScenario, int nbDays, int firstDay,
			State* pStateIni, PPreferences pPreferences);
	LiveNurse(const Nurse& nurse, PScenario pScenario, int nbDays, int firstDay,
			State* pStateIni, PPreferences pPreferences, int nurseId);
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
	const int originalNurseId_;

	// Initial state
	State* pStateIni_;

	// Wishes of days off
	PPreferences pPreferences_;

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
	std::vector<Roster> maxFreeRosters_;

	// vector containing for each day the state of the nurse
	// the size is the number of days of the roster plus one, since the initial
	// and the final states are of importance
	//
	std::vector<State> states_;

	// position of the nurse: this field is deduced from the list of skills
	//
	PPosition pPosition_;

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

	// basic getters
	//
	PPosition pPosition() const {return pPosition_;}
	State state(int day) {return states_[day];}

	// advanced getters
	//
	int totalTimeWorked() {return pStateIni_->totalTimeWorked_;}
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

	const std::map<int,std::vector<Wish> > & wishesOff() const;
  const std::map<int,std::vector<Wish> > & wishesOn() const;

	// returns true if the nurse wishes the day-shift off
	//
	bool wishesOff(int day, int shift) const;
	bool wishesOn(int day, int shift) const;

	// returns level if the nurse wishes the day-shift off : -1 otherwise
	//
	int  wishesOffLevel(int day, int shift) const;
	int  wishesOnLevel(int day, int shift) const;

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
	void checkConstraints(const Roster& roster, const std::vector<State>& states, StatCtNurse& stat);

	// Build States from the roster
	//
	void buildStates();

  // Print the contract type + preferences
  void printContractAndPreferences(PScenario pScenario) const;

};


// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//
bool comparePositions(PPosition p1, PPosition p2);

// Compare two nurses based on their position
// the function is used to sort the nurses in ascending rank of their
// position
// if their positions have the same rank, then the smaller nurse is found
// by a lexicographic comparison of the rarity of the skills of the nurses
//
bool compareNurses(PLiveNurse n1, PLiveNurse n2);


//-----------------------------------------------------------------------------
//
//  C l a s s   Print Solution
//    Allow to display a solution
//
//-----------------------------------------------------------------------------
struct PrintSolution{
	PrintSolution() {}
	virtual ~PrintSolution() {}
	virtual void save(std::vector<int>& weekIndices, std::string outdir) = 0;
	virtual std::string currentSolToString() const = 0;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r P a r a m
//
//  Structure that gather parameters for a solver. Can be given as an input of
//  the solve function of any solver
//
//-----------------------------------------------------------------------------
enum WeightStrategy { MAX, MEAN, RANDOMMEANMAX, BOUNDRATIO, NO_STRAT };
enum OptimalityLevel {UNTIL_FEASIBLE,TWO_DIVES,REPEATED_DIVES,OPTIMALITY};
static std::map<std::string,OptimalityLevel> stringToOptimalityLevel=
	{{"UNTIL_FEASIBLE",UNTIL_FEASIBLE},{"TWO_DIVES",TWO_DIVES},{"REPEATED_DIVES",REPEATED_DIVES},{"OPTIMALITY",OPTIMALITY}};

// Algorithms for the overall solution
//
enum Algorithm{GENCOL, STOCHASTIC_GENCOL, NONE};
static const std::map<std::string,Algorithm> AlgorithmsByName =
    {{"GENCOL",GENCOL},{"STOCHASTIC_GENCOL",STOCHASTIC_GENCOL},{"NONE",NONE}};

enum MySolverType { S_SCIP, S_CLP, S_Gurobi, S_Cplex, S_CBC };
static std::map<std::string,MySolverType> MySolverTypesByName =
    {{"CLP",S_CLP},{"Gurobi",S_Gurobi},{"Cplex",S_Cplex},{"CBC",S_CBC},{"SCIP",S_SCIP}};

// Solution statuses
//
enum Status{UNSOLVED,FEASIBLE,INFEASIBLE,OPTIMAL,TIME_LIMIT};
static const std::map<Status,std::string> statusToString = {{UNSOLVED,"UNSOLVED"},{FEASIBLE,"FEASIBLE"},{INFEASIBLE,"INFEASIBLE"},{OPTIMAL,"OPTIMAL"},{TIME_LIMIT,"TIME_LIMIT"}};


// Subproblem (and also master) type
enum SPType { LONG_ROTATION = 0, ALL_ROTATION = 1, ROSTER = 2};
static const std::map<std::string,SPType> SPTypesByName =
    {{"LONG",LONG_ROTATION},{"ALL",ALL_ROTATION},{"ROSTER",ROSTER}};


class SolverParam{

public:
	SolverParam(){}
	SolverParam(int verbose, OptimalityLevel level=TWO_DIVES,
		std::string outfile="outfiles/SolverResult.txt",
		std::string logfile = "outfiles/SolverLog.txt");

	//maximal solving time in s
	int maxSolvingTimeSeconds_ = LARGE_TIME;

	// tolerance
	double epsilon_ = 1e-5; // precision for the solver

	//print parameters
	int verbose_ = 0;
	bool printEverySolution_ = false;
  std::string outfile_ = "outfiles/";
  std::string logfile_ = "";
	std::vector<int> weekIndices_ = {};
	PrintSolution* saveFunction_ = nullptr;

	bool printRelaxationSol_ = false;
	bool printIntermediarySol_ = false;
	bool printFractionOfInteger_ = false;
	bool printBcpSummary_ = false;
	bool printNodes_ = false;
	bool printFinalSol_ = false;
	bool printBranchStats_ = false;
	bool printRelaxationLp_ = false;

	/* GENERAL SOLUTION PARAMETERS */
	// the solver contains only nurses that belong to the same connex component of positions
	bool isConnexComponentOfPositions_ = true;

	/* PARAMETERS OF THE BRANCH AND BOUND */
	//relative and absolute gap (with the current LB)
	//if sol below absoluteGap_, we stop immediately (the difference between two solution costs is at least 5);
	//if sol below minRelativeGap_, we stop after nbDiveIfMinGap_*dive without new incumbent;
	//if sol over relativeGap_, we stop after nbDiveIfRelGap_*dive without new incumbent.
	// the two last conditions are respected just for certain strategy.
	// If we look for optimality, we also use the last feature maxRelativeLPGapToKeepChild_:
  // if the gap between the tree best lb and the node lb is higher than this gap -> backtrack
	double absoluteGap_ = 5;
	double minRelativeGap_ = .05;
	double relativeGap_ = .1;
  double maxRelativeLPGapToKeepChild_ = .01;

	int nbDiveIfMinGap_ = 1;
	int nbDiveIfRelGap_ = 2;
	int nbDiveIfBranchOnColumns_ = 1; //number of consecutive dive if branching on columns
	int branchColumnDisjoint_ = true; // branch on columns with a disjoint argument
	int branchColumnUntilValue_ = false; // branch on columns a limit on the total rounded value

	bool solveToOptimality_ = false;

	//stop the algorithm after finding X solutions
	//if 0, the algorithm computes the relaxation if the algorithm is a column generation procedure
	int stopAfterXSolution_ = 5;

	//perform the heuristic searching an integer solution after visiting X nodes in the tree
	//if -1, do not perform the heuristic
	int performHeuristicAfterXNode_ = -1;
	double heuristicMinIntegerPercent_ = 50;

	// parameters of the stabilization : initial costs and bounds of the
	// stabilization variables
	bool isStabilization_ = false;
	bool isStabUpdateCost_ = true;
	double stabCostIni_=30.0;
	double stabCostMargin_ = 5.0;
	bool isStabUpdateBounds_ = false;
	double stabBoundIni_=1.0E-4;
	double stabBoundFactor_ = 2.0;

	// other technique against degeneracy: simply stop the solution process when
	// it starts begin degenerate
	int stopAfterXDegenerateIt_ = 9999;

	// fathom a node is upper bound is smaller than the lagrangian bound
	bool isLagrangianFathom_ = true;
	bool isLagrangianFathomRootNode_ = false; // fathom also root node

  // Parameters to remove a column from the master (both should be violated)
  // max number of consecutive iterations that a column can stay outside of the basis
  int max_inactive_iteration_ = 50;
  // min activity rate (frequency of presence in the basis since creation)
  // that a column must have to remain in the master
  double min_activity_rate_ = .1;

	/* PARAMETERS OF THE PRICER */

	// only add disjoint columns if active
	// by default false, but true for roster generation
	bool isColumnDisjoint_ = false;
  // minimum reduced cost that the best solution must have to be considered when disjoint column generation is active
  double minReducedCostDisjoint_=-1e2;
  int nbMaxColumnsToAdd_ = 20;
  int nbSubProblemsToSolve_ = 15;

	//primal-dual strategy
	WeightStrategy weightStrategy_ =  NO_STRAT;

	// Subproblem options
	//
	int sp_default_strategy_ = 0;
	int sp_nbrotationspernurse_ = 20;
	int sp_nbnursestoprice_ = 15;
  double sp_max_reduced_cost_bound_ = 0.0;
  // when the subproblems is looking for all the non dominated solutions,
  // this ratio is used to set the limit on the number of non dominated solutions to find.
  // It's a ratio expressed as a function on the max number of columns to add to the master (nbMaxColumnsToAdd_)
  double sp_columns_ratio_for_number_paths_ = 2;
  // when the number of columns generated by a strategy is lower than the ratio,
  // the pricer will use the next strategy
  double sp_min_columns_ratio_for_increase_ = .1;
  // type of the subproblem used
  SPType sp_type_ = LONG_ROTATION;

public:
	// Initialize all the parameters according to a small number of options that
	// represent the strategies we want to test
	//
	void initialize(int verbose, OptimalityLevel level);

	// Set the parameters relative to the optimality level
	//
	void setOptimalityLevel(OptimalityLevel level);
};


//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r
//
//  Solves the offline problem
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

class Solver{

public:

	// Generic constructor and destructor
	Solver() {}
	virtual ~Solver();

	// Specific constructor
	Solver(PScenario pScenario, PDemand pDemand,
		PPreferences pPreferences, std::vector<State>* pInitState);


	// Main method to solve the rostering problem for a given input and an initial solution
	virtual double solve(std::vector<Roster> solution = {}) { return DBL_MAX;}

	// Main method to solve the rostering problem for a given input and an initial solution and parameters
	virtual double solve(const SolverParam& parameters, std::vector<Roster> solution = {}){
	  param_ = parameters;
		return solve(solution);
	}


	//Resolve the problem with another demand and keep the same preferences
	//
	virtual double resolve(PDemand pDemand, const SolverParam& parameters, std::vector<Roster> solution = {}){
		pDemand_ = pDemand;
		return solve(parameters, solution);
	}

	// Should be protected (and not private) because Solver will have subclasses
protected:

	//-----------------------------------------------------------------------------
	// Inputs of the solver: they are all recorded as pointers
	//-----------------------------------------------------------------------------

	// Recall the "const" attributes as pointers : Scenario informations
	//
	PScenario pScenario_;

	// Minimum and optimum demand for each day, shift and skill
	//
	PDemand pDemand_;

	// Preferences of the nurses (that vector must be of same length and in the
	// same order as the nurses)
	//
	PPreferences pPreferences_;

	// pointer to the state of each nurse at the beginning of the time horizon
	//
	std::vector<State>* pInitState_;

	// Timer started at the creation of the solver and stopped at destruction
	Tools::Timer* pTimerTotal_;

	// current parameters of the solver (change with each solve)
	SolverParam param_;

	//-----------------------------------------------------------------------------
	// Manipulated data
	//-----------------------------------------------------------------------------

	// vector of LiveNurses. Initially a copy of the scenario nurses, they may
	// then be preprocessed and get enw attributes
	//
	std::vector<PLiveNurse> theLiveNurses_;

	// Preprocessed minimum and maximum number of working days on all the weeks
	//
	std::vector<double> minTotalShifts_, maxTotalShifts_, maxTotalWeekends_;

	// Interval inside of which there is no penalty for the total number of
	// working days (for each nurse)
	// This interval is computed from the max/min number of working days averaged
	// over the number of remaining weeks
	std::vector<double> minTotalShiftsAvg_;
	std::vector<double> maxTotalShiftsAvg_;

	// Penalties for values outside of [minTotalShiftsAvg_,maxTotalShiftsAvg_]
	std::vector<double> weightTotalShiftsAvg_;

	// Number of worked week-ends below which there is no penalty for the
	// total number of working week-ends
	// This interval is computed from the max number of working week-ends averaged
	// over the number of remaining weeks
	std::vector<double> maxTotalWeekendsAvg_;

	// Penalties for the number of working weekends on the current period
	// (for each nurse)
	std::vector<double> weightTotalWeekendsAvg_;

	//Number of min, max and weekends allowed for all nurses under the same contract
	std::vector<double> minTotalShiftsContractAvg_, maxTotalShiftsContractAvg_, maxTotalWeekendsContractAvg_;
	//Penalties for exceeding the average number of shifts allowed for all the nurses under the same contract
	std::vector<double> weightTotalShiftsContractAvg_, weightTotalWeekendsContractAvg_;

	//Penalties
	std::vector<double> weightTotalShiftsMin_, weightTotalShiftsMax_, weightTotalWeekendsMax_;

	//-----------------------------------------------------------------------------
	// Outputs of the solver
	//-----------------------------------------------------------------------------

	// Status of the solver
	//
	Status status_ = INFEASIBLE;

	// a solution is a vector of rosters, one for each nurse
	// it is recorded in a vector (roster i in the vector corresponds to nurse i)
	//
	std::vector<Roster> solution_;

	// Objective value of the current solution
	// Warning: this value may not be updated every time it should be
	//
	double objValue_;

	// staffing in the solution : a 3D vector that contains the number of nurses
	//  for each triple (day,shift,skill)
	//
	vector3D<int> satisfiedDemand_;

	// total cost under-staffing cost and under staffing cost for each triple
	// (day,shift,skill)
	//
	int totalCostUnderStaffing_;
	vector3D<int> costUnderStaffing_;

	// vectors of nurses, skills and shifts that shall be sorted before running
	// the greedy algorithms
	//
	std::vector<PLiveNurse> theNursesSorted_;
	std::vector<int> shiftsSorted_;
	std::vector<int> skillsSorted_;

public:
    double epsilon() const { return param_.epsilon_; }

	//------------------------------------------------
	// Solution with rolling horizon process
	//------------------------------------------------

	// list of starting days for which the rotations will be relaxed
	bool isPartialRelaxDays_=false;
	std::vector<bool> isRelaxDay_;
	bool isRelaxDay(int day)  const  {return !isPartialRelaxDays_?false:isRelaxDay_[day];}

	// list of starting days for which the rotations will be fixed
	bool isPartialFixDays_=false;
	std::vector<bool> isFixDay_;
	bool isFixDay(int day) const {return !isPartialFixDays_?false:isFixDay_[day];}

	// list of nurses whose roster is fixed in cirrent resolution
	bool isPartialFixNurses_=false;
	std::vector<bool> isFixNurse_;
	bool isFixNurse(int n) const {return !isPartialFixNurses_?false:isFixNurse_[n];}

	// relax/unrelax the integrality constraints of the variables corresponding to input days
	virtual void relaxDays(std::vector<bool> isRelax) {}
	virtual void unrelaxDays(std::vector<bool> isUnrelax) {}

	// fix/unfix all the variables corresponding to the input vector of days
	virtual void fixDays(std::vector<bool> isFixDay) {}
	virtual void unfixDays(std::vector<bool> isUnfixDay) {}

	// fix/unfix all the variables corresponding to the input vector of nurses
	virtual void fixNurses(std::vector<bool> isFixNurse) {}
	virtual void unfixNurses(std::vector<bool> isUnfixNurse) {}

	// Solve the problem with a method that allows for a warm start
	virtual double rollingSolve(const SolverParam& parameters, int firstDay) {return 0.0;}

	// Special solve function for LNS
	// It is a priori the same as a regular, but it might be modified if needed
	virtual double LNSSolve(const SolverParam& parameters) {return 0.0;}

	// Solve the problem using a decomposition of the set nurses by connex components
	// of the rcspp of positions
	virtual double solveByConnexPositions() {return 0.0;}

	//------------------------------------------------
	// Preprocess the data
	//------------------------------------------------

	// total potential staffing with and without penalty
	//
	int maxTotalStaffNoPenalty_;
	int maxTotalStaffAvgWork_;

	// potential staffing for each skill, with and without penalt
	std::vector<double> maxStaffPerSkillNoPenalty_;
	std::vector<double> maxStaffPerSkillAvgWork_;

	// rarity of the skills
	// it may depend on how many nurses have a skill and what the demand for this
	// skill is
	std::vector<double> skillRarity_;

	// indicators related to the preprocessing
	//
	bool isPreprocessedSkills_;
	bool isPreprocessedNurses_;

	// Load a solution in the solver and build the states of the live nurses
	//
	void loadSolution(std::vector<Roster> &solution);


	//Initialization of the rostering problem with/without solution
	//
	virtual void initialize(const SolverParam& parameters, std::vector<Roster> solution) {}

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
	void setBoundsAndWeights(WeightStrategy strategy);

	void computeWeightsTotalShiftsForStochastic();

	void computeWeightsTotalShiftsForPrimalDual(WeightStrategy strategy);

	void computeBoundsAccordingToDemandSize();

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
	virtual bool checkFeasibility();

//  // Return a solver with the algorithm specified in the options_
//  //
//  virtual Solver* setSolverWithInputAlgorithm(Algorithm algorithm, MySolverType type) const;
//  virtual Solver* setSolverWithInputAlgorithm(PScenario pScenario, std::vector<State> * stateEndOfSchedule, Algorithm algorithm, MySolverType type) const;

	// build the, possibly fractional, roster corresponding to the solution
	// currently stored in the model
	virtual vector3D<double> getFractionalRoster() const {return {};}

	// count the fraction of current solution that is integer
	//
	double computeFractionOfIntegerInCurrentSolution();

	// compute the fractional penalty due to week-end that should be paid if in a
	// model with complete plannings
	double computeFractionalWeekendPenalty();

	// get the total cost of the current solution
	// the solution is simply given by the roster of each nurse
	double computeSolutionCost(int nbDays);

	double computeSolutionCost(){return computeSolutionCost(pDemand_->nbDays_);}

   // get aggregate information on the solution and write them in a string
   //
   std::string solutionStatisticsToString();

   //------------------------------------------------
   // Display functions
   //------------------------------------------------

	// return the status of the solution
	//
	Status getStatus() const {return status_;}
	void setStatus(Status status) {status_ = status;}

	// return/set solution_
	//
	std::vector<Roster> getSolution() { return solution_; }
	void addRosterToSolution(Roster& roster) {solution_.push_back(roster);}
	void setSolution(std::vector<Roster>& solution) {solution_=solution;}

	// get the timer
	//
	Tools::Timer* getTimerTotal() const {return pTimerTotal_;}

	// return the solution, but only for the k first days
	//
	std::vector<Roster> getSolutionAtDay(int k);

	// convert the internal solution of a solver into a interpretable one
	virtual void storeSolution() {}
	virtual std::string costsConstrainstsToString() const {return "";}

	// return the final states of the nurses
	//
	std::vector<State> getFinalStates();

	// Returns the states(k+1) since the states start at 0
	// (hence, the state at the end of day k is state(k+1)
	//
	std::vector<State> getStatesOfDay(int k);

	// display the whole solution in the required format
	//
  std::string solutionToString();

	// display the whole solution week by week for nbWeeks weeks in the required format
	//
	std::vector<std::string> solutionToString(int nbWeeks);

	// display the solution between firstDay and firstDay+nbDays in the required format
	//
  std::string solutionToString(int firstDay, int nbDays,  int firstWeek);

	// display the solution in a more readable format and append advanced
	// information on the solution quality
	//
	std::string solutionToLogString();

	PScenario getScenario() const { return pScenario_; }

  const std::vector<PLiveNurse>& getLiveNurses() const { return theLiveNurses_; }
	const std::vector<PLiveNurse>& getSortedLiveNurses() const { return theNursesSorted_; }

	std::vector<State>* pInitialStates() const { return pInitState_; }

	// Returns the number of days over which the solver solves the problem
	int getFirstDay() const {return pDemand_->firstDay_;}
	int getNbDays() const {return pDemand_->nbDays_;}

	// Returns the number of nurses
	int getNbNurses() const {return theLiveNurses_.size();}

	// Returns the number of shifts
	int getNbShifts() const {return pDemand_->nbShifts_;}

	// Extend the rosters in the solution with the days covered by the input solution
	void extendSolution(std::vector<Roster> solutionExtension);

	// Print the current best solution
	virtual std::string currentSolToString() const { return ""; };

	// When a solution of multiple consecutive weeks is available, display the complete
	// solution in the log and write the solution of the weeks separately
	//
	bool displaySolutionMultipleWeeks(InputPaths inputPaths);

};


#endif /* SOLVER_H_ */
