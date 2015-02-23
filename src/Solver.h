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
	vector<int> costConsDays_; // the same vector also accounts for consecutive days off
	vector<int> costConsShifts_;
	vector<int> costPref_;
	vector<int> costWeekEnd_;

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
	// Informative data
	//----------------------------------------------------------------------------

	// maximum and minimum number of working days for each nurse in the period of
	// the demand without getting any penalty for consecutive shifts
	// RqJO: this neglects the constraint of complete week-ends and the
	// preferences ; they should be added later
	//
	int maxWorkDays_, minWorkDays_;

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


public:
	
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

	// check the soft constraints and record the costs of the violations and the
	// remaining margin for the satisfied ones.
	//
	void checkSoftConstraints();


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
	Solver(Scenario* pScenario, Demand* pDemand,
	Preferences* pPreferences, vector<State>* pInitState);

	// Main method to solve the rostering problem for a given input
	virtual void solve() {}

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
	// then be preprocessed and get new attributes
	//
	vector<LiveNurse*> theLiveNurses_;

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
	vector3D totalStaffing_;

	// total cost under-staffing cost and under staffing cost for each triple
	// (day,shift,skill)
	//
	int totalCostUnderStaffing_;
	vector3D costUnderStaffing_;

public:

	//------------------------------------------------
	// Preprocess the data
	//------------------------------------------------

	// total potential staffing with and without penalty
	//
	int maxTotalStaffNoPenalty_;
	int maxTotalStaff_;

	// potential staffing for each skill, with and without penalt
	vector<int> maxStaffPerSkill_;
	vector<int> maxStaffPerSkillNoPenalty_;

	// rarity of the skills
	// it may depend on how many nurses have a skill and what the demand for this
	// skill is
	vector<double> skillRarity_;


public:

	// go through the nurses to collect data regarding the potential shift and
	// skill coverage of the nurses
	//
	void preprocessTheNurses();

	// compute the rarity indicator for each skill
	//
	void getSkillsRarity();

	// check the feasibility of the demand with these nurses
	//
	bool checkFeasibility() {return true;};

	//------------------------------------------------
	// Display functions
	//------------------------------------------------

	// display the whole solution
	//
	string solutionToString();

};


#endif /* SOLVER_H_ */
