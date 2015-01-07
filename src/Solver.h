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
	void solve(SolverInput input) {}

// Should be protected (and not private) because Solver will have subclasses
protected:

	//-----------------------------------------------------------------------------
	// Inputs of the solver: they are all recorded as pointers
	//-----------------------------------------------------------------------------

	// Recall the "const" attributes as pointers : Scenario informations
	//
	Scenario* pScenario_;
	// vector<Nurse>* pTheNurses_;

	// Minimum and optimum demand for each day, shift and skill
	//
	Demand *pDemand_;

	// Preferences of the nurses (that vector must be of same length and in the
	// same order as the nurses)
	//
	Preferences* pPreferences_;

	// pointer to the state of each nurse at the beginning of the time horizon
	//
	vector<State>* pInitState_;

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

	// maximum and minimum number of working days for each nurse in the period of
	// the demand without getting any penalty for consecutive shifts
	// RqJO: this neglects the constraint of complete week-ends and the
	// preferences ; they should be added later
	//
	vector<int> maxWorkDays_, minWorkDays_;

	// total potential staffing with and without penalty
	//
	int maxTotalStaffNoPenalty_;
	int maxTotalStaff_;

	// potential staffing for each skill, with and without penalt
	vector<int> maxStaffPerSkill_;
	vector<int> maxStaffPerSkillNoPenalty_;

	// go through the nurses to collect data regarding the potential shift and
	// skill coverage of the nurses
	//
	void preprocessTheNurses();

	// check the feasibility of the demand with these nurses
	//
	bool checkFeasibility() {return true;};

};

#endif /* SOLVER_H_ */
