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
#include "SubProblem.h"


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
	LiveNurse(const Nurse& nurse);
	~LiveNurse();

public:

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

	// the current roster assigned to the nurse
	//
	Roster roster_;

	// a vector of rosters with no penalty and a maximum number of worked days
	//
	vector<Roster> maxFreeRosters_;


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
	Demand *pDemand_;

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
	vector<int> skillRarity_;


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

};


#endif /* SOLVER_H_ */
