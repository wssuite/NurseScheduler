/*
 * Solver.h
 *
 *  Created on: 18 déc. 2014
 *      Author: samuel
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
	Solver();
	virtual ~Solver();

	// Specific constructor
	Solver(Scenario* pScenario, vector<Nurse>* pTheNurses) :
		pScenario_(pScenario), pTheNurses_(pTheNurses){
	};

	// Main method to solve the rostering problem
	// minDemand: minimum number of nurses requested per day, per shift, per skill
	// optDemand: optimal number of nurses requested per day, per shift, per skill
	virtual Roster solve(SolverInput input) = 0;

private:

   // pointer to the Scenario under consideration
   //
   Scenario* pScenario_;

   // pointer to the vector of nurses
   //
   vector<Nurse>* pTheNurses_;

   // roster that shall be returned in the end of the solve function.
   // inserted here as an attribute in case it should iteratively be modified during the algorithm (easier to store as an attribute)
   //
   Roster currentRoster_;

};

#endif /* SOLVER_H_ */
