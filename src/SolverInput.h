/*
 * SolverInput.h
 *
 *  Created on: 18 déc. 2014
 *      Author: samuel
 */

#ifndef SOLVERINPUT_H_
#define SOLVERINPUT_H_

#include "MyTools.h"
#include "Nurse.h"
#include "Roster.h"
#include "Scenario.h"
#include "SolverInput.h"



//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r I n p u t
//
//  All the information to make an instance
//
//-----------------------------------------------------------------------------

class SolverInput{

public:

	// Constructor and Destructor
	SolverInput();
	~SolverInput();

	// Specific constructor




private:

	// Recall the "const" attributes as pointers : Nurses and Scenario informations
	//
	Scenario* pScenario_;
	vector<Nurse>* pTheNurses_;

	// Pointers to the minimum and optimum demand for each day, shift and skill
	//
	vector3D<int>* pMinDemand_, pOptDemand_;

	// pointer to the preferences of the nurses nurse (that vector must be of same length and in the same order as the nurses)
	vector<Preference>* pPreferences_;

	// pointer to the state of each nurse at the beginning of the time horizon
	//
	vector<State>* pInitState_;


	// Useful information on the demand vectors to be computed before running the algorithm if necessary
	//
	vector2D<int> totalMinDemandPerSkillPerDay_;			// Minimum aggregate number of nurses requested for a given skill for each day
	vector2D<int> totalOptDemandPerSkillPerDay_;			// Idem for optimal coverage
	int meanMinDemandPerSkill_; 							// Average minimum aggregate number of nurses requested for a given skill for one complete day over the horizon
	int meanOptDemandPerSkill_; 							// Idem for optimal coverage


};




#endif /* SOLVERINPUT_H_ */
