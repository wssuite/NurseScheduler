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



// Attributes are public because they are const
protected:

	// Recall the "const" attributes as pointers : Nurses and Scenario informations
	//
	Scenario* pScenario_;
	vector<Nurse>* pTheNurses_;

	// Pointers to the minimum and optimum demand for each day, shift and skill
	//
	vector3D* pMinDemand_, pOptDemand_;

	// pointer to the preferences of the nurses nurse (that vector must be of same length and in the same order as the nurses)
	Preferences* pPreferences_;

	// pointer to the state of each nurse at the beginning of the time horizon
	//
	vector<State>* pInitState_;

};




#endif /* SOLVERINPUT_H_ */
