/*
 * SolverInput.h
 *
 *  Created on: 18 déc. 2014
 *      Author: samuel
 */

#ifndef SOLVERINPUT_H_
#define SOLVERINPUT_H_

#include "tools/MyTools.h"
#include "data/Nurse.h"
#include "data/Roster.h"
#include "data/Scenario.h"


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
	PScenario pScenario_;
    std::vector<Nurse>* pTheNurses_;

	// Pointers to the minimum and optimum demand for each day, shift and skill
	//
	vector3D<int>* pMinDemand_, pOptDemand_;

	// pointer to the preferences of the nurses nurse (that vector must be of same length and in the same order as the nurses)
	PPreferences pPreferences_;

	// pointer to the state of each nurse at the beginning of the time horizon
	//
  std::vector<State>* pInitState_;

};




#endif /* SOLVERINPUT_H_ */
