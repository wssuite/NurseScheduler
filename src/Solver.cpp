/*
* Solver.cpp
*
*  Created on: 22 déc. 2014
*      Author: jeremy
*/


#include "Solver.h"


// Specific constructor
Solver::Solver(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState):
  pScenario_(pScenario),  pDemand_(pDemand),
  pPreferences_(pPreferences), pInitState_(pInitState) {

  }

// Destructor
Solver::~Solver(){}
