/*
* Solver.cpp
*
*  Created on: 22 déc. 2014
*      Author: jeremy
*/


#include "Solver.h"


// Specific constructor
Solver::Solver(Scenario* pScenario, vector<Nurse>* pTheNurses, vector3D* pMinDemand,
  vector3D* pOptDemand, Preferences* pPreferences, vector<State>* pInitState_):
  pScenario_(pScenario), pTheNurses_(pTheNurses), pMinDemand_(pMinDemand),
  pOptDemand_(pOptDemand), pPreferences_(pPreferences), pInitState_(pInitState) {
    
  }

// Destructor
Solver::~Solver(){}
