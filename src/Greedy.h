/*
* Greedy.h
*
*  Created on: 21 jan. 2015
*      Author: jeremy
*/

#ifndef GREEDY_H_
#define GREEDY_H_

#include "Solver.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   G r e e d y
//
//  Quick solution of the problem with a greedy
//
//-----------------------------------------------------------------------------

class Greedy: public Solver{

public:

  // Generic constructor and destructor
  Greedy() {}
  virtual ~Greedy();

  // Specific constructor
  Greedy(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState):
  Solver(pScenario, pDemand, pPreferences, pInitState) {}

  // Main method to solve the rostering problem for a given input
  void solve();

private:
  // Build the sequence of positions reflecting the order in which the positions
  // will be treated in the greedy
  //
  vector<Position*> sortPositions();

};


#endif
