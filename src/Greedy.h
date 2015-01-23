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

  // Specific constructor and destructor
  Greedy(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState);
  ~Greedy() {}


  // Main method to solve the rostering problem for a given input
  void solve();

protected:
  // vector of live nurses that will be sorted based on the compare function
  //
  vector<LiveNurse*> theNursesSorted_;

  // vector defining the sequence according to which the positions should be
  // treated
  // the vector contains the indices of the positions
  vector<int> sequencePosition_;

private:

  // compare functions that can be used to sort the nurse before assigning them
  // schedules in the greedy algorithm
  //
  bool compareNurses(const LiveNurse  &n1, const LiveNurse &n2);

  // Build the sequence of positions reflecting the order in which the positions
  // will be treated in the greedy
  //
  void sortPositions();

};


#endif
