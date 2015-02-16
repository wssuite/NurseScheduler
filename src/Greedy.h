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

  // Constructive greedy algorithm
  // Goes through the the demands in a chronological order and assign the nurse
  // that seems most appropriate to each task (shift/skill)
  //
  void constructiveGreedy();


protected:
  // vector of live nurses that will be sorted based on the compare function
  //
  vector<LiveNurse*> theNursesSorted_;

  // vector defining the sequence according to which the positions should be
  // treated
  // the vector contains the indices of the positions
  vector<int> sequencePosition_;

  // maximum rank for a nurse
  //
  int rankMax_;

private:

  //----------------------------------------------------------------------------
  // For the constructive greedy
  //----------------------------------------------------------------------------

  // Returns true if the input nurse will respect the hard constraints if she is
  // assigned the input task
  //
  bool isFeasibleTask(const LiveNurse &nurse, int day, int shift, int skill);

  // Method that return the cost for completing the input task with the input
  // nurse
  // The cost depends on the state of the nurse, but the method will not check
  // the feasibility of the task
  //
  double costTask(const LiveNurse &nurse, int day, int shift, int skill);

  // Assign the unassigned nurses with best costs to the demand input tasks
  // nbAssigned is the number of nurses that have actually obtained a new task
  void assignBestNursesToTask(int day, int shift, int skills, int demand,
    vector<LiveNurse*>& pUnassignedNurses, int &nbAssigned);

  //----------------------------------------------------------------------------
  // For the complete roster greedy
  //----------------------------------------------------------------------------

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
