/*
* Greedy.h
*
*  Created on: 21 jan. 2015
*      Author: jeremy
*/

#ifndef GREEDY_H_
#define GREEDY_H_

#include "solvers/Solver.h"

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
  double solve(vector<Roster> solution = {});

  // Constructive greedy algorithm
  // Goes through the the demands in a chronological order and assign the nurse
  // that seems most appropriate to each task (shift/skill)
  // Returns true if the minimum demand could be covered, and false otherwise
  //
  bool constructiveGreedy();


protected:

  // vector defining the sequence according to which the positions should be
  // treated
  // the vector contains the indices of the positions
  vector<int> sequencePosition_;

  // maximum rank for a nurse
  //
  int rankMax_;

  // weights that are used to penalize future violations when computing the cost
  // of an assignment
  //
  double weightNbForbidden_;
  double weightRank_;
  double weightCoverMin_;

  // vector of excess in available nurses for a task with respect to minimum
  // demand
  // necessary to avoid affecting tasks to higher ranked nurses for better cost
  // and then risk to hit an infeasible solution
  vector3D shiftDemand_;

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
  double costTask(LiveNurse &nurse, int day, int shift, int skill,
    vector<State>* states = NULL);

  // Assign the unassigned nurses with best costs to the demand input tasks
  // nbAssigned is the number of nurses that have actually obtained a new task
  //
  void assignBestNursesToTask(int day, int shift, int skills, int demand,
    vector<LiveNurse*>& pUnassignedNurses, int &nbAssigned, bool isMinimum);

  // Necessary actions when assigning a task to a nurse
  //
  void assignTaskToNurse(LiveNurse &nurse, int day, int shift, int skill);

  // When assigning a new task to a nurse that has an unassigned day just
  // just before, find tasks/rest periods to assign in the preceeding block of
  // unassinged days
  // The input day is the last day of the block
  //
  void fillTheGaps(LiveNurse &nurse, int day);

  // Recursive function that tries to add rest or work to the input statesBlock
  // and returns the best result if nbUnassigned days are treated
  // store the assigned shifts and skills in vectors of int
  //
  double bestStatesBlock_rec(LiveNurse &nurse, vector<State> &statesBlock,
    vector<int> &shifts, vector<int> &skills, int dayFirst,int nbUnassigned, double costIni);

  //----------------------------------------------------------------------------
  // For the initialization of the constructive greedy
  //----------------------------------------------------------------------------

  // Preprocess the demand to infer implicit demand on day d that is made
  // necessary by the demand on day d+1
  // For instance, assume that for a given skill sk, there is no demand for the
  // shift Early on day d, and a demand of 2 on day d+1 for the same shift. Then,
  // Due to the list of forbidden successor, it is then necessary that at least
  // two nurses with skill sk are either resting or taking the shift Early on day d
  //
  void computeImplicitDemand();


};


#endif
