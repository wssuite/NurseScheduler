//
//  Roster.h
//  RosterDesNurses
//

#ifndef __Solver__
#define __Solver__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"
#include "Roster.h"

using std::vector;


//-----------------------------------------------------------------------------
// struct task:
// a task is a shift performed on a given day with a given skill
//
//-----------------------------------------------------------------------------
typedef struct {int shift; int skill;} task;

//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r
//
//  Schedule of a single nurse
//
//-----------------------------------------------------------------------------

class Roster{

public:

  // Constructor form no particular planning
  //
  Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map<int,std::set<int>>* pWishesOff, const State& initialState);

  // Constructor: initialize planning from an input set of tasks for the nurse
  //
  Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map<int,std::set<int>>* pWishesOff, const State& initialState,
  vector<task> inputTasks);

  // Destructor
  ~Roster();

private:

  // number of days in the roster and index of the first day
  //
  int nbDays_, firstDay_;

  // pointer to the scenario, the nurse under consideration and her wishes in
  // terms of days off
  // (the key of the map is the day and the value is the set of wishes)
  //
  const Scenario* pScenario_;
  Nurse* pNurse_;
  std::map<int,std::set<int>>* pWishesOff_;

  // vector containing for each day the assignment (shift,skill) of the nurse
  // the size is exactly the number of days of the roster
  //
  vector<task> tasks_;

  // vector containing for each day the state of the nurse
  // the size is the number of days of the roster plus one, since the initial
  // and the final states are of importance
  //
  vector<State> states_;

  // vectors of booleans that for each day, is equal to 1 if the nurse is about
  // to go from a working day to a day off or from a day off to a working day
  //
  vector<bool> switchOff_;

  // vectors of booleans that for each day, is equal to 1 if the nurse is about
  // to take another type of shift
  //
  vector<bool> switchShift_;

  // costs for the violation of soft constraints
  //
  vector<int> costConsDays_; // the same vector also accounts for consecutive days off
  vector<int> costConsShifts_;
  vector<int> costPreferences_;
  vector<int> costCompleteWeekEnd_;

  // vector of booleans equal to true if the shift assigned on each day
  // violates the consecutive shift-type succession constraint
  //
  vector<bool> violationSuccShifts_; 

public:
  // assign a task at on a given day and update the states of the nurse
  //
  void assignTask(task t, int day);

  // check the satisfaction of the hard constraints and record the violations
  //
  void checkHardConstraints();

  // check the soft constraints and record the costs of the violations and the
  // remaining margin for the satisfied ones.
  //
  void checkSoftConstraints();

};


#endif /* defined(__ATCSolver__CftSolver__)*/
