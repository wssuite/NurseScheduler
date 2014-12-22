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
struct task {
  int day;
  int shift;
  int skill;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r
//
//  Schedule of a single nurse
//
//-----------------------------------------------------------------------------

class Roster{

public:

  // Constructor and destructor
  //
  Roster(int nbDays);
  ~Roster();

private:

  // number of days in the roster
  //
  int nbDays_;

  // pointer to the nurse under consideration
  //
  Nurse* pNurse_;

  // vector containing for each day the assignment of the nurse
  // the size is exactly the number of days of the roster
  //
  vector<std::pair<int, int>> tasks_;

  // vector containing for each day the state of the nurse
  // the size is the number of days of the roster plus one, since the initial
  // and the final states are of importance
  //
  vector<State> states_;

public:

  // assign a task and update the states of the nurse

};

//-----------------------------------------------------------------------------
//
//  C l a s s   S o l u t i o n
//
//  Overall schedule for all the nurses
//  Necessary to check the linking constraints on multiple nurses (insufficient
//  staffing for optimal coverage)
//
//-----------------------------------------------------------------------------

class Solution {



private:

  // pointer to the Scenario under consideration
  //
  Scenario* pScenario_;

  // pointer to the vector of nurses
  //
  vector<Nurse>* pTheNurses_;

  // number of days and number of shifts per day
  //
  int nbShifts_, nbDays_;

  // Schedule of each nurse
  //
  vector<task> schedule;

  // staffing in the roster : a 3D vector that contains the number of nurses
  //  for each task (i.e. each triple (day,shift,skill))
  //
  vector3D totalStaffing_;

  // total cost under-staffing cost and under staffing cost for task
  //
  int totalCostUnderStaffing_;
  vector3D costUnderStaffing_;

public:
  // update the roster by assigning a task to a nurse, or removing a task from
  // the schedule of a nurse
  //
  void addTask(int nurse, task t);
  void removeTast(int nurse, task t);

  // compute the constraint violation costs of all the nurses from scratch
  //  the method computes both costNurse_ and totalCostNurse_
  void computeNurseCost();

  // compute the staffing cost of the current planning from scratch
  // the method computes both costUnderStaffing_ and totalCostUnderStaffing_
  //
  void computeStaffingCost();

  // update the cost of the planning after a simple modification in the roster
  // dayShiftAdded is the
  //
  void updateCost();

  // Write the solution corresponding to the current roster
  //
  void writeSolution(std::string strCustomOutputFile);

};


#endif /* defined(__ATCSolver__CftSolver__)*/
