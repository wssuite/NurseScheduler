//
//  Roster.h
//  RosterDesNurses
//

#ifndef __Roster__
#define __Roster__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"

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
//  Overall schedule for all the nurses
//  Necessary to check the linking constraints on multiple nurses (insufficient
//  staffing for optimal coverage)
//
//-----------------------------------------------------------------------------

class Roster {

public:

// Constructor and destructor
//
Roster(Scenario* pScenario, vector<Nurse>* pTheNurses) :
        pScenario_(pScenario), pTheNurses_(pTheNurses)
}
~Roster();

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

   // staffing in the roster : a 3D vector that contains the number of nurses
   //  for each task (i.e. each triple (day,shift,skill))
   //
   vector3D totalStaffing_;

   // total cost under-staffing cost and under staffing cost for task
   //
   int totalCostUnderStaffing_;
   vector3D costUnderStaffing_;

public:
   // update the roster by assigning a task to a nurse, removing a task from
   // the roster, or swapping activity from one removed task to one added task
   inline void addAssignment(task t) {
      totalStaffing[t.day][t.shift][t.skill]++;}
   void removeAssignement(task taskRemoved);
   void swapAssignement(task taskRemoved, task taskAdded);

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
   void updateCost(dayshift dayShiftAdded, dayshift dayShiftRemoved);

   // Write the solution corresponding to the current roster
   //
   void writeSolution(std::string strCustomOutputFile);

};


#endif /* defined(__ATCSolver__CftSolver__)*/
