//
//  Nurse.h
//  RosterDesNurses
//

#ifndef __Nurse__
#define __Nurse__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "MyTools.h"

using std::vector



//-----------------------------------------------------------------------------
//
//  C l a s s   N u r s e
//
//  Class that contains all the attributes describing the characteristics and
//  the planning of each nurse
//
//-----------------------------------------------------------------------------

class Nurse {

public:

// Constructor and destructor
//
Nurse(char* name, int contract, int nbSkills, std::vector<int> skills) :
name_(name), contract_(contract), nbSkills_(nbSkills), skills_(skills)  {}
~Nurse();

private:

// name of the nurse
//
std::string name_;

//-----------------------------------------------------------------------------
// Constant characteristics of the nurses (no set method)
//-----------------------------------------------------------------------------

// total number of weeks and current week being planned
//
int nbWeeks_, thisWeek_;

// number of skills and vector of the skills indices
//
int nbSkills_;
vector<int> skills_;

// soft constraints of the nurse: min and max numbers of total assignments,
// min and max consecutive working days, min and max consectuve days off,
// maximum number of working week-ends and presence of absence of the
// complete week end constraints
//
int minTotalShifts_, maxTotalShifts_;
int minConsDaysWork_, maxConsDaysWork_;
int minConsDaysOff_, maxConsDaysOff_;
int maxTotalWeekEnds_;
int isCompleteWeekEnd_;


public:
   // compute the cost of the current planning from scratch
   //
   void computeCost();

   // update the cost of the planning after simple modification
   //
   void updateCost(int dayShiftAdded,)


};


#endif /* defined(__ATCSolver__CftSolver__) */
