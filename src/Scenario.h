//
//  Scenario.h
//  RosterDesNurses
//

#ifndef __Scenario__
#define __Scenario__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "MyTools.h"
#include "Nurse.h"

using std::vector;

// the penalties for violating the soft constraints on the nurses' schedules
// are in the problem definition
// they are set as static constant values in case they need to be shared with
// other classes (e.g. solvers)
//
static const int costConsShifts_       = 15;
static const int costConsDaysWork_     = 30;
static const int costConsDaysOff_      = 30;
static const int costPreferences_      = 10;
static const int costCompleteWeekEnd_  = 30;
static const int costTotalShifts_      = 20;
static const int costTotalWeekEnds_    = 30;

//-----------------------------------------------------------------------------
//
//  C l a s s   S c e n a r i o
//
//  Class that contains all the attributes describing the scenario
//
//-----------------------------------------------------------------------------

class Scenario {

public:

// Constructor and destructor
//
Scenario() {}
~Scenario();

private:
// name of the scenario
//
std::string name_;

// total number of weeks and current week being planned
//
int nbWeeks_, thisWeek_;

// number of skills, and map matching the name of each skill to an index
//
int nbSkills_;
std::map<char*,int> skillToInt_;

// number of shifts, map matching the name of each shift to an index,
// minimum and maximum number consecutive assignments for each shift,
// and penalty for violations of these bounds
//
int nbShifts_;
std::map<char*,int> shiftToInt_;
std::vector<int> minConsShifts_, maxConsShifts_;

// for each shift, the number of forbidden successors and a table containing
// the indices of these forbidden successors
//
std::vector<int> nbForbiddenSuccessors_;
int** pForbiddenSuccessors_;

// number of contracts, and map matching the name of each contrat to an index
//
int nbContracts_;
std::map<char*,int> contractToInt_;

// descriptions of each contract: min and max numbers of total assignments,
// min and max consecutive working days, min and max consectuve days off,
// maximum number of working week-ends and presence of absence of the
// complete week end constraints
// each set of attributes is followed by the penalty for the violation of
// the corresponding soft constraints
//
std::vector<int> minTotalShifts_, maxTotalShifts_;
std::vector<int> minConsDaysWork_, maxConsDaysWork_;
std::vector<int> minConsDaysOff_, maxConsDaysOff_;
std::vector<int> maxTotalWeekEnds_;
int costTotalWeekEnds_;
std::vector<bool> isCompleteWeekEnd_;


// number of nurses, and vector of all the nurses
//
int nbNurses_;
std::vector<Nurse> theNurses_;


public:
   // Initialize the attributes of the scenario with the content of the input
   // file
   //
   void readScenario(std::string fileName);

};


#endif /* defined(__ATCSolver__CftSolver__) */
