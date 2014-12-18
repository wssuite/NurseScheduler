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
using std::map;
using std::string;

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
Scenario(string name, int nbWeeks,
         int nbSkills, vector<string> intToSkill, map<string,int> skillToInt,
         int nbShifts, vector<string> intToShift, map<string,int> shiftToInt,
         vector<int> minConsShifts, vector<int> maxConsShifts,
         vector<int> nbForbiddenSuccessors, vector2D pForbiddenSuccessors) :
        name_(name), nbWeeks_(nbWeeks),
        nbSkills_(nbSkills), intToSkill_(intToSkill), skillToInt_(skillToInt),
        nbShifts_(nbShifts), intToShift_(intToShift), shiftToInt_(shiftToInt),
    nbForbiddenSuccessors_(nbForbiddenSuccessors), pForbiddenSuccessors_(pForbiddenSuccessors)  {
}
~Scenario();

//constant attributes are public
public:
// name of the scenario
//
const std::string name_;

// total number of weeks and current week being planned
//
const int nbWeeks_;

// number of skills, a map and a vector matching the name of each skill to an
// index and reversely
//
const int nbSkills_;
const vector<string> intToSkill_;
const map<string,int> skillToInt_;

// number of shifts, a map and a vector matching the name of each shift to an
// index and reversely
// minimum and maximum number consecutive assignments for each shift,
// and penalty for violations of these bounds
//
const int nbShifts_;
const vector<string> intToShift_;
const map<string,int> shiftToInt_;
const vector<int> minConsShifts_, maxConsShifts_;

// for each shift, the number of forbidden successors and a table containing
// the indices of these forbidden successors
//
const vector<int> nbForbiddenSuccessors_;
const vector2D pForbiddenSuccessors_;

// commenté ci-dessous pour alléger le code, plutôt mettre ces choses directement
// dans le reader vu que c'est déjà dans les nurses
// // number of contracts, and map matching the name of each contrat to an index
// //
// int nbContracts_;
// std::map<char*,int> contractToInt_;
//
// // descriptions of each contract: min and max numbers of total assignments,
// // min and max consecutive working days, min and max consectuve days off,
// // maximum number of working week-ends and presence of absence of the
// // complete week end constraints
// // each set of attributes is followed by the penalty for the violation of
// // the corresponding soft constraints
// //
// vector<int> minTotalShifts_, maxTotalShifts_;
// vector<int> minConsDaysWork_, maxConsDaysWork_;
// vector<int> minConsDaysOff_, maxConsDaysOff_;
// vector<int> maxTotalWeekEnds_;
// vector<bool> isCompleteWeekEnd_;

private:
// index of the week that is being scheduled
//
int thisWeek_;

// number of nurses, and vector of all the nurses
//
int nbNurses_;
vector<Nurse> theNurses_;


public:
// getters for the class attributes
//
int nbWeeks() {
        return nbWeeks_;
}
int thisWeek() {
        return thisWeek_;
}


// getters for the attributes of the nurses
//
int minTotalShiftsOf(int whichNurse) {
        return theNurses_[whichNurse].minTotalShifts_;
}
int maxTotalShiftsOf(int whichNurse) {
        return theNurses_[whichNurse].maxTotalShifts_;
}
int minConsDaysWorkOf(int whichNurse) {
        return theNurses_[whichNurse].minConsDaysWork_;
}
int maxConsDaysWorkOf(int whichNurse) {
        return theNurses_[whichNurse].maxConsDaysWork_;
}
int minConsDaysOffOf(int whichNurse) {
        return theNurses_[whichNurse].maxConsDaysOff_;
}
int maxConsDaysOffOf(int whichNurse) {
        return theNurses_[whichNurse].maxConsDaysOff_;
}
int maxTotalWeekEndsOf(int whichNurse) {
        return theNurses_[whichNurse].maxTotalWeekEnds_;
}
bool isCompleteWeekEndsOf(int whichNurse) {
        return theNurses_[whichNurse].isCompleteWeekEnds_;
}

// Initialize the attributes of the scenario with the content of the input
// file
//
void readScenario(string fileName);

};


#endif /* defined(__ATCSolver__CftSolver__) */
