//
//  Nurse.h
//  RosterDesNurses
//

#ifndef __Nurse__
#define __Nurse__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "MyTools.h"

using std::vector;
using std::map;
using std::set;

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//  Describes the current (or initial) state of a nurse at D-day
//
//-----------------------------------------------------------------------------
struct State{

	// number of consecutive days worked ending at D, and of consecutive days worked on the same shift ending at D (including RESTSHIFT = 0)
	// and shift worked on D-Day.
	//
	int consDaysWorked_, consShifts_;

	// Type of shift worked on D-Day
	//
	int shift_;

	// Constructor and Destructor
	State();
	~State();

	// Constructor with attributes
	State(int consDaysWorked, int consShifts, int shift) :
		consDaysWorked_(consDaysWorked), consShifts_(consShifts), shift_(shift){};

	// Function that appends a new day worked on a given shift to the previous ones
	//
	void updateWithNewDay(int newShift);
};

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   P r e f e r e n c e s
//
//  Describes the preferences of a nurse for a certain period of time
//  They are given as a vector (entry = nurseId).
//  Each element is a map<int,set<int>> whose keys are the days, and values are the sets of wished shift(s) OFF on that day.
//
//-----------------------------------------------------------------------------
struct Preferences{

	// Number of nurses
	int nNurses_;

	// For each nurse, maps the day to the set of shifts that he/she wants to have off
	//
	vector<map<int,set<int>>> wishesOff_;

	// Constructor and destructor
	Preferences();
	~Preferences();

	// Constructor with initialization to a given number of nurses
	Preferences(int nNurses);

	// For a given day, and a given shift, adds it to the wish-list for OFF-SHIFT
	void addShiftOff(int nurse, int day, int shift);

	// Adds the whole day to the wish-list
	void addDayOff(int nurse, int day);

	// True if the nurses wants that shift off
	bool wantsTheShiftOff(int nurse, int day, int shift);

	// True if the nurses wants the whole day off
	bool wantsTheShiftOff(int nurse, int day);
};



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
      Nurse(char* name, int nbSkills, std::vector<int> skills,
            int minTotalShifts, int maxTotalShifts,
            int minConsDaysWork, int maxConsDaysWork,
            int minConsDaysOff, int maxConsDaysOff,
            int maxTotalWeekEnds, bool isCompleteWeekEnds) :
              name_(name), nbSkills_(nbSkills), skills_(skills),
              minTotalShifts_(minTotalShifts), maxTotalShifts_(maxTotalShifts),
              minConsDaysWork_(minConsDaysWork), maxConsDaysWork_(maxConsDaysWork),
              minConsDaysOff_(minConsDaysOff), maxConsDaysOff_(maxConsDaysOff),
              maxTotalWeekEnds_(maxTotalWeekEnds), isCompleteWeekEnds_(isCompleteWeekEnds) {
      }
      ~Nurse();


// the constant attibutes of the nurses are public
public:

// name of the nurse
//
      const std::string name_;

//-----------------------------------------------------------------------------
// Constant characteristics of the nurses (no set method)
//-----------------------------------------------------------------------------

// number of skills and vector of the skills indices
//
      const int nbSkills_;
      const vector<int> skills_;

// soft constraints of the nurse: min and max numbers of total assignments,
// min and max consecutive working days, min and max consectuve days off,
// maximum number of working week-ends and presence of absence of the
// complete week end constraints
//
      const int minTotalShifts_, maxTotalShifts_;
      const int minConsDaysWork_, maxConsDaysWork_;
      const int minConsDaysOff_, maxConsDaysOff_;
      const int maxTotalWeekEnds_;
      const int isCompleteWeekEnds_;



      };


#endif /* defined(__ATCSolver__CftSolver__) */
