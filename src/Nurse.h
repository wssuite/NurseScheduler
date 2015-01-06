//
//  Nurse.h
//  RosterDesNurses
//

#ifndef __Nurse__
#define __Nurse__

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "MyTools.h"

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------
class Contract{

public:
	// Name of the contract
	//
	const string name_;

	// Minimum and maximum total number of shifts over the time period
	//
	const int minTotalShifts_, maxTotalShifts_;

	// Minimum and maximum number of consecutive days worked
	//
  const int minConsDaysWork_, maxConsDaysWork_;

  // Minimum and maximum number of consecutive days off
  //
  const int minConsDaysOff_, maxConsDaysOff_;

  // Maximum number of weekends worked, and complete weekend constraint
  //
  const int maxTotalWeekends_;
  const int isCompleteWeekends_;

  // Constructor and Destructor
  //
  Contract(string name, int minTotalShifts, int maxTotalShifts,
  		int minConsDaysWork, int maxConsDaysWork,
  		int minConsDaysOff, int maxConsDaysOff,
  		int maxTotalWeekends, int isCompleteWeekends) :
  			name_(name), minTotalShifts_(minTotalShifts), maxTotalShifts_(maxTotalShifts),
  			minConsDaysWork_(minConsDaysWork), maxConsDaysWork_(maxConsDaysWork),
  			minConsDaysOff_(minConsDaysOff), maxConsDaysOff_(maxConsDaysOff),
  			maxTotalWeekends_(maxTotalWeekends), isCompleteWeekends_(isCompleteWeekends) {
  };

  // Display methods: toString + override operator<< (easier)
  //
  string toString();
  friend std::ostream& operator<< (std::ostream& outs, Contract obj) {return outs << obj.toString();}
};


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   H i s t o r y
//
//  Describes the current (or initial) state of a nurse at D-day
//
//-----------------------------------------------------------------------------
class State{

public:
	// Index of the day in the planning horizon
	// WARNING : THE FIRST DAY IS ALWAYS SUPPOSED TO BE A MONDAY !!!!!!!!!!!!!
	//           If it may not be the case, the code should be modified, namely when counting the weekends worked
	//
	int dayId_;

	// Total nummber of days and weekends worked
	//
	int totalDaysWorked_, totalWeekendsWorked_;

	// number of consecutive days worked ending at D, and of consecutive days worked on the same shift ending at D (including RESTSHIFT = 0)
	// and shift worked on D-Day.
	//
	int consDaysWorked_, consShifts_, consDaysOff_;

	// Type of shift worked on D-Day. It can be a rest shift (=0).
	//
	int shift_;

public:
	// Constructor and Destructor
	State() {}
	~State();

	// Constructor with attributes
	State(int dayId, int totalDaysWorked, int totalWeekendsWorked,
			int consDaysWorked, int consShifts, int consDaysOff, int shift) :
				dayId_(dayId), totalDaysWorked_(totalDaysWorked), totalWeekendsWorked_(totalWeekendsWorked),
				consDaysWorked_(consDaysWorked), consShifts_(consShifts),
				consDaysOff_(consDaysOff), shift_(shift){};

	// Function that appends a new day worked on a given shift to the previous ones
	//
	void addNewDay(int newShift);

	// Function that appends a new day worked on a given shift to an input state
	// to update this state
	//
	void addDayToState(const State& prevState, int newShift);



    // Display methods: toString + override operator<< (easier)
    //
    string toString();
    friend std::ostream& operator<< (std::ostream& outs, State obj) {return outs << obj.toString();}

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
class Preferences{

public:
	// Constructor and destructor
	Preferences();
	~Preferences();

	// Constructor with initialization to a given number of nurses
	Preferences(int nbNurses, int nbDays, int nbShifts);

public:
	// Number of nurses
	//
	int nbNurses_;

	// Number of days considered in that case
	//
	int nbDays_;

	// Total number of possible shifts
	//
	int nbShifts_;

	// For each nurse, maps the day to the set of shifts that he/she wants to have off
	//
	vector<map<int,std::set<int> > > wishesOff_;

public:

	// For a given day, and a given shift, adds it to the wish-list for OFF-SHIFT
	void addShiftOff(int nurse, int day, int shift);

	// Adds the whole day to the wish-list
	void addDayOff(int nurse, int day);

	// True if the nurses wants that shift off
	bool wantsTheShiftOff(int nurse, int day, int shift);

	// True if the nurses wants the whole day off
	bool wantsTheDayOff(int nurse, int day);

    // Display methods: toString + override operator<< (easier)
    //
    string toString();
    friend std::ostream& operator<< (std::ostream& outs, Preferences obj) {return outs << obj.toString();}
};


//-----------------------------------------------------------------------------
//
//	C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

class Demand {

public:

	// generic constructor and destructor
	Demand(int nbDays, int nbShifts, int nbSkills,
		vector3D minDemand, vector3D optDemand);
	~Demand();

// constant attributes of the demand
//
public:

	// name of the demand
	//
	std::string name_;

	// number of days covered by the demand
	//
	const int nbDays_, nbShifts_, nbSkills_;

	// minimum and optimal demand for each day, shift and skill
	//
	const vector3D minDemand_;
	const vector3D optDemand_;

// preprocessed attributes aggregating the information of the demand
//
public:
	// total demand in the minimal and optimal demands
	//
	int minTotal_, optTotal_;

	// total demand per skill in the minimal and optimal demands
	//
	vector<int> minPerSkill_, optPerSkill_;

	// total demand per shift in the minimal and optimal demands
	//
	vector<int> minPerShift_, optPerShift_;

	// total demand per day in the minimal and optimal demands
	//
	vector<int> minPerDay_, optPerDay_;

	// highest demands per skill over the considered period
	//
	vector<int> minHighestPerSkill_, optHighestPerSkill_;

public:

	// compute all the potentially helpful attributes of a demand
	// this includes the total demand per skill, per shift,
	void preprocessDemand();

	// write the preprocessed information in the input stream
	//
	void displayPreprocess(Tools::LogOutput* outs);

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
	// Note : need both with const Contract and (non-const) Contract because non-const is used in our code,
	//        and const is needed so that we can override the operator= and have vector<Nurse>. We need to
	//        override it because vector members should have some properties (assignable a.o., which implies
	//        non-const)
	//
	Nurse(int id, string name, int nbSkills, vector<int> skills, Contract* contract) :
				id_(id), name_(name), nbSkills_(nbSkills), skills_(skills), pContract_(contract){}
	Nurse(int id, string name, int nbSkills, vector<int> skills, const Contract* contract) :
				id_(id), name_(name), nbSkills_(nbSkills), skills_(skills), pContract_(contract){}
	~Nurse();


	// the constant attibutes of the nurses are public
public:


	//-----------------------------------------------------------------------------
	// Constant characteristics of the nurses (no set method)
	//-----------------------------------------------------------------------------
	// Id of the nurse (=entry in the vector<Nurse> theNurse of the Scenario)
	//
	const int id_;

	// name of the nurse
	//
	const std::string name_;

	// number of skills and vector of the skills indices
	//
	const int nbSkills_;
	const vector<int> skills_;

	// Her contract type
	//
	const Contract* pContract_;

	// soft constraints of the nurse: min and max numbers of total assignments,
	// min and max consecutive working days, min and max consectuve days off,
	// maximum number of working week-ends and presence of absence of the
	// complete week end constraints
	//
	int minTotalShifts_, maxTotalShifts_;
	int minConsDaysWork_, maxConsDaysWork_;
	int minConsDaysOff_, maxConsDaysOff_;
	int maxTotalWeekEnds_;
	int isCompleteWeekEnds_;

    // Display methods: toString + override operator<< (easier)
    //
    string toString();
    friend std::ostream& operator<< (std::ostream& outs, Nurse obj) {return outs << obj.toString();}

    // Assignment (requested to build a vector<Nurse>)
    //
    Nurse& operator=(const Nurse& n);



};


















#endif /* defined(__ATCSolver__CftSolver__) */
