//
//  Scenario.h
//  RosterDesNurses
//

#ifndef __Scenario__
#define __Scenario__

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "tools/MyTools.h"
#include "data/Nurse.h"

using std::map;
using std::pair;
using std::string;
using std::vector;

// the penalties for violating the soft constraints on the nurses' schedules
// are in the problem definition
// they are set as static constant values in case they need to be shared with
// other classes (e.g. solvers)
//

static const int NB_LEVEL = 3;         //  0: strong; 1: moderate; 2: weak
static const int WEIGHT_OPTIMAL_DEMAND    = 30;
static const int WEIGHT_CONS_SHIFTS       = 15;
static const int WEIGHT_CONS_DAYS_WORK    = 30;
static const int WEIGHT_CONS_DAYS_OFF     = 30;
static const int WEIGHT_PREFERENCES_OFF[NB_LEVEL]       = {10, 10, 10};
static const int WEIGHT_PREFERENCES_ON [NB_LEVEL]       = {-50, -40, -30};
static const int WEIGHT_COMPLETE_WEEKEND  = 30;
static const int WEIGHT_TOTAL_SHIFTS      = 20;
static const int WEIGHT_TOTAL_WEEKENDS    = 30;


class Scenario;

//-----------------------------------------------------------------------------
//
//  C l a s s   S t a t e
//
//  Describes the current (or initial) state of a nurse at D-day
//
//-----------------------------------------------------------------------------
class State{

public:
   // Constructor and Destructor
   State():dayId_(0), totalTimeWorked_(0), totalWeekendsWorked_(0),
    consDaysWorked_(0), consShifts_(0), consDaysOff_(0), shiftType_(0) {}
   ~State();

   // Constructor with attributes
   State(int dayId, int totalTimeWorked, int totalWeekendsWorked,
	 int consDaysWorked, int consShifts, int consDaysOff, int shiftType, int shift) :
         dayId_(dayId), totalTimeWorked_(totalTimeWorked), totalWeekendsWorked_(totalWeekendsWorked),
         consDaysWorked_(consDaysWorked), consShifts_(consShifts),
         consDaysOff_(consDaysOff), shiftType_(shiftType), shift_(shift){};

   // Function that appends a new day worked on a given shiftType to the previous ones
   //
   // void addNewDay(int newShiftType);

   // Function that appends a new day worked on a given shiftType to an input state
   // to update this state
   //
   void addDayToState(const State& prevState, int newShiftType, int newShift, int timeWorked);


   // // Function that appends a new day worked on a given shift to an input state
   // // to update this state
   // //
   // void addDayToState(const State& prevState, int newShift, const Scenario* pScenario);


   // Display methods: toString + override operator<< (easier)
   //
   string toString();
   friend std::ostream& operator<< (std::ostream& outs, State obj) {return outs << obj.toString();}

public:
   // Index of the day in the planning horizon
   // WARNING : THE FIRST DAY IS ALWAYS SUPPOSED TO BE A MONDAY !!!!!!!!!!!!!
   //           If it may not be the case, the code should be modified, namely when counting the weekends worked
   //
   int dayId_;

   // Total nummber of days and weekends worked
   //
   int totalTimeWorked_, totalWeekendsWorked_;

   // number of consecutive days worked ending at D, and of consecutive days worked on the same shiftType ending at D (including RESTSHIFT = 0)
   // and shiftType worked on D-Day.
   //
   int consDaysWorked_, consShifts_, consDaysOff_;

   // Type of shift worked on D-Day. It can be a rest shift (=0).
   // A negative value -d means that the nurse has not been assigned a task for
   // the last d days
   //
   int shiftType_;

   int shift_;
};

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
		 vector<int> hoursToWork, vector<int> shiftIDToShiftTypeID,
		 int nbShiftsType, vector<string> intToShiftType, map<string,int> shiftTypeToInt,
		 vector<vector<int> > shiftTypeIDToShiftID, vector<int> minConsShiftsType, vector<int> maxConsShiftsType,
		 vector<int> nbForbiddenSuccessors, vector2D forbiddenSuccessors,
		 int nbContracts, vector<string> intToContract, map<string,Contract*> contracts,
		 int nbNurses, vector<Nurse>& theNurses, map<string,int> nurseNameToInt);

	// Hybrid copy constructor : this is only called when constructing a new scenario that copies most parameters
	// from the input scenario but for only a subgroup of nurses
	//
	Scenario(Scenario* pScenario,  vector<Nurse>& theNurses, Demand* pDemand, Preferences* pWeekPreferences);

	// copy constructor
	//
	Scenario(Scenario* pScenario);

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

	const int nbShifts_;
	const vector<string> intToShift_;
	const map<string,int> shiftToInt_;
        const vector<int> hoursToWork_, shiftIDToShiftTypeID_;

	// number of typeshifts, a map and a vector matching the name of each type shift to an
	// index and reversely
	// minimum and maximum number consecutive assignments for each shift,
	// and penalty for violations of these bounds
	//
	const int nbShiftsType_;
	const vector<string> intToShiftType_;
	const map<string,int> shiftTypeToInt_;
        const vector<vector<int> > shiftTypeIDToShiftID_;
	// const vector<int> minConsShifts_, maxConsShifts_;

	// Vector of possible contract types
	//
	const int nbContracts_;
	const vector<string> intToContract_;
	const map<string, Contract*> contracts_;

	// number of nurses, and vector of all the nurses
	//
	const int nbNurses_;
	const vector<Nurse> theNurses_;
	map<string,int> nurseNameToInt_;


private:

        // pour forcer l'usage des fonctions (en passant par le type du shift)
  
        const vector<int> minConsShiftType_, maxConsShiftType_;


	// for each shift, the number of forbidden successors and a table containing
	// the indices of these forbidden successors
	//
	const vector<int> nbForbiddenSuccessors_;
	const vector2D forbiddenSuccessors_;
  
	//------------------------------------------------
	// From the Week data file
	//------------------------------------------------
	// Name of the week
	string weekName_;
	// Current week demand for each DAY, SHIFT, and SKILL
	//
	Demand* pWeekDemand_;

	// Shift off requests : Preferences for each nurse : which (day,shift) do they want off ?
	//
	int nbShiftOffRequests_;
	int nbShiftOnRequests_;
	Preferences weekPreferences_;
	//------------------------------------------------


	//------------------------------------------------
	// From the History data file
	//------------------------------------------------
	// Initial historical state of the nurses
	//
	vector<State> initialState_;
	// range of the weeks that are being scheduled
	//
	int thisWeek_;
	int nbWeeksLoaded_;
	//------------------------------------------------


	//------------------------------------------------
	// From the custom file
	//------------------------------------------------
	//------------------------------------------------

	//------------------------------------------------
	// From the preprocessing of the nurses
	//------------------------------------------------
	// Vector of existing positions
	//
	int nbPositions_;
	vector<Position*> pPositions_;
	vector<vector<Nurse> > nursesPerPosition_;
	vector< vector<Position*> > componentsOfConnexPositions_;
	vector<vector<Nurse> > nursesPerConnexComponentOfPositions_;


	//------------------------------------------------

public:

	//------------------------------------------------
	// Getters and setters
	//------------------------------------------------

	// getters for the private class attributes
	//
	int nbWeeks() {return nbWeeks_;}
	int thisWeek() {return thisWeek_;}
	int nbWeeksLoaded() {return nbWeeksLoaded_;}
	string weekName() {return weekName_;}
	Demand* pWeekDemand() {return pWeekDemand_;}
	int nbShifts() {return nbShifts_;}
	int nbShiftOffRequests() {return nbShiftOffRequests_;}
	int nbShiftOnRequests() {return nbShiftOnRequests_;}
	Preferences* pWeekPreferences() {return &weekPreferences_;}
	vector<State>* pInitialState() {return &initialState_;}
	int nbSkills() {return nbSkills_;}
	int nbPositions() {return nbPositions_;}
	vector<Position*> pPositions() {return pPositions_;}
	Position* pPosition(int p) {return pPositions_[p];}
	int nbNurses() {return nbNurses_;}
	int nbOfConnexComponentsOfPositions() {return componentsOfConnexPositions_.size();}
	vector<Position*> componentOfConnexPositions(int c) {return componentsOfConnexPositions_[c];}
	vector<Nurse>& nursesInConnexComponentOfPositions(int c) {return nursesPerConnexComponentOfPositions_[c];}

  int nbForbiddenSuccessorsShift(int shift) {
    int  shiftType = shiftIDToShiftTypeID_[shift];
    return nbForbiddenSuccessors_[shiftType];
  }
  
  int nbForbiddenSuccessorsShiftType(int shiftType) {
    return nbForbiddenSuccessors_[shiftType];
  }

  vector<int> nbForbiddenSuccessors() {
    return nbForbiddenSuccessors_;
  }
  
	// getter for the maximum number of consecutive worked days before the planning horizon
	//
	inline int maxConDaysWorkedInHistory(){
		int ANS = 0;
		for(auto p : initialState_){
			if((p.consDaysWorked_ > ANS) && (p.shiftType_ > 0))
				ANS = p.consDaysWorked_;
		}
		return ANS;
	}

	// getters for the attributes of the nurses
	//
	const int minTotalShiftsOf(int whichNurse) {
		return theNurses_[whichNurse].minTotalShifts();
	}
	int maxTotalShiftsOf(int whichNurse) {
		return theNurses_[whichNurse].maxTotalShifts();
	}
	int minConsDaysWorkOf(int whichNurse) {
		return theNurses_[whichNurse].minConsDaysWork();
	}
	int maxConsDaysWorkOf(int whichNurse) {
		return theNurses_[whichNurse].maxConsDaysWork();
	}
	int minConsDaysOffOf(int whichNurse) {
		return theNurses_[whichNurse].maxConsDaysOff();
	}
	int maxConsDaysOffOf(int whichNurse) {
		return theNurses_[whichNurse].maxConsDaysOff();
	}
	int maxTotalWeekendsOf(int whichNurse) {
		return theNurses_[whichNurse].maxTotalWeekends();
	}
	bool isCompleteWeekendsOf(int whichNurse) {
		return theNurses_[whichNurse].needCompleteWeekends();
	}

  // getters for consecutive type of shifts
  
  int minConsShiftsOfTypeOf(int whichShift);
  int maxConsShiftsOfTypeOf(int whichShift);
  
  int minConsShiftsOf(int whichShiftType);
  int maxConsShiftsOf(int whichShiftType);

	// getters for the attribute of the demand
	//
	int firstDay() {return pWeekDemand_->firstDay_;}
	int nbDays() {return pWeekDemand_->nbDays_;}

	// Setters to class attributes

	// when reading the week file (Demand and preferences)
	//
	inline void setWeekName(string weekName){ weekName_ = weekName;}
	inline void setWeekDemand(Demand* pDemand) {pWeekDemand_ = pDemand;}
	inline void setTNbShiftOffRequests(int nbShiftOffRequests){ nbShiftOffRequests_ = nbShiftOffRequests; }
	inline void setTNbShiftOnRequests(int nbShiftOnRequests){ nbShiftOnRequests_ = nbShiftOnRequests; }
	inline void setWeekPreferences(Preferences weekPreferences){ weekPreferences_ = weekPreferences; }

	// when reading the history file
	//
	inline void setThisWeek(int thisWeek){ thisWeek_ = thisWeek; }
   inline void addAWeek(){ ++nbWeeksLoaded_; }
	inline void setInitialState(vector<State> initialState){ initialState_ = initialState;}

	// return true if the shift shNext is a forbidden successor of shLast
	//
	bool isForbiddenSuccessorShift_Shift(int shNext, int shLast);
	bool isForbiddenSuccessorShift_ShiftType(int shNext, int shTypeLast);
	bool isForbiddenSuccessorShiftType_Shift(int shTypeNext, int shLast);
	bool isForbiddenSuccessorShiftType_ShiftType(int shTypeNext, int shTypeLast);

	// update the scenario to treat a new week
	//
	void updateNewWeek(Demand* pDemand, Preferences &preferences, vector<State> &initialStates);

	// Link the scenario with the Demand and the Preferences
	//
	inline void linkWithDemand(Demand* pDemand){
	   weekName_ = pDemand->name_;
	   pWeekDemand_ = pDemand;
	}

	inline void linkWithPreferences(Preferences preferences){
		int nbShiftOffRequests_ = 0;
		for (Nurse nurse:theNurses_) {
			nbShiftOffRequests_ += preferences.howManyShiftsOff(nurse.id_);
		}
		int nbShiftOnRequests_ = 0;
		for (Nurse nurse:theNurses_) {
			nbShiftOnRequests_ += preferences.howManyShiftsOn(nurse.id_);
		}
		weekPreferences_ = preferences;
	}

public:

	//------------------------------------------------
	// Display functions
	//------------------------------------------------

	// display the whole scenario
	//
	string toString();

	//------------------------------------------------
	// Preprocess functions
	//------------------------------------------------

	// preprocess the nurses to get the types
	//
	void preprocessTheNurses();

	// compute the connex components of the positions graph
	// (one edge between two positions indicate that they share a skill)
	//
	void computeConnexPositions() ;
};

#endif /* defined(__ATCSolver__CftSolver__) */
