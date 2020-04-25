//
//  Scenario.h
//  RosterDesNurses
//

#ifndef __Scenario__
#define __Scenario__


#include "tools/MyTools.h"
#include "Demand.h"

// the penalties for violating the soft constraints on the nurses' schedules
// are in the problem definition
// they are set as static constant values in case they need to be shared with
// other classes (e.g. solvers)
//

static const int NB_LEVEL = 4;         //  0: weak 1: moderate; 2:strong; 3:compulsory;
static const int WEIGHT_OPTIMAL_DEMAND    = 30;
static const int WEIGHT_CONS_SHIFTS       = 15;
static const int WEIGHT_CONS_DAYS_WORK    = 30;
static const int WEIGHT_CONS_DAYS_OFF     = 30;
static const int WEIGHT_PREFERENCES_OFF[NB_LEVEL]       = {10, 20, 50, 1000};
static const int WEIGHT_PREFERENCES_ON [NB_LEVEL]       = {-10, -20, -50, 1000};
static const int WEIGHT_COMPLETE_WEEKEND  = 30;
static const int WEIGHT_TOTAL_SHIFTS      = 20;
static const int WEIGHT_TOTAL_WEEKENDS    = 30;


class Scenario;
class Nurse;
class Contract;
class Position;
class Preferences;

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
   std::string toString();
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
	Scenario(std::string name, int nbWeeks,
		 int nbSkills, std::vector<std::string> intToSkill, std::map<std::string,int> skillToInt,
		 int nbShifts, std::vector<std::string> intToShift, std::map<std::string,int> shiftToInt,
		 std::vector<int> hoursToWork, std::vector<int> shiftIDToShiftTypeID,
		 int nbShiftsType, std::vector<std::string> intToShiftType, std::map<std::string,int> shiftTypeToInt,
		 std::vector<std::vector<int> > shiftTypeIDToShiftID, std::vector<int> minConsShiftsType, std::vector<int> maxConsShiftsType,
		 std::vector<int> nbForbiddenSuccessors, vector2D<int> forbiddenSuccessors,
		 int nbContracts, std::vector<std::string> intToContract, std::map<std::string,Contract*> contracts,
		 int nbNurses, std::vector<Nurse>& theNurses, std::map<std::string,int> nurseNameToInt);

	// Hybrid copy constructor : this is only called when constructing a new scenario that copies most parameters
	// from the input scenario but for only a subgroup of nurses
	//
	Scenario(Scenario* pScenario,  std::vector<Nurse>& theNurses, Demand* pDemand, Preferences* pWeekPreferences);

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

	// number of skills, a std::map and a std::vector matching the name of each skill to an
	// index and reversely
	//
	const int nbSkills_;
	const std::vector<std::string> intToSkill_;
	const std::map<std::string,int> skillToInt_;

	// number of shifts, a std::map and a std::vector matching the name of each shift to an
	// index and reversely

	const int nbShifts_;
	const std::vector<std::string> intToShift_;
	const std::map<std::string,int> shiftToInt_;
  const std::vector<int> timeDurationToWork_, shiftIDToShiftTypeID_;

  bool isRestShift(int shift) const { return shiftIDToShiftTypeID_[shift] == 0; }
  bool isWorkShift(int shift) const { return !isRestShift(shift); }

	// number of typeshifts, a std::map and a std::vector matching the name of each type shift to an
	// index and reversely
	// minimum and maximum number consecutive assignments for each shift,
	// and penalty for violations of these bounds
	//
	const int nbShiftsType_;
	const std::vector<std::string> intToShiftType_;
	const std::map<std::string,int> shiftTypeToInt_;
	const vector2D<int> shiftTypeIDToShiftID_;

	// std::vector of possible contract types
	//
	const int nbContracts_;
	const std::vector<std::string> intToContract_;
	const std::map<std::string, Contract*> contracts_;

	// number of nurses, and std::vector of all the nurses
	//
	const int nbNurses_;
	const std::vector<Nurse> theNurses_;
	std::map<std::string,int> nurseNameToInt_;


private:

        // pour forcer l'usage des fonctions (en passant par le type du shift)

        const std::vector<int> minConsShiftType_, maxConsShiftType_;


	// for each shift, the number of forbidden successors and a table containing
	// the indices of these forbidden successors
	//
	const std::vector<int> nbForbiddenSuccessors_;
	const vector2D<int> forbiddenSuccessors_;

	//------------------------------------------------
	// From the Week data file
	//------------------------------------------------
	// Name of the week
	std::string weekName_;
	// Current week demand for each DAY, SHIFT, and SKILL
	//
	Demand* pWeekDemand_ = nullptr;

	// Shift off requests : Preferences for each nurse : which (day,shift) do they want off ?
	//
	int nbShiftOffRequests_;
	int nbShiftOnRequests_;
	Preferences* pWeekPreferences_ = nullptr;
	//------------------------------------------------


	//------------------------------------------------
	// From the History data file
	//------------------------------------------------
	// Initial historical state of the nurses
	//
	std::vector<State> initialState_;
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
	// std::vector of existing positions
	//
	int nbPositions_;
	std::vector<Position*> pPositions_;
	vector2D<Nurse> nursesPerPosition_;
	vector2D<Position*> componentsOfConnexPositions_;
	vector2D<Nurse> nursesPerConnexComponentOfPositions_;


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
	std::string weekName() {return weekName_;}
	Demand* pWeekDemand() {return pWeekDemand_;}
	int nbShifts() {return nbShifts_;}
	int nbShiftOffRequests() {return nbShiftOffRequests_;}
	int nbShiftOnRequests() {return nbShiftOnRequests_;}
	Preferences* pWeekPreferences() {return pWeekPreferences_;}
	std::vector<State>* pInitialState() {return &initialState_;}
	int nbSkills() {return nbSkills_;}
	int nbPositions() {return nbPositions_;}
	const std::vector<Position*>& pPositions() const {return pPositions_;}
	Position* pPosition(int p) const {return pPositions_[p];}
	int nbNurses() {return nbNurses_;}
	int nbOfConnexComponentsOfPositions() {return componentsOfConnexPositions_.size();}
	const std::vector<Position*>& componentOfConnexPositions(int c) const {return componentsOfConnexPositions_[c];}
	const std::vector<Nurse>& nursesInConnexComponentOfPositions(int c) const {return nursesPerConnexComponentOfPositions_[c];}

  int nbForbiddenSuccessorsShift(int shift) {
    int  shiftType = shiftIDToShiftTypeID_[shift];
    return nbForbiddenSuccessors_[shiftType];
  }

  int nbForbiddenSuccessorsShiftType(int shiftType) {
    return nbForbiddenSuccessors_[shiftType];
  }

  const std::vector<int>& nbForbiddenSuccessors() {
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
	const int minTotalShiftsOf(int whichNurse) const;
	int maxTotalShiftsOf(int whichNurse) const;
	int minConsDaysWorkOf(int whichNurse) const;
	int maxConsDaysWorkOf(int whichNurse) const;
	int minConsDaysOffOf(int whichNurse) const;
	int maxConsDaysOffOf(int whichNurse) const;
	int maxTotalWeekendsOf(int whichNurse) const;
	bool isCompleteWeekendsOf(int whichNurse) const;

  // getters for consecutive type of shifts

  int minConsShiftsOfTypeOf(int whichShift);
  int maxConsShiftsOfTypeOf(int whichShift);

  int minConsShiftsOf(int whichShiftType);
  int maxConsShiftsOf(int whichShiftType);

  // Cost function for consecutive identical shifts
  //
  double consShiftCost(int sh, int n){
    if(minConsShiftsOfTypeOf(sh) - n > 0) return (WEIGHT_CONS_SHIFTS * ( minConsShiftsOfTypeOf(sh) - n ) );
    if(n - maxConsShiftsOfTypeOf(sh) > 0) return (WEIGHT_CONS_SHIFTS * ( n - maxConsShiftsOfTypeOf(sh) ) );
    return 0;
  }

  double consShiftTypeCost(int sh, int n){
    if(minConsShiftsOf(sh) - n > 0) return (WEIGHT_CONS_SHIFTS * ( minConsShiftsOf(sh) - n ) );
    if(n - maxConsShiftsOf(sh) > 0) return (WEIGHT_CONS_SHIFTS * ( n - maxConsShiftsOf(sh) ) );
    return 0;
  }

	// getters for the attribute of the demand
	//
	int firstDay() {return pWeekDemand_->firstDay_;}
	int nbDays() {return pWeekDemand_->nbDays_;}

	// Setters to class attributes

	// when reading the week file (Demand and preferences)
	//
	inline void setWeekName(std::string weekName){ weekName_ = weekName;}
	inline void setWeekDemand(Demand* pDemand) {
    delete pWeekDemand_;
    pWeekDemand_ = pDemand;
  }
	inline void setTNbShiftOffRequests(int nbShiftOffRequests){ nbShiftOffRequests_ = nbShiftOffRequests; }
	inline void setTNbShiftOnRequests(int nbShiftOnRequests){ nbShiftOnRequests_ = nbShiftOnRequests; }
	void setWeekPreferences(Preferences* weekPreferences);

	// when reading the history file
	//
	inline void setThisWeek(int thisWeek){ thisWeek_ = thisWeek; }
   inline void addAWeek(){ ++nbWeeksLoaded_; }
	inline void setInitialState(std::vector<State> initialState){ initialState_ = initialState;}

	// return true if the shift shNext is a forbidden successor of shLast
	//
	bool isForbiddenSuccessorShift_Shift(int shNext, int shLast);
	bool isForbiddenSuccessorShift_ShiftType(int shNext, int shTypeLast);
	bool isForbiddenSuccessorShiftType_Shift(int shTypeNext, int shLast);
	bool isForbiddenSuccessorShiftType_ShiftType(int shTypeNext, int shTypeLast);

	// update the scenario to treat a new week
	//
	void updateNewWeek(Demand* pDemand, Preferences* pPreferences, std::vector<State> &initialStates);

	// Link the scenario with the Demand and the Preferences
	//
	inline void linkWithDemand(Demand* pDemand){
	  delete pWeekDemand_;
	   weekName_ = pDemand->name_;
	   pWeekDemand_ = pDemand;
	}

	void linkWithPreferences(Preferences* pPreferences);

	//------------------------------------------------
	// Display functions
	//------------------------------------------------

	// display the whole scenario
	//
	std::string toString();

	//------------------------------------------------
	// Preprocess functions
	//------------------------------------------------

	// preprocess the nurses to get the types
	//
	void preprocessTheNurses();

	// compute the connex components of the positions rcspp
	// (one edge between two positions indicate that they share a skill)
	//
	void computeConnexPositions() ;
};

#endif /* defined(__ATCSolver__CftSolver__) */
