/*
 * SubProblem.h
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include "RCGraph.h"
#include "Solver.h"
#include "MasterProblem.h"

static const set<pair<int,int> > EMPTY_FORBIDDEN_LIST;
static const set<int> EMPTY_FORBIDDEN_STARTING_DATES_LIST;

// Parameters (called in the solve function)
//
struct SubproblemParam{

	SubproblemParam(){}
	SubproblemParam(int strategy, LiveNurse* pNurse){
		initSubprobemParam(strategy, pNurse);
	}
	~SubproblemParam(){};

	void initSubprobemParam(int strategy, LiveNurse * pNurse){
		maxRotationLength_ = pNurse->maxConsDaysWork();
		switch(strategy){

		// 0 -> [Legal only]
		//		short = day-0 and last-day,
		//		max   = CD_max
		//		sink  = one / last day
		//
		case 0: shortRotationsStrategy_=2; maxRotationLength_+=0; oneSinkNodePerLastDay_ = true; break;

		// 1 -> [Exhaustive search]
		//		short = true,
		// 		max   = LARGE
		//		sink  = one / last day
		//
		case 1:	shortRotationsStrategy_=1;	maxRotationLength_=pNurse->pStateIni_->consDaysWorked_+pNurse->nbDays_; oneSinkNodePerLastDay_ = true; break;

		// 2 -> [Short + not too long]
		//		short = true,
		//		max   = CD_max+3
		//		sink  = one / last day
		//
		case 2: shortRotationsStrategy_=1; maxRotationLength_+=3; oneSinkNodePerLastDay_ = true; break;

		// 3 -> [No short + barely above legal]
		//		short = false,
		//		max   = CD_max+1
		//		sink  = one / last day
		//
		case 3: shortRotationsStrategy_=2; maxRotationLength_+=1; oneSinkNodePerLastDay_ = true; break;

		// UNKNOWN STRATEGY
		default:
			std::cout << "# Unknown strategy for the subproblem (" << strategy << ")" << std::endl;
			break;
		}
	}

	// *** PARAMETERS ***

	// 0 -> no short rotations
	// 1 -> all short rotations
	// 2 -> day-0 and last-day short rotations only
	// 3 -> day-0 short rotations only
	// 4 -> last-day short rotations only
	int shortRotationsStrategy_ = 0;

	// maximal length for a rotation
	int maxRotationLength_ = 0;

	// true  -> one sink node per day
	// false -> one single sink node for the network
	bool oneSinkNodePerLastDay_ = false;

	// Getters for the class fields
	//
	int maxRotationLength(){ return maxRotationLength_; }
	int shortRotationsStrategy(){return shortRotationsStrategy_;}
	bool oneSinkNodePerLastDay(){return oneSinkNodePerLastDay_;}

	// Setters
	//
	void maxRotationLength(int value){maxRotationLength_ = value;}
	void shortRotationsStrategy(int value){shortRotationsStrategy_ = value;}
	void oneSinkNodePerLastDay(bool value){oneSinkNodePerLastDay_ = value;}


};


//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

class SubProblem {

public:

	SubProblem();
	virtual ~SubProblem();

	// Constructor that correctly sets the resource (time + bounds), but NOT THE COST
	//
  SubProblem(Scenario* scenario, int nbDays, const Contract* contract, vector<State>* pInitState);

	// Initialization function for all global variables (not those of the graph)
	//
	virtual void init(vector<State>* pInitState);

	// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
	//
	virtual bool solve(LiveNurse* nurse, DualCosts * costs, SubproblemParam param,
			set<pair<int,int> > forbiddenDayShifts = EMPTY_FORBIDDEN_LIST,
			set<int> forbiddenStartingDays = EMPTY_FORBIDDEN_STARTING_DATES_LIST, bool optimality = false,
			double redCostBound = 0);

	// Returns all rotations saved during the process of solving the SPPRC
	//
	inline vector< Rotation > getRotations(){ return theRotations_; }

	// Returns true if the corresponding shift has no maximum limit of consecutive worked days
	//
	inline bool isUnlimited(int sh){return isUnlimited_[sh];}

	virtual void build();

protected:



	//----------------------------------------------------------------
	//
	// Necessary information: Scenario, contract type, minimum number
	// of paths to return, reduced costs.
	//
	// Optional information: Max rotation length.
	//
	//----------------------------------------------------------------


	// Pointer to the scenario considered
	//
	Scenario * pScenario_;

	// Number of days of the scenario (usually a multiple of 7)
	//
	int nDays_;

	// Contract type
	//
	const Contract * pContract_;

	// (Minimum) number of paths to return to the MP
	//
	int nPathsMin_;

	// Current live nurse considered
	//
	LiveNurse * pLiveNurse_;

	// All costs from Master Problem
	//
	DualCosts * pCosts_;

	// Bound on the reduced cost: if greater than this, the rotation is not added
	//
	double maxReducedCostBound_;

	// Vector that contains a boolean for each shift. TRUE if the maximum consecutive number of these shifts is higher than the maximal rotation length (or number of days); false otherwise
	//
	vector<bool> isUnlimited_;

	// Maximum number of consecutive days already worked by a nurse before the beginning of that period
	//
	int maxOngoingDaysWorked_;

	//----------------------------------------------------------------
	//
	// Answers: rotations, number of paths found
	//
	//----------------------------------------------------------------

	// Saved Rotations
	//
	vector< Rotation > theRotations_;

	// Number of paths found
	//
	int nPaths_;

	// Number of rotations found (that match the bound condition) at that iteration
	//
	int nLongFound_;

	// Best reduced cost found
	//
	double bestReducedCost_;

	//----------------------------------------------------------------
	//
	// Solving options.
	//
	//----------------------------------------------------------------

	SubproblemParam param_;

	//-----------------------
	// THE BASE COSTS
	//-----------------------
	// WARNING : for those that never change, of no use also.

    // For each day k (<= nDays_ - CDMin), contains WEIGHT_COMPLETE_WEEKEND if [it is a Saturday (resp. Sunday) AND the contract requires complete weekends]; 0 otherwise.
	vector<double> startWeekendCosts_, endWeekendCosts_;
	// Costs due to preferences of the nurse: for each day k (<= nDays_ - CDMin), shift s, contains WEIGHT_PREFERENCES if (k,s) is a preference of the nurse; 0 otherwise.
	vector<vector <double> > preferencesCosts_;

	// Cost function for consecutive identical shifts
	double consShiftCost(int sh, int n);
	double consShiftTypeCost(int shType, int n);
	// Cost function for consecutive days
	double consDaysCost(int n);
	// Initializes the startWeekendCost vector
	void initStartWeekendCosts();

	int CDMin_;	// Minimum number of consecutive days worked for free
	int daysMin_;      // principal node network begins at this index-1;  1 if no ShortSucc, CDMin otherwis
  int nLabels_; // Number of labels to use
  int maxRotationLength_;  // MAXIMUM LENGTH OF A ROTATION (in consecutive worked days)


	//-----------------------
	// THE GRAPH
	//-----------------------
	RCGraph g_;

	// Creates all nodes of the graph (including resource window)
	void createNodes();
	// Initiate variables for the nodes structures (vectors, etc.)
	// nDays : length of the scenario
	void initNodesStructures();
	// Add a node to the principal network of the graph, for shift sh, day k, and number of consecutive similar shifts cons
	void addNodeToPrincipalNetwork(int sh, int k, int cons);


	// Data structures to get the arcs id from other data
  //	vector3D arcsFromSource_;		// Index: (shiftType, day, nCons) of destination
	vector4D arcsFromSource_;		// Index: (shiftType, day, nCons, shift) of destination
	// vector3D arcsShiftToNewShift_;	// Index: (shift1, shift2, day1)
	// vector3D arcsShiftToSameShift_;	// Index: (shift, day, nCons) of origin

	// vector2D arcsRepeatShift_;		// Index: (shift, day) of origin
	vector2D arcsPrincipalToRotsizein_;	// Index: (shiftType, day) of origin
	vector<int> arcsSinkDayToSink_;						// Index: (day) of the end of rotation

	  /////
    vector3D principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node id
    std::vector<int> maxvalConsByShift_;						// For each shift, number of levels that the subnetwork contains
    std::map<int,int> principalToShift_;						// For each node of the principal network, maps it ID to the shift it represents
    std::map<int,int> principalToDay_;						// For each node of the principal network, maps it ID to the day it represents
    std::map<int,int> principalToCons_;						// For each node of the principal network, maps it ID to the number of consecutive shifts it represents

    vector4D arcsShiftToNewShift_;		// Index: (shiftType1, shiftType2, day1, shift)
    vector4D arcsShiftToSameShift_;		// Index: (shiftType, day, nCons, shift) of origin
    vector3D arcsShiftToEndsequence_;	// Index: (shiftType, day, nCons) of origin
    vector3D arcsRepeatShift_;		// Index: (shiftType, day, shift) of origin

    // Nodes of the ROTATION_LENGTH subnetwork
    std::vector<int> rotationLengthEntrance_;				// For each day, entrance node to the ROTATION_LENGTH subnetwork
    std::vector<map<int,int> > rotationLengthNodes_;			// For each day, maps the length of the rotation to the corresponding check node
    std::map<int,int> rotationLengthNodesLAT_;				// For each rotation length node, the corresponding EAT
    std::vector<int> rotationLengthExit_; // For each day, exit node to the ROTATION_LENGTH subnetwork

    std::vector<map<int,int> > arcsRotsizeinToRotsizeDay_;	// Index: (day,size) of the rotation [destination]
    std::vector<map<int,int> > arcsRotsizeToRotsizeoutDay_;	// Index: (day,size) of the rotation [origin]
    /////

	// Creates all arcs of the graph
	void createArcs();
	// Initiate variables for the arcs structures (integers, vectors, etc.)
	void initArcsStructures();
	// Create the specific types of arcs
	virtual void createArcsSourceToPrincipal();
	void createArcsPrincipalToPrincipal();
	void createArcsAllRotationSize();

  // Updates the costs depending on the reduced costs given for the nurse
  virtual void updateArcCosts();
    virtual bool feasibleArcSource(int k, int n, int shiftID);
    virtual double costArcSource(int a, int k, int shiftID);
    virtual double costArcPrincipal(int a, int k, int shiftID);
    virtual double costArcEnd(int a, int k);


	//----------------------------------------------------------------
	//
	// Update of the costs / network for solve function
	//
	//----------------------------------------------------------------

  double costOfVeryShortRotation(int firstDay, vector<int> succ);

	// FUNCTIONS -- SOLVE
	//

	virtual bool solveLongRotations(bool optimality);

	// Function called when optimal=true in the arguments of solve
	virtual bool solveLongRotationsOptimal();
	virtual bool solveLongRotationsHeuristic();

	// Initializes some cost vectors that depend on the nurse
	virtual void initStructuresForSolve();
	// Resets all solutions data (rotations, number of solutions, etc.)
	void resetSolutions();
	// Returns the rotation made from the given path
	virtual Rotation buildColumn(const RCSolution& sol);
	// Adds a single rotation to the list of solutions
	void addSingleRotationToListOfSolution();

	void updatedMaxRotationLengthOnNodes();

	// FORBIDDEN ARCS AND NODES
	//
	vector< vector<bool> > dayShiftStatus_;
	vector<bool> startingDayStatus_;

	// Returns true if the succession succ starting on day k does not violate any forbidden day-shift
	bool canSuccStartHere(vector<int> succ, int firstDay);
	// Forbids some days / shifts
	void forbid(set<pair<int,int> > forbiddenDayShifts);
	// Authorizes some days / shifts
	void authorize(set<pair<int,int> > forbiddenDayShifts);
	// Forbids some starting days
	void forbidStartingDays(set<int> forbiddenStartingDays);
	// Authorizes some starting days
	void authorizeStartingDays(set<int> forbiddenStartingDays);
	// Know if node
	inline bool isDayShiftForbidden(int k, int s){return ! dayShiftStatus_[k][s];}
	inline bool isStartingDayforbidden(int k){return startingDayStatus_[k];}
	// Forbid a node / arc
	void forbidDayShift(int k, int s);
	void forbidStartingDay(int k);
	// Authorize a node / arc
	void authorizeDayShift(int k, int s);
	void authorizeStartingDay(int k);
	void resetAuthorizations();

	// Given an arc, returns the normal travel time (i.e. travel time when authorized)
	// Test for random forbidden day-shift
	set< pair<int,int> > randomForbiddenShifts(int nbForbidden);


public:

	// Some getters
	//
	inline int nPaths(){return nPaths_;}
	inline int nLongFound(){return nLongFound_;}

	// Print functions.
	//
	void printRotation(Rotation rot);
	void printAllRotations();
	void printForbiddenDayShift();
	void printContractAndPrefenrences();


	// DBG
	Tools::Timer* timeInS_;
	Tools::Timer* timeInNL_;
};




class SubProblemShort : public SubProblem {

public:

	SubProblemShort();
	virtual ~SubProblemShort();

	// Constructor that correctly sets the resource (time + bounds), but NOT THE COST
	//
    SubProblemShort(Scenario* scenario, int nbDays, const Contract* contract, vector<State>* pInitState);

	// Initialization function for all global variables (not those of the graph)
	//
	void init(vector<State>* pInitState);

	// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
	//
	bool solve(LiveNurse* nurse, DualCosts * costs, SubproblemParam param,
			set<pair<int,int> > forbiddenDayShifts = EMPTY_FORBIDDEN_LIST,
			set<int> forbiddenStartingDays = EMPTY_FORBIDDEN_STARTING_DATES_LIST, bool optimality = false,
			double redCostBound = 0);


protected:

	// Number of rotations found (that match the bound condition) at that iteration
	//
	int nVeryShortFound_;

	//-----------------------
	// THE BASE COSTS
	//-----------------------

	// SHORT SUCCESSIONS (computed when creating them)
	vector< vector<double> > baseArcCostOfShortSucc_;										// For each size c \in [0,CDMin], for each short rotation of size c, contains its base cost (independent from the date)
        
	//-----------------------
	// THE SHORT SUCCESSIONS
	//-----------------------

	// SHORTSUCC -> OBJECTS
	//
	// Short successions (no starting date) -> those of all length
	vector3D allowedShortSuccBySize_;														// For each size c \in [0,CDMin], contains all allowed short successions of that size (satisfies succession constraints)
	vector2D lastShiftOfShortSucc_;															// For each size c \in [0,CDMin], for each short rotation of size c, contains the corresponding last shift performed
	vector2D nLastShiftOfShortSucc_;														// For each size c \in [0,CDMin], for each short rotation of size c, contains the number of consecutive days the last shift has been performed
	// Objects for short successions of maximal size CDMin
  //	int CDMin_;																				// Minimum number of consecutive days worked for free

	vector3D allShortSuccCDMinByLastShiftCons_;												// For each shift s, for each number of days n, contains the list of short successions of size CDMin ending with n consecutive days of shift s
	inline vector<int> shortSuccCDMin(int id){return allowedShortSuccBySize_[CDMin_][id];}	// Returns the short succession of size CDMin from its ID

	// SHORTSUCC -> FUNCTIONS
	//
	// Initializes all short successions, base costs, and corresponding vectors. Should only be called ONCE.
	void initShortSuccessions();

    // Returns the rotation made from the given path
    Rotation rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_res_cont resource);

    int   computeHoursInRotation(int sh, int k, int n);
    int   computeHoursInRotation(int id);

	//----------------------------------------------------------------
	//
	// Cost computation of the "very" short rotations (< CD_min)
	//
	//----------------------------------------------------------------
	bool priceVeryShortRotationsFirstDay();
	bool priceVeryShortRotationsLastDay();
	bool priceVeryShortRotations();
  //	double costOfVeryShortRotation(int firstDay, vector<int> succ);


	//----------------------------------------------------------------
	//
	// Update of the costs / network for solve function
	//
	//----------------------------------------------------------------

	// FUNCTIONS -- SOLVE
	//

	bool solveShortRotations();


	// Initializes some cost vectors that depend on the nurse
	void initStructuresForSolve();

	void createArcsSourceToPrincipal();

	// DATA -- COSTS
	//
	// Data structures that associates an arc to the chosen short succession of lowest cost
	map<int,int> shortSuccCDMinIdFromArc_;						// Maps the arcs to the corresponding short rotation ID
	vector3D idBestShortSuccCDMin_;								// For each day k (<= nDays_ - CDMin), shift s, number n, contains the best short succession of size CDMin that starts on day k, and ends with n consecutive days of shift s
	vector<vector<vector<double> > > arcCostBestShortSuccCDMin_;// For each day k (<= nDays_ - CDMin), shift s, number n, contains the cost of the corresponding arc

	// FUNCTIONS -- COSTS
	//
	// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
	void priceShortSucc();
	// Given a short succession and a start date, returns the cost of the corresponding arc
	double costArcShortSucc(int size, int id, int startDate);
	// Updates the costs depending on the reduced costs given for the nurse
	void updateArcCosts();


public:

	// // Some getters
	// //
	inline int nVeryShortFound(){return nVeryShortFound_;}

	// Print functions.
	//
	void printShortSucc();
	void printShortArcs();

};


#endif /* SUBPROBLEM_H_ */
