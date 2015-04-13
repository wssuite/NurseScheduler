/*
 * SubProblem.h
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include "MyTools.h"
#include "Solver.h"
#include "MasterProblem.h"


#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>


static int MAX_COST = 99999;
static int MAX_TIME = 99999;

// Different node types and their names
//
enum NodeType {
	SOURCE_NODE, PRINCIPAL_NETWORK, ROTATION_LENGTH_ENTRANCE,
	ROTATION_LENGTH, SINK_NODE, NONE_NODE};
static const vector<string> nodeTypeName = {
		"SOURCE_NODE", "PPL_NETWORK", "ROTSIZE__IN",
		"ROTSIZE    ", "SINK_NODE  ", "NONE       "};


// Different arc types and their names
//
enum ArcType{
	SOURCE_TO_PRINCIPAL, SHIFT_TO_NEWSHIFT, SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE,
	REPEATSHIFT, PRINCIPAL_TO_ROTSIZE, ROTSIZEIN_TO_ROTSIZE, ROTSIZE_TO_SINK,
	NONE_ARC};
static const vector<string> arcTypeName = {
		"SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ",
		"REPEATSHIFT    ", "PPL_TO_ROTSIZE ", "ROTSZIN_TO_RTSZ", "ROTSIZE_TO_SINK",
		"NONE           "};

static const set<pair<int,int> > EMPTY_FORBIDDEN_LIST;


//////////////////////////////////////////////////////////////////////////
//
// The following structures are the properties of the nodes and arcs:
//   For the nodes: id, earliest arrival time, and latest arrival time
//   For the arcs : id, cost and time (time is the resource here)
//   For the graph: usual stuff + nodes and and special properties
//
//////////////////////////////////////////////////////////////////////////

// Nodes specific properties for RC
//
struct Vertex_Properties{

	// Constructor
	//
	Vertex_Properties( int n = 0, NodeType t = NONE_NODE, int e = 0, int l = 0 ) : num( n ), type( t ), eat( e ), lat( l ) {}

	// id
	//
	int num;

	// type
	//
	NodeType type;

	// earliest arrival time
	//
	int eat;

	// latest arrival time
	//
	int lat;
};

// Arcs specific properties for RC
//
struct Arc_Properties{

	// Constructor
	//
	Arc_Properties( int n = 0, ArcType ty = NONE_ARC, int c = 0, int t = 0 ) : num( n ), type(ty), cost( c ), time( t ) {}

	// id
	//
	int num;

	// type
	//
	ArcType type;

	// traversal cost
	//
	int cost;

	// traversal time
	//
	int time;
};

// Graph with RC generic structure
//
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, Vertex_Properties, Arc_Properties> Graph;

//////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The following functions / data structures are used for shortest path problem without resource constraints (cost only):
//   Resource container ("resource" = cost)
//   Comparison override --> in SubProblem.cpp
//   Cost extension function
//   Dominance function
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Cost Container
//
struct spp_no_rc_res_cont{

	// Constructor
	//
	spp_no_rc_res_cont( int c = 0 ) : cost( c ) {};

	// Assign
	spp_no_rc_res_cont& operator=( const spp_no_rc_res_cont& other ){
		if( this == &other ) return *this;
		this->~spp_no_rc_res_cont();
		new( this ) spp_no_rc_res_cont( other );
		return *this;
	}

	// Cost of the path
	//
	int cost;
};

// Cost extension
//
class ref_no_res_cont{
public:
	inline bool operator()(const Graph& g, spp_no_rc_res_cont& new_cont, const spp_no_rc_res_cont& old_cont, boost::graph_traits<Graph>::edge_descriptor ed ) const{
		new_cont.cost = old_cont.cost + g[ed].cost;
		return true;
	}
};

// Dominance function model
//
class dominance_no_res_cont{
public:
	inline bool operator()( const spp_no_rc_res_cont& res_cont_1, const spp_no_rc_res_cont& res_cont_2 ) const {
		// must be "<=" here!!!
		// must NOT be "<"!!!
		return res_cont_1.cost <= res_cont_2.cost;
		// this is not a contradiction to the documentation
		// the documentation says:
		// "A label $l_1$ dominates a label $l_2$ if and only if both are resident
		// at the same vertex, and if, for each resource, the resource consumption
		// of $l_1$ is less than or equal to the resource consumption of $l_2$,
		// and if there is at least one resource where $l_1$ has a lower resource
		// consumption than $l_2$."
		// one can think of a new label with a resource consumption equal to that
		// of an old label as being dominated by that old label, because the new
		// one will have a higher number and is created at a later point in time,
		// so one can implicitly use the number or the creation time as a resource
		// for tie-breaking
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// The following functions / data structures are used for the time resource:
//   Resource container ("resource" = cost + time)
//   Comparison override --> in SubProblem.cpp
//   Cost extension function
//   Dominance function
//
/////////////////////////////////////////////////////////////////////////////

// data structures for shortest path problem with time windows (spp_tw)
// ResourceContainer model

struct spp_spptw_res_cont{

	// Constructor
	//
	spp_spptw_res_cont( int c = 0, int t = 0 ) : cost( c ), time( t ) {}

	// Assign
	//
	spp_spptw_res_cont& operator=( const spp_spptw_res_cont& other ){
		if( this == &other ) return *this;
		this->~spp_spptw_res_cont();
		new( this ) spp_spptw_res_cont( other );
		return *this;
	}

	// Current cost
	//
	int cost;

	// Current time consumption
	//
	int time;
};

// Resources extension model (arc has cost + travel time)
class ref_spptw{
public:
	inline bool operator()( const Graph& g, spp_spptw_res_cont& new_cont,	const spp_spptw_res_cont& old_cont,	boost::graph_traits<Graph>::edge_descriptor ed ) const{
		const Arc_Properties& arc_prop = get( boost::edge_bundle, g )[ed];
		const Vertex_Properties& vert_prop = get( boost::vertex_bundle, g )[target( ed, g )];
		new_cont.cost = old_cont.cost + arc_prop.cost;
		int& i_time = new_cont.time;
		i_time = old_cont.time + arc_prop.time;
		i_time < vert_prop.eat ? i_time = vert_prop.eat : 0;
		return i_time <= vert_prop.lat ? true : false;
	}
};

// Dominance function model
class dominance_spptw{
public:
	inline bool operator()( const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2 ) const	{
		// must be "<=" here!!!
		// must NOT be "<"!!!
		return res_cont_1.cost <= res_cont_2.cost and res_cont_1.time <= res_cont_2.time;
		// this is not a contradiction to the documentation
		// the documentation says:
		// "A label $l_1$ dominates a label $l_2$ if and only if both are resident
		// at the same vertex, and if, for each resource, the resource consumption
		// of $l_1$ is less than or equal to the resource consumption of $l_2$,
		// and if there is at least one resource where $l_1$ has a lower resource
		// consumption than $l_2$."
		// one can think of a new label with a resource consumption equal to that
		// of an old label as being dominated by that old label, because the new
		// one will have a higher number and is created at a later point in time,
		// so one can implicitly use the number or the creation time as a resource
		// for tie-breaking
	}
};

/////////////////////////////////////////////////////////////////////////////



struct Costs{
public:

	Costs(){}
	Costs(vector< vector<double> > * workCosts, vector<double> * startWorkCosts, vector<double> * endWorkCosts, double workedWeekendCost):
		pWorkCosts_(workCosts), pStartWorkCosts_(startWorkCosts), pEndWorkCosts_(endWorkCosts), workedWeekendCost_(workedWeekendCost) {
		Tools::initDoubleVector2D(pPreferenceCosts_, 0, 0);
	}

	// CONSTRUCTOR NOT TO BE USED, ONLY FOR RANDOM GENERATED COSTS...
	Costs(vector< vector<double> > workCosts, vector<double> startWorkCosts, vector<double> endWorkCosts, double workedWeekendCost):
		workedWeekendCost_(workedWeekendCost){
		pWorkCosts_ = new vector<vector<double> > (workCosts);
		pStartWorkCosts_ = new vector<double> (startWorkCosts);
		pEndWorkCosts_ = new vector<double> (endWorkCosts);
		pPreferenceCosts_ = new vector<vector<double> > ();
	}

	// GETTERS
	//
	inline double dayShiftWorkCost(int day, int shift){return ((*pWorkCosts_)[day][shift]);}
	inline double startWorkCost(int day){return ((*pStartWorkCosts_)[day]);}
	inline double endWorkCost(int day){return ((*pEndWorkCosts_)[day]);}
	inline double workedWeekendCost(){return workedWeekendCost_;}


	// SETTERS
	//
	inline void setPreferenceCosts(vector< vector<double> > * v){pPreferenceCosts_ = v;}


protected:

	// Indexed by : (day, shift) !! 0 = shift 1 !!
    vector< vector<double> > * pWorkCosts_;

    // Indexed by : day
    vector<double> * pStartWorkCosts_;

    // Indexed by : day
    vector<double> * pEndWorkCosts_;

    // Reduced cost of the weekends
    double workedWeekendCost_;

	// For each day k (<= nDays_ - CDMin), shift s, contains WEIGHT_PREFERENCES if (k,s) is a preference of the nurse; 0 otherwise.
	//
	vector< vector<double> > * pPreferenceCosts_;

    /*
     * Possibility of adding other costs (such as baseCosts, not necessarily dual)
     */


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
	~SubProblem();

	// Constructor that correctly sets the resource (time + bounds), but NOT THE COST
	//
	SubProblem(Scenario* scenario, Demand * demand, const Contract* contract);

	// Initialization function for all global variables (not those of the graph)
	//
	void init();

	// Test function for Shortest Path Problem with Resource Constraint
	//
	void testGraph_spprc();

	// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
	//
	bool solve(LiveNurse* nurse, Costs * costs, set<pair<int,int> > forbiddenDayShifts = EMPTY_FORBIDDEN_LIST, bool optimality = false, int maxRotationLength=MAX_TIME);

	// Returns all rotations saved during the process of solving the SPPRC
	//
	inline vector< Rotation > getRotations(){return theRotations_;}

	// Returns true if the corresponding shift has no maximum limit of consecutive worked days
	//
	inline bool isUnlimited(int sh){return isUnlimited_[sh];}


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

	// Pointer to the demand
	Demand* pDemand_;

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
	Costs * pCosts_;

	// Maximum length of a rotation (in consecutive worked days)
	//
	int maxRotationLength_;

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



	//----------------------------------------------------------------
	//
	// Construction of the network.
	//
	// INDEPENDENT FROM ANY NURSE / REDUCED COST !!!
	//
	//----------------------------------------------------------------

	// THE GRAPH
	Graph g_;

	//-----------------------
	// THE BASE COSTS
	//-----------------------

	// SHORT SUCCESSIONS (computed when creating them)
	vector< vector<double> > baseArcCostOfShortSucc_;										// For each size c \in [0,CDMin], for each short rotation of size c, contains its base cost (independent from the date)

	// All arcs have a base cost
	// WARNING : for short ones, is of no use because must be priced first.
	// WARNING : for those that never change, of no use also.
	vector<double> arcBaseCost_;

    // For each day k (<= nDays_ - CDMin), contains WEIGHT_COMPLETE_WEEKEND if [it is a Saturday (resp. Sunday) AND the contract requires complete weekends]; 0 otherwise.
	vector<double> startWeekendCosts_, endWeekendCosts_;

	// Cost function for consecutive identical shifts
	double consShiftCost(int sh, int n);
	// Cost function for consecutive days
	double consDaysCost(int n);
	// Initializes the startWeekendCost vector
	void initStartWeekendCosts();


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
	int CDMin_;																				// Minimum number of consecutive days worked for free
	vector3D allShortSuccCDMinByLastShiftCons_;												// For each shift s, for each number of days n, contains the list of short successions of size CDMin ending with n consecutive days of shift s
	inline vector<int> shortSuccCDMin(int id){return allowedShortSuccBySize_[CDMin_][id];}	// Returns the short succession of size CDMin from its ID

	// SHORTSUCC -> FUNCTIONS
	//
	// Initializes all short successions, base costs, and corresponding vectors. Should only be called ONCE.
	void initShortSuccessions();




	//-----------------------
	// THE NODES
	//-----------------------

	// NODES -> OBJECTS
	//
	int nNodes_;										// Total number of nodes in the graph
	vector<NodeType> allNodesTypes_;					// vector of their types
	// Source
	int sourceNode_;
	// Nodes of the PRINCIPAL_NETWORK subnetwork
	vector3D principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node id
	vector<int> maxvalConsByShift_;						// For each shift, number of levels that the subnetwork contains
	map<int,int> principalToShift_;						// For each node of the principal network, maps it ID to the shift it represents
	map<int,int> principalToDay_;						// For each node of the principal network, maps it ID to the day it represents
	map<int,int> principalToCons_;						// For each node of the principal network, maps it ID to the number of consecutive shifts it represents
	// Nodes of the ROTATION_LENGTH subnetwork
	int rotationLengthEntrance_;						// Entrance node to the ROTATION_LENGTH subnetwork
	map<int,int> rotationLengthNodes_;					// Maps the length of the rotation to the corresponding check node
	// Sink Node
	int sinkNode_;

	// NODES -> FUNCTIONS
	//
	// Creates all nodes of the graph (including resource window)
	void createNodes();
	// Basic function for adding a node
	void addSingleNode(NodeType type, int eat, int lat);
	// Initiate variables for the nodes structures (vectors, etc.)
	// nDays : length of the scenario
	void initNodesStructures();
	// Add a node to the principal network of the graph, for shift sh, day k, and number of consecutive similar shifts cons
	void addNodeToPrincipalNetwork(int sh, int k, int cons);
	// Get info from the node ID
	inline NodeType nodeType(int v){return get( &Vertex_Properties::type, g_)[v];}
	inline int nodeEat(int v){return get( &Vertex_Properties::eat, g_)[v];}
	inline int nodeLat(int v){return get( &Vertex_Properties::lat, g_)[v];}




	//-----------------------
	// THE ARCS
	//-----------------------

	// ARCS -> OBJECTS
	//
	int nArcs_;											// Total number of arcs in the graph
	vector<ArcType> allArcsTypes_;						// Vector of their types
	vector< boost::graph_traits< Graph>::edge_descriptor > arcsDescriptors_;
	// Data structures to get the arcs id from other data
	vector3D arcsFromSource_;							// Index: (shift, day, nCons) of destination
	vector3D arcsShiftToNewShift_;						// Index: (shift1, shift2, day1)
	vector3D arcsShiftToSameShift_;						// Index: (shift, day, nCons) of origin
	vector3D arcsShiftToEndsequence_;					// Index: (shift, day, nCons) of origin
	vector2D arcsRepeatShift_;							// Index: (shift, day) of origin
	vector2D arcsPrincipalToRotsizein_;					// Index: (shift, day) of origin
	map<int,int> arcsRotsizeinToRotsize_;				// Index: (size) of the rotation [destination]
	map<int,int> arcsRotsizeToRotsizeout_;				// Index: (size) of the rotation [origin]

	// ARCS -> FUNCTIONS
	//
	// Creates all arcs of the graph
	void createArcs();
	// Basic function for adding an arc
	void addSingleArc(int origin, int destination, int baseCost, int travelTime, ArcType type);
	// Initiate variables for the arcs structures (integers, vectors, etc.)
	void initArcsStructures();
	// Create the specific types of arcs
	void createArcsSourceToPrincipal();
	void createArcsPrincipalToPrincipal();
	void createArcsAllRotationSize();

	// Initialize and update the costs of the arcs.
	// Some arcs always have the same cost (0). Hence, their cost may not be changed
	void initBaseCostArcs();
	void updateCostArcs();

	// Get info with the arc ID
	inline ArcType arcType(int a) {return allArcsTypes_[a];}
	inline int arcOrigin(int a) {return source(arcsDescriptors_[a], g_);}
	inline int arcDestination(int a) {return target(arcsDescriptors_[a], g_);}
	inline int arcLength(int a) {return get( &Arc_Properties::time, g_, arcsDescriptors_[a]);}
	inline int arcCost(int a) {return get( &Arc_Properties::cost, g_, arcsDescriptors_[a]);}







	//----------------------------------------------------------------
	//
	// Update of the costs / network for solve function
	//
	// INDEPENDENT FROM ANY NURSE / REDUCED COST !!!
	//
	//----------------------------------------------------------------

	// FUNCTIONS -- SOLVE
	//
	// Initializes some cost vectors that depend on the nurse
	void initStructuresForSolve(LiveNurse* nurse, Costs * costs, set<pair<int,int> > forbiddenDayShifts, int maxRotationLength);
	// Resets all solutions data (rotations, number of solutions, etc.)
	void resetSolutions();
	// Transforms the solutions found into proper rotations. Returns true if at least one has been added (can be false if none was found or if all were nonnegative and negativeOnly=true)
	bool addRotationsFromPaths(vector< vector< boost::graph_traits<Graph>::edge_descriptor > > paths, vector<spp_spptw_res_cont> resources, bool negativeOnly);
	// Returns the rotation made from the given path
	Rotation rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource);

	// DATA -- COSTS
	//
	// Data structures that associates an arc to the chosen short succession of lowest cost
	map<int,int> shortSuccCDMinIdFromArc_;						// Maps the arcs to the corresponding short rotation ID
	vector3D idBestShortSuccCDMin_;								// For each day k (<= nDays_ - CDMin), shift s, number n, contains the best short succession of size CDMin that starts on day k, and ends with n consecutive days of shift s
	vector<vector <double> > preferencesCosts_;					// Costs due to preferences of the nurse
	vector<vector<vector<double> > > arcCostBestShortSuccCDMin_;// For each day k (<= nDays_ - CDMin), shift s, number n, contains the cost of the corresponding arc

	// FUNCTIONS -- COSTS
	//
	// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
	void priceShortSucc();
	// Given a short succession and a start date, returns the cost of the corresponding arc
	double costArcShortSucc(int size, int id, int startDate);
	// Single cost/time change
	inline void updateCost(int arc, double cost){boost::put( &Arc_Properties::cost, g_, arcsDescriptors_[arc], cost );}
	// Updates the costs depending on the reduced costs given for the nurse
	void updateArcCosts();
	// For tests, must be able to randomly generate costs
	void generateRandomCosts(double minVal, double maxVal);

	// DATA -- FORBIDDEN ARCS AND NODES
	//
	vector< vector<bool> > dayShiftStatus_;
	vector<bool> arcStatus_;
	vector<bool> nodeStatus_;

	// FUNCTIONS -- FORBIDDEN ARCS AND NODES
	//
	// Returns true if the succession succ starting on day k does not violate any forbidden day-shift
	bool canSuccStartHere(vector<int> succ, int firstDay);
	// Forbids some days / shifts
	void forbid(set<pair<int,int> > forbiddenDayShifts);
	// Know if node / arc is forbidden
	inline bool isArcForbidden(int a){return ! arcStatus_[a];}
	inline bool isNodeForbidden(int v){return ! nodeStatus_[v];}
	inline bool isDayShiftForbidden(int k, int s){return ! dayShiftStatus_[k][s];}
	// Forbid a node / arc
	void forbidArc(int a);
	void forbidNode(int v);
	void forbidDayShift(int k, int s);
	// Authorize a node / arc
	void authorizeArc(int a);
	void authorizeNode(int v);
	void authorizeDayShift(int k, int s);
	void resetAuthorizations();
	// Updates the travel time of an arc / node
	inline void updateTime(int a, double time){boost::put( &Arc_Properties::time, g_, arcsDescriptors_[a], time );}
	inline void updateLat(int v, int time){boost::put( &Vertex_Properties::lat, g_, v, time);}
	// Given an arc, returns the normal travel time (i.e. travel time when authorized)
	int normalTravelTime(int a);
	// Test for random forbidden day-shift
	set< pair<int,int> > randomForbiddenShifts(int nbForbidden);





	//----------------------------------------------------------------
	//
	// Utilities functions
	//
	//----------------------------------------------------------------
	bool succContainsDayShift(int size, int succId, int startDate, int thatDay, int thatShift);
	int mapAntecedent(map<int,int> m, int val);


public:

	// Print functions.
	//
	void printGraph();
	string printNode(int v);
	void printAllNodes();
	string printArc(int a);
	void printAllArcs();
	string shortNameNode(int v);
	string printSummaryOfGraph();
	void printShortSucc();
	void printPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont ressource);
	void printRotation(Rotation rot);
	void printAllRotations();
	void printForbiddenDayShift();
	void printShortArcs();


};


#endif /* SUBPROBLEM_H_ */
