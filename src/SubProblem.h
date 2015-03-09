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

// Different node types and their names
//
enum NodeType {
	SOURCE_NODE, SHORT_ROTATION, PRINCIPAL_NETWORK, ROTATION_LENGTH_ENTRANCE,
	ROTATION_LENGTH, ROTATION_LENGTH_EXIT, SINK_NODE, NONE_NODE};
static const vector<string> nodeTypeName = {
		"SOURCE_NODE", "SHORT_ROTAT", "PPL_NETWORK", "ROTSIZE__IN",
		"ROTSIZE    ", "ROTSIZE_OUT", "SINK_NODE  ", "NONE       "};


// Different arc types and their names
//
enum ArcType{
	SOURCE_TO_SHORT, SHORT_TO_SINK, SHORT_TO_PRINCIPAL,	SHIFT_TO_NEWSHIFT,
	SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE, REPEATSHIFT, PRINCIPAL_TO_ROTSIZE,
	ROTSIZEIN_TO_ROTSIZE, ROTSIZE_TO_ROTSIZEOUT, ROTSIZEOUT_TO_SINK, NONE_ARC};
static const vector<string> arcTypeName = {
		"SOURCE_TO_SHORT", "SHORT_TO_SINK  ", "SHORT_TO_PRPAL ", "SHIFT_TO_NEWSH ",
		"SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ", "REPEATSHIFT    ", "PPL_TO_ROTSIZE ",
		"ROTSZIN_TO_RTSZ", "ROTSZ_TO_RSZOUT", "RTSZOUT_TO_SINK", "NONE           "};

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
	SubProblem(Scenario* scenario, Demand * demand, Contract* contract);

	// Initialization function for all global variables (not those of the graph)
	//
	void init();

	// Test function for Shortest Path Problem with Resource Constraint
	//
	void testGraph_spprc();

	// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
	//
	bool solve(LiveNurse* nurse, vector<vector<double> > * dualCosts, set<pair<int,int> > forbiddenDayShifts = EMPTY_FORBIDDEN_LIST, bool optimality = false, int maxRotationLength=-1);

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
	Contract * pContract_;

	// (Minimum) number of paths to return to the MP
	//
	int nPathsMin_;

	// 1 cost / day / worked shift
	//
	vector<vector<double> > * costs_;

	// Maximum length of a rotation (in consecutive worked days)
	//
	int maxRotationLength_;

	// Vector that contains a boolean for each shift. TRUE if the maximum consecutive number of these shifts is higher than the maximal rotation length (or number of days); false otherwise
	//
	vector<bool> isUnlimited_;


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


	// LE GRAPHE
	Graph g_;



	//----------------------------------------------------------------
	//
	// Data and functions for the NODES of the network
	//
	//----------------------------------------------------------------

	//-----------------//
	// DATA STRUCTURES //
	//-----------------//

	// Total number of nodes in the graph (is also the id of the node to add if a new node is to be added) and the vector of their types
	//
	int nNodes_;
	vector<NodeType> allNodesTypes_;

	// Source
	//
	int sourceNode_;

	// Nodes of the SHORT_ROTATION subnetwork
	//
	vector< vector<Rotation*> > shortRotations_;		// List of all short rotations (contains their sequence of tasks)
	vector< vector<int> > shortRotationsNodes_;			// For each length (#days), the list of all nodes that correspond to short rotations of this length
	map<int,int> lastShiftOfShort_;						// For each short rotation, the id of the last shift worked
	map<int,int> nLastShiftOfShort_;					// The number of consecutive similar shifts that ends the short rotation
	map<int,int> lastDayOfShort_;
	map<int,Rotation> nodeToShortRotation_;				// Maps the node ID to the corresponding short rotation


	// Nodes of the PRINCIPAL_NETWORK subnetwork
	//
	vector3D principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node
	vector<int> maxvalConsByShift_;						// For each shift, number of levels that the subnetwork contains
	map<int,int> principalToShift_;						// For each node of the principal network, maps it ID to the shift it represents
	map<int,int> principalToDay_;						// For each node of the principal network, maps it ID to the day it represents
	map<int,int> principalToCons_;						// For each node of the principal network, maps it ID to the number of consecutive shifts it represents

	// Nodes of the ROTATION_LENGTH subnetwork
	int rotationLengthEntrance_;						// Entrance node to the ROTATION_LENGTH subnetwork
	map<int,int> rotationLengthNodes_;					// Maps the length of the rotation to the corresponding check node
	int rotationLengthExit_;							// Exit node to the ROTATION_LENGTH subnetwork

	// Sink Node
	//
	int sinkNode_;

	//-----------//
	// FUNCTIONS //
	//-----------//

	// Creates all nodes of the graph (including resource window)
	//
	void createNodes();

	// Basic function for adding a node
	//
	void addSingleNode(NodeType type, int eat, int lat);

	// Initiate variables for the nodes structures (vectors, etc.)
	// nDays : length of the scenario
	//
	void initNodesStructures();

	// Returns a vector3D: For each duration from 0 (empty) to CD_min, returns the list of LEGID shift successions.
	// They do not depend on the starting date, only allowed successions w.r.t. the forbidden successors.
	//
	vector3D allowedShortSuccessions();

	// Add a short rotation to the graph, that starts at k0 and contains the given succession of tasks of duration length
	//
	void addShortRotationToGraph(int k0, vector<int> shiftSuccession, int length);

	// Add a node to the principal network of the graph, for shift sh, day k, and number of consecutive similar shifts cons
	//
	void addNodeToPrincipalNetwork(int sh, int k, int cons);

	// Get info from the node ID
	//
	inline NodeType nodeType(int v){return get( &Vertex_Properties::type, g_)[v];}
	inline int nodeEat(int v){return get( &Vertex_Properties::eat, g_)[v];}
	inline int nodeLat(int v){return get( &Vertex_Properties::lat, g_)[v];}



	//----------------------------------------------------------------
	//
	// Data and functions for the ARCS of the network
	//
	//----------------------------------------------------------------


	//-----------------//
	// DATA STRUCTURES //
	//-----------------//

	// Total number of arcs in the graph (is also the id of the arc to add if a new node is to be added) and the vector of their types
	//
	int nArcs_;
	vector< boost::graph_traits< Graph>::edge_descriptor > arcsDescriptors_;
	vector<ArcType> allArcsTypes_;

	// Vector of all arcs whose origin is the source node
	//
	vector<int> arcsFromSource_;



	//-----------//
	// FUNCTIONS //
	//-----------//

	// Creates all arcs of the graph
	//
	void createArcs();

	// Basic function for adding an arc
	//
	void addSingleArc(int origin, int destination, int cost, int travelTime, ArcType type);

	// Initiate variables for the arcs structures (integers, vectors, etc.)
	//
	void initArcsStructures();

	// Create the specific types of arcs
	//
	void createArcsSourceToShort();
	void createArcsShortToSink();
	void createArcsShortToPrincipal();
	void createArcsPrincipalToPrincipal();
	void createArcsAllRotationSize();

	// Initialize and update the costs of the arcs.
	// Some arcs always have the same cost (0). Hence, their cost may not be changed
	//
	void initCostArcs();
	void updateCostArcs();

	// Get info with the arc ID
	inline ArcType arcType(int a) {return allArcsTypes_[a];}
	inline int arcOrigin(int a) {return source(arcsDescriptors_[a], g_);}
	inline int arcDestination(int a) {return target(arcsDescriptors_[a], g_);}
	inline int arcLength(int a) {return get( &Arc_Properties::time, g_, arcsDescriptors_[a]);}
	inline int arcCost(int a) {return get( &Arc_Properties::cost, g_, arcsDescriptors_[a]);}


public:

	// Print functions.
	//
	void printGraph();
	string printNode(int v);
	string printArc(int a);
	string shortNameNode(int v);
	string printSummaryOfGraph();


};


#endif /* SUBPROBLEM_H_ */
