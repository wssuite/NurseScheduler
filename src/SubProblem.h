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
enum NodeType {SOURCE_NODE, SHORT_ROTATION, PRINCIPAL_NETWORK, ROTATION_LENGTH_ENTRANCE, ROTATION_LENGTH, ROTATION_LENGTH_EXIT, SINK_NODE, NONE_NODE};
static const vector<string> nodeTypeName = {"SOURCE_NODE", "SHORT_ROTAT", "PPL_NETWORK", "ROTSIZE__IN", "ROTSIZE     ", "ROTSIZE_OUT", "SINK_NODE  ", "NONE       "};


// Different arc types and their names
//
enum ArcType{SOURCE_TO_SHORT, SHORT_TO_SINK, SHORT_TO_PRINCIPAL, PRINCIPAL_TO_SINK,
	SHIFT_TO_NEWSHIFT, SHIFT_TO_SAMESHIFT, SHIFT_TO_ENDSEQUENCE, REPEATSHIFT,
	ROTSIZEIN_TO_ROTSIZE, ROTSIZE_TO_ROTSIZEOUT, ROTSIZEOUT_TO_SINK, NONE_ARC};
static const vector<string> arcTypeName = {"SOURCE_TO_SHORT", "SHORT_TO_SINK  ", "SHORT_TO_PRPAL ", "PRPAL_TO_SINK  ",
		"SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ", "REPEATSHIFT    ",
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

	// Constructeur du graphe (seulement le graphe; pas de couts / bornes / temps de trajet)
	//
	SubProblem(Scenario* scenario, Contract* contract);

	// Fonction de test de l'algor de plus courts chemins
	//
	void testGraph_spprc();

	// Solve : Retourne true si des chemins de couts reduit negatif ont ete trouves; false sinon.
	//
	bool solve(LiveNurse* nurse, vector<vector<double> > * dualCosts, set<pair<int,int> > forbiddenDayShifts = EMPTY_FORBIDDEN_LIST, bool optimality = false);

	// Retourne l'ensemble des rotations mises en memoire par le SP
	//
	inline vector< Rotation > getRotations(){return theRotations_;}


protected:



	//----------------------------------------------------------------
	//
	// Necessary information: Scenario, contract type, minimum number
	// of paths to return, reduced costs.
	//
	// Optional information: Max rotation length.
	//
	//----------------------------------------------------------------


	// Pointeur vers le scenario a resoudre
	//
	Scenario * pScenario_;

	// Contract type
	//
	Contract * pContract_;

	// Nombre de chemins (minimum) a retourner au MP
	//
	int nPathsMin_;

	// 1 cout / jour / shift travaille
	//
	vector<vector<double> > * costs_;

	// Maximum length of a rotation (in consecutive worked days)
	//
	int maxRotationLength_;


	//----------------------------------------------------------------
	//
	// Answers: rotations, number of paths found
	//
	//----------------------------------------------------------------

	// Les Rotations sauvegardees
	//
	vector< Rotation > theRotations_;

	// Nombre de chemins trouves
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
	map<int,Rotation> nodeToShortRotation_;				// Maps the node ID to the corresponding short rotation


	// Nodes of the PRINCIPAL_NETWORK subnetwork
	//
	vector3D principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node
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
	void initNodesStructures(int nDays);

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

	// Creates the specific types of arcs
	void createArcsFromSource();

	// Update the arcs (their cost in particular)
	//
	void updateArcs();

	// Get info with the arc ID
	inline ArcType arcType(int a) {return allArcsTypes_[a];}
	//inline int arcOrigin(int a) {boost::graph_traits<Graph>::edge_descriptor id = a; return source(id, g_);}
	//inline int arcDestination(int a) {boost::graph_traits<Graph>::edge_descriptor id = a; return target(id, g_);}
	//inline int arcLength(int a) {boost::graph_traits<Graph>::edge_descriptor id = a; return get( &Arc_Properties::time, g_, id);}
	//inline int arcCost(int a) {boost::graph_traits<Graph>::edge_descriptor id = a; return get( &Arc_Properties::cost, g_, id);}


public:

	// Print functions.
	//
	void printGraph();


};


#endif /* SUBPROBLEM_H_ */
