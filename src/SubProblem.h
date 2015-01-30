/*
 * SubProblem.h
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// Pour l'instant cette classe n'est pas reliee au reste du
// projet. Ne sert qu'a tester Boost et les algos de plus courts
// chemins avec contraintes de ressources.
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//
// The following structures are the properties of the nodes and arcs:
//   For the nodes: id, earliest arrival time, and latest arrival time
//   For the arcs : id, cost and time (time is the ressource here)
//   For the graph: usual stuff + nodes and ard special properties
//
//////////////////////////////////////////////////////////////////////////

// Nodes specific properties for RC
//
struct SPPRC_Example_Graph_Vert_Prop{

	// Constructor
	//
	SPPRC_Example_Graph_Vert_Prop( int n = 0, int e = 0, int l = 0 ) : num( n ), eat( e ), lat( l ) {}

	// id
	//
	int num;

	// earliest arrival time
	//
	int eat;

	// latest arrival time
	//
	int lat;
};

// Arcs specific properties for RC
//
struct SPPRC_Example_Graph_Arc_Prop{

	// Constructor
	//
	SPPRC_Example_Graph_Arc_Prop( int n = 0, int c = 0, int t = 0 ) : num( n ), cost( c ), time( t ) {}

	// id
	//
	int num;

	// traversal cost
	//
	int cost;

	// traversal time
	//
	int time;
};

// Graph with RC generic structure
//
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, SPPRC_Example_Graph_Vert_Prop, SPPRC_Example_Graph_Arc_Prop> SPPRC_Example_Graph;

//////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The following functions / data structures are used for spp without resource constraints (cost only):
//   Resource container ("resource" = cost)
//   Comparison override --> in SubProblem.cpp
//   Cost extension function
//   Dominance function
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

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
	inline bool operator()(const SPPRC_Example_Graph& g, spp_no_rc_res_cont& new_cont, const spp_no_rc_res_cont& old_cont, boost::graph_traits<SPPRC_Example_Graph>::edge_descriptor ed ) const{
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

///////////////////////////////////////////////////////////////////////////////////////////////////////

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
	inline bool operator()( const SPPRC_Example_Graph& g, spp_spptw_res_cont& new_cont,	const spp_spptw_res_cont& old_cont,	boost::graph_traits<SPPRC_Example_Graph>::edge_descriptor ed ) const{
		const SPPRC_Example_Graph_Arc_Prop& arc_prop = get( boost::edge_bundle, g )[ed];
		const SPPRC_Example_Graph_Vert_Prop& vert_prop = get( boost::vertex_bundle, g )[target( ed, g )];
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

class SubProblem {

	// En travaux...

public:
	void testGraph_spprc();



};


#endif /* SUBPROBLEM_H_ */
