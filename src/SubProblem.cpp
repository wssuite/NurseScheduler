/*
 * SubProblem.cpp
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#include "SubProblem.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/r_c_shortest_paths.hpp>

using std::stringstream;
using std::vector;



//////////////////////////////////////////////////
//
// Comparison override for 1 and 2 ressource paths
//
//////////////////////////////////////////////////

// 1 resource comparisons (== and <)
//
bool operator==( const spp_no_rc_res_cont& res_cont_1, const spp_no_rc_res_cont& res_cont_2 ){
	return ( res_cont_1.cost == res_cont_2.cost );
}

bool operator<( const spp_no_rc_res_cont& res_cont_1, const spp_no_rc_res_cont& res_cont_2 ){
	return ( res_cont_1.cost < res_cont_2.cost );
}

// 1 cost + 1 time resource comparisons (== and <)
bool operator==( const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2 ){
	return ( res_cont_1.cost == res_cont_2.cost	and res_cont_1.time == res_cont_2.time );
}

bool operator<( const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2 ){
	if( res_cont_1.cost > res_cont_2.cost )	return false;
	else if( res_cont_1.cost == res_cont_2.cost ) return res_cont_1.time < res_cont_2.time;
	else return true;
}
//////////////////////////////////////////////////



//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
//
SubProblem::SubProblem() {}

SubProblem::SubProblem(Scenario * scenario, Contract * contract):
	pScenario_(scenario), pContract_ (contract){

	maxRotationLength_ = pScenario_->nbWeeks_*7;

	createNodes();
	createArcs();

	nPathsMin_ = 0;



	printGraph();


}

SubProblem::~SubProblem(){}


//--------------------------------------------
//
// Functions for the NODES of the graph
//
//--------------------------------------------

// Function that creates the nodes of the network
void SubProblem::createNodes(){

	// Primary information needed
	int CD_min = pContract_->minConsDaysWork_;									// Minimum consecutive days worked for free
	int CD_max = pContract_->maxConsDaysWork_;									// Maximum consecutive days worked for free
	int nDays  = pScenario_->nbWeeks_ * 7;										// Number of days in the time period
	int nShifts= pScenario_->nbShifts_;											// Number of different shifts

	// INITIALIZATION
	nNodes_ = 0;
	initNodesStructures(nDays);

	// 1. SOURCE NODE
	//
	sourceNode_ = nNodes_;
	addSingleNode(SOURCE_NODE, 0, maxRotationLength_);

	// 2. SUBNETWORK FOR SHORT ROTATIONS (DEPENDS ON CD_min)
	//
	vector3D allowedShortShiftSuccessionsBySize = allowedShortSuccessions();	// A- Compute all allowed successions of shifts of length <= CD_min
	for(int c=1; c<=CD_min; c++){												// B- Duplication for different starting days: For each short succession length (c)
		vector2D successions = allowedShortShiftSuccessionsBySize[c];
		for(int k0=0; k0<nDays-c+1; k0++){										// For each possible starting date (k0)
			for(int m=0; m<allowedShortShiftSuccessionsBySize[c].size(); m++){	// For each allowed succession of length c
				addShortRotationToGraph(k0, successions[m], c);					// Creation of the node for that short rotation with starting date k0
			}
		}
	}


	// 3. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
	//
	for(int sh=1; sh<nShifts; sh++){											// For each possible worked shift
		int maxConsForShift = pScenario_->maxConsShifts_.at(sh);
		for(int k=0; k<nDays; k++){												// For each date
			for(int cons=0; cons<maxConsForShift; cons++){						// For each level of network
				addNodeToPrincipalNetwork(sh, k, cons);							// Add a node to the principal network
			}
		}
	}

	// 4. ROTATION LENGTH CHECK
	//
	// Entrance
	rotationLengthEntrance_ = nNodes_;											// Single node for the entrance in subnetwork
	addSingleNode(ROTATION_LENGTH_ENTRANCE, 0, maxRotationLength_);
	// Check nodes
	for(int l=CD_max; l<maxRotationLength_; l++){								// Check nodes: from CD_max (longest free) to maximum rotation length
		rotationLengthNodes_.insert(pair<int,int>(l,nNodes_));
		addSingleNode(ROTATION_LENGTH, l, maxRotationLength_);
	}
	// Exit
	rotationLengthExit_ = nNodes_;												// Single node for the exit (could probably be replaced by the sink)
	addSingleNode(ROTATION_LENGTH_EXIT, 0, maxRotationLength_);

	// 5. SINK NODE
	//
	sinkNode_ = nNodes_;
	addSingleNode(SINK_NODE, 0, maxRotationLength_);
}

// Addition of a single node (101)
void SubProblem::addSingleNode(NodeType type, int eat, int lat){
	add_vertex( Vertex_Properties( nNodes_, type, eat, lat ), g_ );
	allNodesTypes_.push_back(type);
	nNodes_++;

}

// Function that returns the vector of all allowed successions (by length from 0 [empty] to CD_min)
vector3D SubProblem::allowedShortSuccessions(){
	vector3D ans;

	// Primary information needed
	//
	int CD_min = pContract_->minConsDaysWork_;											// Minimum consecutive days worked for free
	int nShifts= pScenario_->nbShifts_;													// Number of different shifts

	vector2D v;
	ans.push_back(v);																	// Put an empty list of size 0

	// For all size from 1 to CD_min, compute all allowed shift successions.
	//
	for(int c=1; c<=CD_min; c++){
		vector< vector<int> > allowedShortShiftSuccessionsOfSizeC;
		if(c==1){																		// Size 1 -> special case, initialization -> add all single shift rotations
			for(int s=1; s<nShifts; s++){
				vector<int> shiftSuccession;
				shiftSuccession.push_back(s);
				allowedShortShiftSuccessionsOfSizeC.push_back(shiftSuccession);
			}
		}
		else {																			// Larger -> extend the previous one
			vector< vector<int> > legidShortShiftSuccessionsOfSizeCMinusOne = ans[c-1];
			for(int i=0; i<legidShortShiftSuccessionsOfSizeCMinusOne.size(); i++){		// For each short rotation of size c-1
				vector<int> succ = legidShortShiftSuccessionsOfSizeCMinusOne[i];
				int sh = succ[succ.size()-1];
				for(int newSh=1; newSh<nShifts; newSh++){								// For each possible shift at date c
					if(! pScenario_->isForbiddenSuccessor(newSh,sh)){					// IF the succession is allowed, then add it
						vector<int> newSucc;
						for(int s=0; s<succ.size(); s++) newSucc.push_back(succ[s]);
						newSucc.push_back(newSh);
						allowedShortShiftSuccessionsOfSizeC.push_back(newSucc);
					}
				}
			}
		}
		ans.push_back(allowedShortShiftSuccessionsOfSizeC);								// For all values of c, add the vector of possible successions of size c to the list
	}


	// Print for debug
	//
	for(int i=1; i<ans.size(); i++){
		std::cout << "#   +-----------------------+" << std::endl;
		std::cout << "#   | SUCCESSIONS OF SIZE " << i << " |" << std::endl;
		std::cout << "#   +-----------------------+" << std::endl;
		vector2D succs = ans[i];
		for(int j=0; j<succs.size(); j++){
			vector<int> succ = succs[j];
			std::cout << "#    | ";
			for(int k=0; k<succ.size(); k++){
				std::cout << " " << pScenario_->intToShift_[succ[k]];
			}
			std::cout << std::endl;
		}
	}

	return ans;
}

// Initiate variables for the nodes structures (vectors, etc.)
void SubProblem::initNodesStructures(int nDays){

	// Data
	//
	int CD_min = pContract_->minConsDaysWork_;		// Minimum consecutive days worked for free
	int nShifts= pScenario_->nbShifts_;				// Number of different shifts

	// All nodes
	//
	nNodes_ = 0;
	vector<NodeType> v; allNodesTypes_ = v;

	// Short rotations
	//
	for(int c=0; c<=CD_min; c++){
		vector<Rotation*> vr; shortRotations_.push_back(vr);
		vector<int> vi; shortRotationsNodes_.push_back(vi);
	}

	// Principal networks
	//
	for(int sh=0; sh<nShifts; sh++){
		vector2D v2D;
		if(sh!=0){
			for(int k=0; k<nDays; k++){
				vector<int> v1D;
				for(int cons=0; cons<pScenario_->maxConsShifts_.at(sh); cons++){
					v1D.push_back(-1);
				}
				v2D.push_back(v1D);
			}
		}
		principalNetworkNodes_.push_back(v2D);
	}
}

// Function that adds a short succession to the graph, that starts at k0, and consists in performing the given shift succession
void SubProblem::addShortRotationToGraph(int k0, vector<int> shiftSuccession, int length){

	// Create the effective rotation depending on starting date
	//
	Rotation rot (k0, shiftSuccession);

	// Add it to the list and create the node simultaneously
	//
	nodeToShortRotation_.insert(pair<int,Rotation> (nNodes_, rot));
	shortRotations_[length].push_back( & nodeToShortRotation_.at(nNodes_) );

	// Store the last shift
	//
	int lastShift = shiftSuccession[length-1];
	lastShiftOfShort_.insert(pair<int,int>( nNodes_, lastShift ));

	// Store the number of successive similar shifts in the end of that rotation
	//
	int n=0;
	while(n<length and shiftSuccession[length-1-n]==lastShift)n++;
	nLastShiftOfShort_.insert(pair<int,int>( nNodes_ , n ));

	// Store the id
	//
	shortRotationsNodes_[length].push_back(nNodes_);

	// Create the node
	//
	addSingleNode(SHORT_ROTATION, 0, maxRotationLength_);
}

// Add a node to the principal network of the graph, for shift sh, day k, and number of consecutive similar shifts cons
void SubProblem::addNodeToPrincipalNetwork(int sh, int k, int cons){

	// Store its ID in the vector3D
	//
	allNodesTypes_.push_back(PRINCIPAL_NETWORK);
	principalNetworkNodes_[sh][k][cons] = nNodes_;

	// Store the information backwards
	//
	principalToShift_.insert(pair<int,int>(nNodes_, sh));
	principalToDay_.insert(pair<int,int>(nNodes_, k));
	principalToCons_.insert(pair<int,int>(nNodes_, cons));

	// Create the node
	//
	addSingleNode(PRINCIPAL_NETWORK, 0, maxRotationLength_);
}


//--------------------------------------------
//
// Functions for the ARCS of the graph
//
//--------------------------------------------

// Function that creates the arcs of the network
void SubProblem::createArcs(){
	initArcsStructures();
	createArcsFromSource();
}

// Adds a single arc (origin, destination, cost, travel time, type)
void SubProblem::addSingleArc(int o, int d, int c, int t, ArcType type){
	add_edge( o, d, Arc_Properties( nArcs_, type, c, t ), g_ );
	allArcsTypes_.push_back(type);
	nArcs_++;
}

// Initializes the data structures used for the arcs
void SubProblem::initArcsStructures(){
	nArcs_ = 0;
	vector<int> v; arcsFromSource_ = v;
}

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void SubProblem::createArcsFromSource(){

	int CD_min = pContract_->minConsDaysWork_;									// Minimum consecutive days worked for free
	int CD_max = pContract_->maxConsDaysWork_;									// Maximum consecutive days worked for free
	int nDays  = pScenario_->nbWeeks_ * 7;										// Number of days in the time period
	int nShifts= pScenario_->nbShifts_;											// Number of different shifts

	for(int length=1; length<CD_min; length++){
		vector<int> shortRotLengthC = shortRotationsNodes_[length];
		for(int m=0; m<shortRotLengthC.size(); m++){
			int dest = shortRotLengthC[m];
			addSingleArc(sourceNode_, dest, MAX_COST, length, SOURCE_TO_SHORT);
		}
	}



}

// To change the costs of the arcs
void SubProblem::updateArcs(){}



// PRINT FUNCTIONS

// Print the graph
void SubProblem::printGraph(){
	stringstream rep;

	rep << "# " << std::endl;
	rep << "# GRAPH OF THE SUBPROBLEM " << std::endl;
	rep << "# " << std::endl;
	rep << "#   NODES " << std::endl;

	// THE NODES

	for(int v=0; v<nNodes_; v++){

		NodeType type_v = get( &Vertex_Properties::type, g_)[v];
		string nameType_v = nodeTypeName[type_v];
		int eat_v = get( &Vertex_Properties::eat, g_)[v];
		int lat_v = get( &Vertex_Properties::lat, g_)[v];


		rep << "# NODE  " << v << " \t" << nameType_v << " \t[" << eat_v << " " << lat_v << "] \t";

		if(type_v == SHORT_ROTATION){
			string lastShiftName = pScenario_->intToShift_.at(lastShiftOfShort_.at(v));
			int nLastShiftCons = nLastShiftOfShort_.at(v);
			Rotation r = nodeToShortRotation_.at(v);
			rep << "  k0= " << r.firstDay_ << "  ";
			for(int k=0; k<r.length_; k++) rep << " " << pScenario_->intToShift_[r.shifts_.at(r.firstDay_+k)];
			rep << "    --    FIN= " << lastShiftName << "-" << nLastShiftCons;

		}

		else if (type_v == PRINCIPAL_NETWORK){
			string shiftName = pScenario_->intToShift_[principalToShift_.at(v)];
			int k = principalToDay_.at(v);
			int cons = principalToCons_.at(v);

			rep << "  " << shiftName << "_" << k << "-" << cons;
		}

		rep << std::endl;

	}

	rep << std::endl;

	// THE ARCS

	for(int a=0; a<nArcs_; a++){
		rep << "# ARC   " << a << " \t" << arcTypeName[arcType(a)];
		//rep << " \t(" << arcOrigin(a) << "," << arcDestination(a) << ") \t";
		//rep << "c= " << arcCost(a) << " \tt=" << arcLength(a);
		rep << std::endl;
	}


	std::cout << rep.str();
}

























const int num_nodes = 5;
enum nodes { A, B, C, D, E };
char name[] = "ABCDE";

void SubProblem::testGraph_spprc(){


	std::cout << "# " << std::endl;
	std::cout << "# ====================================================================================" << std::endl;
	std::cout << "# = Fonction de test du probleme de plus court chemin avec contraintes de ressources =" << std::endl;
	std::cout << "# ====================================================================================" << std::endl;
	std::cout << "# " << std::endl;

	Graph g;

	add_vertex( Vertex_Properties( A, NONE_NODE, 0, 0 ), g );
	add_vertex( Vertex_Properties( B, NONE_NODE, 5, 20 ), g );
	add_vertex( Vertex_Properties( C, NONE_NODE, 6, 10 ), g );
	add_vertex( Vertex_Properties( D, NONE_NODE, 3, 12 ), g );
	add_vertex( Vertex_Properties( E, NONE_NODE, 0, 100 ), g );

	add_edge( A, C, Arc_Properties( 0, NONE_ARC, 1, 5 ), g );
	add_edge( B, B, Arc_Properties( 1, NONE_ARC, 2, 5 ), g );
	add_edge( B, D, Arc_Properties( 2, NONE_ARC, 1, 2 ), g );
	add_edge( B, E, Arc_Properties( 3, NONE_ARC, 2, 7 ), g );
	add_edge( C, B, Arc_Properties( 4, NONE_ARC, 7, 3 ), g );
	add_edge( C, D, Arc_Properties( 5, NONE_ARC, 3, 8 ), g );
	add_edge( D, E, Arc_Properties( 6, NONE_ARC, 1, 3 ), g );
	add_edge( E, A, Arc_Properties( 7, NONE_ARC, 1, 5 ), g );
	add_edge( E, B, Arc_Properties( 8, NONE_ARC, 1, 4 ), g );


	// the unique shortest path from A to E in the dijkstra-example.cpp is
	// A -> C -> D -> E
	// its length is 5
	// the following code also yields this result

	// with the above time windows, this path is infeasible
	// now, there are two shortest paths that are also feasible with respect to
	// the vertex time windows:
	// A -> C -> B -> D -> E and
	// A -> C -> B -> E
	// however, the latter has a longer total travel time and is therefore not
	// pareto-optimal, i.e., it is dominated by the former path
	// therefore, the code below returns only the former path

	// spp without resource constraints
	boost::graph_traits<Graph>::vertex_descriptor s = A;
	boost::graph_traits<Graph>::vertex_descriptor t = E;

	vector< vector< boost::graph_traits< Graph>::edge_descriptor> > opt_solutions;
	vector<spp_no_rc_res_cont> pareto_opt_rcs_no_rc;

	boost::r_c_shortest_paths(
			g,
			get( &Vertex_Properties::num, g ),
			get( &Arc_Properties::num, g ),
			s,
			t,
			opt_solutions,
			pareto_opt_rcs_no_rc,
			spp_no_rc_res_cont( 0 ),
			ref_no_res_cont(),
			dominance_no_res_cont(),
			std::allocator< boost::r_c_shortest_paths_label<Graph, spp_no_rc_res_cont> >(), boost::default_r_c_shortest_paths_visitor() );

	std::cout << "SPP without resource constraints:" << std::endl;
	std::cout << "Number of optimal solutions: ";
	std::cout << static_cast<int>( opt_solutions.size() ) << std::endl;
	for( int i = 0; i < static_cast<int>( opt_solutions.size() ); ++i )
	{
		std::cout << "The " << i << "th shortest path from A to E is: ";
		std::cout << std::endl;
		for( int j = static_cast<int>( opt_solutions[i].size() ) - 1; j >= 0; --j )
			std::cout << name[source( opt_solutions[i][j], g )] << std::endl;
		std::cout << "E" << std::endl;
		std::cout << "Length: " << pareto_opt_rcs_no_rc[i].cost << std::endl;
	}
	std::cout << std::endl;

	// spptw
	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;

	r_c_shortest_paths(
			g,
			get( &Vertex_Properties::num, g ),
			get( &Arc_Properties::num, g ),
			s,
			t,
			opt_solutions_spptw,
			pareto_opt_rcs_spptw,
			spp_spptw_res_cont( 0, 0 ),
			ref_spptw(),
			dominance_spptw(),
			std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
			boost::default_r_c_shortest_paths_visitor() );

	std::cout << "SPP with time windows:" << std::endl;
	std::cout << "Number of optimal solutions: ";
	std::cout << static_cast<int>( opt_solutions.size() ) << std::endl;
	for( int i = 0; i < static_cast<int>( opt_solutions.size() ); ++i )
	{
		std::cout << "The " << i << "th shortest path from A to E is: ";
		std::cout << std::endl;
		for( int j = static_cast<int>( opt_solutions_spptw[i].size() ) - 1;
				j >= 0;
				--j )
			std::cout << name[source( opt_solutions_spptw[i][j], g )] << std::endl;
		std::cout << "E" << std::endl;
		std::cout << "Length: " << pareto_opt_rcs_spptw[i].cost << std::endl;
		std::cout << "Time: " << pareto_opt_rcs_spptw[i].time << std::endl;
	}

	// utility function check_r_c_path example
	std::cout << std::endl;
	bool b_is_a_path_at_all = false;
	bool b_feasible = false;
	bool b_correctly_extended = false;
	spp_spptw_res_cont actual_final_resource_levels( 0, 0 );
	boost::graph_traits<Graph>::edge_descriptor ed_last_extended_arc;
	check_r_c_path( g,
			opt_solutions_spptw[0],
			spp_spptw_res_cont( 0, 0 ),
			true,
			pareto_opt_rcs_spptw[0],
			actual_final_resource_levels,
			ref_spptw(),
			b_is_a_path_at_all,
			b_feasible,
			b_correctly_extended,
			ed_last_extended_arc );
	if( !b_is_a_path_at_all )
		std::cout << "Not a path." << std::endl;
	if( !b_feasible )
		std::cout << "Not a feasible path." << std::endl;
	if( !b_correctly_extended )
		std::cout << "Not correctly extended." << std::endl;
	if( b_is_a_path_at_all && b_feasible && b_correctly_extended )
	{
		std::cout << "Actual final resource levels:" << std::endl;
		std::cout << "Length: " << actual_final_resource_levels.cost << std::endl;
		std::cout << "Time: " << actual_final_resource_levels.time << std::endl;
		std::cout << "OK." << std::endl;
	}
}
