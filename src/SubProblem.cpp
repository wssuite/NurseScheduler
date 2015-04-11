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



/////////////////////////////////////////////////
//
// Comparison override for 1 and 2 resource paths
//
/////////////////////////////////////////////////

// 1 resource comparisons (== and <)
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
SubProblem::SubProblem() {}

SubProblem::SubProblem(Scenario * scenario, Demand * demand, const Contract * contract):
	pScenario_(scenario), pDemand_(demand), pContract_ (contract),
	CDMin_(contract->minConsDaysWork_), maxRotationLength_(demand->nbDays_), nDays_(demand->nbDays_){


	init();

	initShortSuccessions();

	createNodes();
	createArcs();

	// Set all arc and node status to authorized
	for(int v=0; v<nNodes_; v++) nodeStatus_.push_back(true);
	for(int a=0; a<nArcs_; a++) arcStatus_.push_back(true);
	for(int k=0; k<nDays_; k++){
		vector<bool> v;
		for(int s=0; s<pScenario_->nbShifts_; s++)
			v.push_back(true);
		dayShiftStatus_.push_back(v);
	}

	nPathsMin_ = 0;

	std::cout << "# A new subproblem has been created for contract " << contract->name_ << std::endl;

	//printGraph();
	//std::cout << printSummaryOfGraph();
	//printShortSucc();

}

SubProblem::~SubProblem(){}

// Initialization function
void SubProblem::init(){

	// Maximum number of consecutive days worked by a nurse ending at day -1
	//

	// TODO : Déterminer le nombre de jours maximum des séries en cours pour les nurses
	maxOngoingDaysWorked_ = 2;

	// Initialization of isUnlimited_ and nLevelsByShift_
	//
	vector<bool> v; isUnlimited_ = v; isUnlimited_.push_back(false);
	vector<int> w; maxvalConsByShift_ = w; maxvalConsByShift_.push_back(0);

	for(int sh=1; sh<pScenario_->nbShifts_; sh++){
		isUnlimited_.push_back( pScenario_->maxConsShifts_[sh] >= nDays_ + maxOngoingDaysWorked_ );
		int nl = isUnlimited_[sh] ? pScenario_->minConsShifts_[sh] : pScenario_->maxConsShifts_[sh];
		maxvalConsByShift_.push_back( nl );
	}



	// id and arcCost of best succession (given a triplet s,k,n)
	idBestShortSuccCDMin_.clear();
	arcCostBestShortSuccCDMin_.clear();
	for(int s=0; s<pScenario_->nbShifts_; s++){
		vector2D v2; vector<vector<double> > w2;
		int n = maxvalConsByShift_[s]+1;
		Tools::initVector2D(&v2, nDays_, n);
		idBestShortSuccCDMin_.push_back(v2);
		Tools::initDoubleVector2D(&w2, nDays_, n);
		arcCostBestShortSuccCDMin_.push_back(w2);
	}
}

// Initializes the short successions. Should only be used ONCE (when creating the SubProblem).
void SubProblem::initShortSuccessions(){

	// Primary information needed
	//
	int nShifts= pScenario_->nbShifts_;

	// Put an empty list of size 0 for all data because there exists no succession of length 0/
	//
	vector2D v2; vector<int> v1,v1bis; vector<double> vd1;
	allowedShortSuccBySize_.push_back(v2);
	lastShiftOfShortSucc_.push_back(v1);
	nLastShiftOfShortSucc_.push_back(v1bis);
	baseArcCostOfShortSucc_.push_back(vd1);

	// Initialize the other way round
	for(int s=0; s<nShifts; s++){
		vector2D allShortSuccCDMinByLastShiftCons2;
		for(int n=0; n<=maxvalConsByShift_[s]; n++){
			vector<int> allShortSuccCDMinByLastShiftCons1;
			allShortSuccCDMinByLastShiftCons2.push_back(allShortSuccCDMinByLastShiftCons1);
		}
		allShortSuccCDMinByLastShiftCons_.push_back(allShortSuccCDMinByLastShiftCons2);
	}

	// For all size from 1 to CD_min, compute all allowed shift successions.
	//
	for(int c=1; c<=CDMin_; c++){

		// Blank data
		vector2D allSuccSizeC;
		vector<int> lastShiftSucc;
		vector<int> nLastShiftSucc;
		vector<double> arcCostSucc;

		// Compute new succession
		//

		// Size 1 -> special case, initialization -> add all single shift rotations
		//
		if(c==1){
			for(int s=1; s<nShifts; s++){
				vector<int> shiftSuccession; shiftSuccession.push_back(s);		// Create Succession
				allSuccSizeC.push_back(shiftSuccession);						// Add it to the possibilities
				lastShiftSucc.push_back(s);										// Record its last shift
				nLastShiftSucc.push_back(1);									// Only 1 successive performed so far
				arcCostSucc.push_back(0);										// No succession ended yet
			}
		}

		// Larger but not last -> Extend the previous one by extending each of size c-1 in all possible ways
		//
		else if(c<CDMin_){
			// For each short rotation of size c-1
			for(int i=0; i<allowedShortSuccBySize_[c-1].size(); i++){
				vector<int> succ (allowedShortSuccBySize_[c-1][i]);
				int lastSh = succ[succ.size()-1];
				int nLast = nLastShiftOfShortSucc_[c-1][i];
				double cost = baseArcCostOfShortSucc_[c-1][i];
				// For each possible new shift s.t. the succession is allowed
				for(int newSh=1; newSh<nShifts; newSh++){
					if(! pScenario_->isForbiddenSuccessor(newSh,lastSh)){

						vector<int> newSucc (succ); newSucc.push_back(newSh);		// Create Succession
						allSuccSizeC.push_back(newSucc);							// Add it to the possibilities
						lastShiftSucc.push_back(newSh);								// Record its last shift
						if (newSh == lastSh){										// Depending on the previous one, update number of consecutive and cost
							nLastShiftSucc.push_back(nLast+1);
							arcCostSucc.push_back(cost);
						} else {
							nLastShiftSucc.push_back(1) ;
							arcCostSucc.push_back(cost + consShiftCost(lastSh, nLast));
						}

					}
				}
			}
		}

		// Maximum allowed size -> more things to consider
		//
		else {
			// For each short rotation of size c-1
			for(int i=0; i<allowedShortSuccBySize_[c-1].size(); i++){
				vector<int> succ (allowedShortSuccBySize_[c-1][i]);
				int lastSh = succ[succ.size()-1];
				int nLast = nLastShiftOfShortSucc_[c-1][i];
				double cost = baseArcCostOfShortSucc_[c-1][i];
				// For each possible new shift s.t. the succession is allowed
				for(int newSh=1; newSh<nShifts; newSh++){
					if(! pScenario_->isForbiddenSuccessor(newSh,lastSh)){

						vector<int> newSucc (succ); newSucc.push_back(newSh);		// Create Succession
						allSuccSizeC.push_back(newSucc);							// Add it to the possibilities
						lastShiftSucc.push_back(newSh);								// Record its last shift
						int newNLast = 1;
						double newCost = cost;
						if (newSh == lastSh){	// BUT : add the cost if longer than the maximum allowed
							newNLast += nLast;
						} else {
							newCost += consShiftCost(lastSh, nLast) ;
						}
						nLastShiftSucc.push_back(newNLast);
						arcCostSucc.push_back(newCost);

						// Since it is the longest one, record the tables the other way round
						//
						int n = min(newNLast, maxvalConsByShift_[newSh]);
						allShortSuccCDMinByLastShiftCons_[newSh][n].push_back(allSuccSizeC.size()-1);

					}
				}
			}
		}

		// Store all vectors
		//
		allowedShortSuccBySize_.push_back(allSuccSizeC);
		lastShiftOfShortSucc_.push_back(lastShiftSucc);
		nLastShiftOfShortSucc_.push_back(nLastShiftSucc);
		baseArcCostOfShortSucc_.push_back(arcCostSucc);

	}

	// Print for debug
	//
	/*
	for(int i=1; i<allowedShortSuccBySize_.size(); i++){
		std::cout << "#   +-----------------------+" << std::endl;
		std::cout << "#   | SUCCESSIONS OF SIZE " << i << " |" << std::endl;
		std::cout << "#   +-----------------------+" << std::endl;
		for(int j=0; j<allowedShortSuccBySize_[i].size(); j++){
			vector<int> succ = allowedShortSuccBySize_[i][j];
			std::cout << "#    | ";
			std::cout << "c=" << baseArcCostOfShortSucc_[i][j];
			baseArcCostOfShortSucc_[i][j] == 0 ? cout<< " " : cout << "";
			cout << " ; ";
			//std::cout << "lastSh=" << pScenario_->intToShift_[lastShiftOfShortSucc_[i][j]] << " ; ";
			//std::cout << "nLast=" << nLastShiftOfShortSucc_[i][j] << " ; ";
			for(int k=0; k<succ.size(); k++){
				std::cout << " " << pScenario_->intToShift_[succ[k]];
			}
			std::cout << std::endl;
		}
	}
	for(int s=0; s<nShifts; s++){
		for(int n=0; n<=maxvalConsByShift_[s]; n++){
			cout << "# Succession ending with " << n << " days of " << pScenario_->intToShift_[s] << " :";
			for(int i : allShortSuccCDMinByLastShiftCons_[s][n]){
				cout << " ";
				vector<int> v = allowedShortSuccBySize_[CDMin_][i];
				for(int sh : v)
				cout << ((pScenario_->intToShift_[sh]).at(0));
			}
			cout << endl;
		}
	}
	getchar();
	*/

}

// Cost function for consecutive identical shifts
//
double SubProblem::consShiftCost(int sh, int n){
	if(pScenario_->minConsShifts_[sh] - n > 0) return (WEIGHT_CONS_SHIFTS * ( pScenario_->minConsShifts_[sh] - n ) );
	if(n - pScenario_->maxConsShifts_[sh] > 0) return (WEIGHT_CONS_SHIFTS * ( n - pScenario_->maxConsShifts_[sh] ) );
	return 0;
}

// Cost function for consecutive identical shifts
//
double SubProblem::consDaysCost(int n){
	if(pContract_->minConsDaysWork_ - n > 0) return (WEIGHT_CONS_DAYS_WORK * ( pContract_->minConsDaysWork_ - n ) );
	if(n - pContract_->maxConsDaysWork_ > 0) return (WEIGHT_CONS_DAYS_WORK * ( n - pContract_->maxConsDaysWork_ ) );
	return 0;
}



//--------------------------------------------
//
// Solve function
//
//--------------------------------------------

// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
bool SubProblem::solve(LiveNurse* nurse, Costs * costs, set<pair<int,int> > forbiddenDayShifts, bool optimality, int maxRotationLength){

	std::cout << "# " << std::endl;
	std::cout << "# Solving subproblem for nurse " << nurse->name_ << " (id:" <<  nurse->id_ << "), " << pContract_->name_ << std::endl;
	//std::cout << printSummaryOfGraph();

	// TODO : modify if needed
	//
	resetAuthorizations();
	resetSolutions();

	// Forbid arcs or nodes if needed
	//
	// TODO : delete following line and replace fDayShifts by forbiddenDayShifts in the rest of the function
	//set< pair<int,int> > fDayShifts = randomForbiddenShifts(25);
	//printForbiddenDayShift();

	// Basic data (nurse, reduced costs, maximum rotation length)
	//
	pLiveNurse_ = nurse;
	pCosts_ = costs;
	maxRotationLength_ = maxRotationLength;

	initStructuresForSolve(nurse, costs, forbiddenDayShifts, maxRotationLength);
	//generateRandomCosts(-50, 50);


	forbid(forbiddenDayShifts);
	printForbiddenDayShift();

	updateArcCosts();

	//printShortArcs();



	// Solving the problem
	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;
	spp_spptw_res_cont rc (0,0);

	r_c_shortest_paths(
			g_,
			get( &Vertex_Properties::num, g_ ),
			get( &Arc_Properties::num, g_ ),
			sourceNode_,
			sinkNode_,
			opt_solutions_spptw,
			pareto_opt_rcs_spptw,
			rc,
			ref_spptw(),
			dominance_spptw(),
			std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
			boost::default_r_c_shortest_paths_visitor() );

	return addRotationsFromPaths(opt_solutions_spptw, pareto_opt_rcs_spptw, true);
}

// Transforms the solutions found into proper rotations.
//
bool SubProblem::addRotationsFromPaths(vector< vector< boost::graph_traits<Graph>::edge_descriptor> > paths, vector<spp_spptw_res_cont> resources,
		bool negativeOnly){
	bool oneFound = false;
	// For each path of the list, record the corresponding rotation (if negativeOnly=true, do it only if the dualCost < 0)
	for(int p=0; p < paths.size(); ++p){
		Rotation rot = rotationFromPath(paths[p], resources[p]);
		if( (! negativeOnly) or (rot.dualCost_ < 0)){
			//printPath(paths[p], resources[p]);
			theRotations_.push_back(rot);
			nPaths_ ++;
			oneFound = true;
		}
	}
	printAllRotations();
	return oneFound;
}

// Adds a rotation made from the given path to the current list of answers and increases their counter
//
Rotation SubProblem::rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource){

	int firstDay = -1;
	int dualCost = 0;
	vector<int> shiftSuccession;

	// All arcs are consecutively considered
	//
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		ArcType aType = arcType(a);
		int origin = boost::source( path[j], g_ );
		int destin = boost::target( path[j], g_ );

		// A. Arc from source (equivalent to short rotation
		if(aType == SOURCE_TO_PRINCIPAL){
			firstDay =  principalToDay_[destin] - CDMin_ + 1;
			// TODO : *** Utiliser un swap à la place par exemple ***
			for(int s: static_cast<vector <int> >( allowedShortSuccBySize_[CDMin_][ shortSuccCDMinIdFromArc_.at(a) ] )){
				shiftSuccession.push_back(s);
			}
		}

		// B. Arc to a new day
		else if(aType == SHIFT_TO_NEWSHIFT or aType == SHIFT_TO_SAMESHIFT or aType == REPEATSHIFT){
			shiftSuccession.push_back( principalToShift_[destin] );
		}

		dualCost += arcCost(a);
	}

	Rotation rot (firstDay, shiftSuccession, pLiveNurse_);
	//rot.computeCost(pScenario_, pLiveNurse_->pWishesOff_, nDays_);
	rot.dualCost_ = dualCost;
	return rot;
}

// Resets all solutions data (rotations, number of solutions, etc.)
//
void SubProblem::resetSolutions(){
	theRotations_.clear();
	nPaths_ = 0;

}



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
	int nShifts= pScenario_->nbShifts_;											// Number of different shifts

	// INITIALIZATION
	nNodes_ = 0;
	initNodesStructures();

	// 1. SOURCE NODE
	//
	sourceNode_ = nNodes_;
	addSingleNode(SOURCE_NODE, 0, maxRotationLength_);

	// 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
	//
	for(int sh=1; sh<nShifts; sh++){											// For each possible worked shift
		for(int k=CDMin_-1; k<nDays_; k++){											// For each date
			for(int cons=1; cons<=maxvalConsByShift_[sh]; cons++){					// For each level of network
				addNodeToPrincipalNetwork(sh, k, cons);								// Add a node to the principal network
			}
		}
	}

	// 3. ROTATION LENGTH CHECK
	//
	// Entrance
	rotationLengthEntrance_ = nNodes_;											// Single node for the entrance in subnetwork
	addSingleNode(ROTATION_LENGTH_ENTRANCE, 0, maxRotationLength_);
	// Check nodes
	for(int l=CD_max; l<maxRotationLength_; l++){								// Check nodes: from CD_max (longest free) to maximum rotation length
		rotationLengthNodes_.insert(pair<int,int>(l,nNodes_));
		addSingleNode(ROTATION_LENGTH, 0, l);
	}

	// 4. SINK NODE
	//
	sinkNode_ = nNodes_;
	addSingleNode(SINK_NODE, 0, maxRotationLength_);
}

// Initiate variables for the nodes structures (vectors, etc.)
void SubProblem::initNodesStructures(){

	// Data
	//
	int nShifts= pScenario_->nbShifts_;				// Number of different shifts

	// All nodes
	//
	nNodes_ = 0;
	vector<NodeType> v; allNodesTypes_ = v;

	// Principal networks
	//
	for(int sh=0; sh<nShifts; sh++){
		vector2D v2D;
		if(sh!=0){
			for(int k=0; k<nDays_; k++){
				vector<int> v1D;
				for(int cons=0; cons<=maxvalConsByShift_[sh]; cons++){
					v1D.push_back(-1);
				}
				v2D.push_back(v1D);
			}
		}
		principalNetworkNodes_.push_back(v2D);
	}
}

// Addition of a single node (101)
void SubProblem::addSingleNode(NodeType type, int eat, int lat){
	add_vertex( Vertex_Properties( nNodes_, type, eat, lat ), g_ );
	allNodesTypes_.push_back(type);
	nNodes_++;

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
	createArcsSourceToPrincipal();
	createArcsPrincipalToPrincipal();
	createArcsAllRotationSize();
}

// Adds a single arc (origin, destination, cost, travel time, type)
void SubProblem::addSingleArc(int o, int d, int baseCost, int t, ArcType type){
	boost::graph_traits< Graph>::edge_descriptor e = (add_edge( o, d, Arc_Properties( nArcs_, type, baseCost, t ), g_ )).first;
	arcsDescriptors_.push_back(e);
	allArcsTypes_.push_back(type);
	arcBaseCost_.push_back(baseCost);
	nArcs_++;
}

// Initializes the data structures used for the arcs
void SubProblem::initArcsStructures(){
	nArcs_ = 0;
	vector< boost::graph_traits< Graph>::edge_descriptor > ve; arcsDescriptors_ = ve;

	// Initialization of info -> arcId data structures
	arcsFromSource_.clear();
	arcsShiftToNewShift_.clear();
	arcsShiftToSameShift_.clear();
	arcsShiftToEndsequence_.clear();
	arcsRepeatShift_.clear();
	arcsPrincipalToRotsizein_.clear();
	arcsRotsizeinToRotsize_.clear();
	arcsRotsizeToRotsizeout_.clear();
	// VECTORS 3 D
	for(int s=0; s<pScenario_->nbShifts_; s++){
		vector2D v2, w2, x2;
		Tools::initVector2D(&v2, nDays_, maxvalConsByShift_[s]+1); arcsFromSource_.push_back(v2);
		Tools::initVector2D(&w2, nDays_, maxvalConsByShift_[s]+1); arcsShiftToSameShift_.push_back(w2);
		Tools::initVector2D(&x2, nDays_, maxvalConsByShift_[s]+1); arcsShiftToEndsequence_.push_back(x2);
	}
	Tools::initVector3D(&arcsShiftToNewShift_, pScenario_->nbShifts_, pScenario_->nbShifts_, nDays_);
	// VECTORS 2 D
	Tools::initVector2D(&arcsRepeatShift_, pScenario_->nbShifts_, nDays_);
	Tools::initVector2D(&arcsPrincipalToRotsizein_, pScenario_->nbShifts_, nDays_);
	// MAPS
	arcsRotsizeinToRotsize_.clear();
	arcsRotsizeToRotsizeout_.clear();

}

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void SubProblem::createArcsSourceToPrincipal(){

	// DATA
	int nShifts = pScenario_->nbShifts_;
	int origin, destin;

	for(int sh=1; sh<nShifts; sh++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int nCons=1; nCons<=maxvalConsByShift_[sh]; nCons ++){
				origin = sourceNode_;
				destin = principalNetworkNodes_[sh][k][nCons];
				arcsFromSource_[sh][k][nCons] = nArcs_;
				addSingleArc(origin, destin, 0, CDMin_, SOURCE_TO_PRINCIPAL);
			}
		}
	}

}

// Create all arcs within the principal network
void SubProblem::createArcsPrincipalToPrincipal(){

	// DATA
	int nShifts = pScenario_->nbShifts_;
	int origin, destin;

	// FOR EACH OF THE SUBNETWORKS AND EACH OF THE DAYS
	//
	for (int sh=1; sh<nShifts; sh++){
		for(int k=CDMin_-1; k<nDays_-1; k++){

			//   1. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS NOT REACHED YET
			//
			for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons ++){
				origin = principalNetworkNodes_[sh][k][nCons];
				destin = principalNetworkNodes_[sh][k+1][nCons+1];
				arcsShiftToSameShift_[sh][k][nCons] = nArcs_;
				addSingleArc(origin, destin, 0, 1, SHIFT_TO_SAMESHIFT);
			}

			//   2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
			//
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = principalNetworkNodes_[sh][k+1][maxvalConsByShift_[sh]];
			arcsRepeatShift_[sh][k] = nArcs_;
			double cost = isUnlimited(sh) ? 0 : WEIGHT_CONS_SHIFTS;
			addSingleArc(origin, destin, cost, 1, REPEATSHIFT);

			// 3. WORK ONE MORE DAY ON A DIFFERENT SHIFT (CHANGE SUBNETWORK)
			//
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			for(int newSh=1; newSh<nShifts; newSh++){
				if(newSh != sh and ! pScenario_->isForbiddenSuccessor(newSh,sh)){
					int destin = principalNetworkNodes_[newSh][k+1][1];
					arcsShiftToNewShift_[sh][newSh][k] = nArcs_;
					addSingleArc(origin, destin, 0, 1, SHIFT_TO_NEWSHIFT);
				}
			}

			// 4. END THE CONSECUTIVE SHIFT SEQUENCE
			//
			for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons++){
				origin = principalNetworkNodes_[sh][k][nCons];
				destin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
				arcsShiftToEndsequence_[sh][k][nCons] = nArcs_;
				addSingleArc(origin, destin, consShiftCost(sh, nCons), 0, SHIFT_TO_ENDSEQUENCE);
			}
		}

		// SPECIAL CASE: LAST DAY
		//
		for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons++){
			origin = principalNetworkNodes_[sh][nDays_-1][nCons];
			destin = principalNetworkNodes_[sh][nDays_-1][maxvalConsByShift_[sh]];
			arcsShiftToEndsequence_[sh][nDays_-1][nCons] = nArcs_;
			addSingleArc(origin, destin, consShiftCost(sh, nCons), 0, SHIFT_TO_ENDSEQUENCE);
		}
	}
}

// Create all arcs that involve the rotation size checking subnetwork (incoming, internal, and exiting that subnetwork)
void SubProblem::createArcsAllRotationSize(){

	int nShifts = pScenario_->nbShifts_;
	int origin, destin;

	// 1. ALL INCOMING ARCS
	//
	for(int sh=1; sh<nShifts; sh++){											// For all shifts
		for(int k=CDMin_-1; k<nDays_; k++){										// For all days
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = rotationLengthEntrance_;
			arcsPrincipalToRotsizein_[sh][k] = nArcs_;
			addSingleArc(origin, destin, 0, 0, PRINCIPAL_TO_ROTSIZE);			// Allow to stop rotation that day
		}
	}

	// 2. ALL INTERNAL ARCS
	//
	for(map<int,int>::iterator itRLN = rotationLengthNodes_.begin(); itRLN != rotationLengthNodes_.end(); ++itRLN){
		// From entrance to checknode
		origin = rotationLengthEntrance_;
		destin = itRLN->second;
		arcsRotsizeinToRotsize_.insert(pair<int,int>(itRLN->first, nArcs_));
		addSingleArc(origin, destin, consDaysCost(itRLN->first), 0, ROTSIZEIN_TO_ROTSIZE);
		// From checknode to exit
		origin = itRLN->second;
		destin = sinkNode_;
		arcsRotsizeToRotsizeout_.insert(pair<int,int>( itRLN->first, nArcs_));
		addSingleArc(origin, destin, 0, 0, ROTSIZE_TO_SINK);
	}
}



//--------------------------------------------
//
// Functions for the pricing of the short rotations
//
//--------------------------------------------

// Initializes some cost vectors that depend on the nurse
void SubProblem::initStructuresForSolve(LiveNurse* nurse, Costs * costs, set<pair<int,int> > forbiddenDayShifts, int maxRotationLength){

	// Start and End weekend costs
	//
	Tools::initDoubleVector(&startWeekendCosts_,nDays_);
	Tools::initDoubleVector(&endWeekendCosts_,nDays_);
	if(pLiveNurse_->needCompleteWeekends()){
		for(int k=0; k<nDays_; k++){
			if(Tools::isSaturday(k)) endWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
			else if(Tools::isSunday(k)) startWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
		}
	}

	// Preference costs. WARNING: STARTING DATE OF THE SCENARIO IS NOT THAT OF THE PREFERENCE ??
	// TODO: Check if shifting is necessary
	//
	Tools::initDoubleVector2D(&preferencesCosts_, nDays_, pScenario_->nbShifts_);
	for(map<int,set<int> >::iterator prefList = pLiveNurse_->pWishesOff_->begin(); prefList != pLiveNurse_->pWishesOff_->end(); ++ prefList){
		for(int sh : prefList->second){
			preferencesCosts_[prefList->first][sh] = WEIGHT_PREFERENCES;
		}
	}

	// id and arcCost of best succession (given a triplet s,k,n)
	idBestShortSuccCDMin_.clear();
	arcCostBestShortSuccCDMin_.clear();
	for(int s=0; s<pScenario_->nbShifts_; s++){
		vector2D v2; vector<vector<double> > w2;
		int n = maxvalConsByShift_[s]+1;
		Tools::initVector2D(&v2, nDays_, n);
		idBestShortSuccCDMin_.push_back(v2);
		Tools::initDoubleVector2D(&w2, nDays_, n);
		arcCostBestShortSuccCDMin_.push_back(w2);
	}

}

// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
//
void SubProblem::priceShortSucc(){

	//cout << "# " << endl;
	//cout << "# PRICE SHORT SUCCESSIONS:" << endl;

	for(int s=1; s<pScenario_->nbShifts_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){

				//cout << "# (" << pScenario_->intToShift_[s].at(0) << "," << k << "," << n << ") <- ";

				int bestSuccId = -1;
				double bestCost = MAX_COST;
				idBestShortSuccCDMin_[s][k][n] = -1;
				arcCostBestShortSuccCDMin_[s][k][n] = MAX_COST;
				for(int i=0; i<(allShortSuccCDMinByLastShiftCons_[s][n]).size(); i++){
					int curSuccId = allShortSuccCDMinByLastShiftCons_[s][n][i];
					vector<int> succ = allowedShortSuccBySize_[CDMin_][curSuccId];

					//cout << " ";
					//for(int ss=0; ss<succ.size(); ss++) cout << pScenario_->intToShift_[succ[ss]].at(0);
					//cout << " (";
					//canSuccStartHere( succ, k-CDMin_+1 ) ? cout << "OK" : cout << "NO";
					//cout << ") ";

					// SUCCESSION IS TAKEN INTO ACCOUNT ONLY IF IT DOES NOT VIOLATE ANY FORBIDDEN DAY-SHIFT COUPLE
					if(canSuccStartHere( succ, k-CDMin_+1 )){
						double curCost = costArcShortSucc(CDMin_, curSuccId, k-CDMin_+1);
						if(curCost < bestCost){
							idBestShortSuccCDMin_[s][k][n] = curSuccId;
							arcCostBestShortSuccCDMin_[s][k][n] = curCost;
							bestSuccId = curSuccId;
							bestCost = curCost;
						}
					}


				}
				// IF NO VALID SUCCESSION, THEN FORBID THE ARC
				if(bestCost == MAX_COST){
					forbidArc( arcsFromSource_[s][k][n] );
				} else {
					//cout << "  \t\t  [chosen : ";
					//vector<int> chosenSucc = allowedShortSuccBySize_[CDMin_][idBestShortSuccCDMin_[s][k][n]];
					//for(int ss=0; ss<chosenSucc.size(); ss++) cout << pScenario_->intToShift_[chosenSucc[ss]].at(0);
					//cout << "]";
				}

				//cout << endl;

			}
		}
	}
}

// Given a short succession and a start date, returns the cost of the corresponding arc
//
double SubProblem::costArcShortSucc(int size, int succId, int startDate){
	double ANS = 0;

	// If the rotation starts on the first day and if the nurse was already working before
	//
	if(startDate == 0 and pLiveNurse_->pStateIni_->consDaysWorked_ > 0){
		int nConsIni = pLiveNurse_->pStateIni_->consDaysWorked_;
		int shiftIni = pLiveNurse_->pStateIni_->shift_;
		int nConsecCurrent = 1;

		// IF the nurse continues working on same shift -> subtract the cost of the current sequence (later replaced by the one for the longer sequence)
		if(shiftIni == allowedShortSuccBySize_[size][succId][0]){
			ANS -= consShiftCost(shiftIni, nConsIni);
			nConsecCurrent += nConsIni;
		}

		// SUCCESSIVE SHIFTS
		for(int i=1; i<size; i++){
			int prevShift = allowedShortSuccBySize_[size][succId][i-1];
			int shift = allowedShortSuccBySize_[size][succId][i];
			if(shift != prevShift){
				ANS += consShiftCost(prevShift, nConsecCurrent);
				nConsecCurrent = 1;
			}
		}
		// IF MORE THAN THE MAXIMUM OF CONSECUTIVE ALLOWED
		if(nConsecCurrent > pScenario_->maxConsShifts_[allowedShortSuccBySize_[size][succId][size-1]])
			ANS += consShiftCost(allowedShortSuccBySize_[size][succId][size-1], nConsecCurrent);
		// WEEKEND REDUCED COST (does not count if day 0 is a Sunday because already taken into account in the MP)
		if(Tools::containsWeekend( 1 , size-1 )) ANS += pCosts_->workedWeekendCost();
	}

	// If the rotation does not start on the first day
	//
	else {
		// SUCCESSIVE SHIFTS
		ANS += baseArcCostOfShortSucc_[size][succId];
		// COMPLETE WEEKEND
		ANS += startWeekendCosts_[startDate];
		// WEEKEND REDUCED COST
		if(Tools::containsWeekend(startDate, startDate + size - 1))	ANS += pCosts_->workedWeekendCost();
		// FIRST DAY (BACK TO WORK)
		ANS += pCosts_->startWorkCost(startDate);
	}

	// COSTS THAT ALWAYS APPLY
	// PREFERENCES + REDCOST OF EACH SHIFT
	for(int i=0; i<size; i++){
		int day = startDate + i, shift = allowedShortSuccBySize_[size][succId][i];
		ANS += preferencesCosts_[ day ][ shift ];
		ANS += pCosts_->dayShiftWorkCost( day, shift-1 );
	}

	return ANS;
}



//--------------------------------------------
//
// Functions for the costs
//
//--------------------------------------------

// Updates the costs depending on the reduced costs given for the nurse
//
void SubProblem::updateArcCosts(){

	priceShortSucc();

	// A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
	//
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=CDMin_-1; k<nDays_; k++)
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				int a = arcsFromSource_[s][k][n];
				double c = arcCostBestShortSuccCDMin_[s][k][n];
				updateCost( a , c );
				shortSuccCDMinIdFromArc_.insert(pair<int,int>( arcsFromSource_[s][k][n], idBestShortSuccCDMin_[s][k][n]));

				if(idBestShortSuccCDMin_[s][k][n] > -1){
					//cout << "# -X-  (" << pScenario_->intToShift_[s].at(0) << "," << k << "," << n << ")  <-  (id=" << idBestShortSuccCDMin_[s][k][n] << ") ";
					//vector<int> succ = allowedShortSuccBySize_[CDMin_][ idBestShortSuccCDMin_[s][k][n] ];
					//for(int i=0; i<succ.size(); i++) cout << pScenario_->intToShift_[succ[i]].at(0);
					//cout << endl;
				}

			}

	// B. ARCS : SHIFT_TO_NEWSHIFT [baseCost = 0]
	//
	for(int s1=1; s1<pScenario_->nbShifts_; s1++)
		for(int s2=1; s2<pScenario_->nbShifts_; s2++)
			for(int k=CDMin_-1; k<nDays_-1; k++){
				int a = arcsShiftToNewShift_[s1][s2][k];
				if(a > 0){
					double c = arcBaseCost_[a];
					c += preferencesCosts_[k+1][s2] ;
					c += pCosts_->dayShiftWorkCost(k+1,s2-1);
					c += Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
					updateCost( a , c );
				}
			}

	// C. ARCS : SHIFT_TO_SAMESHIFT [baseCost = 0]
	//
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=CDMin_-1; k<nDays_-1; k++)
			for(int n=1; n<maxvalConsByShift_[s]; n++){
				int a = arcsShiftToSameShift_[s][k][n];
				double c = arcBaseCost_[a];
				c += preferencesCosts_[k+1][s] ;
				c += pCosts_->dayShiftWorkCost(k+1,s-1);
				c += Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
				updateCost( a , c );
			}

	// D. ARCS : SHIFT_TO_ENDSEQUENCE [They never change]

	// E. ARCS : REPEATSHIFT [baseCost contains consecutive shift cost]
	//
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=CDMin_-1; k<nDays_-1; k++){
			int a = arcsRepeatShift_[s][k];
			double c = arcBaseCost_[a];
			c += preferencesCosts_[k+1][s];
			c += pCosts_->dayShiftWorkCost(k+1,s-1);
			c += Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
			updateCost( a , c );
		}

	// F. ARCS : PRINCIPAL_TO_ROTSIZE [baseCost contains complete weekend constraint]
	//
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=CDMin_-1; k<nDays_; k++){
			int a = arcsPrincipalToRotsizein_[s][k];
			double c = arcBaseCost_[a];
			c += pCosts_->endWorkCost(k);
			updateCost( a , c );
		}

	// G. ARCS : ROTSIZEIN_TO_ROTSIZE [baseCost contains rotation length cost. They never change]

	// H. ARCS : ROTSIZE_TO_ROTSIZEOUT [They never change]

	// I. ARCS : ROTSIZEOUT_TO_SINK [Never changes]
}

// Returns true if the given successions contains the given shift
//
bool SubProblem::succContainsDayShift(int size, int succId, int startDate, int thatDay, int thatShift){
	return	thatDay >= startDate
			and thatDay < startDate + size
			and allowedShortSuccBySize_[ size ][ succId ][ thatDay-startDate ] == thatShift;
}

// For tests, must be able to randomly generate costs
//
void SubProblem::generateRandomCosts(double minVal, double maxVal){

	vector<vector<double> > randomWorkCosts = Tools::randomDoubleVector2D(nDays_, pScenario_->nbShifts_, minVal, maxVal);
	vector<double> randomStartWorkCosts = Tools::randomDoubleVector(nDays_, minVal, maxVal);
	vector<double> randomEndWorkCosts = Tools::randomDoubleVector(nDays_, minVal, maxVal);
	double randomWorkedWeekendCost = (maxVal - minVal) * ( (double)rand() / (double)RAND_MAX ) + minVal;

	pCosts_ = new Costs(randomWorkCosts, randomStartWorkCosts, randomEndWorkCosts, randomWorkedWeekendCost);

}



//--------------------------------------------
//
// Functions to forbid / authorize arcs and nodes
//
//--------------------------------------------


// Returns true if the succession succ starting on day k does not violate any forbidden day-shift
//
bool SubProblem::canSuccStartHere(vector<int> succ, int firstDay){
	//std::cout << "# Checking " << firstDay << "-";
	//for(int i=0; i<succ.size(); i++) std::cout << pScenario_->intToShift_[succ[i]].at(0);

	for(int i=0; i<succ.size(); i++){
		if( ! dayShiftStatus_[firstDay+i][succ[i]] ){
			//std::cout << "   impossible because of " << (firstDay+i) << "-" << pScenario_->intToShift_[succ[i]].at(0) << endl;
			return false;
		}
	}
	//std::cout << "   OK" << endl;
	return true;
}

// Forbids the nodes that correspond to forbidden shifts
//
void SubProblem::forbid(set<pair<int,int> > forbiddenDayShifts){
	for(pair<int,int> p : forbiddenDayShifts){
		//std::cout << "# Trying to forbid " << p.first << "-" << pScenario_->intToShift_[p.second].at(0) << endl;
		forbidDayShift(p.first,p.second);
	}
}

// Forbid an arc
//
void SubProblem::forbidArc(int a){
	if(!isArcForbidden(a)){
		arcStatus_[a] = false;
		updateTime(a,MAX_TIME);
	}
}

// Forbid a node
//
void SubProblem::forbidNode(int v){
	if(!isNodeForbidden(v)){
		nodeStatus_[v] = false;
		updateLat(v,0);
	}

}

// Authorize an arc
//
void SubProblem::authorizeArc(int a){
	if(isArcForbidden(a)){
		arcStatus_[a] = true;
		updateTime(a,normalTravelTime(a));
	}
}

// Authorize a node
//
void SubProblem::authorizeNode(int v){
	nodeStatus_[v] = true;
	int lat = maxRotationLength_;
	if(nodeType(v) == ROTATION_LENGTH) lat = mapAntecedent(rotationLengthNodes_, v);
	updateLat(v,lat);
}

// Given the arc type, returns the normal travel time (when authorized)
//
int SubProblem::normalTravelTime(int a){
	ArcType atype = arcType(a);
	if(atype == SOURCE_TO_PRINCIPAL) return CDMin_;
	else if(atype == SHIFT_TO_NEWSHIFT or atype == SHIFT_TO_SAMESHIFT or atype == REPEATSHIFT) return 1;
	else return 0;
}

// Forbids a day-shift couple : Forbid the nodes + mark (day,shift) as forbidden for the short rotation pricer
//
void SubProblem::forbidDayShift(int k, int s){
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = false;
	// Forbid arcs from principal network corresponding to that day-shift only if k >= CDMin_
	if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[s]; n++){
			forbidNode( principalNetworkNodes_[s][k][n] );
		}
	}
}

// (re)Authorizes the day-shift couple BUT does not take it into account in the short rotation pricer (too complicated, will be called in the next solve() anyway)
void SubProblem::authorizeDayShift(int k, int s){
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = true;
	// Authorize arcs from principal network corresponding to that day-shift
	if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[s]; n++)
			authorizeNode( principalNetworkNodes_[s][k][n] );
	}
}

// Reset all authorizations to true
//
void SubProblem::resetAuthorizations(){
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=0; k<nDays_; k++)
			authorizeDayShift(k,s);
}

// Generate random forbidden shifts
set< pair<int,int> > SubProblem::randomForbiddenShifts(int nbForbidden){
	set< pair<int,int> > ans;
	for(int f=0; f<nbForbidden; f++){
		int k = nDays_ * ( (double)rand() / (double)RAND_MAX );
		int s = (pScenario_->nbShifts_ - 1) * ( (double)rand() / (double)RAND_MAX ) + 1;
		ans.insert(pair<int,int>(k,s));
	}
	return ans;
}


//----------------------------------------------------------------
//
// Utilities functions
//
//----------------------------------------------------------------

// Returns the key that corresponds to the given value; -1 otherwise.
int SubProblem::mapAntecedent(map<int,int> m, int val){
	for(map<int,int>::iterator it = m.begin(); it != m.end(); ++it)
		if(it->second == val) return it->first;
	return -1;
}



//--------------------------------------------
//
// PRINT FUNCTIONS
//
//--------------------------------------------

// Print the graph
void SubProblem::printGraph(){

	// TITLE
	std::cout << "# " << std::endl;
	std::cout << "# GRAPH OF THE SUBPROBLEM " << std::endl;
	std::cout << "# " << std::endl;

	// THE NODES
	//
	std::cout << "#   NODES (" << nNodes_ << ")" << std::endl;
	for(int v=0; v<nNodes_; v++) std::cout << printNode(v) << std::endl;
	std::cout << "# " << std::endl;

	// THE ARCS
	//
	std::cout << "#   ARCS (" << nArcs_ << "]" << std::endl;
	for(int a=0; a<nArcs_; a++) std::cout << printArc(a) << std::endl;
	std::cout << "# " << std::endl;

	// SUMMARY
	//
	std::cout << printSummaryOfGraph();



}

// Prints the line of a node
string SubProblem::printNode(int v){
	stringstream rep;
	rep << "# NODE  " << v << " \t" << nodeTypeName[nodeType(v)] << " \t[" << nodeEat(v) << " " << nodeLat(v) << "] \t" << shortNameNode(v);
	return rep.str();
}

// Prints the line of an arc
string SubProblem::printArc(int a){
	stringstream rep;
	rep << "# ARC   " << a << " \t" << arcTypeName[arcType(a)] << " \t";
	rep << "(" << arcOrigin(a) << "," << arcDestination(a) << ") \t" << "c= " << arcCost(a) ;
	arcCost(a) < 10000 ? rep << "     " : rep << " ";
	rep << "\tt=" << arcLength(a);
	rep << " \t[" << shortNameNode(arcOrigin(a)) << "] -> [" << shortNameNode(arcDestination(a)) << "]";
	return rep.str();

}

// Short name for a node
string SubProblem::shortNameNode(int v){

	stringstream rep;
	NodeType type_v = get( &Vertex_Properties::type, g_)[v];

	if(type_v == SOURCE_NODE){
		rep << "SOURCE";
	}

	else if (type_v == PRINCIPAL_NETWORK){
		int k = principalToDay_.at(v);
		int cons = principalToCons_.at(v);
		rep << (pScenario_->intToShift_[principalToShift_.at(v)])[0] << "-" << k << "-" << cons;
	}

	else if (type_v == ROTATION_LENGTH_ENTRANCE){
		rep << "LEN_IN";
	}

	else if (type_v == ROTATION_LENGTH){
		rep << "LEN_<=" << nodeLat(v);
	}

	else if (type_v == SINK_NODE){
		rep << "SINK";
	}

	else{
		rep << "NONE";
	}

	return rep.str();
}

// Summary of the graph
string SubProblem::printSummaryOfGraph(){
	stringstream rep;
	map<NodeType,int> nNodesPerType;
	map<ArcType,int> nArcsPerType;
	rep << "# +------------------+" << std::endl;
	rep << "# | SUBPROBLEM GRAPH |" << std::endl;
	rep << "# +------------------+" << std::endl;
	rep << "# " << std::endl;
	rep << "#     [ " << nDays_ << " days, " << (pScenario_->nbShifts_-1) << " shifts ]" << std::endl;
	rep << "# " << std::endl;
	// COUNT THE NODES
	for(int t = SOURCE_NODE; t!=NONE_NODE; t++){
	   NodeType ty = static_cast<NodeType>(t);
	   nNodesPerType.insert(pair<NodeType,int>(ty,0));
	}
	for(int v=0; v<nNodes_; v++) nNodesPerType.at(nodeType(v))++;
	// DISPLAY NODES
	rep << "#     -------------------------" << std::endl;
	rep << "#   > NODES                    " << std::endl;
	rep << "#     -------------------------" << std::endl;
	for(int t = SOURCE_NODE; t!=NONE_NODE; t++){
	   NodeType ty = static_cast<NodeType>(t);
	   rep << "#        " << nodeTypeName[ty] << "      " <<nNodesPerType.at(ty) << std::endl;
	}
	rep << "#     -------------------------" << std::endl;
	rep << "#        TOTAL            " << nNodes_ << std::endl;
	rep << "#     -------------------------" << std::endl;
	rep << "# " << std::endl;
	rep << "# " << std::endl;

	// COUNT THE ARCS
	for(int t = SOURCE_TO_PRINCIPAL; t!=NONE_ARC; t++){
	   ArcType ty = static_cast<ArcType>(t);
	   nArcsPerType.insert(pair<ArcType,int>(ty,0));
	}
	for(int a=0; a<nArcs_; a++) nArcsPerType.at(arcType(a))++;
	// DISPLAY ARCS
	rep << "#     -------------------------" << std::endl;
	rep << "#   > ARCS                     " << std::endl;
	rep << "#     -------------------------" << std::endl;
	for(int t = SOURCE_TO_PRINCIPAL; t!=NONE_ARC; t++){
		ArcType ty = static_cast<ArcType>(t);
		rep << "#        " << arcTypeName[ty] << "  " <<nArcsPerType.at(ty) << std::endl;
	}
	rep << "#     -------------------------" << std::endl;
	rep << "#        TOTAL            " << nArcs_ << std::endl;
	rep << "#     -------------------------" << std::endl;

	return rep.str();
}

// Summary of the short successions generated
void SubProblem::printShortSucc(){
	std::cout << "#   +------------+" << std::endl;
	std::cout << "#   | CD_min = " << CDMin_ << std::endl;
	std::cout << "#   +------------+" << std::endl;
	int nShortSucc = 0;
	int nShortRot = 0;
	for(int i=0; i<allowedShortSuccBySize_.size(); i++){
		vector2D v2 = allowedShortSuccBySize_[i];
		std::cout << "#   | " << v2.size()  << " short successions of size " << i << std::endl;
		nShortSucc += v2.size();
		nShortRot += v2.size() * (nDays_-i+1);

	}
	std::cout << "#   +------------+" << std::endl;
	std::cout << "#   | " << nShortSucc << std::endl;
	std::cout << "#   +------------+" << std::endl;
}

// Print the path (arcs, nodes, cost of each arc in the current network, etc.)
//
void SubProblem::printPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont ressource){

	// The successive nodes, and corresponding arc costs / time
	//
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		std::cout << "# \t| [ " << shortNameNode(source( path[j], g_ )) << " ]";
		std::cout << "\t\tLength:" << arcCost(a) << "\t\tTime:" << arcLength(a);
		std::cout << std::endl;
	}

	// Last node and total
	//
	std::cout << "# \t| [" << shortNameNode(sinkNode_) << "]" << std::endl;
	std::cout << "# \t| ~TOTAL~   \t\tLength: " << ressource.cost << "\t\tTime: " << ressource.time << std::endl;
	std::cout << "# \t| " << std::endl;
	std::cout << "# \t| Rotation: |";

	// Print it "as a rotation"
	//
	int k=0;
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		int origin = boost::source( path[j], g_ );
		int destin = boost::target( path[j], g_ );
		if(origin == sourceNode_){
			int firstDay =  principalToDay_[destin] - CDMin_ + 1;
			while(k<firstDay){
				std::cout << " |";
				k++;
			}
			int succId = shortSuccCDMinIdFromArc_[a];
			for(int s: allowedShortSuccBySize_[CDMin_][ shortSuccCDMinIdFromArc_[a] ]){
				std::cout << pScenario_->intToShift_[s].at(0) << "|";
				k++;
			}
		}
		else if(allNodesTypes_[origin] == PRINCIPAL_NETWORK and k == principalToDay_[origin]) {
			std::cout << pScenario_->intToShift_[principalToShift_[origin]].at(0) << "|";
			k++;
		}
	}
	while(k < nDays_){
		std::cout << " |";
		k++;
	}
	std::cout << std::endl;
}

// Print a given rotation
void SubProblem::printRotation(Rotation rot){

	std::cout << "# \t| ROTATION:" << "  cost=" << rot.cost_ << "  dualCost=" << rot.dualCost_ << "  firstDay=" << rot.firstDay_ << "  length=" << rot.length_ << std::endl;
	std::cout << "# \t            |";
	vector<int> allTasks (nDays_);
	for(map<int,int>::iterator itTask = rot.shifts_.begin(); itTask != rot.shifts_.end(); ++itTask)
		allTasks[itTask->first] = itTask->second;
	for(int i=0; i<allTasks.size(); i++){
		if(allTasks[i] < 1) std::cout << " |";
		else std::cout << pScenario_->intToShift_[allTasks[i]].at(0) << "|";
	}
	std::cout << std::endl;
}

// Prints all rotations in the current list
void SubProblem::printAllRotations(){
	std::cout << "# HERE ARE ALL " << nPaths_ << " ROTATIONS OF THE CURRENT SOLUTION LIST :" << std::endl;
	for(Rotation r : theRotations_){
		printRotation(r);
	}
	std::cout << "# " << endl;
}

// Print the list of currently forbidden day and shifts
void SubProblem::printForbiddenDayShift(){
	std::cout << "# List of currently forbidden day-shift pairs :";
	bool anyForbidden = false;
	for(int k=0; k<nDays_; k++){
		bool alreadyStarted=false;
		for(int s=1; s<pScenario_->nbShifts_; s++){
			if(isDayShiftForbidden(k,s)){
				anyForbidden = true;
				if(!alreadyStarted){
					std::cout << std::endl << "#      | Day " << k << " :";
					alreadyStarted = true;
				}
				std::cout << " " << pScenario_->intToShift_[s].at(0);
			}
		}
	}
	if(!anyForbidden) std::cout << " NONE";
	std::cout << std::endl;
}

// Prints all active pairs ( arcFromSource - corresponding short successions )
void SubProblem::printShortArcs(){
	for(int s=1; s<pScenario_->nbShifts_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				if(!isArcForbidden(arcsFromSource_[s][k][n])){
					int v = principalNetworkNodes_[s][k][n];
					int succId = idBestShortSuccCDMin_[s][k][n];
					vector<int> succ = allowedShortSuccBySize_[CDMin_][succId];
					std::cout << "# " << shortNameNode(v) << " <- (id=" << succId << ") ";
					for(int i=0; i<succ.size(); i++){
						std::cout << pScenario_->intToShift_[succ[i]].at(0);
					}
					std::cout << std::endl;
				}
			}
		}
	}
}









/*************************************************************************
 * A GARDER AU CAS OU COMME EXEMPLE POUR CERTAINES FONCTIONS / SYNTAXES. *
 *************************************************************************/

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
			std::allocator< boost::r_c_shortest_paths_label<Graph, spp_no_rc_res_cont> >(),
			boost::default_r_c_shortest_paths_visitor() );

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
