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


static int MAX_COST = 99999;
static int MAX_TIME = 99999;


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

SubProblem::SubProblem(Scenario * scenario, int nbDays, const Contract * contract, vector<State>* pInitState):
					pScenario_(scenario), nDays_(nbDays), pContract_ (contract),
					CDMin_(contract->minConsDaysWork_), maxRotationLength_(nbDays) {

	init(pInitState);

	initShortSuccessions();

	createNodes();
	createArcs();

	// Set all arc and node status to authorized
	for(int v=0; v<nNodes_; v++) nodeStatus_.push_back(true);
	for(int a=0; a<nArcs_; a++) arcStatus_.push_back(true);
	for(int k=0; k<nDays_; k++){
		vector<bool> v;
		for(int s=0; s<pScenario_->nbShiftsType_; s++)
			v.push_back(true);
		dayShiftStatus_.push_back(v);
	}
	for(int k=0; k<nDays_; k++) startingDayStatus_.push_back(true);

	nPathsMin_ = 0;

	//std::cout << "# A new subproblem has been created for contract " << contract->name_ << std::endl;

	//printGraph();
	//printShortSucc();

	timeInS_ = new Tools::Timer(); timeInS_->init();
	timeInNL_ = new Tools::Timer(); timeInNL_->init();

}

SubProblem::~SubProblem(){}

// Initialization function
void SubProblem::init(vector<State>* pInitState){

	// Maximum number of consecutive days worked by a nurse ending at day -1
	//
	maxOngoingDaysWorked_ = 0;
	for(unsigned int i=0; i<pInitState->size(); i++){
		maxOngoingDaysWorked_ = max( (pInitState->at(i)).consDaysWorked_, maxOngoingDaysWorked_ );
	}

	// Initialization of isUnlimited_ and nLevelsByShift_
	//
	vector<bool> v; isUnlimited_ = v; isUnlimited_.push_back(false);
	vector<int> w; maxvalConsByShift_ = w; maxvalConsByShift_.push_back(0);

	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++){
		isUnlimited_.push_back(
				       pScenario_->maxConsShiftsOfTypeOf(sh) >= nDays_ + maxOngoingDaysWorked_
				       or pScenario_->maxConsShiftsOfTypeOf(sh) >= NB_SHIFT_UNLIMITED
		);
		int nl = isUnlimited_[sh] ? pScenario_->minConsShiftsOfTypeOf(sh) : pScenario_->maxConsShiftsOfTypeOf(sh);
		maxvalConsByShift_.push_back( nl );
	}



	// id and arcCost of best succession (given a triplet s,k,n)
	idBestShortSuccCDMin_.clear();
	arcCostBestShortSuccCDMin_.clear();
	for(int s=0; s<pScenario_->nbShiftsType_; s++){
		vector2D v2; vector<vector<double> > w2;
		int n = maxvalConsByShift_[s]+1;
		Tools::initVector2D(&v2, nDays_, n);
		idBestShortSuccCDMin_.push_back(v2);
		Tools::initDoubleVector2D(&w2, nDays_, n);
		arcCostBestShortSuccCDMin_.push_back(w2);
	}

	preferencesCosts_.clear();
	for(int k=0; k<nDays_; k++){
		vector<double> v;
		// for(int s=0; s<pScenario_->nbShiftsType_; s++) v.push_back(0);
		for(int s=0; s<pScenario_->nbShifts_; s++) v.push_back(0);
		preferencesCosts_.push_back(v);
	}

	nLongFound_=0;
	nVeryShortFound_=0;
}

// Initializes the short successions. Should only be used ONCE (when creating the SubProblem).
void SubProblem::initShortSuccessions(){

	// Primary information needed
	//
	int nShiftsType= pScenario_->nbShiftsType_;

	// Put an empty list of size 0 for all data because there exists no succession of length 0/
	//
	vector2D v2; vector<int> v1,v1bis; vector<double> vd1;
	allowedShortSuccBySize_.push_back(v2);
	lastShiftOfShortSucc_.push_back(v1);
	nLastShiftOfShortSucc_.push_back(v1bis);
	baseArcCostOfShortSucc_.push_back(vd1);

	// Initialize the other way round
	for(int s=0; s<nShiftsType; s++){
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
			for(int s=1; s<nShiftsType; s++){
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
			for(unsigned int i=0; i<allowedShortSuccBySize_[c-1].size(); i++){
				vector<int> succ (allowedShortSuccBySize_[c-1][i]);
				int lastSh = succ[succ.size()-1];
				int nLast = nLastShiftOfShortSucc_[c-1][i];
				double cost = baseArcCostOfShortSucc_[c-1][i];
				// For each possible new shift s.t. the succession is allowed
				for(int newSh=1; newSh<nShiftsType; newSh++){
					if(! pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh)){

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
			for(unsigned int i=0; i<allowedShortSuccBySize_[c-1].size(); i++){
				vector<int> succ (allowedShortSuccBySize_[c-1][i]);
				int lastSh = succ[succ.size()-1];
				int nLast = nLastShiftOfShortSucc_[c-1][i];
				double cost = baseArcCostOfShortSucc_[c-1][i];
				// For each possible new shift s.t. the succession is allowed
				for(int newSh=1; newSh<nShiftsType; newSh++){
					if(! pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh)){

						vector<int> newSucc (succ); newSucc.push_back(newSh);		// Create Succession
						allSuccSizeC.push_back(newSucc);							// Add it to the possibilities
						lastShiftSucc.push_back(newSh);								// Record its last shift
						int newNLast = 1;
						double newCost = cost;
						if(newSh == lastSh){	// BUT : add the cost if longer than the maximum allowed
							newNLast += nLast;
							if(newNLast >= pScenario_->maxConsShiftsOfTypeOf(newSh)){
								newCost += consShiftCost(lastSh, nLast) ;
							}
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

}

// Cost function for consecutive identical shifts
//
double SubProblem::consShiftCost(int sh, int n){
  if(pScenario_->minConsShiftsOfTypeOf(sh) - n > 0) return (WEIGHT_CONS_SHIFTS * ( pScenario_->minConsShiftsOfTypeOf(sh) - n ) );
  if(n - pScenario_->maxConsShiftsOfTypeOf(sh) > 0) return (WEIGHT_CONS_SHIFTS * ( n - pScenario_->maxConsShiftsOfTypeOf(sh) ) );
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
bool SubProblem::solve(LiveNurse* nurse, DualCosts * costs, SubproblemParam param, set<pair<int,int> > forbiddenDayShifts,
		set<int> forbiddenStartingDays, bool optimality, double redCostBound){


	bestReducedCost_ = 0;
	param_ = param;

	// Maximum rotations length: update the bounds on the nodes if needed
	//
	int newMaxRotationLength = min(nDays_+maxOngoingDaysWorked_, max(pContract_->maxConsDaysWork_, param_.maxRotationLength_));
//	int newMaxRotationLength = min(nDays_, max(pContract_->maxConsDaysWork_, maxRotationLength));
	if(newMaxRotationLength != maxRotationLength_){
		maxRotationLength_ = newMaxRotationLength;
		updatedMaxRotationLengthOnNodes();
	}

	maxReducedCostBound_ = redCostBound - EPSILON;			// Cost bound
	pLiveNurse_ = nurse;									// Reset the nurse
	pCosts_ = costs;										// Reset the costs
	resetAuthorizations();									// Reset authorizations
	resetSolutions();										// Delete all already existing solutions
	initStructuresForSolve();								// Initialize structures
	nLongFound_=0;											// Initialize number of solutions found at 0 (long rotations)
	nVeryShortFound_=0;										// Initialize number of solutions found at 0 (short rotations)
	forbid(forbiddenDayShifts);								// Forbid arcs
	forbidStartingDays(forbiddenStartingDays);				// Forbid starting days

	if(false) printContractAndPrefenrences();				// Set to true if you want to display contract + preferences (for debug)

	timeInS_->start();
	bool ANS_short = solveShortRotations();
	timeInS_->stop();
	timeInNL_->start();
	bool ANS_long = solveLongRotations(optimality);
	timeInNL_->stop();

	// Reset authorizations
	authorize(forbiddenDayShifts);
	authorizeStartingDays(forbiddenStartingDays);

	return ANS_short or ANS_long;
}

// For the short rotations, depends on the chosen option + on wether we want optimality (more important)
bool SubProblem::solveShortRotations(){
	bool ANS = false;

	// SOLVE_SHORT_NONE
	if(param_.shortRotationsStrategy_==0) {}

	// SOLVE_SHORT_ALL
	else if(param_.shortRotationsStrategy_==1){
		bool tmp_ANS = priceVeryShortRotations();
		ANS = ANS or tmp_ANS;
	}

	// SOLVE_SHORT_DAY_0_AND_LAST_ONLY
	else if(param_.shortRotationsStrategy_==2){
		bool tmp_ANS = priceVeryShortRotationsFirstDay();
		bool tmp_ANS2 = priceVeryShortRotationsLastDay();
		ANS = tmp_ANS2 or ANS or tmp_ANS;
	}

	// SOLVE_SHORT_DAY_0_ONLY
	else if(param_.shortRotationsStrategy_==3){
		bool tmp_ANS = priceVeryShortRotationsFirstDay();
		ANS = ANS or tmp_ANS;
	}

	// SOLVE_SHORT_DAY_0_ONLY
	else if(param_.shortRotationsStrategy_==4){
		bool tmp_ANS = priceVeryShortRotationsLastDay();
		ANS = ANS or tmp_ANS;
	}

	else {
		cout << "# INVALID / OBSOLETE OPTION FOR SHORT ROTATIONS" << endl;
		getchar();
		return false;
	}
	return ANS;
}

// For the long rotations, depends wether we want optimality or not
bool SubProblem::solveLongRotations(bool optimality){
	if(optimality){
		updateArcCosts();						// Update costs
		return solveLongRotationsOptimal();		// Solve shortest path problem
	}
	else
		return solveLongRotationsHeuristic();	// Generate new rotations with greedy
}

// Function called when optimal=true in the arguments of solve -> shortest path problem is to be solved
bool SubProblem::solveLongRotationsOptimal(){

	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;

	// ONE SINGLE SINK FOR ALL DAYS
	//
	if(!(param_.oneSinkNodePerLastDay_)){
		r_c_shortest_paths(
				g_,
				get( &Vertex_Properties::num, g_ ),
				get( &Arc_Properties::num, g_ ),
				sourceNode_,
				sinkNode_,
				opt_solutions_spptw,
				pareto_opt_rcs_spptw,
				spp_spptw_res_cont (0,0),
				ref_spptw(),
				dominance_spptw(),
				std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
				boost::default_r_c_shortest_paths_visitor() );
		return addRotationsFromPaths(opt_solutions_spptw, pareto_opt_rcs_spptw);
	}

	// ONE SINK FOR EACH DAY
	//
	else {
		std::vector<boost::graph_traits<Graph>::vertex_descriptor> allSinks;
		for(int k=CDMin_-1; k<nDays_; k++){
			allSinks.push_back( sinkNodesByDay_[k] );
		}
		r_c_shortest_paths_several_sinks(
				g_,
				get( &Vertex_Properties::num, g_ ),
				get( &Arc_Properties::num, g_ ),
				sourceNode_,
				allSinks,
				opt_solutions_spptw,
				pareto_opt_rcs_spptw,
				spp_spptw_res_cont (0,0),
				ref_spptw(),
				dominance_spptw(),
				std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
				boost::default_r_c_shortest_paths_visitor() );
		return addRotationsFromPaths(opt_solutions_spptw, pareto_opt_rcs_spptw);
	}
}

// Transforms the solutions found into proper rotations.
//
bool SubProblem::addRotationsFromPaths(vector< vector< boost::graph_traits<Graph>::edge_descriptor> > paths, vector<spp_spptw_res_cont> resources){
	int nFound = 0;
	// For each path of the list, record the corresponding rotation (if negativeOnly=true, do it only if the dualCost < 0)
	for(unsigned int p=0; p < paths.size(); ++p){

		// 1. Check if it is valid
		bool b_is_a_path_at_all = false;
		bool b_feasible = false;
		bool b_correctly_extended = false;
		spp_spptw_res_cont actual_final_resource_levels( 0, 0 );
		boost::graph_traits<Graph>::edge_descriptor ed_last_extended_arc;
		check_r_c_path( g_,
				paths[p],
				spp_spptw_res_cont( 0, 0 ),
				true,
				resources[p],
				actual_final_resource_levels,
				ref_spptw(),
				b_is_a_path_at_all,
				b_feasible,
				b_correctly_extended,
				ed_last_extended_arc );
		if( b_is_a_path_at_all && b_feasible && b_correctly_extended )
		{
			if(resources[p].cost < maxReducedCostBound_){
				Rotation rot = rotationFromPath(paths[p], resources[p]);
				theRotations_.push_back(rot);
				nPaths_ ++;
				nLongFound_++;
				nFound ++;
				bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
			}
		} else {
			if( !b_is_a_path_at_all )
				std::cout << "Not a path." << std::endl;
			if( !b_feasible )
				std::cout << "Not a feasible path." << std::endl;
			if( !b_correctly_extended )
				std::cout << "Not correctly extended." << std::endl;
		}

//		if(resources[p].cost < maxReducedCostBound_){
//			Rotation rot = rotationFromPath(paths[p], resources[p]);
//			theRotations_.push_back(rot);
//			nPaths_ ++;
//			nLongFound_++;
//			nFound ++;
//		}
	}
	//printAllRotations();
	//	std::cout << "# -> " << nFound << std::endl;
	return (nFound > 0);
}

// Adds a rotation made from the given path to the current list of answers and increases their counter
//
Rotation SubProblem::rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource){

	int firstDay = -1;
	vector<int> shiftSuccession;

	// All arcs are consecutively considered
	//
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		ArcType aType = arcType(a);
		int destin = boost::target( path[j], g_ );

		// A. Arc from source (equivalent to short rotation
		if(aType == SOURCE_TO_PRINCIPAL){
			firstDay =  principalToDay_[destin] - CDMin_ + 1;
			for(int s: static_cast<vector <int> >( allowedShortSuccBySize_[CDMin_][ shortSuccCDMinIdFromArc_.at(a) ] )){
				shiftSuccession.push_back(s);
			}
		}

		// B. Arc to a new day
		else if(aType == SHIFT_TO_NEWSHIFT or aType == SHIFT_TO_SAMESHIFT or aType == REPEATSHIFT){
			shiftSuccession.push_back( principalToShift_[destin] );
		}
	}

	Rotation rot (firstDay, shiftSuccession, pLiveNurse_->id_, MAX_COST, resource.cost);
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
	int CD_max = pContract_->maxConsDaysWork_;									// Maximum consecutive days worked for free
	int nShiftsType= pScenario_->nbShiftsType_;											// Number of different shifts

	// INITIALIZATION
	nNodes_ = 0;
	initNodesStructures();

	// 1. SOURCE NODE
	//
	sourceNode_ = nNodes_;
	addSingleNode(SOURCE_NODE, 0, maxRotationLength_);

	// 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
	//
	for(int sh=1; sh<nShiftsType; sh++){											// For each possible worked shift
		for(int k=CDMin_-1; k<nDays_; k++){											// For each date
			for(int cons=1; cons<=maxvalConsByShift_[sh]; cons++){					// For each level of network
				addNodeToPrincipalNetwork(sh, k, cons);								// Add a node to the principal network
			}
		}
	}

	// 3. ROTATION LENGTH CHECK
	//
	// For each of the days, do a rotation-length-checker
	for(int k=0; k<nDays_; k++){
		rotationLengthEntrance_.push_back(nNodes_);									// One node for the entrance in subnetwork per day
		addSingleNode(ROTATION_LENGTH_ENTRANCE, 0, maxRotationLength_);
		map<int,int> checkNodesForThatDay;
		// Check nodes
		for(int l=CD_max; l<=maxRotationLength_; l++){								// Check nodes: from CD_max (longest free) to maximum rotation length, for each day
			checkNodesForThatDay.insert(pair<int,int>(l,nNodes_));
			rotationLengthNodesLAT_.insert(pair<int,int>(nNodes_,l));
			addSingleNode(ROTATION_LENGTH, 0, l);
		}
		rotationLengthNodes_.push_back(checkNodesForThatDay);
		// Sink day
		sinkNodesByDay_.push_back(nNodes_);											// Daily sink node
		addSingleNode(SINK_DAY, 0, maxRotationLength_);

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
	int nShiftsType= pScenario_->nbShiftsType_;				// Number of different shifts

	allNodesTypes_.clear();
	principalNetworkNodes_.clear();
	principalToShift_.clear();
	principalToDay_.clear();
	principalToCons_.clear();
	rotationLengthEntrance_.clear();
	rotationLengthNodes_.clear();
	rotationLengthNodesLAT_.clear();
	sinkNodesByDay_.clear();

	// All nodes
	//
	nNodes_ = 0;
	vector<NodeType> v; allNodesTypes_ = v;

	// Principal networks
	//
	for(int sh=0; sh<nShiftsType; sh++){
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
void SubProblem::addSingleArc(int o, int d, double baseCost, int t, ArcType type){
	boost::graph_traits< Graph>::edge_descriptor e = (add_edge( o, d, Arc_Properties( nArcs_, type, baseCost, t ), g_ )).first;
	arcsDescriptors_.push_back(e);
	allArcsTypes_.push_back(type);
	arcBaseCost_.push_back(baseCost);
	nArcs_++;

	if(nodeType(o) == PRINCIPAL_NETWORK
			and nodeType(d) == PRINCIPAL_NETWORK
			and principalToDay_[o] < principalToDay_[d]-1){
		printArc(nArcs_-1);
		getchar();
	}

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
	arcsRotsizeinToRotsizeDay_.clear();
	arcsRotsizeToRotsizeoutDay_.clear();
	// // VECTORS 3 D
	// for(int s=0; s<pScenario_->nbShiftsType_; s++){
	// 	vector2D v2, w2, x2;
	// 	Tools::initVector2D(&v2, nDays_, maxvalConsByShift_[s]+1); arcsFromSource_.push_back(v2);
	// 	Tools::initVector2D(&w2, nDays_, maxvalConsByShift_[s]+1); arcsShiftToSameShift_.push_back(w2);
	// 	Tools::initVector2D(&x2, nDays_, maxvalConsByShift_[s]+1); arcsShiftToEndsequence_.push_back(x2);
	// }
	// Tools::initVector3D(&arcsShiftToNewShift_, pScenario_->nbShiftsType_, pScenario_->nbShiftsType_, nDays_);
	// // VECTORS 2 D
	// Tools::initVector2D(&arcsRepeatShift_, pScenario_->nbShiftsType_, nDays_);
	// Tools::initVector2D(&arcsPrincipalToRotsizein_, pScenario_->nbShiftsType_, nDays_, -1);

	// VECTORS 3 D
	for(int sh=0; sh<pScenario_->nbShiftsType_; sh++){
	  vector2D v2, w2, x2;
	  vector3D y2, z2;
	  unsigned int nShifts = pScenario_->shiftTypeIDToShiftID_[sh].size();
	  
	  Tools::initVector2D(&v2, nDays_, maxvalConsByShift_[sh]+1);           arcsFromSource_.push_back(v2);
	  Tools::initVector2D(&w2, nDays_, maxvalConsByShift_[sh]+1);           arcsShiftToEndsequence_.push_back(w2);
	  Tools::initVector2D(&x2, nDays_, nShifts);                            arcsRepeatShift_.push_back(x2);
	  Tools::initVector3D(&y2, nDays_, maxvalConsByShift_[sh]+1, nShifts);  arcsShiftToSameShift_.push_back(y2);
	  Tools::initVector3D(&z2, pScenario_->nbShiftsType_, nDays_, nShifts); arcsShiftToNewShift_.push_back(z2);
	}
	// VECTORS 2 D
	Tools::initVector2D(&arcsPrincipalToRotsizein_, pScenario_->nbShiftsType_, nDays_, -1);
}

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void SubProblem::createArcsSourceToPrincipal(){

	// DATA
	int nShiftsType = pScenario_->nbShiftsType_;
	int origin, destin;

	for(int sh=1; sh<nShiftsType; sh++){
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
	int nShiftsType = pScenario_->nbShiftsType_;
	int origin, destin;

	// FOR EACH OF THE SUBNETWORKS AND EACH OF THE DAYS
	//
	for (int sh=1; sh<nShiftsType; sh++){
	  unsigned int nShifts = pScenario_->shiftTypeIDToShiftID_[sh].size();
		for(int k=CDMin_-1; k<nDays_-1; k++){

			//   1. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS NOT REACHED YET
			//
			for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons ++){
				origin = principalNetworkNodes_[sh][k][nCons];
				destin = principalNetworkNodes_[sh][k+1][nCons+1];

				for (unsigned int s = 0; s < nShifts; s++) {
				  int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
				  arcsShiftToSameShift_[sh][k][nCons][s] = nArcs_;
				  // addSingleArc(origin, destin, 0, 1, SHIFT_TO_SAMESHIFT);
				  addSingleArc(origin, destin, 0, pScenario_->hoursToWork_[shiftID], SHIFT_TO_SAMESHIFT);
				}
			}

			//   2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
			//
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = principalNetworkNodes_[sh][k+1][maxvalConsByShift_[sh]];
			double cost = isUnlimited(sh) ? 0 : WEIGHT_CONS_SHIFTS;

			for (unsigned int s = 0; s < nShifts; s++) {
			  int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
			  arcsRepeatShift_[sh][k][s] = nArcs_;
			// addSingleArc(origin, destin, cost, 1, REPEATSHIFT);
			  addSingleArc(origin, destin, cost, pScenario_->hoursToWork_[shiftID], REPEATSHIFT);
			}
			
			// 3. WORK ONE MORE DAY ON A DIFFERENT SHIFT (CHANGE SUBNETWORK)
			//
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			for(int newSh=1; newSh<nShiftsType; newSh++){
			  if(newSh != sh and ! pScenario_->isForbiddenSuccessorShiftType_ShiftType(newSh,sh)){
			    unsigned int nNewShifts = pScenario_->shiftTypeIDToShiftID_[newSh].size();
			    int destin = principalNetworkNodes_[newSh][k+1][1];

			    for (unsigned int s = 0; s < nNewShifts; s++) {
			      int  newShiftID = pScenario_-> shiftTypeIDToShiftID_[newSh][s];
			      arcsShiftToNewShift_[sh][newSh][k][s] = nArcs_;
			      // addSingleArc(origin, destin, 0, 1, SHIFT_TO_NEWSHIFT);
			      addSingleArc(origin, destin, 0, pScenario_->hoursToWork_[newShiftID], SHIFT_TO_NEWSHIFT);
			    }
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
			addSingleArc(origin, destin, 0, 0, SHIFT_TO_ENDSEQUENCE);
		}
	}
}

// Create all arcs that involve the rotation size checking subnetwork (incoming, internal, and exiting that subnetwork)
void SubProblem::createArcsAllRotationSize(){

	int nShiftsType = pScenario_->nbShiftsType_;
	int origin, destin;

	// 1. ALL INCOMING ARCS
	//
	for(int sh=1; sh<nShiftsType; sh++){											// For all shifts
		for(int k=CDMin_-1; k<nDays_; k++){										// For all days
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = rotationLengthEntrance_[k];
			arcsPrincipalToRotsizein_[sh][k] = nArcs_;
			addSingleArc(origin, destin, 0, 0, PRINCIPAL_TO_ROTSIZE);			// Allow to stop rotation that day
		}
	}

	// 2. ALL INTERNAL ARCS
	//
	for(int k=CDMin_-1; k<nDays_; k++){

		map<int,int> rotLengthNodesForDay = rotationLengthNodes_[k];
		map<int,int> arcsRotsizeinToRotsize;
		map<int,int> arcsRotsizeToRotsizeout;
		for(map<int,int>::iterator itRLN = rotLengthNodesForDay.begin(); itRLN != rotLengthNodesForDay.end(); ++itRLN){
			// From entrance of that day to checknode
			origin = rotationLengthEntrance_[k];
			destin = itRLN->second;
			arcsRotsizeinToRotsize.insert(pair<int,int>(itRLN->first, nArcs_));
			addSingleArc(origin, destin, consDaysCost(itRLN->first), 0, ROTSIZEIN_TO_ROTSIZE);
			// From checknode to exit of that day
			origin = itRLN->second;
			destin = sinkNodesByDay_[k];
			arcsRotsizeToRotsizeout.insert(pair<int,int>( itRLN->first, nArcs_));
			addSingleArc(origin, destin, 0, 0, ROTSIZE_TO_SINK);
		}

		arcsRotsizeinToRotsizeDay_.push_back(arcsRotsizeinToRotsize);
		arcsRotsizeToRotsizeoutDay_.push_back(arcsRotsizeToRotsizeout);

		// link all sink nodes to the main sink node
		origin = sinkNodesByDay_[k];
		destin = sinkNode_;
		addSingleArc(origin, destin, 0,0, SINKDAY_TO_SINK);

	}


}



//--------------------------------------------
//
// Functions for the pricing of the short rotations
//
//--------------------------------------------

// Initializes some cost vectors that depend on the nurse
void SubProblem::initStructuresForSolve(){

	// Start and End weekend costs
	//
	startWeekendCosts_.clear();
	endWeekendCosts_.clear();
	Tools::initDoubleVector(&startWeekendCosts_,nDays_);
	Tools::initDoubleVector(&endWeekendCosts_,nDays_);
	if(pLiveNurse_->needCompleteWeekends()){
		for(int k=0; k<nDays_; k++){
			if(Tools::isSaturday(k)) endWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
			else if(Tools::isSunday(k)) startWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
		}
	}

	// Preference costs.
	//
	for(int k=0; k<nDays_; k++)
		// for(int s=0; s<pScenario_->nbShiftsType_; s++)
		for(int s=0; s<pScenario_->nbShifts_; s++)
			preferencesCosts_[k][s] = 0;
	for(map<int,set<int> >::iterator it = pLiveNurse_->pWishesOff_->begin(); it != pLiveNurse_->pWishesOff_->end(); ++it){
		for(int s : it->second){
			preferencesCosts_[it->first][s] = WEIGHT_PREFERENCES;

		}
	}

	// id and arcCost of best succession (given a triplet s,k,n)
	idBestShortSuccCDMin_.clear();
	arcCostBestShortSuccCDMin_.clear();
	for(int s=0; s<pScenario_->nbShiftsType_; s++){
		vector2D v2; vector<vector<double> > w2;
		int n = maxvalConsByShift_[s]+1;
		Tools::initVector2D(&v2, nDays_, n, -1);
		idBestShortSuccCDMin_.push_back(v2);
		Tools::initDoubleVector2D(&w2, nDays_, n, MAX_COST);
		arcCostBestShortSuccCDMin_.push_back(w2);
	}

}

// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
//
void SubProblem::priceShortSucc(){

	map<int,int> specialArcsSuccId;
	map<int,double> specialArcsCost;

	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){

				idBestShortSuccCDMin_[s][k][n] = -1;
				arcCostBestShortSuccCDMin_[s][k][n] = MAX_COST;

				// CHECK THE ROTATIONS ONLY IF THE FIRST DAY IS ALLOWED
				if(startingDayStatus_[k-CDMin_+1]){

					for(unsigned int i=0; i<(allShortSuccCDMinByLastShiftCons_[s][n]).size(); i++){
						int curSuccId = allShortSuccCDMinByLastShiftCons_[s][n][i];
						vector<int> succ = allowedShortSuccBySize_[CDMin_][curSuccId];

						// SUCCESSION IS TAKEN INTO ACCOUNT ONLY IF IT DOES NOT VIOLATE ANY FORBIDDEN DAY-SHIFT COUPLE
						if(canSuccStartHere( succ, k-CDMin_+1 )){
							double curCost = costArcShortSucc(CDMin_, curSuccId, k-CDMin_+1);

							// ONLY CASE WHEN THE DESTINATION NODE MAY HAVE TO CHANGE:
							// 1. Start date is 0
							// 2. Size of short succession is < than the number of levels maxValByShift[s]
							// 3. Number of last shifts cons in succession is CDMin_
							// 4. The shift is the same as the last one worked by the nurse at initial state
							if(k==CDMin_-1 and CDMin_<maxvalConsByShift_[s] and n==CDMin_ and s==pLiveNurse_->pStateIni_->shift_){
								// a. Determine the destination of that arc
								int nConsWithPrev = CDMin_ + pLiveNurse_->pStateIni_->consShifts_;
								int nDestination = min( nConsWithPrev , maxvalConsByShift_[s] );
								int a = arcsFromSource_[s][k][nDestination];
								// b. Store the succession ID + the special cost for that arc
								specialArcsSuccId.insert(pair<int,int>(a,curSuccId));
								specialArcsCost.insert(pair<int,double>(a,curCost));
							}

							// OTHER CASES ("REGULAR ONES")
							else if(curCost < arcCostBestShortSuccCDMin_[s][k][n]){
								idBestShortSuccCDMin_[s][k][n] = curSuccId;
								arcCostBestShortSuccCDMin_[s][k][n] = curCost;
							}
						}
					}
				}

				// IF NO VALID SUCCESSION OR IF THE FIRST DAY IS FORBIDDEN AS A STARTING DAY, THEN FORBID THE ARC
				if(!startingDayStatus_[k-CDMin_+1] || arcCostBestShortSuccCDMin_[s][k][n] >= MAX_COST-1){
					forbidArc( arcsFromSource_[s][k][n] );
				}
			}
		}
	}

	// FOR THE SHIFTS ON THE FIRST DAY THAT EXTEND THE ONGOING WORK AT INITIAL STATE
	//
	for(map<int,int>::iterator itId = specialArcsSuccId.begin(); itId != specialArcsSuccId.end(); ++itId){
		int a = itId->first;
		int d = arcDestination(a);
		int s = principalToShift_[d];
		int k = principalToDay_[d];
		int n = principalToCons_[d];
		if(specialArcsCost.find(a) == specialArcsCost.end()){
			cout << "# Problem within pricing of some short rotations (press Enter to go on)" << endl;
			getchar();
		}
		double cost = specialArcsCost.at(a);
		authorizeArc(a);
		idBestShortSuccCDMin_[s][k][n] = itId->second;
		arcCostBestShortSuccCDMin_[s][k][n] = cost;
	}
}

// Given a short succession and a start date, returns the cost of the corresponding arc
//
double SubProblem::costArcShortSucc(int size, int succId, int startDate){
	double ANS = 0;
	vector<int> succ = allowedShortSuccBySize_[size][succId];

	// A. COST: BASE COST
	//
	ANS += baseArcCostOfShortSucc_[size][succId];


	// B. COST: SPECIAL CASE FOR THE FIRST DAY
	//
	if(startDate ==0){

		int shiftIni = pLiveNurse_->pStateIni_->shift_;
		int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
		int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

		int firstShift = succ[0];
		int nConsFirstShift = 0;
		unsigned int ii=0;
		while(ii<succ.size() and succ[ii]==firstShift){
			nConsFirstShift ++;
			ii++;
		}


		// 1. The nurse was resting: pay more only if the rest is too short
		if(shiftIni == 0){
			int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
			ANS += max(0, diffRest*WEIGHT_CONS_DAYS_OFF);
		}

		// 2. The nurse was working
		else {

			// a. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be read later
			int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
			ANS -= max(0, diffWork*WEIGHT_CONS_DAYS_WORK);

			// b. (i)   The nurse was working on a different shift: if too short, add the corresponding cost
			if(shiftIni != firstShift){
			  int diff = pScenario_->minConsShiftsOfTypeOf(shiftIni) - nConsShiftIni;
				ANS += max(0, diff*(WEIGHT_CONS_SHIFTS));
			}

			// b. (ii)  The nurse was working on the same shift AND the short rotation contains other shifts (easy case for add/subtract)
			//            - Subtract the cost due to the consecutive beginning
			//            - Subtract the cost due to the consecutive end of the initial state
			//            - Add the consecutive cost for all shifts
			else if(nConsFirstShift < CDMin_) {
			  int diffShift = nConsShiftIni - pScenario_->maxConsShiftsOfTypeOf(shiftIni);
				ANS -= max(0, diffShift*WEIGHT_CONS_SHIFTS);
				ANS -= consShiftCost(firstShift, nConsFirstShift);
				ANS += consShiftCost(firstShift, (nConsFirstShift + nConsShiftIni));
			}

			// b. (iii) The nurse was working on the same shift AND the short rotation only contains that shift (recompute the cost -> easier)
			else {
				ANS -= baseArcCostOfShortSucc_[size][succId];
				int shiftTypeIni = pScenario_->shiftIDToShiftTypeID_[shiftIni];
				if( (nConsFirstShift + nConsShiftIni) >= maxvalConsByShift_[shiftTypeIni] )
					ANS += consShiftCost(shiftIni, (nConsFirstShift + nConsShiftIni));
			}
		}
	}

	// C. COST: COMPLETE WEEKEND
	//
	ANS += startWeekendCosts_[startDate];

	// D. COST: PREFERENCES
	//
	for(int i=0; i<size; i++) ANS += preferencesCosts_[ startDate + i ][ allowedShortSuccBySize_[size][succId][i] ];



	// E. REDCOST: WEEKENDS
	//
	int nbWeekends = Tools::containsWeekend(startDate, startDate + size - 1);
	ANS -= nbWeekends * pCosts_->workedWeekendCost();

	// F. REDCOST: FIRST DAY (BACK TO WORK)
	//
	ANS -= pCosts_->startWorkCost(startDate);

	// G. REDCOST: EACH DAY/SHIFT REDUCED COST
	//
	for(int i=0; i<size; i++) ANS -= pCosts_->dayShiftWorkCost( startDate+i , allowedShortSuccBySize_[size][succId][i]-1 );



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
	shortSuccCDMinIdFromArc_.clear();
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				int a = arcsFromSource_[s][k][n];
				double c = arcCostBestShortSuccCDMin_[s][k][n];
				updateCost( a , c );
				shortSuccCDMinIdFromArc_.insert(pair<int,int>( arcsFromSource_[s][k][n], idBestShortSuccCDMin_[s][k][n]));
			}
		}

		// For all those that start on the first day, must update the travel time
		//
		for(int n=1; n<=maxvalConsByShift_[s]; n++){
			if(idBestShortSuccCDMin_[s][CDMin_-1][n] > -1){
				int a = arcsFromSource_[s][CDMin_-1][n];
				updateTime(a, (CDMin_+ pLiveNurse_->pStateIni_->consDaysWorked_ ));
			}
		}
	}

	// B. ARCS : SHIFT_TO_NEWSHIFT [baseCost = 0]
	//
	for(int s1=1; s1<pScenario_->nbShiftsType_; s1++)
		for(int s2=1; s2<pScenario_->nbShiftsType_; s2++)
		  if(s2 != s1 and ! pScenario_->isForbiddenSuccessorShiftType_ShiftType(s2,s1)){
			for(int k=CDMin_-1; k<nDays_-1; k++){
			  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[s2].size(); s++) {
				int a = arcsShiftToNewShift_[s1][s2][k][s];
				if(a > 0){
					double c = arcBaseCost_[a];
					int  shiftID = pScenario_-> shiftTypeIDToShiftID_[s2][s];
					// c += preferencesCosts_[k+1][s2] ;
				        c += preferencesCosts_[k+1][shiftID] ;
					c -= pCosts_->dayShiftWorkCost(k+1,s2-1);
					c -= Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
					updateCost( a , c );
				}
			  }
			}
		  }

	// C. ARCS : SHIFT_TO_SAMESHIFT [baseCost = 0]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=CDMin_-1; k<nDays_-1; k++)
			for(int n=1; n<maxvalConsByShift_[sh]; n++){
			    for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
				int a = arcsShiftToSameShift_[sh][k][n][s];
				double c = arcBaseCost_[a];
				int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
				// c += preferencesCosts_[k+1][sh] ;
				c += preferencesCosts_[k+1][shiftID] ;
				c -= pCosts_->dayShiftWorkCost(k+1,sh-1);
				if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
				updateCost( a , c );
			    }
			}

	// D. ARCS : SHIFT_TO_ENDSEQUENCE [They never change]

	// E. ARCS : REPEATSHIFT [baseCost contains consecutive shift cost]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=CDMin_-1; k<nDays_-1; k++){
		  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
			int a = arcsRepeatShift_[sh][k][s];
			double c = arcBaseCost_[a];
			int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
			// c += preferencesCosts_[k+1][sh];
			c += preferencesCosts_[k+1][shiftID];
			c -= pCosts_->dayShiftWorkCost(k+1,sh-1);
			if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
			updateCost( a , c );
		  }
		}

	// F. ARCS : PRINCIPAL_TO_ROTSIZE [baseCost contains complete weekend constraint]
	//
	for(int s=1; s<pScenario_->nbShiftsType_; s++)
		for(int k=CDMin_-1; k<nDays_; k++){
			int a = arcsPrincipalToRotsizein_[s][k];
			double c = arcBaseCost_[a];
			c += endWeekendCosts_[k];
			c -= pCosts_->endWorkCost(k);
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



//--------------------------------------------
//
// Functions to update the maximum length of a rotation
//
//--------------------------------------------

// Updates the maximum arrival time at all nodes so that the maximum length of a rotation is updated.
//   + Updates only
void SubProblem::updatedMaxRotationLengthOnNodes(){
	for(int v=0; v<nNodes_; v++){
		if(nodeStatus_[v]){
			if(nodeType(v) != ROTATION_LENGTH){
				updateLat(v,maxRotationLength_);
			}
		}
	}
}



//--------------------------------------------
//
// Functions to forbid / authorize arcs and nodes
//
//--------------------------------------------


// Returns true if the succession succ starting on day k does not violate any forbidden day-shift
//
bool SubProblem::canSuccStartHere(vector<int> succ, int firstDay){
	// If the starting date is forbidden, return false
	if(!(startingDayStatus_[firstDay]))
		return false;
	// If the succession with the previous shift (day -1) is not allowed
	if(firstDay==0 and pScenario_->isForbiddenSuccessorShift_Shift(succ[0],pLiveNurse_->pStateIni_->shift_))
		return false;
	// If some day-shift is forbidden...
	for(unsigned int i=0; i<succ.size(); i++){
		if( ! dayShiftStatus_[firstDay+i][succ[i]]){
			return false;
		}
	}
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

// Authorize the nodes that correspond to forbidden shifts
//
void SubProblem::authorize(set<pair<int,int> > forbiddenDayShifts){
	for(pair<int,int> p : forbiddenDayShifts){
		//std::cout << "# Trying to forbid " << p.first << "-" << pScenario_->intToShift_[p.second].at(0) << endl;
		authorizeDayShift(p.first,p.second);
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
	if(isNodeForbidden(v)){
		nodeStatus_[v] = true;
		int lat = maxRotationLength_;
		if(nodeType(v) == ROTATION_LENGTH) lat = rotationLengthNodesLAT_.at(v);
		updateLat(v,lat);
	}
}

// Given the arc type, returns the normal travel time (when authorized)
//
int SubProblem::normalTravelTime(int a){
	ArcType atype = arcType(a);
	if(atype == SOURCE_TO_PRINCIPAL){
		if(principalToDay_[arcDestination(a)] == CDMin_-1){
			return (CDMin_ + pLiveNurse_->pStateIni_->consDaysWorked_);
		} else {
			return CDMin_;
		}
	}
	else if(atype == SHIFT_TO_NEWSHIFT or atype == SHIFT_TO_SAMESHIFT or atype == REPEATSHIFT) return 1;
	else return 0;
}

// Forbids a day-shift couple : Forbid the nodes + mark (day,shift) as forbidden for the short rotation pricer
//
void SubProblem::forbidDayShift(int k, int s){
	//if rest, do nothing
	if(s==0) return;
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = false;

	int shiftType = pScenario_->shiftIDToShiftTypeID_[s];
	// Forbid arcs from principal network corresponding to that day-shift only if k >= CDMin_
	if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[shiftType]; n++){
			forbidNode( principalNetworkNodes_[shiftType][k][n] );
		}
	}
}

// (re)Authorizes the day-shift couple BUT does not take it into account in the short rotation pricer (too complicated, will be called in the next solve() anyway)
void SubProblem::authorizeDayShift(int k, int s){
	//if rest, do nothing
	if(s==0) return;
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = true;
	int shiftType = pScenario_->shiftIDToShiftTypeID_[s];
	// Authorize arcs from principal network corresponding to that day-shift
	if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[shiftType]; n++)
			authorizeNode( principalNetworkNodes_[shiftType][k][n] );
	}
}

// Forbids some starting days
//
void SubProblem::forbidStartingDays(set<int> forbiddenStartingDays){
	for(int k: forbiddenStartingDays)
		forbidStartingDay(k);
}

// Authorizes some starting days
//
void SubProblem::authorizeStartingDays(set<int> authorizedStartingDays){
	for(int k: authorizedStartingDays)
		authorizeStartingDay(k);
}

// Forbids a starting date: no rotation can now start on that day. Gives a prohibitive resource consumption on all short
// rotation arcs that correspond to rotations starting on that day.
//
void SubProblem::forbidStartingDay(int k){
	// Mark the starting day as forbidden
	startingDayStatus_[k] = false;
	// IF EVERYTHING WENT RIGHT, IT SHOULD BE ENOUGH BECAUSE IT WILL BE INDICATED DURING THE PRINCING OF THE SHORT
	// ROTATIONS !!!
}

// Authorizes a starting date: rotations may now start on that day.
//
void SubProblem::authorizeStartingDay(int k){
	// Mark the starting day as allowed
	startingDayStatus_[k] = true;
}

// Reset all authorizations to true
//
void SubProblem::resetAuthorizations(){
	for(int s=1; s<pScenario_->nbShiftsType_; s++)
		for(int k=0; k<nDays_; k++)
			authorizeDayShift(k,s);

	for(int a=0; a<nArcs_; a++)
		authorizeArc(a);
}

// Generate random forbidden shifts
set< pair<int,int> > SubProblem::randomForbiddenShifts(int nbForbidden){
	set< pair<int,int> > ans;
	for(int f=0; f<nbForbidden; f++){
		int k = Tools::randomInt(0, nDays_ - 1);
		int s = Tools::randomInt(1, pScenario_->nbShiftsType_ - 1);
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





//----------------------------------------------------------------
//
// Cost computation of the "very" short rotations (< CD_min)
//
//----------------------------------------------------------------

// Brutally try all possible short rotations from the very first day
bool SubProblem::priceVeryShortRotationsFirstDay(){
	int nFound = 0;
	if(startingDayStatus_[0]){
		for(int c=1; c<CDMin_; c++){
			vector2D succs = allowedShortSuccBySize_[c];
			for(unsigned int i=0; i<succs.size(); i++){
				vector<int> succ = allowedShortSuccBySize_[c][i];
				double redCost = costOfVeryShortRotation(0,succ);
				int rotationLength = succ.size() + (pLiveNurse_->pStateIni_->shift_ > 0 ? pLiveNurse_->pStateIni_->consDaysWorked_ : 0);
				if(redCost < maxReducedCostBound_ and rotationLength <= maxRotationLength_){
					Rotation rot (0, succ, pLiveNurse_->id_, MAX_COST, redCost);
					theRotations_.push_back(rot);
					nPaths_ ++;
					nVeryShortFound_++;
					nFound ++;
					bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
				}
			}
		}
	}
	return nFound > 0;
}

// Brutally try all possible short rotations that end on the last day
bool SubProblem::priceVeryShortRotationsLastDay(){
	int nFound = 0;
	for(int c=1; c<CDMin_; c++){
		if(startingDayStatus_[nDays_-c]){
			vector2D succs = allowedShortSuccBySize_[c];
			for(unsigned int i=0; i<succs.size(); i++){
				vector<int> succ = allowedShortSuccBySize_[c][i];
				double redCost = costOfVeryShortRotation( nDays_-c ,succ);
				if(redCost < maxReducedCostBound_){
					Rotation rot (nDays_-c, succ, pLiveNurse_->id_, MAX_COST, redCost);
					theRotations_.push_back(rot);
					nPaths_ ++;
					nVeryShortFound_++;
					nFound ++;
					bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
				}
			}
		}
	}
	return nFound > 0;
}

// Brutally try all possible short rotations from every first day
bool SubProblem::priceVeryShortRotations(){
	int nFound = 0;
	for(int c=1; c<CDMin_; c++){
		vector2D succs = allowedShortSuccBySize_[c];
		for(unsigned int i = 0; i<succs.size(); i++){
			vector<int> succ = allowedShortSuccBySize_[c][i];
			for(int k=0; k <= nDays_ - c; k++){
				if(startingDayStatus_[k]){
					double redCost = costOfVeryShortRotation(k,succ);
					if(redCost < maxReducedCostBound_){
						Rotation rot (k, succ, pLiveNurse_->id_, MAX_COST, redCost);
						theRotations_.push_back(rot);
						nPaths_ ++;
						nVeryShortFound_++;
						nFound ++;
						bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
					}
				}
			}
		}
	}
	return nFound > 0;
}

// Compute the cost of a single short rotation
double SubProblem::costOfVeryShortRotation(int startDate, vector<int> succ){

	int endDate = startDate + succ.size() - 1;

	// check if the first day is forbidden
	if(!startingDayStatus_[startDate])
		return MAX_COST;

	//check if any shift is forbidden
	for(int k=startDate; k<=endDate; k++)
		if(isDayShiftForbidden(k, succ[k-startDate]))
			return MAX_COST;

	// Regular costs
	double consDaysRegCost=0, consShiftsRegCost=0, completeWeekendRegCost=0, preferencesRegCost=0, shortRestRegCost=0;
	// Reduced costs
	double dayShiftsRedCost=0, startRedCost=0, endRedCost=0, weekendRedCost=0;

	// Initialize values
	int shift=0, consShifts=0, consDays=succ.size();

	// A. SPECIAL CASE OF THE FIRST DAY
	//
	if(startDate==0){
		// The nurse was working
		if(pLiveNurse_->pStateIni_->shift_ > 0){

			if(pScenario_->isForbiddenSuccessorShift_Shift( succ[0], pLiveNurse_->pStateIni_->shift_)){
				return MAX_COST;
			}

			// Change initial values
			shift = pLiveNurse_->pStateIni_->shift_;
			consShifts = pLiveNurse_->pStateIni_->consShifts_;
			consDays += pLiveNurse_->pStateIni_->consDaysWorked_;
			// If worked too much, subtract the already counted surplus
			consDaysRegCost -= max(0, pLiveNurse_->pStateIni_->consDaysWorked_ - pContract_->maxConsDaysWork_) * WEIGHT_CONS_DAYS_WORK;
			consShiftsRegCost -= max(0, pLiveNurse_->pStateIni_->consShifts_ - pScenario_->maxConsShiftsOfTypeOf(shift)) * WEIGHT_CONS_SHIFTS;
		}
		// The nurse was resting
		else {
			// Cost of a too short rest
			shortRestRegCost += max(0, pContract_->minConsDaysOff_ - pLiveNurse_->pStateIni_->consDaysOff_) * WEIGHT_CONS_DAYS_OFF;
		}
	}

	// B. REGULAR COST: CONSECUTIVE NUMBER OF DAYS (only if it does not end on last day)
	//
	if((int) (startDate+succ.size()) < nDays_ and consDays){
		consDaysRegCost += consDaysCost(consDays);
	}

	// C. REGULAR COST: CONSECUTIVE SHIFTS
	//
	for(int k=startDate; k<=endDate; k++){
		int newShift = succ[k-startDate];
		if(newShift == shift){
			consShifts ++;
		} else {
			consShiftsRegCost += consShiftCost(shift,consShifts);
			consShifts = 1;
			shift = newShift;
		}
		if(k==endDate and (k<nDays_-1 or consShifts > pScenario_->maxConsShiftsOfTypeOf(shift))) consShiftsRegCost += consShiftCost(shift, consShifts);
	}

	// D. REGULAR COST: COMPLETE WEEKENDS
	//
	completeWeekendRegCost = startWeekendCosts_[startDate] + endWeekendCosts_[endDate];

	// E. REGULAR COST: PREFERENCES
	//
	for(int k=startDate; k<=endDate; k++) preferencesRegCost += preferencesCosts_[k][ succ[k-startDate] ];



	// F. REDUCED COST: WEEKENDS
	//
	weekendRedCost -= Tools::containsWeekend(startDate, endDate) * pCosts_->workedWeekendCost();

	// F. REDUCED COST: FIRST DAY (BACK TO WORK)
	//
	startRedCost -= pCosts_->startWorkCost(startDate);

	// G. REDUCED COST: LAST DAY (BACK TO WORK)
	//
	endRedCost -= pCosts_->endWorkCost(endDate);

	// H. REDUCED COST: EACH DAY/SHIFT REDUCED COST
	//
	for(int k=startDate; k<=endDate; k++) dayShiftsRedCost -= pCosts_->dayShiftWorkCost( k, succ[k-startDate] - 1 );


	// I. RETURN THE TOTAL COST
	//
	double regCost = consDaysRegCost + consShiftsRegCost + completeWeekendRegCost + preferencesRegCost + shortRestRegCost;
	double redCost = dayShiftsRedCost + startRedCost + endRedCost + weekendRedCost;
	double ANS = regCost + redCost;

	if(false){
		cout << "# " << endl;
		cout << "#+---------------------------------------------+ (" << pLiveNurse_->name_ << ")" << endl;
		cout << "# " << startDate << "-";
		for(unsigned int i=0; i<succ.size(); i++) cout << pScenario_->intToShift_[succ[i]].at(0);
		cout << endl;
		cout << "# length " << consDays;
		if(startDate == 0 and pLiveNurse_->pStateIni_->shift_ > 0){
			cout << " (incl. " << pLiveNurse_->pStateIni_->consDaysWorked_ << " before planing horizon)";
		}
		cout << endl;

		cout << "# REG- Consecutive days cost   : " << consDaysRegCost << endl;
		cout << "# REG- Consecutive shifts cost : " << consShiftsRegCost << endl;
		cout << "# REG- Complete weekends cost  : " << completeWeekendRegCost << endl;
		cout << "# REG- Preferences cost        : " << preferencesRegCost << endl;
		cout << "# REG- Short rest before cost  : " << shortRestRegCost << endl;
		cout << "# REG-                 ~TOTAL~ : " << regCost << endl;
		cout << "# " << endl;
		cout << "# RED- Day-shifts cost         : " << dayShiftsRedCost << endl;
		cout << "# RED- Start work cost         : " << startRedCost << endl;
		cout << "# RED- End work cost           : " << endRedCost << endl;
		cout << "# RED- Weekend dual cost       : " << weekendRedCost << endl;
		cout << "# RED-                 ~TOTAL~ : " << redCost << endl;
		cout << "#+---------------------------------------------+" << endl;
		cout << "#                      ~TOTAL~ : " << ANS << endl;
		cout << "#+---------------------------------------------+" << endl;
		cout << "# " << endl;
		//getchar();
	}

	return ANS;
}





//----------------------------------------------------------------
//
// Shortest path function with several sinks
// (modified from boost so that we can give several sink nodes)
//
//----------------------------------------------------------------

// r_c_shortest_paths_several_sinks function -> calls r_c_shortest_paths_dispatch
template<class Graph,
class VertexIndexMap,
class EdgeIndexMap,
class Resource_Container,
class Resource_Extension_Function,
class Dominance_Function,
class Label_Allocator,
class Visitor>
void SubProblem::r_c_shortest_paths_several_sinks( const Graph& g,
		const VertexIndexMap& vertex_index_map,
		const EdgeIndexMap& edge_index_map,
		typename boost::graph_traits<Graph>::vertex_descriptor s,
		std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> t,
		// each inner vector corresponds to a pareto-optimal path
		std::vector<std::vector<typename boost::graph_traits<Graph>::edge_descriptor> >&
		pareto_optimal_solutions,
		std::vector<Resource_Container>& pareto_optimal_resource_containers,
		// to initialize the first label/resource container
		// and to carry the type information
		const Resource_Container& rc,
		const Resource_Extension_Function& ref,
		const Dominance_Function& dominance,
		// to specify the memory management strategy for the labels
		Label_Allocator la,
		Visitor vis ){
	r_c_shortest_paths_dispatch_several_sinks( g,
			vertex_index_map,
			edge_index_map,
			s,
			t,
			pareto_optimal_solutions,
			pareto_optimal_resource_containers,
			true,
			rc,
			ref,
			dominance,
			la,
			vis );
}

// r_c_shortest_paths_dispatch function (body/implementation)
template<class Graph,
class VertexIndexMap,
class EdgeIndexMap,
class Resource_Container,
class Resource_Extension_Function,
class Dominance_Function,
class Label_Allocator,
class Visitor>
void SubProblem::r_c_shortest_paths_dispatch_several_sinks( const Graph& g,
		const VertexIndexMap& vertex_index_map,
		const EdgeIndexMap& /*edge_index_map*/,
		typename boost::graph_traits<Graph>::vertex_descriptor s,
		std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> t,
		// each inner vector corresponds to a pareto-optimal path
		std::vector
		<std::vector
		<typename boost::graph_traits
		<Graph>::edge_descriptor> >& pareto_optimal_solutions,
		std::vector
		<Resource_Container>& pareto_optimal_resource_containers,
		bool b_all_pareto_optimal_solutions,
		// to initialize the first label/resource container
		// and to carry the type information
		const Resource_Container& rc,
		Resource_Extension_Function& ref,
		Dominance_Function& dominance,
		// to specify the memory management strategy for the labels
		Label_Allocator /*la*/,
		Visitor vis ){
	pareto_optimal_resource_containers.clear();
	pareto_optimal_solutions.clear();

	size_t i_label_num = 0;
	typedef
			typename
			Label_Allocator::template rebind
			<boost::r_c_shortest_paths_label
			<Graph, Resource_Container> >::other LAlloc;
	LAlloc l_alloc;
	typedef
			ks_smart_pointer
			<boost::r_c_shortest_paths_label<Graph, Resource_Container> > Splabel;
	std::priority_queue<Splabel, std::vector<Splabel>, std::greater<Splabel> >
	unprocessed_labels;

	int nbCreatedLabels = 0;
	int nbDeletedLabels = 0;
	bool b_feasible = true;
	boost::r_c_shortest_paths_label<Graph, Resource_Container>* first_label =
			l_alloc.allocate( 1 );
	l_alloc.construct
	( first_label,
			boost::r_c_shortest_paths_label
			<Graph, Resource_Container>( i_label_num++,
					rc,
					0,
					typename boost::graph_traits<Graph>::
					edge_descriptor(),
					s ) );
	nbCreatedLabels++;

	Splabel splabel_first_label = Splabel( first_label );
	unprocessed_labels.push( splabel_first_label );
	std::vector<std::list<Splabel> > vec_vertex_labels_data( num_vertices( g ) );
	boost::iterator_property_map<typename std::vector<std::list<Splabel> >::iterator,
	VertexIndexMap>
	vec_vertex_labels(vec_vertex_labels_data.begin(), vertex_index_map);
	vec_vertex_labels[s].push_back( splabel_first_label );
	typedef
			std::vector<typename std::list<Splabel>::iterator>
	vec_last_valid_positions_for_dominance_data_type;
	vec_last_valid_positions_for_dominance_data_type
	vec_last_valid_positions_for_dominance_data( num_vertices( g ) );
	boost::iterator_property_map<
	typename vec_last_valid_positions_for_dominance_data_type::iterator,
	VertexIndexMap>
	vec_last_valid_positions_for_dominance
	(vec_last_valid_positions_for_dominance_data.begin(),
			vertex_index_map);
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		put(vec_last_valid_positions_for_dominance, v, vec_vertex_labels[v].begin());
	}
	std::vector<size_t> vec_last_valid_index_for_dominance_data( num_vertices( g ), 0 );
	boost::iterator_property_map<std::vector<size_t>::iterator, VertexIndexMap>
	vec_last_valid_index_for_dominance
	(vec_last_valid_index_for_dominance_data.begin(), vertex_index_map);
	std::vector<bool>
	b_vec_vertex_already_checked_for_dominance_data( num_vertices( g ), false );
	boost::iterator_property_map<std::vector<bool>::iterator, VertexIndexMap>
	b_vec_vertex_already_checked_for_dominance
	(b_vec_vertex_already_checked_for_dominance_data.begin(),
			vertex_index_map);

	std::vector<Splabel> label_trash;

	while( !unprocessed_labels.empty()  && vis.on_enter_loop(unprocessed_labels, g) )
	{
		Splabel cur_label = unprocessed_labels.top();
		assert (cur_label->b_is_valid);
		unprocessed_labels.pop();
		vis.on_label_popped( *cur_label, g );
		// an Splabel object in unprocessed_labels and the respective Splabel
		// object in the respective list<Splabel> of vec_vertex_labels share their
		// embedded r_c_shortest_paths_label object
		// to avoid memory leaks, dominated
		// r_c_shortest_paths_label objects are marked and deleted when popped
		// from unprocessed_labels, as they can no longer be deleted at the end of
		// the function; only the Splabel object in unprocessed_labels still
		// references the r_c_shortest_paths_label object
		// this is also for efficiency, because the else branch is executed only
		// if there is a chance that extending the
		// label leads to new undominated labels, which in turn is possible only
		// if the label to be extended is undominated
		assert (cur_label->b_is_valid);
		if( !cur_label->b_is_dominated )
		{
			typename boost::graph_traits<Graph>::vertex_descriptor
			i_cur_resident_vertex = cur_label->resident_vertex;
			std::list<Splabel>& list_labels_cur_vertex =
					get(vec_vertex_labels, i_cur_resident_vertex);
			if( list_labels_cur_vertex.size() >= 2
					&& vec_last_valid_index_for_dominance[i_cur_resident_vertex]
														  < list_labels_cur_vertex.size() )
			{
				typename std::list<Splabel>::iterator outer_iter =
						list_labels_cur_vertex.begin();
				bool b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = false;
				while( outer_iter != list_labels_cur_vertex.end() )
				{
					Splabel cur_outer_splabel = *outer_iter;
					assert (cur_outer_splabel->b_is_valid);
					typename std::list<Splabel>::iterator inner_iter = outer_iter;
					if( !b_outer_iter_at_or_beyond_last_valid_pos_for_dominance
							&& outer_iter ==
									get(vec_last_valid_positions_for_dominance,
											i_cur_resident_vertex) )
						b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = true;
					if( !get(b_vec_vertex_already_checked_for_dominance, i_cur_resident_vertex)
							|| b_outer_iter_at_or_beyond_last_valid_pos_for_dominance )
					{
						++inner_iter;
					}
					else
					{
						inner_iter =
								get(vec_last_valid_positions_for_dominance,
										i_cur_resident_vertex);
						++inner_iter;
					}
					bool b_outer_iter_erased = false;
					while( inner_iter != list_labels_cur_vertex.end() )
					{
						Splabel cur_inner_splabel = *inner_iter;
						assert (cur_inner_splabel->b_is_valid);
						if( dominance( cur_outer_splabel->
								cumulated_resource_consumption,
								cur_inner_splabel->
								cumulated_resource_consumption ) )
						{
							typename std::list<Splabel>::iterator buf = inner_iter;
							++inner_iter;
							list_labels_cur_vertex.erase( buf );
							if( cur_inner_splabel->b_is_processed )
							{
								cur_inner_splabel->b_is_valid = false;
								l_alloc.destroy( cur_inner_splabel.get() );
								l_alloc.deallocate( cur_inner_splabel.get(), 1 );
								nbDeletedLabels++;
							}
							else
								cur_inner_splabel->b_is_dominated = true;
							continue;
						}
						else
							++inner_iter;
						if( dominance( cur_inner_splabel->
								cumulated_resource_consumption,
								cur_outer_splabel->
								cumulated_resource_consumption ) )
						{
							typename std::list<Splabel>::iterator buf = outer_iter;
							++outer_iter;
							list_labels_cur_vertex.erase( buf );
							b_outer_iter_erased = true;
							assert (cur_outer_splabel->b_is_valid);
							if( cur_outer_splabel->b_is_processed )
							{
								label_trash.push_back(cur_outer_splabel);
								// DBG : Attention !!!! j'ai du decommenter ces lignes pour
								// eviter d'enormes fuites memoire
								// cur_outer_splabel->b_is_valid = false;
								// l_alloc.destroy( cur_outer_splabel.get() );
								// l_alloc.deallocate( cur_outer_splabel.get(), 1 );
								nbDeletedLabels++;
							}
							else
								cur_outer_splabel->b_is_dominated = true;
							break;
						}
					}
					if( !b_outer_iter_erased )
						++outer_iter;
				}
				if( list_labels_cur_vertex.size() > 1 )
					put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex,
							(--(list_labels_cur_vertex.end())));
				else
					put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex,
							list_labels_cur_vertex.begin());
				put(b_vec_vertex_already_checked_for_dominance,
						i_cur_resident_vertex, true);
				put(vec_last_valid_index_for_dominance, i_cur_resident_vertex,
						list_labels_cur_vertex.size() - 1);
			}
		}
		assert (b_all_pareto_optimal_solutions || cur_label->b_is_valid);

		// ------------------------------------------------------------------------- START SAMUEL
		//if( !b_all_pareto_optimal_solutions && cur_label->resident_vertex == t )
		if( !b_all_pareto_optimal_solutions && std::find(t.begin(), t.end(), cur_label->resident_vertex) !=t.end())
		{
			// ------------------------------------------------------------------------- END SAMUEL

			// the devil don't sleep
			if( cur_label->b_is_dominated )
			{
				cur_label->b_is_valid = false;
				l_alloc.destroy( cur_label.get() );
				l_alloc.deallocate( cur_label.get(), 1 );
				nbDeletedLabels++;
			}
			while( unprocessed_labels.size() )
			{
				Splabel l = unprocessed_labels.top();
				assert (l->b_is_valid);
				unprocessed_labels.pop();
				// delete only dominated labels, because nondominated labels are
				// deleted at the end of the function
				if( l->b_is_dominated )
				{
					l->b_is_valid = false;
					l_alloc.destroy( l.get() );
					l_alloc.deallocate( l.get(), 1 );
					nbDeletedLabels++;
				}
			}
			break;
		}
		if( !cur_label->b_is_dominated )
		{
			cur_label->b_is_processed = true;
			vis.on_label_not_dominated( *cur_label, g );
			typename boost::graph_traits<Graph>::vertex_descriptor cur_vertex =
					cur_label->resident_vertex;
			typename boost::graph_traits<Graph>::out_edge_iterator oei, oei_end;
			for( boost::tie( oei, oei_end ) = out_edges( cur_vertex, g );
					oei != oei_end;
					++oei )
			{
				b_feasible = true;
				boost::r_c_shortest_paths_label<Graph, Resource_Container>* new_label =
						l_alloc.allocate( 1 );
				l_alloc.construct( new_label,
						boost::r_c_shortest_paths_label
						<Graph, Resource_Container>
				( i_label_num++,
						cur_label->cumulated_resource_consumption,
						cur_label.get(),
						*oei,
						target( *oei, g ) ) );
				nbCreatedLabels++;

				b_feasible =
						ref( g,
								new_label->cumulated_resource_consumption,
								new_label->p_pred_label->cumulated_resource_consumption,
								new_label->pred_edge );

				if( !b_feasible )
				{
					vis.on_label_not_feasible( *new_label, g );
					new_label->b_is_valid = false;
					l_alloc.destroy( new_label );
					l_alloc.deallocate( new_label, 1 );
					nbDeletedLabels++;
				}
				else
				{
					const boost::r_c_shortest_paths_label<Graph, Resource_Container>&
					ref_new_label = *new_label;
					vis.on_label_feasible( ref_new_label, g );
					Splabel new_sp_label( new_label );
					vec_vertex_labels[new_sp_label->resident_vertex].
					push_back( new_sp_label );
					unprocessed_labels.push( new_sp_label );
				}
			}
		}
		else
		{
			assert (cur_label->b_is_valid);
			vis.on_label_dominated( *cur_label, g );
			cur_label->b_is_valid = false;
			l_alloc.destroy( cur_label.get() );
			l_alloc.deallocate( cur_label.get(), 1 );
			nbDeletedLabels++;
		}
	}

	// ------------------------------------------------------------------------- START SAMUEL
	typename std::list<Splabel>::const_iterator csi;
	typename std::list<Splabel>::const_iterator csi_end;
	for(int sink=0; sink<t.size(); sink++){
		std::list<Splabel> dsplabels = get(vec_vertex_labels, t[sink]);
		csi = dsplabels.begin();
		csi_end = dsplabels.end();
		// if d could be reached from o
		if( !dsplabels.empty() )
		{
			for( ; csi != csi_end; ++csi )
			{
				std::vector<typename boost::graph_traits<Graph>::edge_descriptor>
				cur_pareto_optimal_path;
				const boost::r_c_shortest_paths_label<Graph, Resource_Container>* p_cur_label =
						(*csi).get();
				assert (p_cur_label->b_is_valid);
				pareto_optimal_resource_containers.
				push_back( p_cur_label->cumulated_resource_consumption );
				while( p_cur_label->num != 0 )
				{
					cur_pareto_optimal_path.push_back( p_cur_label->pred_edge );
					p_cur_label = p_cur_label->p_pred_label;
					assert (p_cur_label->b_is_valid);
				}
				pareto_optimal_solutions.push_back( cur_pareto_optimal_path );
				if( !b_all_pareto_optimal_solutions )
					break;
			}
		}
	}
	/*
    std::list<Splabel> dsplabels = get(vec_vertex_labels, t);
      typename std::list<Splabel>::const_iterator csi = dsplabels.begin();
      typename std::list<Splabel>::const_iterator csi_end = dsplabels.end();
      // if d could be reached from o
      if( !dsplabels.empty() )
      {
        for( ; csi != csi_end; ++csi )
        {
          std::vector<typename graph_traits<Graph>::edge_descriptor>
            cur_pareto_optimal_path;
          const r_c_shortest_paths_label<Graph, Resource_Container>* p_cur_label =
            (*csi).get();
          assert (p_cur_label->b_is_valid);
          pareto_optimal_resource_containers.
            push_back( p_cur_label->cumulated_resource_consumption );
          while( p_cur_label->num != 0 )
          {
            cur_pareto_optimal_path.push_back( p_cur_label->pred_edge );
            p_cur_label = p_cur_label->p_pred_label;
            assert (p_cur_label->b_is_valid);
          }
          pareto_optimal_solutions.push_back( cur_pareto_optimal_path );
          if( !b_all_pareto_optimal_solutions )
            break;
        }
      }
	 */
	// ------------------------------------------------------------------------- END SAMUEL 

	BGL_FORALL_VERTICES_T(i, g, Graph) {
		const std::list<Splabel>& list_labels_cur_vertex = vec_vertex_labels[i];
		csi_end = list_labels_cur_vertex.end();
		for( csi = list_labels_cur_vertex.begin(); csi != csi_end; ++csi )
		{
			assert ((*csi)->b_is_valid);
			(*csi)->b_is_valid = false;
			l_alloc.destroy( (*csi).get() );
			l_alloc.deallocate( (*csi).get(), 1 );
			nbDeletedLabels++;
		}
	}
	for (Splabel label: label_trash) {
		assert(label->b_is_valid );
		label->b_is_valid = false;
		l_alloc.destroy( label.get() );
		l_alloc.deallocate( label.get(), 1 );
	}
	// std::cout << "Nb created labels " << nbCreatedLabels << std::endl;
	// std::cout << "Nb deleted labels " << nbDeletedLabels << std::endl;
} // r_c_shortest_paths_dispatch






//----------------------------------------------------------------
//
// Greedy heuristic for the shortest path problem with resource
// constraints.
//
//----------------------------------------------------------------
bool SubProblem::solveLongRotationsHeuristic(){

	// ALL OTHER DAYS
	//
	int nFound = 0;
	for(int startDate=1; startDate<nDays_-CDMin_; startDate ++){

		vector<int> bestSucc;
		double bestCost = 0;
		int currentDate = startDate -1;
		int length = 0;
		int lastSh = 0;
		int nConsShift = 0;
		bool hasFoundBetter = true;

		// Try to extend the rotation as long as it is of negative reduced cost AND does not exceed the maximum length allowed
		//
		while(hasFoundBetter and length < maxRotationLength_-1 and currentDate < nDays_-1){
			hasFoundBetter = false;

			// Initialization, rotations of length CDMin_ (so that all rotations generated are long)
			//
			if(length == 0){
				double bestCostCDMinDays = 0;
				vector<int> bestSuccCDMinDays;
				for(int newSh=1; newSh<pScenario_->nbShiftsType_; newSh++){
					for(unsigned int i=0; i<allowedShortSuccBySize_[CDMin_].size(); i++){
						vector<int> succ = allowedShortSuccBySize_[CDMin_][i];
						double potentialCost = costOfVeryShortRotation(startDate,succ);
						// Store the rotation if it is good
						if(potentialCost < maxReducedCostBound_){
							vector<int> succToStore = succ;
							Rotation rot (startDate, succToStore, pLiveNurse_->id_, MAX_COST, potentialCost);
							nPaths_ ++;
							theRotations_.push_back(rot);
							nLongFound_ ++;
							nFound ++;
							bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
						}
						// Chose to extend that one if it is the best
						if(potentialCost < bestCostCDMinDays){
							bestSuccCDMinDays = succ;
							bestCostCDMinDays = potentialCost;
							hasFoundBetter = true;
						}

					}
				}

				// If a good one has been found -> make it the basis for the next extension
				if(hasFoundBetter){
					bestSucc = bestSuccCDMinDays;
					bestCost = bestCostCDMinDays;
					lastSh = bestSucc[CDMin_-1];
					length = CDMin_;
					currentDate = currentDate + CDMin_;
					nConsShift = 0;
					while(bestSucc[bestSucc.size()-nConsShift-1] == lastSh) nConsShift ++;
				}
				continue;
			}

			// OTHER CASES: EXTEND THE ROTATION
			//
			else{

				int bestNewShift;
				double bestNewCost = bestCost;
				hasFoundBetter = false;

				for(int newSh=1; newSh<pScenario_->nbShiftsType_; newSh++){
					double potentialCost = bestCost;

					// If the succession is allowed -> try to extend the rotation
					if(!pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh) and !isDayShiftForbidden(currentDate+1,newSh)){

						// REGULAR COSTS
						//

						// REG COST: CONSECUTIVE SHIFTS -> if last day + too short, totally cancel that cost !!
						if(currentDate < nDays_-2){
							if(newSh==lastSh){
								potentialCost -= consShiftCost(lastSh,nConsShift);
								potentialCost += consShiftCost(lastSh,nConsShift+1);
							} else{
								potentialCost += consShiftCost(newSh,1);
							}
						}
						// Last day:
						else {
							if(newSh == lastSh){
								potentialCost -= consShiftCost(lastSh,nConsShift);
								if(newSh == lastSh and nConsShift+1 > pScenario_->maxConsShiftsOfTypeOf(newSh))
									potentialCost += consShiftCost(lastSh,nConsShift+1);
							}
						}

						// if last day + too short -> cancel the previous cost:


						// REG COST: CONSECUTIVE DAYS
						potentialCost -= consDaysCost(length);
						if(currentDate < nDays_-2 or length + 1 >= pContract_->minConsDaysWork_){
							potentialCost += consDaysCost(length+1);
						}

						// REG COST: COMPLETE WEEKENDS ?
						if(pContract_->needCompleteWeekends_){
							if(Tools::isSaturday(currentDate+1)) potentialCost += WEIGHT_COMPLETE_WEEKEND;
							if(Tools::isSunday(currentDate+1)) potentialCost -= WEIGHT_COMPLETE_WEEKEND;
						}

						// REG COST: PREFERENCES
						if(pLiveNurse_->wishesOff(currentDate+1,newSh)) potentialCost += WEIGHT_PREFERENCES;


						// DUAL COSTS
						//

						// RED COST: WEEKEND
						if(Tools::isSaturday(currentDate+1)) potentialCost -= pCosts_->workedWeekendCost();

						// RED COST: FIRST DAY -> no need here because already in the beginning

						// RED COST: LAST DAY
						potentialCost += pCosts_->endWorkCost(currentDate);
						potentialCost -= pCosts_->endWorkCost(currentDate+1);

						// RED COST: DAY-SHIFT
						potentialCost -= pCosts_->dayShiftWorkCost(currentDate+1,newSh-1);


						if(false){
							cout << "# COST COMPUTATION for extension to " << startDate << "-";
							for(unsigned int i=0; i<bestSucc.size(); i++) cout << (pScenario_->intToShift_[bestSucc[i]])[0];
							cout << (pScenario_->intToShift_[newSh])[0] << endl;
							cout << "#    Previous: " << bestCost << endl;
							cout << "#          -=: " << consShiftCost(lastSh,nConsShift) << " + " << consShiftCost(lastSh,nConsShift+1) << endl;
							cout << "#          +=: " << consShiftCost(newSh,1) << endl;
							cout << "#          -=: " << consDaysCost(length) << " + " << consDaysCost(length+1) << endl;
							cout << "#          +=: ";
							if(pContract_->needCompleteWeekends_){
								if(Tools::isSaturday(currentDate+1)) cout << WEIGHT_COMPLETE_WEEKEND;
								if(Tools::isSunday(currentDate+1)) cout << " - " << WEIGHT_COMPLETE_WEEKEND;
							}
							cout << endl;
							cout << "#          +=: ";
							if(pLiveNurse_->wishesOff(currentDate+1,newSh)) cout << WEIGHT_PREFERENCES;
							cout << endl;
							cout << "#          -=: ";
							if(Tools::isSaturday(currentDate+1)) cout << pCosts_->workedWeekendCost();
							cout << endl;
							cout << "#          +=: " << pCosts_->endWorkCost(currentDate) << " - " << pCosts_->endWorkCost(currentDate+1) << endl;
							cout << "#          -=: " << pCosts_->dayShiftWorkCost(currentDate+1,newSh-1) << endl;
						}

						// CHECK IF CAN BE ADDED
						//
						if(potentialCost < maxReducedCostBound_){
							vector<int> succToStore = bestSucc;
							succToStore.push_back(newSh);
							Rotation rot (startDate, succToStore, pLiveNurse_->id_, MAX_COST, potentialCost);
							nPaths_ ++;
							theRotations_.push_back(rot);
							nLongFound_ ++;
							nFound ++;
							bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
						}

						// CHECK IF BETTER THAN PREVIOUS ONE
						//
						if(potentialCost < bestNewCost){
							bestNewShift = newSh;
							bestNewCost = potentialCost;
							hasFoundBetter = true;
						}
					}
				}

				// Update the rotation for the next step
				//
				if(hasFoundBetter){
					bestSucc.push_back(bestNewShift);
					bestCost = bestNewCost;
					currentDate ++;
					length ++;
					nConsShift = (lastSh == bestNewShift) ? (nConsShift+1) : 1;
					lastSh = bestNewShift;
				}
			}
		}
	}
	return nFound > 0;
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
	printAllNodes();

	// THE ARCS
	//
	printAllArcs();

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

// Prints all nodes
void SubProblem::printAllNodes(){
	std::cout << "#   NODES (" << nNodes_ << ")" << std::endl;
	for(int v=0; v<nNodes_; v++) std::cout << printNode(v) << std::endl;
	std::cout << "# " << std::endl;

}
// Prints the line of an arc

string SubProblem::printArc(int a){
	stringstream rep;
	rep << "# ARC   " << a << " \t" << arcTypeName[arcType(a)] << " \t";
	rep << "(" << arcOrigin(a) << "," << arcDestination(a) << ") \t" << "c= " << arcCost(a) ;
	arcCost(a) < 10000 ? rep << "     " : rep << "";
	rep << "\tt=" << arcLength(a);
	rep << " \t[" << shortNameNode(arcOrigin(a)) << "] -> [" << shortNameNode(arcDestination(a)) << "]";
	return rep.str();

}

// Prints all arcs
void SubProblem::printAllArcs(){
	std::cout << "#   ARCS (" << nArcs_ << "]" << std::endl;
	for(int a=0; a<nArcs_; a++) std::cout << printArc(a) << std::endl;
	std::cout << "# " << std::endl;
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

	else if (type_v == SINK_DAY){
		rep << "SINK DAY";
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
	rep << "#     [ " << nDays_ << " days, " << (pScenario_->nbShiftsType_-1) << " shifts ]" << std::endl;
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
	for(unsigned int i=0; i<allowedShortSuccBySize_.size(); i++){
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
void SubProblem::printPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource){

	// The successive nodes, and corresponding arc costs / time
	//
	std::cout << "# " << std::endl;
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		std::cout << "# \t| [ " << shortNameNode(source( path[j], g_ )) << " ]";
		std::cout << "\t\tCost:  " << arcCost(a) << "\t\tTime:" << arcLength(a);
		std::cout << "\t\t[" << (arcStatus_[a] ? " allowed " : "forbidden") << "]" << std::endl;
	}
	std::cout << "# " << std::endl;
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		std::cout << printArc(a) << endl;
	}



	// Last node and total
	//
	std::cout << "# \t| [" << shortNameNode(sinkNode_) << "]" << std::endl;
	std::cout << "# \t| ~TOTAL~   \t\tCost:   " << resource.cost << "\t\tTime: " << resource.time << std::endl;
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
			for(int s: allowedShortSuccBySize_[CDMin_][ shortSuccCDMinIdFromArc_.at(a) ]){
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
	for(unsigned int i=0; i<allTasks.size(); i++){
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
		for(int s=1; s<pScenario_->nbShiftsType_; s++){
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
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				if(!isArcForbidden(arcsFromSource_[s][k][n])){
					int v = principalNetworkNodes_[s][k][n];
					int succId = idBestShortSuccCDMin_[s][k][n];
					vector<int> succ = allowedShortSuccBySize_[CDMin_][succId];
					std::cout << "# " << shortNameNode(v) << " <- (id=" << succId << ") ";
					for(unsigned int i=0; i<succ.size(); i++){
						std::cout << pScenario_->intToShift_[succ[i]].at(0);
					}
					std::cout << std::endl;
				}
			}
		}
	}
}

// Print the contract type + preferences
void SubProblem::printContractAndPrefenrences(){
	std::cout << "# Preferences:" << endl;
	for(map<int,set<int> >::iterator it = pLiveNurse_->pWishesOff_->begin(); it != pLiveNurse_->pWishesOff_->end(); ++it){
		cout <<  "      | " << it->first << "  ->  ";
		for(int s : it->second) cout << pScenario_->intToShift_[s];
		cout << endl;
	}
	std::cout << "# Contract :   ";
	std::cout << "Days [" << pContract_->minConsDaysOff_ << "<" << pContract_->maxConsDaysOff_ << "]   ";
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
	  std::cout << pScenario_->intToShift_[s] << " [" << pScenario_->minConsShiftsOfTypeOf(s) << "<" << pScenario_->maxConsShiftsOfTypeOf(s) << "]   ";
	}
	std::cout << std::endl;
	std::cout << "# " << std::endl;
	std::cout << "# " << std::endl;
}









/*************************************************************************
 * A GARDER AU CAS OU COMME EXEMPLE POUR CERTAINES FONCTIONS / SYNTAXES. *
 *************************************************************************/
enum nodes { A, B, C, D, E, F, G };
char name[] = "ABCDEFG";

void SubProblem::testGraph_spprc(){


	std::cout << "# " << std::endl;
	std::cout << "# ====================================================================================" << std::endl;
	std::cout << "# = Fonction de test du probleme de plus court chemin avec contraintes de ressources =" << std::endl;
	std::cout << "# ====================================================================================" << std::endl;
	std::cout << "# " << std::endl;

	Graph g;

	add_vertex( Vertex_Properties( A, NONE_NODE, 0, 0 ), g );
	add_vertex( Vertex_Properties( B, NONE_NODE, 3, 20 ), g );
	add_vertex( Vertex_Properties( C, NONE_NODE, 5, 20 ), g );
	add_vertex( Vertex_Properties( D, NONE_NODE, 0, 20 ), g );
	add_vertex( Vertex_Properties( E, NONE_NODE, 0, 20 ), g );
	add_vertex( Vertex_Properties( F, NONE_NODE, 1, 2 ), g );
	add_vertex( Vertex_Properties( G, NONE_NODE, 10, 20 ), g );

	add_edge( A, B, Arc_Properties( 0, NONE_ARC, 1, 1 ), g );
	add_edge( A, D, Arc_Properties( 1, NONE_ARC, 1, 1 ), g );
	add_edge( B, C, Arc_Properties( 2, NONE_ARC, 1, 2 ), g );
	add_edge( B, G, Arc_Properties( 3, NONE_ARC, 1, 12 ), g );
	add_edge( D, C, Arc_Properties( 4, NONE_ARC, 1, 30 ), g );
	add_edge( D, G, Arc_Properties( 5, NONE_ARC, 5, 2 ), g );
	add_edge( D, E, Arc_Properties( 6, NONE_ARC, 2, 5 ), g );
	add_edge( D, F, Arc_Properties( 7, NONE_ARC, 5, 2 ), g );
	add_edge( E, G, Arc_Properties( 8, NONE_ARC, 1, 4 ), g );
	add_edge( E, F, Arc_Properties( 9, NONE_ARC, 0, 2 ), g );


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

	boost::graph_traits<Graph>::vertex_descriptor s = A;
	boost::graph_traits<Graph>::vertex_descriptor t = G;

	std::vector<boost::graph_traits<Graph>::vertex_descriptor> tVec;
	tVec.push_back(C);
	tVec.push_back(G);
	tVec.push_back(F);

	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw_single;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw_single;

	r_c_shortest_paths(
			g,
			get( &Vertex_Properties::num, g ),
			get( &Arc_Properties::num, g ),
			s,
			t,
			opt_solutions_spptw_single,
			pareto_opt_rcs_spptw_single,
			spp_spptw_res_cont( 0, 0 ),
			ref_spptw(),
			dominance_spptw(),
			std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
			boost::default_r_c_shortest_paths_visitor() );

	std::cout << "#" << endl;
	std::cout << "# ----------------------------------" << endl;
	std::cout << "#" << endl;
	std::cout << "SPP with time windows:" << std::endl;
	std::cout << "Number of optimal solutions: ";
	std::cout << static_cast<int>( opt_solutions_spptw_single.size() ) << std::endl;
	for( int i = 0; i < static_cast<int>( opt_solutions_spptw_single.size() ); ++i )
	{
		std::cout << "The " << i << "th shortest path from " << name[s] << " to " << name[t] << " is: ";
		std::cout << std::endl;
		for( int j = static_cast<int>( opt_solutions_spptw_single[i].size() ) - 1;
				j >= 0;
				--j )
			std::cout << name[source( opt_solutions_spptw_single[i][j], g )] << std::endl;
		std::cout << name[t] << std::endl;
		std::cout << "Length: " << pareto_opt_rcs_spptw_single[i].cost << std::endl;
		std::cout << "Time: " << pareto_opt_rcs_spptw_single[i].time << std::endl;
	}
	std::cout << "#" << endl;
	std::cout << "# ----------------------------------" << endl;
	std::cout << "#" << endl;



	// spptw
	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;

	r_c_shortest_paths_several_sinks(
			g,
			get( &Vertex_Properties::num, g ),
			get( &Arc_Properties::num, g ),
			s,
			tVec,
			opt_solutions_spptw,
			pareto_opt_rcs_spptw,
			spp_spptw_res_cont( 0, 0 ),
			ref_spptw(),
			dominance_spptw(),
			std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
			boost::default_r_c_shortest_paths_visitor() );

	std::cout << "#" << endl;
	std::cout << "# ----------------------------------" << endl;
	std::cout << "#" << endl;
	std::cout << "SPP with time windows:" << std::endl;
	std::cout << "Number of optimal solutions: ";
	std::cout << static_cast<int>( opt_solutions_spptw.size() ) << std::endl;
	for( int i = 0; i < static_cast<int>( opt_solutions_spptw.size() ); ++i )
	{
		std::cout << "The " << i << "th shortest path from A to sink " << name[target( opt_solutions_spptw[i][0], g )]<< " is: ";
		std::cout << std::endl;
		for( int j = static_cast<int>( opt_solutions_spptw[i].size() ) - 1;
				j >= 0;
				--j ){
			std::cout << name[source( opt_solutions_spptw[i][j], g )] << std::endl;
			if(j==0){
				std::cout << name[target(opt_solutions_spptw[i][j], g)] << std::endl;
			}

		}
		std::cout << "Length: " << pareto_opt_rcs_spptw[i].cost << std::endl;
		std::cout << "Time: " << pareto_opt_rcs_spptw[i].time << std::endl;
	}
	std::cout << "#" << endl;
	std::cout << "# ----------------------------------" << endl;
	std::cout << "#" << endl;

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

	getchar();

}
