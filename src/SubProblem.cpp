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

SubProblem::SubProblem(Scenario * scenario, Demand * demand, const Contract * contract, vector<State>* pInitState):
	pScenario_(scenario), pDemand_(demand), pContract_ (contract),
	CDMin_(contract->minConsDaysWork_), maxRotationLength_(demand->nbDays_), nDays_(demand->nbDays_){

	init(pInitState);

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
	//printShortSucc();

}

SubProblem::~SubProblem(){}

// Initialization function
void SubProblem::init(vector<State>* pInitState){

	// Maximum number of consecutive days worked by a nurse ending at day -1
	//
	maxOngoingDaysWorked_ = 0;
	for(int i=0; i<pInitState->size(); i++){
		maxOngoingDaysWorked_ = max( (pInitState->at(i)).consDaysWorked_, maxOngoingDaysWorked_ );
	}

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
						/* Antoine + Samuel ( re-modif) */
						if(newSh == lastSh){	// BUT : add the cost if longer than the maximum allowed
							newNLast += nLast;
							if(newNLast >= pScenario_->maxConsShifts_[newSh]){
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
bool SubProblem::solve(LiveNurse* nurse, Costs * costs, vector<SolveOption> options, set<pair<int,int> > forbiddenDayShifts, bool optimality, int maxRotationLength){


	std::cout << "# Solving subproblem for nurse " << nurse->name_ << " (id:" <<  nurse->id_ << "), " << pContract_->name_ << " [completeWeekends=";
	if(pContract_->needCompleteWeekends_ == 1) std::cout << "YES"; else std::cout << "NO";
	std::cout << "]" << std::endl;

	// Set to true if you want to display contract + preferences (for debug)
	if(false){
		std::cout << "# Preferences:" << endl;
		for(map<int,set<int> >::iterator it = nurse->pWishesOff_->begin(); it != nurse->pWishesOff_->end(); ++it){
			cout <<  "      | " << it->first << "  ->  ";
			for(int s : it->second) cout << pScenario_->intToShift_[s];
			cout << endl;
		}
		std::cout << "# Contract :   ";
		for(int s=1; s<pScenario_->nbShifts_; s++){
			std::cout << pScenario_->intToShift_[s] << " [" << pScenario_->minConsShifts_[s] << "<" << pScenario_->maxConsShifts_[s] << "]   ";
		}
		std::cout << std::endl;
		std::cout << "# " << std::endl;
		std::cout << "# " << std::endl;
	}

	// Get the parameters informations
	//
	setSolveOptions(options);

	// Maximum rotation length
	//
	maxRotationLength_ = min(nDays_, maxRotationLength);

	// Reset authorizations if needed
	//
	if(isOptionActive(SOLVE_FORBIDDEN_RESET)) resetAuthorizations();

	// Delete all previous solutions if needed
	//
	if(isOptionActive(SOLVE_SOLUTIONS_RESET)) resetSolutions();

	// Basic data (nurse, reduced costs, maximum rotation length)
	//
	pLiveNurse_ = nurse;
	pCosts_ = costs;


	// If needed, generate other forbidden day-shift list
	//
	set<pair<int,int> > forbiddenDayShiftsUsed = forbiddenDayShifts;
	if(isOptionActive(SOLVE_FORBIDDEN_RANDOM)) forbiddenDayShiftsUsed = randomForbiddenShifts(25);

	// Initialize structures
	//
	initStructuresForSolve(nurse, costs, forbiddenDayShiftsUsed, maxRotationLength);

	// If needed, generate other costs
	if(isOptionActive(SOLVE_COST_RANDOM)) generateRandomCosts(0, 0);

	// Forbid new arcs, and update the costs
	//
	forbid(forbiddenDayShifts);
	updateArcCosts();

	// SOLUTION OF THE PROBLEM -> DEPENDS ON THE CHOSEN OPTION

	vector< vector< boost::graph_traits<Graph>::edge_descriptor> > opt_solutions_spptw;
	vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;
	spp_spptw_res_cont rc (0,0);

	//cout << "# Serie en cours : " << pLiveNurse_->pStateIni_->consShifts_ << " de " << pScenario_->intToShift_[pLiveNurse_->pStateIni_->shift_];
	//cout << " (total " << pLiveNurse_->pStateIni_->consDaysWorked_ << " jours)" << endl;

	/*
	for(int a=0; a<nArcs_; a++){
		if(!isArcForbidden(a)){
			double c = arcCost(a) + .5;
			cout << "# " << c << endl;
			updateCost(a, c);
		}
	}
	*/

	// "Original way of doing" = pareto optimal for the whole month
	//
	if(isOptionActive(SOLVE_SINGLE_SINKNODE)){
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
	}

	// Gather all the pareto-fronts that correspond to the different last worked days
	//
	else if(isOptionActive(SOLVE_ONE_SINK_PER_LAST_DAY)){

		// 1. Forbid all ending arcs
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int k=CDMin_-1; k<nDays_; k++){
				forbidArc( arcsPrincipalToRotsizein_[s][k] );
			}
		}
		// 2. For every end-day, get the pareto-front
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int s=1; s<pScenario_->nbShifts_; s++) authorizeArc( arcsPrincipalToRotsizein_[s][k] );

			vector< vector< boost::graph_traits<Graph>::edge_descriptor> > last_day_opt_solutions_spptw;
			vector<spp_spptw_res_cont> last_day_pareto_opt_rcs_spptw;
			spp_spptw_res_cont last_day_rc (0,0);
			r_c_shortest_paths(
					g_,
					get( &Vertex_Properties::num, g_ ),
					get( &Arc_Properties::num, g_ ),
					sourceNode_,
					sinkNode_,
					last_day_opt_solutions_spptw,
					last_day_pareto_opt_rcs_spptw,
					last_day_rc,
					ref_spptw(),
					dominance_spptw(),
					std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
					boost::default_r_c_shortest_paths_visitor()
			);
			Tools::push_several_back(& opt_solutions_spptw, last_day_opt_solutions_spptw);
			Tools::push_several_back(& pareto_opt_rcs_spptw, last_day_pareto_opt_rcs_spptw);
			for(int s=1; s<pScenario_->nbShifts_; s++) forbidArc( arcsPrincipalToRotsizein_[s][k] );

		}
		// 3. Re-authorize all ending arcs
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int k=CDMin_-1; k<nDays_; k++){
				authorizeArc( arcsPrincipalToRotsizein_[s][k] );
			}
		}
	}

	// Gather all the pareto-fronts that correspond to the different first worked days
	//
	else if(isOptionActive(SOLVE_ONE_SINK_PER_FIRST_DAY)){

		// 1. Store a vector with all initially allowed short arcs AND forbid all of them
		vector3D initiallyAuthorized;
		for(int s=0; s<pScenario_->nbShifts_; s++){
			vector2D initiallyAuthorized2;
			for(int k=0; k<nDays_; k++){
				vector<int> initiallyAuthorized1;
				for(int n=1; n<=maxvalConsByShift_[s]; n++){
					int a = arcsFromSource_[s][k][n];
					if(!isArcForbidden(a)){
						initiallyAuthorized1.push_back( a );
						forbidArc(a);
					}
				}
				initiallyAuthorized2.push_back(initiallyAuthorized1);
			}
			initiallyAuthorized.push_back(initiallyAuthorized2);
		}

		// 2. For every first day, get the pareto-front
		for(int k=CDMin_-1; k<nDays_; k++){

			// Re-authorize the arcs starting on day (k-CDMin_)
			for(int s=1; s<pScenario_->nbShifts_; s++)
				for(int a: initiallyAuthorized[s][k])
					authorizeArc(a);

			// Get the pareto-front and add it to the solution list
			vector< vector< boost::graph_traits<Graph>::edge_descriptor> > first_day_opt_solutions_spptw;
			vector<spp_spptw_res_cont> first_day_pareto_opt_rcs_spptw;
			spp_spptw_res_cont first_day_rc (0,0);
			r_c_shortest_paths(
					g_,
					get( &Vertex_Properties::num, g_ ),
					get( &Arc_Properties::num, g_ ),
					sourceNode_,
					sinkNode_,
					first_day_opt_solutions_spptw,
					first_day_pareto_opt_rcs_spptw,
					first_day_rc,
					ref_spptw(),
					dominance_spptw(),
					std::allocator< boost::r_c_shortest_paths_label< Graph, spp_spptw_res_cont> >(),
					boost::default_r_c_shortest_paths_visitor()
			);
			Tools::push_several_back(& opt_solutions_spptw, first_day_opt_solutions_spptw);
			Tools::push_several_back(& pareto_opt_rcs_spptw, first_day_pareto_opt_rcs_spptw);

			// Forbid the arcs starting on that day
			for(int s=1; s<pScenario_->nbShifts_; s++)
				for(int a: initiallyAuthorized[s][k])
					forbidArc(a);
		}

		// 3. Re-set all authorizations as they previsouly were
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int k=CDMin_-1; k<nDays_; k++){
				for(int a: initiallyAuthorized[s][k]) authorizeArc(a);
			}
		}
	}

	// Return TRUE if a rotation was added. Last argument is true IF only negative reduced cost rotations are added
	//
	return addRotationsFromPaths(opt_solutions_spptw, pareto_opt_rcs_spptw);
}

// Store the options in a readable way
void SubProblem::setSolveOptions(vector<SolveOption> options){

	activeOptions_.assign(100, false);
	for(SolveOption o : options) activeOptions_[o] = true;

	// Check for incompatible options...
	//std::cout << "# SOLVE OPTIONS:" << std::endl;
	for(vector<SolveOption> cluster : incompatibilityClusters){
		int n=0;
		SolveOption chosenOption;
		for(SolveOption o : cluster){
			if(isOptionActive(o)){
				chosenOption = o;
				n++;
			}
		}
		if(n==0){
			activeOptions_[cluster[0]] = true;
			//std::cout << "#     [Default]: " << solveOptionName[cluster[0]] << std::endl;
		} else if (n==1){
			//std::cout << "#     [Chosen] : " << solveOptionName[chosenOption] << std::endl;
		} else {
			std::cout << "# WARNING !! TOO MANY OPTIONS -> [Default] " << solveOptionName[cluster[0]] << std::endl;
			activeOptions_[cluster[0]] = true;
			for(int i=1; i<cluster.size(); i++) activeOptions_[cluster[i]] = false;
		}
	}

	//printActiveSolveOptions();
}

// Transforms the solutions found into proper rotations.
//
bool SubProblem::addRotationsFromPaths(vector< vector< boost::graph_traits<Graph>::edge_descriptor> > paths, vector<spp_spptw_res_cont> resources){
	bool oneFound = false;
	// For each path of the list, record the corresponding rotation (if negativeOnly=true, do it only if the dualCost < 0)
	for(int p=0; p < paths.size(); ++p){
		Rotation rot = rotationFromPath(paths[p], resources[p]);
		if( isOptionActive(SOLVE_NEGATIVE_ALLVALUES)
				or ( isOptionActive(SOLVE_NEGATIVE_ONLY) and (rot.dualCost_ < -EPSILON) ) ){
			theRotations_.push_back(rot);
			nPaths_ ++;
			oneFound = true;
			//printPath(paths[p], resources[p]);
			printRotation(rot);
		}
	}
	//printAllRotations();
	return oneFound;
}

// Adds a rotation made from the given path to the current list of answers and increases their counter
//
Rotation SubProblem::rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource){

	int firstDay = -1;
	double dualCost = 0;
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
void SubProblem::addSingleArc(int o, int d, double baseCost, int t, ArcType type){
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
         /*
          * Antoine modif: on paie juste les jours en plus
          */
			addSingleArc(origin, destin, 0, 0, SHIFT_TO_ENDSEQUENCE);
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
	preferencesCosts_.clear();
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

	map<int,int> specialArcsSuccId;
	map<int,double> specialArcsCost;

	for(int s=1; s<pScenario_->nbShifts_; s++){
		for(int k=CDMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){

				idBestShortSuccCDMin_[s][k][n] = -1;
				arcCostBestShortSuccCDMin_[s][k][n] = MAX_COST;
				for(int i=0; i<(allShortSuccCDMinByLastShiftCons_[s][n]).size(); i++){
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

				// IF NO VALID SUCCESSION, THEN FORBID THE ARC
				int a = arcsFromSource_[s][k][n];
				if(arcCostBestShortSuccCDMin_[s][k][n] == MAX_COST){
					forbidArc( a );
				}
			}
		}
	}

	// FOR THE SHIFTS ON THE FIRST DAY THAT EXTEND THE ONGOING WORK AT INITIAL STATE
	//
	map<int,double>::iterator itCost = specialArcsCost.begin();
	for(map<int,int>::iterator itId = specialArcsSuccId.begin(); itId != specialArcsSuccId.end(); ++itId){
		int a = itId->first;
		int d = arcDestination(a);
		int s = principalToShift_[d];
		int k = principalToDay_[d];
		int n = principalToCons_[d];
		authorizeArc(a);
		idBestShortSuccCDMin_[s][k][n] = itId->second;
		arcCostBestShortSuccCDMin_[s][k][n] = itCost->second;
		++ itCost;
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
		int nConsFirstShift = 0; int ii=0; while(succ[ii]==firstShift){nConsFirstShift ++; ii++;}


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
				int diff = pScenario_->minConsShifts_[shiftIni] - nConsShiftIni;
				ANS += max(0, diff*(WEIGHT_CONS_SHIFTS));
			}

			// b. (ii)  The nurse was working on the same shift AND the short rotation contains other shifts (easy case for add/subtract)
			//            - Subtract the cost due to the consecutive beginning
			//            - Subtract the cost due to the consecutive end of the initial state
			//            - Add the consecutive cost for all shifts
			else if(nConsFirstShift < CDMin_) {
				int diffShift = nConsShiftIni - pScenario_->maxConsShifts_[shiftIni];
				ANS -= max(0, diffShift*WEIGHT_CONS_SHIFTS);
				ANS -= consShiftCost(firstShift, nConsFirstShift);
				ANS += consShiftCost(firstShift, (nConsFirstShift + nConsShiftIni));
			}

			// b. (iii) The nurse was working on the same shift AND the short rotation only contains that shift (recompute the cost -> easier)
			else {
				ANS -= baseArcCostOfShortSucc_[size][succId];
				if( (nConsFirstShift + nConsShiftIni) >= maxvalConsByShift_[shiftIni] )
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
	for(int s=1; s<pScenario_->nbShifts_; s++){
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
	for(int s1=1; s1<pScenario_->nbShifts_; s1++)
		for(int s2=1; s2<pScenario_->nbShifts_; s2++)
			for(int k=CDMin_-1; k<nDays_-1; k++){
				int a = arcsShiftToNewShift_[s1][s2][k];
				if(a > 0){
					double c = arcBaseCost_[a];
					c += preferencesCosts_[k+1][s2] ;
					c -= pCosts_->dayShiftWorkCost(k+1,s2-1);
					c -= Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
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
				c -= pCosts_->dayShiftWorkCost(k+1,s-1);
				if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
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
			c -= pCosts_->dayShiftWorkCost(k+1,s-1);
			if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
			updateCost( a , c );
		}

	// F. ARCS : PRINCIPAL_TO_ROTSIZE [baseCost contains complete weekend constraint]
	//
	for(int s=1; s<pScenario_->nbShifts_; s++)
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
	// If the succession with the previous shift (day -1) is not allowed
	if(firstDay==0 and pScenario_->isForbiddenSuccessor(succ[0],pLiveNurse_->pStateIni_->shift_))
		return false;
	// If some day-shift is forbidden...
	for(int i=0; i<succ.size(); i++){
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
	arcCost(a) < 10000 ? rep << "     " : rep << " ";
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
void SubProblem::printPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_spptw_res_cont resource){

	// The successive nodes, and corresponding arc costs / time
	//
	for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
		int a = boost::get(&Arc_Properties::num, g_, path[j]);
		std::cout << "# \t| [ " << shortNameNode(source( path[j], g_ )) << " ]";
		std::cout << "\t\tCost:  " << arcCost(a) << "\t\tTime:" << arcLength(a);
		std::cout << "\t\t[" << (arcStatus_[a] ? " allowed " : "forbidden") << "]" << std::endl;
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
			int succId = shortSuccCDMinIdFromArc_.at(a);
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

// Prints all active solving options
void SubProblem::printActiveSolveOptions(){
	std::cout << "# CHOSEN SOLVING OPTIONS:" << std::endl;
	for(int i=0; i<activeOptions_.size(); i++){
		if(activeOptions_[i])
			std::cout << "#     | " << solveOptionName[i] << std::endl;
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
