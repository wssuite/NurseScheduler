/*
 * SubProblem.cpp
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#include "solvers/mp/SubProblem.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


using std::stringstream;
using std::vector;


static int MAX_COST = 99999;


//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
SubProblem::SubProblem() {}

SubProblem::SubProblem(Scenario * scenario, int nDays, const Contract * contract, vector<State>* pInitState):
					pScenario_(scenario), nDays_(nDays), pContract_ (contract),
					CDMin_(contract->minConsDaysWork_), daysMin_(1), nLabels_(2),
					maxRotationLength_(nDays) {
  init(pInitState);
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
				       pScenario_->maxConsShiftsOf(sh) >= nDays_ + maxOngoingDaysWorked_
				       or pScenario_->maxConsShiftsOf(sh) >= NB_SHIFT_UNLIMITED
		);
		int nl = isUnlimited_[sh] ? pScenario_->minConsShiftsOf(sh) : pScenario_->maxConsShiftsOf(sh);
		maxvalConsByShift_.push_back( nl );
	}

	// // id and arcCost of best succession (given a triplet s,k,n)
	// idBestShortSuccCDMin_.clear();
	// arcCostBestShortSuccCDMin_.clear();
	// for(int s=0; s<pScenario_->nbShiftsType_; s++){
	// 	vector2D v2; vector<vector<double> > w2;
	// 	int n = maxvalConsByShift_[s]+1;
	// 	Tools::initVector2D(&v2, nDays_, n);
	// 	idBestShortSuccCDMin_.push_back(v2);
	// 	Tools::initDoubleVector2D(&w2, nDays_, n);
	// 	arcCostBestShortSuccCDMin_.push_back(w2);
	// }

	preferencesCosts_.clear();
	for(int k=0; k<nDays_; k++){
		vector<double> v;
		// for(int s=0; s<pScenario_->nbShiftsType_; s++) v.push_back(0);
		for(int s=0; s<pScenario_->nbShifts_; s++) v.push_back(0);
		preferencesCosts_.push_back(v);
	}

	nLongFound_=0;
	// nVeryShortFound_=0;
}

void SubProblem::build() {

  g_ = RCGraph(nDays_);

  //	initShortSuccessions();

  createNodes();
  createArcs();

  // Set all status to authorized
  for(int k=0; k<nDays_; k++){
    vector<bool> v;
    for(int s=0; s<pScenario_->nbShifts_; s++)
      v.push_back(true);
    dayShiftStatus_.push_back(v);
  }
  for(int k=0; k<nDays_; k++) startingDayStatus_.push_back(true);

  nPathsMin_ = 0;

  //std::cout << "# A new subproblem has been created for contract " << contract->name_ << std::endl;

  //printGraph();

  timeInS_ = new Tools::Timer(); timeInS_->init();
  timeInNL_ = new Tools::Timer(); timeInNL_->init();
}

// Cost function for consecutive identical shifts
//
double SubProblem::consShiftCost(int sh, int n){
  if(pScenario_->minConsShiftsOfTypeOf(sh) - n > 0) return (WEIGHT_CONS_SHIFTS * ( pScenario_->minConsShiftsOfTypeOf(sh) - n ) );
  if(n - pScenario_->maxConsShiftsOfTypeOf(sh) > 0) return (WEIGHT_CONS_SHIFTS * ( n - pScenario_->maxConsShiftsOfTypeOf(sh) ) );
	return 0;
}

double SubProblem::consShiftTypeCost(int sh, int n){
  if(pScenario_->minConsShiftsOf(sh) - n > 0) return (WEIGHT_CONS_SHIFTS * ( pScenario_->minConsShiftsOf(sh) - n ) );
  if(n - pScenario_->maxConsShiftsOf(sh) > 0) return (WEIGHT_CONS_SHIFTS * ( n - pScenario_->maxConsShiftsOf(sh) ) );
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
	//	nVeryShortFound_=0;										// Initialize number of solutions found at 0 (short rotations)
	forbid(forbiddenDayShifts);								// Forbid arcs
	forbidStartingDays(forbiddenStartingDays);				// Forbid starting days

	if(false) printContractAndPrefenrences();				// Set to true if you want to display contract + preferences (for debug)

	// cout << " 1   ----------------------------------------------" << endl;
	// printGraph();

	// timeInS_->start();
	// bool ANS_short = solveShortRotations();
	// timeInS_->stop();

	// cout << " 2   ----------------------------------------------" << endl;
	// printGraph();
	
	timeInNL_->start();
	bool ANS_long = solveLongRotations(optimality);
	timeInNL_->stop();

	// cout << " 3   ----------------------------------------------" << endl;
	// printGraph();

	// printAllRotations();

	// Reset authorizations
	authorize(forbiddenDayShifts);
	authorizeStartingDays(forbiddenStartingDays);

	//	return ANS_short or ANS_long;
	return ANS_long;
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
	std::vector<RCSolution> solutions = g_.solve(nLabels_, maxReducedCostBound_);
//	if(!(param_.oneSinkNodePerLastDay_)){

  for(const RCSolution& sol: solutions) {
    Rotation rot = buildColumn(sol);
    theRotations_.push_back(rot);
    nPaths_ ++;
    nLongFound_++;
    bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
  }

		return !solutions.empty();
}

// Transforms a solution found into a proper column.
//
Rotation SubProblem::buildColumn(const RCSolution& sol){
  Rotation rot (sol.firstDay, sol.shifts, pLiveNurse_->id_, MAX_COST, sol.cost);
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
	int nShiftsType= pScenario_->nbShiftsType_;			// Number of different shifts

	// INITIALIZATION
	initNodesStructures();

	// 1. SOURCE NODE
	//
	int v = g_.addSingleNode(SOURCE_NODE, {0,0}, {maxRotationLength_, CDMin_});
	g_.setSource(v);

	// 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
	//
	for(int sh=1; sh<nShiftsType; sh++){		// For each possible worked shift
	  for(int k=daysMin_-1; k<nDays_; k++){		// For each date
	    //	for(int k=pContract_->minConsDaysWork_-1; k<nDays_; k++){	// For each date
			for(int cons=1; cons<=maxvalConsByShift_[sh]; cons++){		// For each level of network
				addNodeToPrincipalNetwork(sh, k, cons);		// Add a node to the principal network
			}
		}
	}

	// 3. ROTATION LENGTH CHECK
	//
	// For each of the days, do a rotation-length-checker
	for(int k=0; k<nDays_; k++){
	  // One node for the entrance in subnetwork per day
		v = g_.addSingleNode(ROTATION_LENGTH_ENTRANCE, {0,0}, {maxRotationLength_, CDMin_});
    rotationLengthEntrance_.push_back(v);
		map<int,int> checkNodesForThatDay;
		// Check nodes
		for(int l=0; l<=maxRotationLength_; l++){		// Check nodes: from CD_max (longest free) to maximum rotation length, for each day
			v = g_.addSingleNode(ROTATION_LENGTH, {0, 0}, {l, max(0, CDMin_-l)});
      checkNodesForThatDay.insert(pair<int,int>(l,v));
      rotationLengthNodesLAT_.insert(pair<int,int>(v,l));
		}
		rotationLengthNodes_.push_back(checkNodesForThatDay);
		// Sink day
    v = g_.addSingleNode(ROTATION_LENGTH_EXIT, {0,0}, {maxRotationLength_, CDMin_});
    rotationLengthExit_.push_back(v);
    // Daily sink node
    if(param_.oneSinkNodePerLastDay_)
      g_.addSink(v);
  }

	// 4. SINK NODE
	//
	v = g_.addSingleNode(SINK_NODE, {0,0}, {maxRotationLength_, CDMin_});
  if(!param_.oneSinkNodePerLastDay_)
    g_.addSink(v);

	
//	for(int sh=1; sh<nShiftsType; sh++){		// For each possible worked shift
//	  for(int k=daysMin_-1; k<nDays_; k++){		// For each date
//	    //	for(int k = daysMin_-1; k<pContract_->minConsDaysWork_-1; k++){	// For each date
//			for(int cons=1; cons<=maxvalConsByShift_[sh]; cons++){		// For each level of network
//				addNodeToPrincipalNetwork(sh, k, cons);		// Add a node to the principal network
//			}
//		}
//	}
}

// Initiate variables for the nodes structures (vectors, etc.)
void SubProblem::initNodesStructures(){

	// Data
	//
	int nShiftsType= pScenario_->nbShiftsType_;				// Number of different shifts
	principalNetworkNodes_.clear();
	principalToShift_.clear();
	principalToDay_.clear();
	principalToCons_.clear();
	rotationLengthEntrance_.clear();
	rotationLengthNodes_.clear();
	rotationLengthNodesLAT_.clear();
  rotationLengthExit_.clear();

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

// Add a node to the principal network of the graph, for shift sh, day k, and number of consecutive similar shifts cons
void SubProblem::addNodeToPrincipalNetwork(int sh, int k, int cons){
	// Create the node
	//
	int v = g_.addSingleNode(PRINCIPAL_NETWORK, {0,0}, {maxRotationLength_, CDMin_});

  // Store its ID in the vector3D
  //
  principalNetworkNodes_[sh][k][cons] = v;

  // Store the information backwards
  //
  principalToShift_.insert(pair<int,int>(v, sh));
  principalToDay_.insert(pair<int,int>(v, k));
  principalToCons_.insert(pair<int,int>(v, cons));
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

// Initializes the data structures used for the arcs
void SubProblem::initArcsStructures(){
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
	  vector2D vv2, w2, x2, zz2;
	  vector3D v2, y2, z2;
	  // vector2D v2, w2, x2, zz2;
	  // vector3D y2, z2;
	  unsigned int nShifts = pScenario_->shiftTypeIDToShiftID_[sh].size();

	  //	  Tools::initVector2D(&v2, nDays_, maxvalConsByShift_[sh]+1);           arcsFromSource_.push_back(v2);
	  Tools::initVector2D(&w2, nDays_, maxvalConsByShift_[sh]+1);           arcsShiftToEndsequence_.push_back(w2);
	  Tools::initVector2D(&x2, nDays_, nShifts);                            arcsRepeatShift_.push_back(x2);
	  
	  Tools::initVector3D(&y2, nDays_, maxvalConsByShift_[sh]+1, nShifts);  arcsShiftToSameShift_.push_back(y2);
	  Tools::initVector3D(&v2, nDays_, maxvalConsByShift_[sh]+1, nShifts);  arcsFromSource_.push_back(v2);

	  // int nNewShifts = pScenario_->shiftTypeIDToShiftID_[sh].size();
	  // Tools::initVector3D(&z2, pScenario_->nbShiftsType_, nDays_, nNewShifts); arcsShiftToNewShift_.push_back(z2);
	  // // Tools::initVector3D(&z2, pScenario_->nbShiftsType_, nDays_, nShifts); arcsShiftToNewShift_.push_back(z2);
	  
	  for(int sh2=0; sh2<pScenario_->nbShiftsType_; sh2++){
	    int nNewShifts = pScenario_->shiftTypeIDToShiftID_[sh2].size();
	    Tools::initVector2D(&zz2, nDays_, nNewShifts);     z2.push_back(zz2);
	  }
	  arcsShiftToNewShift_.push_back(z2);
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
    for(int k=daysMin_-1; k<nDays_; k++){
      for(int nCons=1; nCons<=maxvalConsByShift_[sh]; nCons ++){
	origin = g_.source();
	destin = principalNetworkNodes_[sh][k][nCons];

	for (int s=0; s < pScenario_->shiftTypeIDToShiftID_[sh].size(); s++) {
	  int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
	  int a = g_.addSingleArc(origin, destin, 0, {daysMin_, CDMin_-daysMin_}, SOURCE_TO_PRINCIPAL, k, {shiftID});
    arcsFromSource_[sh][k][nCons][s] = a;
  }
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
		for(int k=daysMin_-1; k<nDays_-1; k++){

			//   1. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS NOT REACHED YET
			//
			for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons ++){
				origin = principalNetworkNodes_[sh][k][nCons];
				destin = principalNetworkNodes_[sh][k+1][nCons+1];

				for (unsigned int s = 0; s < nShifts; s++) {
				  int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
				  int a = g_.addSingleArc(origin, destin, 0, {1,-1}, SHIFT_TO_SAMESHIFT, k+1, {shiftID});
          arcsShiftToSameShift_[sh][k][nCons][s] = a;
        }
			}

			//   2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
			//
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = principalNetworkNodes_[sh][k+1][maxvalConsByShift_[sh]];
			double cost = isUnlimited(sh) ? 0 : WEIGHT_CONS_SHIFTS;

			for (unsigned int s = 0; s < nShifts; s++) {
			  int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
			  int a = g_.addSingleArc(origin, destin, cost, {1,-1}, REPEATSHIFT, k+1, {shiftID});
        arcsRepeatShift_[sh][k][s] = a;
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
            arcsShiftToNewShift_[sh][newSh][k][s] =
                g_.addSingleArc(origin, destin, 0, {1,-1}, SHIFT_TO_NEWSHIFT, k+1, {newShiftID});
			    }
			  }
			}

			// 4. END THE CONSECUTIVE SHIFT SEQUENCE
			//
			for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons++){
				origin = principalNetworkNodes_[sh][k][nCons];
				destin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
				arcsShiftToEndsequence_[sh][k][nCons] =
				    g_.addSingleArc(origin, destin, consShiftTypeCost(sh, nCons), {0,0}, SHIFT_TO_ENDSEQUENCE, k);
			}
		}

		// SPECIAL CASE: LAST DAY
		//
		for(int nCons=1; nCons<maxvalConsByShift_[sh]; nCons++){
			origin = principalNetworkNodes_[sh][nDays_-1][nCons];
			destin = principalNetworkNodes_[sh][nDays_-1][maxvalConsByShift_[sh]];
			arcsShiftToEndsequence_[sh][nDays_-1][nCons] =
			    g_.addSingleArc(origin, destin, 0, {0,0}, SHIFT_TO_ENDSEQUENCE, nDays_-1);
		}
	}
}

// Create all arcs that involve the rotation size checking subnetwork (incoming, internal, and exiting that subnetwork)
void SubProblem::createArcsAllRotationSize(){

	int nShiftsType = pScenario_->nbShiftsType_;
	int origin, destin;

	// 1. ALL INCOMING ARCS
	//
	for(int sh=1; sh<nShiftsType; sh++){		// For all shifts
		for(int k=daysMin_-1; k<nDays_; k++){				// For all days
			origin = principalNetworkNodes_[sh][k][maxvalConsByShift_[sh]];
			destin = rotationLengthEntrance_[k];
			arcsPrincipalToRotsizein_[sh][k] =
			    g_.addSingleArc(origin, destin, 0, {0,0}, PRINCIPAL_TO_ROTSIZE, k);	// Allow to stop rotation that day
		}
	}

	// 2. ALL INTERNAL ARCS
	//
	for(int k=daysMin_-1; k<nDays_; k++){

		map<int,int> rotLengthNodesForDay = rotationLengthNodes_[k];
		map<int,int> arcsRotsizeinToRotsize;
		map<int,int> arcsRotsizeToRotsizeout;
		for(map<int,int>::iterator itRLN = rotLengthNodesForDay.begin(); itRLN != rotLengthNodesForDay.end(); ++itRLN){
			// From entrance of that day to checknode
			origin = rotationLengthEntrance_[k];
			destin = itRLN->second;
			int c = consDaysCost(itRLN->first);
			int a = g_.addSingleArc(origin, destin, c, {0,0}, ROTSIZEIN_TO_ROTSIZE, k);
      arcsRotsizeinToRotsize.insert(pair<int,int>(itRLN->first, a));
			// From checknode to exit of that day
			origin = itRLN->second;
			destin = rotationLengthExit_[k];
			a = g_.addSingleArc(origin, destin, 0, {0,0}, ROTSIZE_TO_ROTSIZEOUT, k);
      arcsRotsizeToRotsizeout.insert(pair<int,int>( itRLN->first, a));
		}

		arcsRotsizeinToRotsizeDay_.push_back(arcsRotsizeinToRotsize);
		arcsRotsizeToRotsizeoutDay_.push_back(arcsRotsizeToRotsizeout);

		// link all end rotation nodes to the main sink node
    if(!param_.oneSinkNodePerLastDay_) {
      origin = rotationLengthExit_[k];
      destin = g_.sink();
      g_.addSingleArc(origin, destin, 0, {0,0}, ROTSIZEOUT_TO_SINK, k);
    }

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
	for(map<int,vector<Wish> >::iterator it = pLiveNurse_->pWishesOff_->begin(); it != pLiveNurse_->pWishesOff_->end(); ++it){
		for(Wish s : it->second){
			preferencesCosts_[it->first][s.shift] = WEIGHT_PREFERENCES_OFF[s.level];

		}
	}

	for(map<int,vector<Wish> >::iterator it = pLiveNurse_->pWishesOn_->begin(); it != pLiveNurse_->pWishesOn_->end(); ++it){
		for(Wish s : it->second){
			preferencesCosts_[it->first][s.shift] = WEIGHT_PREFERENCES_ON[s.level];

		}
	}

	// // id and arcCost of best succession (given a triplet s,k,n)
	// idBestShortSuccCDMin_.clear();
	// arcCostBestShortSuccCDMin_.clear();
	// for(int s=0; s<pScenario_->nbShiftsType_; s++){
	// 	vector2D v2; vector<vector<double> > w2;
	// 	int n = maxvalConsByShift_[s]+1;
	// 	Tools::initVector2D(&v2, nDays_, n, -1);
	// 	idBestShortSuccCDMin_.push_back(v2);
	// 	Tools::initDoubleVector2D(&w2, nDays_, n, MAX_COST);
	// 	arcCostBestShortSuccCDMin_.push_back(w2);
	// }
}

bool SubProblem::feasibleArcSource(int k, int n, int shiftID) {
    // check if the nurse can perform this shift on this day
    if(!canSuccStartHere({shiftID}, k))
        return false;

    // check if the level n can be reached
    if(k==0) {
        int firstShiftType = pScenario_->shiftIDToShiftTypeID_[shiftID];
        int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_,
            maxConsShift = pScenario_->maxConsShiftsOf(firstShiftType);
        // if the number of consecutive shift of same type is lower than the max
        if(nConsShiftIni < maxConsShift)
            return (nConsShiftIni+1 == n);
        // otherwise n should be equal to maxCons
        return (n == maxConsShift);
    }
    // else, can just go to the first level n=1
    return (n == 1);
}

double SubProblem::costArcSource(int a, int k, int shiftID){
  double  cost = g_.arcInitialCost(a);
  
  cost += preferencesCosts_[k][shiftID] ;
  cost -= pCosts_->dayShiftWorkCost(k, shiftID);
  if(Tools::isWeekend(k)) cost -= pCosts_->workedWeekendCost();
  cost -= pCosts_->startWorkCost(k);
  cost += startWeekendCosts_[k];

  // if first day, take into account historical state
  if (k == 0) {
    int shiftIni = pLiveNurse_->pStateIni_->shift_;
    int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
    int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

    // check if the nurse is changing of shift type
    int shiftTypeIni = pScenario_->shiftIDToShiftTypeID_[shiftIni];
    int firstShiftType = pScenario_->shiftIDToShiftTypeID_[shiftID];

    // 1. The nurse was resting: pay more only if the rest is too short
    if(shiftTypeIni == 0){
      int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      cost += max(0, diffRest*WEIGHT_CONS_DAYS_OFF);
    }
    // 2. The nurse was working
    else {
      // a. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be added later
      int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
      cost -= max(0, diffWork*WEIGHT_CONS_DAYS_WORK);

      // b. (i)   The nurse was working on a different shift: if too short, add the corresponding cost
      if(shiftTypeIni != firstShiftType){
        int diff = pScenario_->minConsShiftsOfTypeOf(shiftIni) - nConsShiftIni;
        cost += max(0, diff*(WEIGHT_CONS_SHIFTS));
      }

//      // // b. (ii)  The nurse was working on the same shift AND the short rotation contains other shifts (easy case for add/subtract)
//      // //            - Subtract the cost due to the consecutive beginning
//      // //            - Subtract the cost due to the consecutive end of the initial state
//      // //            - Add the consecutive cost for all shifts
//      // else if(nConsFirstShift < pContract_->minConsDaysWork_) {
//      // 	//			else if(nConsFirstShift < CDMin_) {
//      // 	int diffShift = nConsShiftIni - pScenario_->maxConsShiftsOfTypeOf(shiftIni);
//      // 	ANS -= max(0, diffShift*WEIGHT_CONS_SHIFTS);
//      // 	ANS -= consShiftCost(firstShift, nConsFirstShift);
//      // 	ANS += consShiftCost(firstShift, (nConsFirstShift + nConsShiftIni));
//      // }
//
//      // b. (iii) The nurse was working on the same shift AND the short rotation only contains that shift (recompute the cost -> easier)
//      else {
//	// ANS -= baseArcCostOfShortSucc_[size][succId];
//	// if( (nConsFirstShift + nConsShiftIni) >= maxvalConsByShift_[shiftTypeIni] )
//	//   ANS += consShiftCost(shiftIni, (nConsFirstShift + nConsShiftIni));
//
//	cost -= arcBaseCost_[a];
//
//	if ( nConsShiftIni + 1  >= maxvalConsByShift_[shiftTypeIni] )
//	  cost += consShiftCost(shiftIni, ( nConsShiftIni + 1 ));
//      }
    }
  }
  
  return cost;
}


double SubProblem::costArcPrincipal(int a, int k, int shiftID) {
  double  cost = g_.arcInitialCost(a);
  cost += preferencesCosts_[k][shiftID] ;
  cost -= pCosts_->dayShiftWorkCost(k,shiftID);
  if(Tools::isSaturday(k)) cost -= pCosts_->workedWeekendCost();
  return cost;
}

double SubProblem::costArcEnd(int a, int k) {
    double cost = g_.arcInitialCost(a);
    cost += endWeekendCosts_[k];
    cost -= pCosts_->endWorkCost(k);
    return cost;
}

//--------------------------------------------
//
// Functions for the costs
//
//--------------------------------------------

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
		int shiftType = pScenario_->shiftIDToShiftTypeID_[shift];
		int newShiftType = pScenario_->shiftIDToShiftTypeID_[newShift];
		
		if(newShiftType == shiftType){
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
	for(int k=startDate; k<=endDate; k++) dayShiftsRedCost -= pCosts_->dayShiftWorkCost( k, succ[k-startDate] );


	// I. RETURN THE TOTAL COST
	//
	double regCost = consDaysRegCost + consShiftsRegCost + completeWeekendRegCost + preferencesRegCost + shortRestRegCost;
	double redCost = dayShiftsRedCost + startRedCost + endRedCost + weekendRedCost;
	double ANS = regCost + redCost;

	if(false)
	  {
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

// Updates the costs depending on the reduced costs given for the nurse
//
void SubProblem::updateArcCosts(){

  //	A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
	
  for(int sh=1; sh<pScenario_->nbShiftsType_; sh++){
      for(int k=daysMin_-1; k<nDays_-1; k++){
        for(int n=1; n<=maxvalConsByShift_[sh]; n++){
            for (int s=0; s < pScenario_->shiftTypeIDToShiftID_[sh].size(); s++) {
              int a = arcsFromSource_[sh][k][n][s];
              int shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
              if(feasibleArcSource(k, n, shiftID)) {
                  double c = costArcSource(a, k, shiftID);
                  g_.updateCost( a , c );
                  // For an arc that starts on the first day, must update the consumption based on the historical state
                  if(k==0) {
                      int  a = arcsFromSource_[sh][0][n][s];
                      std::vector<int> consumptions = { daysMin_ + pLiveNurse_->pStateIni_->consDaysWorked_,
                                                        CDMin_ - daysMin_ - pLiveNurse_->pStateIni_->consDaysWorked_ };
                      g_.updateConsumptions(a, consumptions);
                  }
              } else g_.forbidArc(a);
            }
      }
    }
  }

	// B. ARCS : SHIFT_TO_NEWSHIFT [baseCost = 0]
	//
	for(int s1=1; s1<pScenario_->nbShiftsType_; s1++)
		for(int s2=1; s2<pScenario_->nbShiftsType_; s2++)
		  if(s2 != s1 and ! pScenario_->isForbiddenSuccessorShiftType_ShiftType(s2,s1)){
			for(int k=daysMin_-1; k<nDays_-1; k++){
			  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[s2].size(); s++) {
				int a = arcsShiftToNewShift_[s1][s2][k][s];
				if(a > 0){
				        int  shiftID = pScenario_-> shiftTypeIDToShiftID_[s2][s];
				        double c = costArcPrincipal(a, k+1, shiftID);
					g_.updateCost( a , c );
				}
			  }
			}
		  }

	// C. ARCS : SHIFT_TO_SAMESHIFT [baseCost = 0]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=daysMin_-1; k<nDays_-1; k++)
			for(int n=1; n<maxvalConsByShift_[sh]; n++){
			    for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
				int  a = arcsShiftToSameShift_[sh][k][n][s];
				int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
				double c = costArcPrincipal(a, k+1, shiftID);
				g_.updateCost( a , c );
			    }
			}

	// D. ARCS : SHIFT_TO_ENDSEQUENCE [They never change]

	// E. ARCS : REPEATSHIFT [baseCost contains consecutive shift cost]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=daysMin_-1; k<nDays_-1; k++){
		  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
			int  a = arcsRepeatShift_[sh][k][s];
			int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
			double c = costArcPrincipal(a, k+1, shiftID);
        g_.updateCost( a , c );
		  }
		}

	// F. ARCS : PRINCIPAL_TO_ROTSIZE [baseCost contains complete weekend constraint]
	//
	for(int s=1; s<pScenario_->nbShiftsType_; s++)
		for(int k=daysMin_-1; k<nDays_; k++){
			int a = arcsPrincipalToRotsizein_[s][k];
            double c = costArcEnd(a, k);
      g_.updateCost( a , c );
		}

	// G. ARCS : ROTSIZEIN_TO_ROTSIZE [baseCost contains rotation length cost. They never change]

	// H. ARCS : ROTSIZE_TO_ROTSIZEOUT [They never change]

	// I. ARCS : ROTSIZEOUT_TO_SINK [Never changes]
}


//--------------------------------------------
//
// Functions to update the maximum length of a rotation
//
//--------------------------------------------

// Updates the maximum arrival time at all nodes so that the maximum length of a rotation is updated.
//   + Updates only
void SubProblem::updatedMaxRotationLengthOnNodes(){
	for(int v=0; v<g_.nodesSize(); v++){
		if(!g_.nodeForbidden(v)){
			if(g_.nodeType(v) != ROTATION_LENGTH){
			    std::vector<int> ubs = g_.nodeUBs(v);
                ubs[MAX_CONS_DAYS] = maxRotationLength_;
				g_.updateUBs(v, ubs);
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

// Authorize the nodes that correspond to forbidden shifts
//
void SubProblem::authorize(set<pair<int,int> > forbiddenDayShifts){
	for(pair<int,int> p : forbiddenDayShifts){
		//std::cout << "# Trying to forbid " << p.first << "-" << pScenario_->intToShift_[p.second].at(0) << endl;
		authorizeDayShift(p.first,p.second);
	}
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
	if(k >= pContract_->minConsDaysWork_-1){
	  //	if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[shiftType]; n++){
			g_.forbidNode( principalNetworkNodes_[shiftType][k][n] );
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
	if(k >= pContract_->minConsDaysWork_-1){
	  //		if(k >= CDMin_-1){
		for(int n=1; n<=maxvalConsByShift_[shiftType]; n++)
      g_.authorizeNode( principalNetworkNodes_[shiftType][k][n] );
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
	for(int s=1; s<pScenario_->nbShifts_; s++)
		for(int k=0; k<nDays_; k++)
			authorizeDayShift(k,s);

		g_.resetAuthorizations();
}

// Generate random forbidden shifts
set< pair<int,int> > SubProblem::randomForbiddenShifts(int nbForbidden){
	set< pair<int,int> > ans;
	for(int f=0; f<nbForbidden; f++){
		int k = Tools::randomInt(0, nDays_ - 1);
		int s = Tools::randomInt(1, pScenario_->nbShifts_ - 1);
		ans.insert(pair<int,int>(k,s));
	}
	return ans;
}

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
	for(int startDate=1; startDate<nDays_-daysMin_; startDate ++){

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


	  ///////////////////////////////////////

	  ////  SERGEB --   A MODIFIER

	  ///////////////////////////////////////


			// // Initialization, rotations of length CDMin_ (so that all rotations generated are long)
			// //
			//	if(length == 0){
			// 	double bestCostCDMinDays = 0;
			// 	vector<int> bestSuccCDMinDays;
			// 	// for(int newSh=1; newSh<pScenario_->nbShiftsType_; newSh++){
			// 		for(unsigned int i=0; i<allowedShortSuccBySize_[CDMin_].size(); i++){
			// 			vector<int> succ = allowedShortSuccBySize_[CDMin_][i];
			// 			double potentialCost = costOfVeryShortRotation(startDate,succ);
			// 			// Store the rotation if it is good
			// 			if(potentialCost < maxReducedCostBound_){
			// 				vector<int> succToStore = succ;
			// 				Rotation rot (startDate, succToStore, pLiveNurse_->id_, MAX_COST, potentialCost);
			// 				nPaths_ ++;
			// 				theRotations_.push_back(rot);
			// 				nLongFound_ ++;
			// 				nFound ++;
			// 				bestReducedCost_ = min(bestReducedCost_, rot.dualCost_);
			// 			}
			// 			// Chose to extend that one if it is the best
			// 			if(potentialCost < bestCostCDMinDays){
			// 				bestSuccCDMinDays = succ;
			// 				bestCostCDMinDays = potentialCost;
			// 				hasFoundBetter = true;
			// 			}

			// 		}
			// 	// }

			// 	// If a good one has been found -> make it the basis for the next extension
			// 	if(hasFoundBetter){
			// 		bestSucc = bestSuccCDMinDays;
			// 		bestCost = bestCostCDMinDays;
			// 		lastSh = bestSucc[CDMin_-1];
			// 		length = CDMin_;
			// 		currentDate = currentDate + CDMin_;
			// 		nConsShift = 0;
					
			// 		int lastShType = pScenario_->shiftIDToShiftTypeID_[lastSh];
			// 		int bestSuccShType = pScenario_->shiftIDToShiftTypeID_[bestSucc[bestSucc.size()-nConsShift-1]];
					
			// 		while(bestSuccShType == lastShType) nConsShift ++;
			// 	}
			// 	continue;
			//	}

			// OTHER CASES: EXTEND THE ROTATION
			//
			//	else{

				int bestNewShift;
				double bestNewCost = bestCost;
				hasFoundBetter = false;
				
				int lastShType = pScenario_->shiftIDToShiftTypeID_[lastSh];

				for(int newSh=1; newSh<pScenario_->nbShifts_; newSh++){
					double potentialCost = bestCost;
					
					int newShType = pScenario_->shiftIDToShiftTypeID_[newSh];

					// If the succession is allowed -> try to extend the rotation
					if(!pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh) and !isDayShiftForbidden(currentDate+1,newSh)){

						// REGULAR COSTS
						//

						// REG COST: CONSECUTIVE SHIFTS -> if last day + too short, totally cancel that cost !!
						if(currentDate < nDays_-2){
							if(newShType==lastShType){
								potentialCost -= consShiftCost(lastSh,nConsShift);
								potentialCost += consShiftCost(lastSh,nConsShift+1);
							} else{
								potentialCost += consShiftCost(newSh,1);
							}
						}
						// Last day:
						else {
							if(newShType == lastShType){
								potentialCost -= consShiftCost(lastSh,nConsShift);
								if(newShType == lastShType and nConsShift+1 > pScenario_->maxConsShiftsOfTypeOf(newSh))
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
						int level = pLiveNurse_->wishesOffLevel(currentDate+1,newSh);

						if (level != -1)
						  potentialCost += WEIGHT_PREFERENCES_OFF[level];

						level = pLiveNurse_->wishesOnLevel(currentDate+1,newSh);

						if (level != -1)
						  potentialCost += WEIGHT_PREFERENCES_ON[level];


						// DUAL COSTS
						//

						// RED COST: WEEKEND
						if(Tools::isSaturday(currentDate+1)) potentialCost -= pCosts_->workedWeekendCost();

						// RED COST: FIRST DAY -> no need here because already in the beginning

						// RED COST: LAST DAY
						potentialCost += pCosts_->endWorkCost(currentDate);
						potentialCost -= pCosts_->endWorkCost(currentDate+1);

						// RED COST: DAY-SHIFT
						potentialCost -= pCosts_->dayShiftWorkCost(currentDate+1,newSh);


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
							if(pLiveNurse_->wishesOff(currentDate+1,newSh)) cout << WEIGHT_PREFERENCES_OFF;
							if(pLiveNurse_->wishesOn(currentDate+1,newSh)) cout << WEIGHT_PREFERENCES_ON;
							cout << endl;
							cout << "#          -=: ";
							if(Tools::isSaturday(currentDate+1)) cout << pCosts_->workedWeekendCost();
							cout << endl;
							cout << "#          +=: " << pCosts_->endWorkCost(currentDate) << " - " << pCosts_->endWorkCost(currentDate+1) << endl;
							cout << "#          -=: " << pCosts_->dayShiftWorkCost(currentDate+1,newSh) << endl;
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

					int bestNewShiftType = pScenario_->shiftIDToShiftTypeID_[bestNewShift];
					
					nConsShift = (lastShType == bestNewShiftType) ? (nConsShift+1) : 1;
					lastSh = bestNewShift;
				}
				//		}
		}
	}
	return nFound > 0;
}





//--------------------------------------------
//
// PRINT FUNCTIONS
//
//--------------------------------------------

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

// Print the contract type + preferences
void SubProblem::printContractAndPrefenrences(){
	std::cout << "# Preferences:" << endl;
	for(map<int,vector<Wish> >::iterator it = pLiveNurse_->pWishesOff_->begin(); it != pLiveNurse_->pWishesOff_->end(); ++it){
		cout <<  "      | " << it->first << "  ->  ";
		for(Wish s : it->second) cout << pScenario_->intToShift_[s.shift];
		cout << endl;
	}
	for(map<int,vector<Wish> >::iterator it = pLiveNurse_->pWishesOn_->begin(); it != pLiveNurse_->pWishesOn_->end(); ++it){
		cout <<  "      | " << it->first << "  ->  ";
		for(Wish s : it->second) cout << pScenario_->intToShift_[s.shift];
		cout << endl;
	}
	std::cout << "# Contract :   ";
	std::cout << "Days [" << pContract_->minConsDaysOff_ << "<" << pContract_->maxConsDaysOff_ << "]   ";
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
	  std::cout << pScenario_->intToShiftType_[s] << " [" << pScenario_->minConsShiftsOf(s) << "<" << pScenario_->maxConsShiftsOf(s) << "]   ";
	}
	std::cout << std::endl;
	std::cout << "# " << std::endl;
	std::cout << "# " << std::endl;
}


///////////////////////////////////////////


//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
SubProblemShort::SubProblemShort() {}

SubProblemShort::SubProblemShort(Scenario * scenario, int nbDays, const Contract * contract, vector<State>* pInitState):
  SubProblem(scenario, nbDays, contract,  pInitState) {

  daysMin_ = contract->minConsDaysWork_;
  nLabels_ = 1;

  initShortSuccessions();
}


SubProblemShort::~SubProblemShort(){}


// Initialization function
void SubProblemShort::init(vector<State>* pInitState){

  SubProblem::init(pInitState);

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

  nVeryShortFound_=0;
}


// Initializes the short successions. Should only be used ONCE (when creating the SubProblem).
void SubProblemShort::initShortSuccessions(){

	// Primary information needed
	//
        int nShiftsType = pScenario_->nbShiftsType_;
	int nShifts     = pScenario_->nbShifts_;

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
			for(int s=1; s<nShifts; s++){
				vector<int> shiftSuccession; shiftSuccession.push_back(s);		// Create Succession
				allSuccSizeC.push_back(shiftSuccession);	// Add it to the possibilities
				lastShiftSucc.push_back(s);			// Record its last shift
				nLastShiftSucc.push_back(1);			// Only 1 successive performed so far
				arcCostSucc.push_back(0);		       	// No succession ended yet
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
				for(int newSh=1; newSh<nShifts; newSh++){
					if(! pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh)){

						vector<int> newSucc (succ); newSucc.push_back(newSh);		// Create Succession
						allSuccSizeC.push_back(newSucc);							// Add it to the possibilities
						lastShiftSucc.push_back(newSh);								// Record its last shift
						
						int newShTypeID = pScenario_->shiftIDToShiftTypeID_[newSh];
						int lastShTypeID = pScenario_->shiftIDToShiftTypeID_[lastSh];

						if (newShTypeID == lastShTypeID){										// Depending on the previous one, update number of consecutive and cost
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
				for(int newSh=1; newSh<nShifts; newSh++){
					if(! pScenario_->isForbiddenSuccessorShift_Shift(newSh,lastSh)){

						vector<int> newSucc (succ); newSucc.push_back(newSh);		// Create Succession
						allSuccSizeC.push_back(newSucc);							// Add it to the possibilities
						lastShiftSucc.push_back(newSh);								// Record its last shift
						int newNLast = 1;
						double newCost = cost;

						
						int newShTypeID = pScenario_->shiftIDToShiftTypeID_[newSh];
						int lastShTypeID = pScenario_->shiftIDToShiftTypeID_[lastSh];

						
						if(newShTypeID == lastShTypeID){	// BUT : add the cost if longer than the maximum allowed
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
						int n = min(newNLast, maxvalConsByShift_[newShTypeID]);
						allShortSuccCDMinByLastShiftCons_[newShTypeID][n].push_back(allSuccSizeC.size()-1);

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


//--------------------------------------------
//
// Solve function
//
//--------------------------------------------

// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
bool SubProblemShort::solve(LiveNurse* nurse, DualCosts * costs, SubproblemParam param, set<pair<int,int> > forbiddenDayShifts,
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
	forbid(forbiddenDayShifts);				 				// Forbid arcs
	forbidStartingDays(forbiddenStartingDays);				// Forbid starting days

	if(false) printContractAndPrefenrences();				// Set to true if you want to display contract + preferences (for debug)

	// cout << " 1   ----------------------------------------------" << endl;
	// printGraph();

	timeInS_->start();
	bool ANS_short = solveShortRotations();
	timeInS_->stop();

	// cout << " 2   ----------------------------------------------" << endl;
	// printGraph();
	
	timeInNL_->start();
	bool ANS_long = solveLongRotations(optimality);
	timeInNL_->stop();

	// cout << " 3   ----------------------------------------------" << endl;
	// printGraph();

	// printAllRotations();

	// Reset authorizations
	authorize(forbiddenDayShifts);
	authorizeStartingDays(forbiddenStartingDays);

	return ANS_short or ANS_long;
}


// For the short rotations, depends on the chosen option + on wether we want optimality (more important)
bool SubProblemShort::solveShortRotations(){
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


// Initializes some cost vectors that depend on the nurse
//
void SubProblemShort::initStructuresForSolve(){

  SubProblem::initStructuresForSolve();

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

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void SubProblemShort::createArcsSourceToPrincipal(){

  // DATA
  int nShiftsType = pScenario_->nbShiftsType_;
  int origin, destin;

  for(int sh=1; sh<nShiftsType; sh++){
    for(int k=daysMin_-1; k<nDays_; k++){
      for(int nCons=1; nCons<=maxvalConsByShift_[sh]; nCons ++){
        origin = g_.source();
        destin = principalNetworkNodes_[sh][k][nCons];
        // shifts will be filled when updating the graph
        int a = g_.addSingleArc(origin, destin, 0, {daysMin_, CDMin_-daysMin_}, SOURCE_TO_PRINCIPAL, k-daysMin_+1, {});
        arcsFromSource_[sh][k][nCons][0] = a;
      }
    }
  }
}


// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
//
void SubProblemShort::priceShortSucc(){

	map<int,int> specialArcsSuccId;
	map<int,double> specialArcsCost;

	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=daysMin_-1; k<nDays_; k++){
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
							if(k==CDMin_-1 and CDMin_<maxvalConsByShift_[s] and n==CDMin_ and s==pLiveNurse_->pStateIni_->shiftType_){
								// a. Determine the destination of that arc
								int nConsWithPrev = CDMin_ + pLiveNurse_->pStateIni_->consShifts_;
								int nDestination = min( nConsWithPrev , maxvalConsByShift_[s] );
								int a = arcsFromSource_[s][k][nDestination][0];
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
					g_.forbidArc( arcsFromSource_[s][k][n][0] );
				}
			}
		}
	}

	// FOR THE SHIFTS ON THE FIRST DAY THAT EXTEND THE ONGOING WORK AT INITIAL STATE
	//
	for(map<int,int>::iterator itId = specialArcsSuccId.begin(); itId != specialArcsSuccId.end(); ++itId){
		int a = itId->first;
		int d = g_.arcDestination(a);
		int s = principalToShift_[d];
		int k = principalToDay_[d];
		int n = principalToCons_[d];
		if(specialArcsCost.find(a) == specialArcsCost.end()){
			cout << "# Problem within pricing of some short rotations (press Enter to go on)" << endl;
			getchar();
		}
		double cost = specialArcsCost.at(a);
		g_.authorizeArc(a);
		idBestShortSuccCDMin_[s][k][n] = itId->second;
		arcCostBestShortSuccCDMin_[s][k][n] = cost;
	}
}

// Given a short succession and a start date, returns the cost of the corresponding arc
//
double SubProblemShort::costArcShortSucc(int size, int succId, int startDate){
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

		int shiftTypeIni = pScenario_->shiftIDToShiftTypeID_[shiftIni];
		int firstShiftType = pScenario_->shiftIDToShiftTypeID_[firstShift];

		// 1. The nurse was resting: pay more only if the rest is too short
		if(shiftTypeIni == 0){
			int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
			ANS += max(0, diffRest*WEIGHT_CONS_DAYS_OFF);
		}

		// 2. The nurse was working
		else {

			// a. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be read later
			int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
			ANS -= max(0, diffWork*WEIGHT_CONS_DAYS_WORK);

			// b. (i)   The nurse was working on a different shift: if too short, add the corresponding cost
			if(shiftTypeIni != firstShiftType){
			  int diff = pScenario_->minConsShiftsOfTypeOf(shiftIni) - nConsShiftIni;
				ANS += max(0, diff*(WEIGHT_CONS_SHIFTS));
			}

			// b. (ii)  The nurse was working on the same shift AND the short rotation contains other shifts (easy case for add/subtract)
			//            - Subtract the cost due to the consecutive beginning
			//            - Subtract the cost due to the consecutive end of the initial state
			//            - Add the consecutive cost for all shifts
			else if(nConsFirstShift < pContract_->minConsDaysWork_) {
			  //			else if(nConsFirstShift < CDMin_) {
			  int diffShift = nConsShiftIni - pScenario_->maxConsShiftsOfTypeOf(shiftIni);
				ANS -= max(0, diffShift*WEIGHT_CONS_SHIFTS);
				ANS -= consShiftCost(firstShift, nConsFirstShift);
				ANS += consShiftCost(firstShift, (nConsFirstShift + nConsShiftIni));
			}

			// b. (iii) The nurse was working on the same shift AND the short rotation only contains that shift (recompute the cost -> easier)
			else {
				ANS -= baseArcCostOfShortSucc_[size][succId];
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
	for(int i=0; i<size; i++) ANS -= pCosts_->dayShiftWorkCost( startDate+i , allowedShortSuccBySize_[size][succId][i] );



	return ANS;
}

//--------------------------------------------
//
// Functions for the costs
//
//--------------------------------------------

// Updates the costs depending on the reduced costs given for the nurse
//
void SubProblemShort::updateArcCosts(){

	priceShortSucc();

	// A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
	//
	shortSuccCDMinIdFromArc_.clear();
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=daysMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				int a = arcsFromSource_[s][k][n][0];
				double c = arcCostBestShortSuccCDMin_[s][k][n];
				g_.updateCost( a , c );
				int best = idBestShortSuccCDMin_[s][k][n];
				shortSuccCDMinIdFromArc_.insert(pair<int,int>( arcsFromSource_[s][k][n][0], best));
        g_.updateShifts(a, allowedShortSuccBySize_[CDMin_][best]);
			}
		}

		// For all those that start on the first day, must update the travel time
		//
		for(int n=1; n<=maxvalConsByShift_[s]; n++){
			if(idBestShortSuccCDMin_[s][CDMin_-1][n] > -1){
			  int a = arcsFromSource_[s][CDMin_-1][n][0];

                std::vector<int> consumptions = { daysMin_ + pLiveNurse_->pStateIni_->consDaysWorked_,
                                                  CDMin_ - daysMin_ - pLiveNurse_->pStateIni_->consDaysWorked_ };
        g_.updateConsumptions(a, consumptions);
			}
		}
	}

	// B. ARCS : SHIFT_TO_NEWSHIFT [baseCost = 0]
	//
	for(int s1=1; s1<pScenario_->nbShiftsType_; s1++)
		for(int s2=1; s2<pScenario_->nbShiftsType_; s2++)
		  if(s2 != s1 and ! pScenario_->isForbiddenSuccessorShiftType_ShiftType(s2,s1)){
			for(int k=daysMin_-1; k<nDays_-1; k++){
			  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[s2].size(); s++) {
				int a = arcsShiftToNewShift_[s1][s2][k][s];
				if(a > 0){
					double c = g_.arcInitialCost(a);
					int  shiftID = pScenario_-> shiftTypeIDToShiftID_[s2][s];
					c += preferencesCosts_[k+1][shiftID] ;
					c -= pCosts_->dayShiftWorkCost(k+1,shiftID);
					c -= Tools::isSaturday(k+1) ? pCosts_->workedWeekendCost() : 0 ;
          g_.updateCost( a , c );
				}
			  }
			}
		  }

	// C. ARCS : SHIFT_TO_SAMESHIFT [baseCost = 0]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=daysMin_-1; k<nDays_-1; k++)
			for(int n=1; n<maxvalConsByShift_[sh]; n++){
			    for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
				int a = arcsShiftToSameShift_[sh][k][n][s];
				double c = g_.arcInitialCost(a);
				int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
				// c += preferencesCosts_[k+1][sh] ;
				c += preferencesCosts_[k+1][shiftID] ;
				c -= pCosts_->dayShiftWorkCost(k+1,shiftID);
				if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
				g_.updateCost( a , c );
			    }
			}

	// D. ARCS : SHIFT_TO_ENDSEQUENCE [They never change]

	// E. ARCS : REPEATSHIFT [baseCost contains consecutive shift cost]
	//
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)
		for(int k=daysMin_-1; k<nDays_-1; k++){
		  for (unsigned int s = 0; s < pScenario_-> shiftTypeIDToShiftID_[sh].size(); s++) {
			int a = arcsRepeatShift_[sh][k][s];
			double c = g_.arcInitialCost(a);
			int  shiftID = pScenario_-> shiftTypeIDToShiftID_[sh][s];
			// c += preferencesCosts_[k+1][sh];
			c += preferencesCosts_[k+1][shiftID];
			c -= pCosts_->dayShiftWorkCost(k+1,shiftID);
			if(Tools::isSaturday(k+1)) c-= pCosts_->workedWeekendCost();
			g_.updateCost( a , c );
		  }
		}

	// F. ARCS : PRINCIPAL_TO_ROTSIZE [baseCost contains complete weekend constraint]
	//
	for(int s=1; s<pScenario_->nbShiftsType_; s++)
		for(int k=daysMin_-1; k<nDays_; k++){
			int a = arcsPrincipalToRotsizein_[s][k];
			double c = g_.arcInitialCost(a);
			c += endWeekendCosts_[k];
			c -= pCosts_->endWorkCost(k);
			g_.updateCost( a , c );
		}

	// G. ARCS : ROTSIZEIN_TO_ROTSIZE [baseCost contains rotation length cost. They never change]

	// H. ARCS : ROTSIZE_TO_ROTSIZEOUT [They never change]

	// I. ARCS : ROTSIZE_TO_ROTSIZEOUT [Never changes]
}

// Adds a rotation made from the given path to the current list of answers and increases their counter
//
Rotation SubProblemShort::rotationFromPath(vector< boost::graph_traits<Graph>::edge_descriptor > path, spp_res_cont resource){

    int firstDay = -1;
    vector<int> shiftSuccession;

    // All arcs are consecutively considered
    //
    for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
        int a = 0; //boost::get(&Arc_Properties::num, g_, path[j]);
        ArcType aType = g_.arcType(a);
        int destin = 0; //boost::target( path[j], g_ );

        // A. Arc from source (equivalent to short rotation
        if(aType == SOURCE_TO_PRINCIPAL){
            firstDay =  principalToDay_[destin] - CDMin_ + 1;
            for(int s: static_cast<vector <int> >( allowedShortSuccBySize_[CDMin_][ shortSuccCDMinIdFromArc_.at(a) ] )){
                shiftSuccession.push_back(s);
            }
        }

            // B. Arc to a new day
        else if(aType == SHIFT_TO_NEWSHIFT or aType == SHIFT_TO_SAMESHIFT or aType == REPEATSHIFT){
            shiftSuccession.push_back( g_.arcShifts(a).front() );
        }
    }

    Rotation rot (firstDay, shiftSuccession, pLiveNurse_->id_, MAX_COST, resource.cost);
    return rot;
}


int SubProblemShort::computeHoursInRotation(int id){
  int  timeWorked = 0;
			
  for (int i = 0; i < allowedShortSuccBySize_[daysMin_][id].size(); i++) {
    int  shift = allowedShortSuccBySize_[daysMin_][id][i];
    timeWorked += pScenario_->hoursToWork_[shift];
  }

  return timeWorked;
}


int SubProblemShort::computeHoursInRotation(int sh, int k, int n){
  return computeHoursInRotation(idBestShortSuccCDMin_[sh][k][n]);
}

//----------------------------------------------------------------
//
// Cost computation of the "very" short rotations (< CD_min)
//
//----------------------------------------------------------------

// Brutally try all possible short rotations from the very first day
bool SubProblemShort::priceVeryShortRotationsFirstDay(){
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
bool SubProblemShort::priceVeryShortRotationsLastDay(){
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
bool SubProblemShort::priceVeryShortRotations(){
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

// Summary of the short successions generated
void SubProblemShort::printShortSucc(){
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

// Prints all active pairs ( arcFromSource - corresponding short successions )
void SubProblemShort::printShortArcs(){
	for(int s=1; s<pScenario_->nbShiftsType_; s++){
		for(int k=daysMin_-1; k<nDays_; k++){
			for(int n=1; n<=maxvalConsByShift_[s]; n++){
				if(!g_.arcForbidden(arcsFromSource_[s][k][n][0])){
					int v = principalNetworkNodes_[s][k][n];
					int succId = idBestShortSuccCDMin_[s][k][n];
					vector<int> succ = allowedShortSuccBySize_[CDMin_][succId];
					std::cout << "# " << g_.shortNameNode(v) << " <- (id=" << succId << ") ";
					for(unsigned int i=0; i<succ.size(); i++){
						std::cout << pScenario_->intToShift_[succ[i]].at(0);
					}
					std::cout << std::endl;
				}
			}
		}
	}
}


