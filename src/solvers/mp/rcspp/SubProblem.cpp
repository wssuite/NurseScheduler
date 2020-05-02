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


using std::stringstream;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;


//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
SubProblem::SubProblem() {}

SubProblem::SubProblem(PScenario scenario, int nDays, PConstContract contract, vector<State>* pInitState):
					pScenario_(scenario), nDays_(nDays), pContract_ (contract),
					CDMin_(contract->minConsDaysWork_), daysMin_(1), nLabels_(2),
					maxRotationLength_(nDays) {
  if(pInitState) init(*pInitState);
}

SubProblem::~SubProblem(){}

// Initialization function
void SubProblem::init(const vector<State>& initStates){

	// Maximum number of consecutive days worked by a nurse ending at day -1
	//
	maxOngoingDaysWorked_ = 0;
	for(const State& state: initStates)
		maxOngoingDaysWorked_ = std::max( state.consDaysWorked_, maxOngoingDaysWorked_ );

	nFound_=0;
}

void SubProblem::build() {

  g_ = RCGraph(nDays_);

  //	initShortSuccessions();

  createNodes();
  createArcs();

  // Set all status to authorized
  Tools::initVector2D(dayShiftStatus_, nDays_, pScenario_->nbShifts_, true);
  Tools::initVector(startingDayStatus_, nDays_, true);
  nPathsMin_ = 0;

  //std::cout << "# A new subproblem has been created for contract " << contract->name_ << std::endl;

  //printGraph();

  timeInS_ = new Tools::Timer(); timeInS_->init();
  timeInNL_ = new Tools::Timer(); timeInNL_->init();
}



//--------------------------------------------
//
// Solve function
//
//--------------------------------------------

// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
bool SubProblem::solve(PLiveNurse nurse, DualCosts * costs, SubproblemParam param, set<pair<int,int> > forbiddenDayShifts,
		set<int> forbiddenStartingDays, bool optimality, double redCostBound) {


	bestReducedCost_ = 0;
  nFound_=0;
	param_ = param;
  maxReducedCostBound_ = redCostBound - EPSILON;			// Cost bound
  pLiveNurse_ = nurse;									// Store the new nurse
  pCosts_ = costs;										// Store the new cost

  // Maximum rotations length: update the bounds on the nodes if needed
  updatedMaxRotationLengthOnNodes( std::min(nDays_+maxOngoingDaysWorked_,
                                            std::max(pContract_->maxConsDaysWork_, param_.maxRotationLength_)) );

	resetAuthorizations();									// Reset authorizations
	resetSolutions();										// Delete all already existing solutions
	initStructuresForSolve();								// Initialize structures
	forbid(forbiddenDayShifts);								// Forbid arcs
	forbidStartingDays(forbiddenStartingDays);				// Forbid starting days

	if(false) pLiveNurse_->printContractAndPreferences(pScenario_);				// Set to true if you want to display contract + preferences (for debug)

	// cout << " 1   ----------------------------------------------" << endl;
	// printGraph();

	 timeInS_->start();
	 preprocess();
	 timeInS_->stop();

	// cout << " 2   ----------------------------------------------" << endl;
	// printGraph();
	
	timeInNL_->start();
	bool ANS = solveRCGraph(optimality);
	timeInNL_->stop();

	// cout << " 3   ----------------------------------------------" << endl;
//	g_.printGraph();

	// printAllSolutions();

#ifdef DBG
  // check that all solution respects forbidden shifts
  for(const RCSolution& sol: theSolutions_)
    checkForbiddenDaysAndShifts(sol);
#endif

	return ANS;
}

bool SubProblem::preprocess() {
  return true;
}

// For the long rotations, depends whether we want optimality or not
bool SubProblem::solveRCGraph(bool optimality){
  updateArcCosts();
#ifdef DBG
//  g_.printGraph();
#endif
	return solveRCGraphOptimal();		// Solve shortest path problem
}

// Function called when optimal=true in the arguments of solve -> shortest path problem is to be solved
bool SubProblem::solveRCGraphOptimal(){
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> sinks = g_.sinks();

  if(param_.oneSinkNodePerLastDay_)
    sinks.resize(sinks.size()-1); // remove last sink (it's the main one)
  else sinks = {sinks.back()}; // keep just the main one

	std::vector<RCSolution> solutions = g_.solve(nLabels_, maxReducedCostBound_, sinks);

  for(const RCSolution& sol: solutions) {
    theSolutions_.push_back(sol);
    nPaths_ ++;
    nFound_++;
    bestReducedCost_ = std::min(bestReducedCost_, sol.cost);
  }

#ifdef DBG
//  printAllSolutions();
#endif

  return !solutions.empty();
}

// Resets all solutions data (rotations, number of solutions, etc.)
//
void SubProblem::resetSolutions(){
	theSolutions_.clear();
	nPaths_ = 0;
}


//--------------------------------------------
//
// Functions for the NODES of the rcspp
//
//--------------------------------------------

// Function that creates the nodes of the network
void SubProblem::createNodes(){
	// INITIALIZATION
  principalGraphs_.clear();
  priceLabelsGraphs_.clear();

	// 1. SOURCE NODE
	//
	int v = addSingleNode(SOURCE_NODE);
	g_.setSource(v);

	// 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
	//
  principalGraphs_.emplace_back(PrincipalGraph(0,  nullptr)); // just add a dummy rcspp to have the right indices
	for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)		// For each possible worked shift
    principalGraphs_.emplace_back(PrincipalGraph(sh, this));

	// 3. ROTATION LENGTH CHECK
	//
	// For each of the days, do a rotation-length-checker
	// Deactivate the min cost for the last day
	for(int k=0; k<nDays_; k++){
    priceLabelsGraphs_.emplace_back(PriceLabelsGraph(CDMin_, maxRotationLength_, this, k<nDays_-1)); // k<nDays_-1 is false for last day
    // Daily sink node
    g_.addSink(priceLabelsGraphs_.back().exit());
  }

	// 4. SINK NODE
	//
	v = addSingleNode(SINK_NODE);
  g_.addSink(v);
}


//--------------------------------------------
//
// Functions for the ARCS of the rcspp
//
//--------------------------------------------

// Function that creates the arcs of the network
void SubProblem::createArcs(){
  // Initialization
  Tools::initVector4D(arcsFromSource_, pScenario_->nbShiftsType_, nDays_, 0, 0, -1);
  Tools::initVector2D(arcsPrincipalToPriceLabelsIn_, pScenario_->nbShiftsType_, nDays_, -1);

  // create arcs
	createArcsSourceToPrincipal();
	createArcsPrincipalToPrincipal();
  createArcsAllPriceLabels();
}

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void SubProblem::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (PrincipalGraph &pg: principalGraphs_)
    for (int k = daysMin_ - 1; k < nDays_; k++)
      for(int dest : pg.getDayNodes(k)) {
        std::vector<int> vec;
        for (int s: pScenario_->shiftTypeIDToShiftID_[pg.shiftType()])
          vec.emplace_back(addSingleArc(origin, dest, 0, {daysMin_, CDMin_ - daysMin_}, SOURCE_TO_PRINCIPAL, k, {s}));
        arcsFromSource_[pg.shiftType()][k].push_back(vec);
      }
}

// Create all arcs within the principal network
void SubProblem::createArcsPrincipalToPrincipal(){
	// CHANGE SUBNETWORK FOR EACH OF THE SUBNETWORKS AND EACH OF THE DAYS
	//
  int nShiftsType = pScenario_->nbShiftsType_;
  for (int sh=1; sh<nShiftsType; sh++){
		for(int k=0; k<nDays_-1; k++){
			int origin = principalGraphs_[sh].getDayNodes(k).back(); // last level for day k
			for(int newSh=1; newSh<nShiftsType; newSh++){
			  // check if succession is allowed
			  if(newSh != sh and ! pScenario_->isForbiddenSuccessorShiftType_ShiftType(newSh, sh)){
			    int destin = principalGraphs_[newSh].getNode(k+1, 0); // entrance level for day k+1
			    g_.addSingleArc(origin, destin, 0, {0,0}, SHIFT_TO_NEWSHIFT, k+1, {});
			  }
			}
		}
	}
}

// Create all arcs that involve the rotation size checking subnetwork (incoming, internal, and exiting that subnetwork)
void SubProblem::createArcsAllPriceLabels(){
  for(int k=0; k<nDays_; k++){				// For all days
    for(int sh=1; sh<pScenario_->nbShiftsType_; sh++){		// For all shifts
      // incoming  arc
      int origin = principalGraphs_[sh].getDayNodes(k).back();
			int destin = priceLabelsGraphs_[k].entrance();
			arcsPrincipalToPriceLabelsIn_[sh][k] =
			    g_.addSingleArc(origin, destin, 0, {0,0}, PRINCIPAL_TO_ROTSIZE, k);	// Allow to stop rotation that day
		}

    // outgoing  arcs
    int origin = priceLabelsGraphs_[k].exit();
    int destin = g_.lastSink();
    g_.addSingleArc(origin, destin, 0, {0,0}, ROTSIZEOUT_TO_SINK, k);
	}
}

//--------------------------------------------
//
// Functions for the pricing of the columns
//
//--------------------------------------------

// Initializes some cost vectors that depend on the nurse
void SubProblem::initStructuresForSolve(){
	// Start and End weekend costs
	//
	Tools::initVector(startWeekendCosts_,nDays_,.0);
	Tools::initVector(endWeekendCosts_,nDays_,.0);
	if(pLiveNurse_->needCompleteWeekends()){
		for(int k=0; k<nDays_; k++){
			if(Tools::isSaturday(k)) endWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
			else if(Tools::isSunday(k)) startWeekendCosts_[k] = WEIGHT_COMPLETE_WEEKEND;
		}
	}

	// Preference costs.
	//
  Tools::initVector2D(preferencesCosts_, nDays_, pScenario_->nbShifts_, .0);
  for(const auto &p: pLiveNurse_->wishesOff())
		for(const Wish& s : p.second)
      preferencesCosts_[p.first][s.shift] = WEIGHT_PREFERENCES_OFF[s.level];
  for(const auto &p: pLiveNurse_->wishesOn())
    for(const Wish& s: p.second)
      preferencesCosts_[p.first][s.shift] = WEIGHT_PREFERENCES_ON[s.level];
}

double SubProblem::startWorkCost(int a) const {
  // retrieve the work cost
  const Arc_Properties& arc_prop = g_.arc(a);
  double cost = workCost(arc_prop, true);
  cost -= pCosts_->startWorkCost(arc_prop.day);
  cost += startWeekendCosts_[arc_prop.day];

  // if first day, take into account historical state. (shift ID cannot be 0)
  if (arc_prop.day == 0) {
    int shiftIni = pLiveNurse_->pStateIni_->shift_;
    int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
    int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

    // check if the nurse is changing of shift type
    int shiftTypeIni = pScenario_->shiftIDToShiftTypeID_[shiftIni];

    // 1. The nurse was resting: pay more only if the rest is too short
    if(shiftTypeIni == 0) {
      int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      cost += std::max(0, diffRest*WEIGHT_CONS_DAYS_OFF);
    }
    // 2. The nurse was working
    else {
      // a. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be added later
      int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
      cost -= std::max(0, diffWork*WEIGHT_CONS_DAYS_WORK);

      // b. (i)   The nurse was working on a different shift: if too short, add the corresponding cost
      int shiftType = pScenario_->shiftIDToShiftTypeID_[arc_prop.shifts.front()];
      if(shiftTypeIni != shiftType){
        int diff = pScenario_->minConsShiftsOfTypeOf(shiftIni) - nConsShiftIni;
        cost += std::max(0, diff*(WEIGHT_CONS_SHIFTS));
      }
      // b. (ii) If working on the same shift type, need to update the consecutive shift cost
      else {
        int nOldConsShift = 0;
        for(int s: arc_prop.shifts) {
          ++ nConsShiftIni;
          ++ nOldConsShift;
          if(pScenario_->shiftIDToShiftTypeID_[s] != shiftType) {
            // remove old cost and add new one
            cost -= pScenario_->consShiftCost(shiftTypeIni, nOldConsShift);
            cost += pScenario_->consShiftCost(shiftTypeIni, nConsShiftIni);
            break;
          }
        }
      }
    }
  }
  
  return cost;
}


double SubProblem::workCost(int a, bool first_day) const {
  return workCost(g_.arc(a));
}

double SubProblem::workCost(const Arc_Properties& arc_prop, bool first_day) const {
  double cost = arc_prop.initialCost;
  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for(int s: arc_prop.shifts) {
    cost += preferencesCosts_[k][s] ;
    cost -= pCosts_->workedDayShiftCost(k, s);
    // add weekend cost on saturday, except when first day is sunday
    if(Tools::isSaturday(k) || (first_day && Tools::isSunday(k)))
      cost -= pCosts_->workedWeekendCost();
    ++k;
  }
  return cost;
}

double SubProblem::endWorkCost(int a) const {
  const Arc_Properties& arc_prop = g_.arc(a);
  double cost = workCost(arc_prop);
  int length = arc_prop.shifts.size(), end = arc_prop.day + (length>1 ? length-1 : 0);
  cost += endWeekendCosts_[end];
  cost -= pCosts_->endWorkCost(end);
  return cost;
}

//--------------------------------------------
//
// Functions for the costs
//
//--------------------------------------------

// Updates the costs depending on the reduced costs given for the nurse
//
void SubProblem::updateArcCosts() {

  //	A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
  for (PrincipalGraph &pg: principalGraphs_)

    for (int k = daysMin_ - 1; k < nDays_; k++)
      for (int n = 0; n <= pg.maxCons(); ++n)
        for (int a: arcsFromSource_[pg.shiftType()][k][n]) {
          const Arc_Properties &arc_prop = g_.arc(a);
          if (!arc_prop.forbidden && canSuccStartHere(arc_prop) && pg.checkFeasibilityEntranceArc(arc_prop, n)) {
            double c = startWorkCost(a);
            g_.updateCost(a, c);
            // For an arc that starts on the first day, must update the consumption based on the historical state
            if (k == daysMin_ - 1) {
              std::vector<int> consumptions = {daysMin_ + pLiveNurse_->pStateIni_->consDaysWorked_,
                                               CDMin_ - daysMin_ - pLiveNurse_->pStateIni_->consDaysWorked_};
              g_.updateConsumptions(a, consumptions);
            }
          } else g_.forbidArc(a);
        }

  // B. ARCS : PRINCIPAL GRAPH
  //
  for (PrincipalGraph &pg: principalGraphs_)
    pg.updateArcCosts();

  // C. ARCS : PRINCIPAL_TO_ROTSIZE
  //
  for (int s = 1; s < pScenario_->nbShiftsType_; s++)
    for (int k = daysMin_ - 1; k < nDays_; k++) {
      int a = arcsPrincipalToPriceLabelsIn_[s][k];
      double c = endWorkCost(a);
      g_.updateCost(a, c);
    }

  //D. ARCS : PRICE LABELS
  for (PriceLabelsGraph &plg: priceLabelsGraphs_)
    plg.updateArcCosts(); // nothing for the moment
}


//--------------------------------------------
//
// Functions to update the maximum length of a rotation
//
//--------------------------------------------

// Updates the maximum arrival time at all nodes so that the maximum length of a rotation is updated.
//   + Updates only
void SubProblem::updatedMaxRotationLengthOnNodes(int maxRotationLength){
  if(maxRotationLength != maxRotationLength_){
    maxRotationLength_ = maxRotationLength;
    for(int v=0; v<g_.nodesSize(); v++) {
      if (g_.nodeType(v) != ROTATION_LENGTH) {
        std::vector<int> ubs = g_.nodeUBs(v);
        ubs[MAX_CONS_DAYS] = maxRotationLength;
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
bool SubProblem::canSuccStartHere(int a) const{
	return canSuccStartHere(g_.arc(a));
}

bool SubProblem::canSuccStartHere(const Arc_Properties& arc_prop) const{
  return canSuccStartHere(arc_prop.day, arc_prop.shifts);
}

bool SubProblem::canSuccStartHere(int k, const std::vector<int>& shifts) const{
  // If the starting date is forbidden, return false
  if(!(startingDayStatus_[k]))
    return false;
  // If the succession with the previous shift (day -1) is not allowed
  if(k==0 and
     pScenario_->isForbiddenSuccessorShift_Shift(shifts.front(),pLiveNurse_->pStateIni_->shift_))
    return false;
  // If some day-shift is forbidden
  for(int s: shifts)
    if( ! dayShiftStatus_[k++][s])
      return false;
  return true;
}

// Forbids the nodes that correspond to forbidden shifts
//
void SubProblem::forbid(const set<pair<int,int> >& forbiddenDayShifts){
	for(const pair<int,int>& p : forbiddenDayShifts){
		//std::cout << "# Trying to forbid " << p.first << "-" << pScenario_->intToShift_[p.second].at(0) << endl;
		forbidDayShift(p.first,p.second);
	}
}

// Authorize the nodes that correspond to forbidden shifts
//
void SubProblem::authorize(const set<pair<int,int> >& forbiddenDayShifts){
	for(const pair<int,int> &p : forbiddenDayShifts){
		//std::cout << "# Trying to forbid " << p.first << "-" << pScenario_->intToShift_[p.second].at(0) << endl;
		authorizeDayShift(p.first,p.second);
	}
}

// Forbids a day-shift couple : Forbid the nodes + mark (day,shift) as forbidden for the short rotation pricer
//
void SubProblem::forbidDayShift(int k, int s){
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = false;

	int sh = pScenario_->shiftIDToShiftTypeID_[s];
  principalGraphs_[sh].forbidDayShift(k,s);
}

// (re)Authorizes the day-shift couple BUT does not take it into account in the short rotation pricer (too complicated, will be called in the next solve() anyway)
void SubProblem::authorizeDayShift(int k, int s){
	// Mark the day-shift as forbidden
	dayShiftStatus_[k][s] = true;

  int sh = pScenario_->shiftIDToShiftTypeID_[s];
  principalGraphs_[sh].authorizeDayShift(k,s);
}

// Forbids some starting days
//
void SubProblem::forbidStartingDays(const set<int>& forbiddenStartingDays){
	for(int k: forbiddenStartingDays)
		forbidStartingDay(k);
}

// Authorizes some starting days
//
void SubProblem::authorizeStartingDays(const set<int>& authorizedStartingDays){
	for(int k: authorizedStartingDays)
		authorizeStartingDay(k);
}

// Forbids a starting date: no rotation can now start on that day. Gives a prohibitive resource consumption on all short
// rotation arcs that correspond to rotations starting on that day.
//
void SubProblem::forbidStartingDay(int k){
	// Mark the starting day as forbidden
	startingDayStatus_[k] = false;
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
	// set all value to true
  Tools::initVector2D(dayShiftStatus_, nDays_, pScenario_->nbShifts_, true);
  Tools::initVector(startingDayStatus_, nDays_, true);
  // reset authorizations for all arcs and nodes
	g_.resetAuthorizations();
}

// Generate random forbidden shifts
set< pair<int,int> > SubProblem::randomForbiddenShifts(int nbForbidden){
	set< pair<int,int> > ans;
	for(int f=0; f<nbForbidden; f++){
		int k = Tools::randomInt(0, nDays_ - 1);
		int s = Tools::randomInt(1, pScenario_->nbShifts_ - 1);
		ans.insert(std::pair<int,int>(k,s));
	}
	return ans;
}

//--------------------------------------------
//
// PRINT FUNCTIONS
//
//--------------------------------------------

// Prints all rotations in the current list
void SubProblem::printAllSolutions() const {
	std::cout << "# HERE ARE ALL " << nPaths_ << " ROTATIONS OF THE CURRENT SOLUTION LIST :" << std::endl;
	for(const RCSolution& sol : theSolutions_){
		std::cout << sol.toString(pScenario_->shiftIDToShiftTypeID_);
	}
	std::cout << "# " << std::endl;
}

// Print the list of currently forbidden day and shifts
void SubProblem::printForbiddenDayShift() const {
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

void SubProblem::checkForbiddenDaysAndShifts(const RCSolution& sol) const {
  if (isStartingDayforbidden(sol.firstDay))
    Tools::throwError("A RC solution starts on a forbidden day " + std::to_string(sol.firstDay) + ": "
                      + sol.toString(pScenario_->shiftIDToShiftTypeID_));
  int k = sol.firstDay;
  for (int s: sol.shifts)
    if (isDayShiftForbidden(k++, s))
      Tools::throwError("A RC solution uses the forbidden day " + std::to_string(k - 1)
                        + " and shift " + std::to_string(s) + ": "
                        + sol.toString(pScenario_->shiftIDToShiftTypeID_));
}