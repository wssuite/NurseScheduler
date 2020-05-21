/*
 * SubProblem.cpp
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#include "SubProblem.h"
#include "rcspp/BoostRCSPP.h"

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

const int SubproblemParam::maxSubproblemStrategyLevel_ = 3;

void SubproblemParam::initSubprobemParam(int strategy, PLiveNurse pNurse, MasterProblem* pMaster) {
  epsilon = pMaster->getModel()->epsilon();
  maxRotationLength_ = pNurse->maxConsDaysWork();
  int nb_max_path = (int) round(pMaster->getModel()->getParameters().sp_columns_ratio_for_number_paths_ *
                                pMaster->getModel()->getParameters().nbMaxColumnsToAdd_);
  switch (strategy) {

    // 0 -> [Heuristic large search]
    //		short = all,
    //    min   = 0,
    //		max   = CD_max+3
    //
    case 0:
      search_strategy_ = SP_BEST_FIRST;
      nb_max_paths_ = nb_max_path;
      violateConsecutiveConstraints_ = true;
      shortRotationsStrategy_ = 3;
      maxRotationLength_ += 3;
      break;

    // 1 -> [Exact legal only]
    //		short = first and last day,
    //    min   = CD_min,
    //		max   = CD_max
    //
    case 1:
      search_strategy_ = SP_BREADTH_FIRST;
      nb_max_paths_ = -1;
      violateConsecutiveConstraints_ = false;
      shortRotationsStrategy_ = 2;
      maxRotationLength_ += 0;
      break;

      // 2 -> [Exact above legal]
      //		short = all,
      //    min   = 0,
      //		max   = CD_max+2
      //
    case 2:
      search_strategy_ = SP_BREADTH_FIRST;
      nb_max_paths_ = -1;
      violateConsecutiveConstraints_ = true;
      shortRotationsStrategy_ = 3;
      maxRotationLength_ += 2;
      break;

      // 3 -> [Exact exhaustive search]
      //		short = all,
      //    min   = 0,
      // 		max   = LARGE
      //
    case 3:
      search_strategy_ = SP_BREADTH_FIRST;
      nb_max_paths_ = -1;
      violateConsecutiveConstraints_ = true;
      shortRotationsStrategy_ = 3;
      maxRotationLength_ = pNurse->pStateIni_->consDaysWorked_ + pNurse->nbDays_;
      break;


      // UNKNOWN STRATEGY
    default:
      std::cout << "# Unknown strategy for the subproblem (" << strategy << ")" << std::endl;
      break;
  }
}

// Constructors and destructor
SubProblem::SubProblem():
    pScenario_(nullptr), nDays_(0), pContract_(nullptr), pLiveNurse_(nullptr), pCosts_(nullptr) {}

SubProblem::SubProblem(PScenario scenario, int nDays, PConstContract contract, vector<State>* pInitState):
					pScenario_(scenario), nDays_(nDays), pContract_ (contract),
					CDMin_(contract->minConsDaysWork_), minConsDays_(1), nLabels_(2),
					maxRotationLength_(nDays) {
  int max = 0;
  for(int t: pScenario_->timeDurationToWork_)
    if(t>max) max = t;
  maxTotalDuration_ = max * nDays; // working everyday on the longest shift
  if(pInitState) init(*pInitState);
}

SubProblem::~SubProblem(){
  delete timeInS_;
  delete timeInNL_;
}

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
		set<int> forbiddenStartingDays, double redCostBound) {


	bestReducedCost_ = 0;
  nFound_=0;
	param_ = param;
  maxReducedCostBound_ = - redCostBound - param.epsilon;			// Cost bound
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
	if(!param.violateConsecutiveConstraints_)
    forbidViolationConsecutiveConstraints();

	if(false) pLiveNurse_->printContractAndPreferences(pScenario_);				// Set to true if you want to display contract + preferences (for debug)

	 timeInS_->start();
	 preprocess();
	 timeInS_->stop();
	
	timeInNL_->start();
	bool ANS = solveRCGraph();
	timeInNL_->stop();

  // and check that all solution respects forbidden shifts
#ifdef DBG
  for(RCSolution& sol: theSolutions_) {
    checkForbiddenDaysAndShifts(sol);
  }
#endif

	return ANS;
}

bool SubProblem::preprocess() {
  return true;
}

// return a function that will post process any label found by the BoostRCSPPSolver
std::function<void (spp_res_cont&)> postProcessResCont(const DualCosts &costs) {
  double constant = costs.constant();
  return [constant](spp_res_cont& res_cont) {
      res_cont.cost -= constant;
  };
}
RCSPPSolver* SubProblem::initRCSSPSolver() {
  return new BoostRCSPPSolver(&g_, maxReducedCostBound_,
      param_.epsilon, param_.search_strategy_, param_.nb_max_paths_, nullptr);
}

// shortest path problem is to be solved
bool SubProblem::solveRCGraph(){
  // update arcs costs
  updateArcCosts();

#ifdef DBG
//  g_.printGraph(nLabels_, minConsDays_);
#endif

  // solve the RC SPP
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> sinks = g_.sinks();

  if(param_.oneSinkNodePerLastDay_ && sinks.size()>1) {
    sinks.resize(sinks.size() - 1); // remove last sink (it's the main one)
    for(int a: arcsTosink_) g_.forbidArc(a);
  }
  else sinks = {sinks.back()}; // keep just the main one

  std::vector<int> labelsMinLevel = {
      pLiveNurse_->maxConsDaysWork(),
      pLiveNurse_->minConsDaysWork(), // cannot be reach (UB)
      pLiveNurse_->maxTotalShifts(),
      pLiveNurse_->minTotalShifts(), // cannot be reach (UB)
      pLiveNurse_->maxTotalWeekends()
  };

  RCSPPSolver* solver = initRCSSPSolver();
	std::vector<RCSolution> solutions = g_.solve(solver, nLabels_, labelsMinLevel, sinks);
	delete solver;

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
// Functions for the ARCS of the rcspp
//
//--------------------------------------------

// Function that creates the arcs of the network
void SubProblem::createArcs(){
  // Initialization
  Tools::initVector4D(arcsFromSource_, pScenario_->nbShiftsType_, nDays_, 0, 0, -1);
  Tools::initVector3D(principalToPrincipal_, pScenario_->nbShiftsType_, pScenario_->nbShiftsType_, nDays_, -1);
  Tools::initVector2D(arcsPrincipalToPriceLabelsIn_, pScenario_->nbShiftsType_, nDays_, -1);

  // create arcs
	createArcsSourceToPrincipal();
	createArcsPrincipalToPrincipal();
  createArcsAllPriceLabels();
}

// Create all arcs within the principal network
void SubProblem::createArcsPrincipalToPrincipal() {
  // CHANGE SUBNETWORK FOR EACH OF THE SUBNETWORKS AND EACH OF THE DAYS
  //
  int nShiftsType = pScenario_->nbShiftsType_;
  for (int sh = 0; sh < nShiftsType; sh++)
    for (int newSh = 0; newSh < nShiftsType; newSh++)
      // check if succession is allowed
      if (newSh != sh and !pScenario_->isForbiddenSuccessorShiftType_ShiftType(newSh, sh))
        for (int k = 0; k < nDays_ - 1; k++) {
          int origin = principalGraphs_[sh].exit(k); // last level for day k
          if (origin == -1) continue; // if undefined, continue
          // entrance level for day k
          int destin = principalGraphs_[newSh].entrance(k);
          if (destin == -1) continue; // if undefined, continue
          // if start work on sunday (rest to start shift)
          bool weekend = (sh==0 && Tools::isSunday(k+1));
          principalToPrincipal_[sh][newSh][k] =
              g_.addSingleArc(origin, destin, 0, {0, 0, 0, 0, weekend}, SHIFT_TO_NEWSHIFT, k);
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
			if(Tools::isSaturday(k)) endWeekendCosts_[k] = pScenario_->weights().WEIGHT_COMPLETE_WEEKEND;
			else if(Tools::isSunday(k)) startWeekendCosts_[k] =  pScenario_->weights().WEIGHT_COMPLETE_WEEKEND;
		}
	}

	// Preference costs.
	//
  Tools::initVector2D(preferencesCosts_, nDays_, pScenario_->nbShifts_, .0);
  for(const auto &p: pLiveNurse_->wishesOff())
		for(const Wish& s : p.second)
      preferencesCosts_[p.first][s.shift] = pScenario_->weights().WEIGHT_PREFERENCES_OFF[s.level];
  for(const auto &p: pLiveNurse_->wishesOn())
    for(const Wish& s: p.second)
      preferencesCosts_[p.first][s.shift] = pScenario_->weights().WEIGHT_PREFERENCES_ON[s.level];
}

double SubProblem::startWorkCost(int a) const {
  // retrieve the work cost
  const Arc_Properties &arc_prop = g_.arc(a);
  double cost = shiftCost(arc_prop, true);
  int start = arc_prop.day;
  if (arc_prop.shifts.empty())
    ++start; // start work the next day as no shift today
  cost -= pCosts_->startWorkCost(start);
  cost += startWeekendCosts_[start];

  return cost;
}


double SubProblem::shiftCost(int a, bool first_day) const {
  return shiftCost(g_.arc(a), first_day);
}

double SubProblem::shiftCost(const Arc_Properties &arc_prop, bool first_day) const {
  double cost = arc_prop.initialCost;
  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for(int s: arc_prop.shifts) {
    cost += preferencesCosts_[k][s];
    if(pScenario_->isWorkShift(s)) {
      cost -= pCosts_->workedDayShiftCost(k, s);
      // add weekend cost on saturday, except when first day is sunday
      if (Tools::isSaturday(k) || (first_day && Tools::isSunday(k)))
        cost -= pCosts_->workedWeekendCost();
    }
    ++k;
  }
  return cost;
}

double SubProblem::endWorkCost(int a) const {
  const Arc_Properties &arc_prop = g_.arc(a);
  double cost = shiftCost(arc_prop);
  int length = arc_prop.shifts.size(), end = arc_prop.day;
  if (length > 1) end += length - 1; // compute the end of the sequence of shifts
  if(arc_prop.shifts.empty() || arc_prop.shifts.back())
  cost += endWeekendCosts_[end];
  cost -= pCosts_->endWorkCost(end);

  return cost;
}

// take into account historical state depending on current shift (should be called wisely)
double SubProblem::historicalCost(int a) const {
  const Arc_Properties &arc_prop = g_.arc(a);

  // WARNING: the following logic is based on the fact that shifts contains only one element
  if (arc_prop.shifts.size() != 1)
    Tools::throwError("The initial state handling in startWorkCost "
                      "is implemented only for a sequence of one shift.");

  double cost = 0;
  int currentShift = arc_prop.shifts.front();
  int shiftTypeIni = pLiveNurse_->pStateIni_->shiftType_;
  int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
  int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

  // if resting
  if (pScenario_->isRestShift(currentShift)) {
    // 1. The nurse was resting
    if (shiftTypeIni == 0) {
      // if the nurse has already exceeded its max amount of rest,
      // add one penalty as current shift is already over the max
      if(pLiveNurse_->pStateIni_->consDaysOff_ >= pLiveNurse_->maxConsDaysOff())
        cost += pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
    }
    // 2. The nurse was working
    else {
      // pay just penalty for min
      int diff = pLiveNurse_->minConsDaysWork() - nConsWorkIni;
      cost += std::max(.0, diff * pScenario_->weights().WEIGHT_CONS_DAYS_WORK);

      int diff2 = pScenario_->minConsShiftsOf(shiftTypeIni) - nConsShiftIni;
      cost += std::max(.0, diff2 * pScenario_->weights().WEIGHT_CONS_SHIFTS);
    }
  }
  // otherwise, currently working
  else {
    // 1. The nurse was resting: pay more only if the rest is too short
    if (shiftTypeIni == 0) {
      int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      cost += std::max(.0, diffRest * pScenario_->weights().WEIGHT_CONS_DAYS_OFF);
    }
    // 2. The nurse was working
    else {
      // a. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be added later
      int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
      cost -= std::max(.0, diffWork * pScenario_->weights().WEIGHT_CONS_DAYS_WORK);

      // b.   The nurse was working on a different shift: if too short, add the corresponding cost
      int shiftType = pScenario_->shiftIDToShiftTypeID_[currentShift];
      if (shiftTypeIni != shiftType) {
        int diff = pScenario_->minConsShiftsOf(shiftTypeIni) - nConsShiftIni;
        cost += std::max(.0, diff * pScenario_->weights().WEIGHT_CONS_SHIFTS);
      }
      // c. If working on the same shift type, need to update the consecutive shift cost just if exceeding the max
      else if(nConsShiftIni >= pScenario_->maxConsShiftsOf(shiftTypeIni))
        cost += pScenario_->weights().WEIGHT_CONS_SHIFTS;
    }
  }

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
    for (unsigned int k = minConsDays_ - 1; k < arcsFromSource_[pg.shiftType()].size(); k++)
      for (int n = 0; n <= pg.maxCons(); ++n)
        for (int a: arcsFromSource_[pg.shiftType()][k][n]) {
          const Arc_Properties &arc_prop = g_.arc(a);
          if (!arc_prop.forbidden && canSuccStartHere(arc_prop) &&
              pg.checkFeasibilityEntranceArc(arc_prop, n)) {
            double c = 0;
            // if first day -> add the historical costs
            if(k==minConsDays_-1)
              c += historicalCost(a);
            // if rest shift, just add shift cost
            if(pg.shiftType()==0)
              c += shiftCost(a);
            // otherwise, call startWorkCost method
            else
              c += startWorkCost(a);
            g_.updateCost(a, c);
            // For an arc that starts on the first day, must update the consumption based on the historical state
            if (k == minConsDays_ - 1)
              g_.updateConsumptions(a, startConsumption(k, arc_prop.shifts));
          } else g_.forbidArc(a);
        }

  // B. ARCS : PRINCIPAL GRAPH
  //
  for (PrincipalGraph &pg: principalGraphs_)
    pg.updateArcCosts();

  // C. ARCS : PRINCIPAL GRAPH TO PRINCIPAL GRAPH
  // usefull if rest principal graph active
  for (int s = 1; s < pScenario_->nbShiftsType_; s++) {
    for (int a: principalToPrincipal_[0][s])
      if(a != -1) g_.updateCost(a, startWorkCost(a));
    for (int a: principalToPrincipal_[s][0])
      if(a != -1) g_.updateCost(a, endWorkCost(a));
  }

  // D. ARCS : PRINCIPAL_TO_PRICE_LABEL
  // starts at 1, as on 0 we rest (so no end work)
  for (int s = 1; s < pScenario_->nbShiftsType_; s++)
    for (int a: arcsPrincipalToPriceLabelsIn_[s])
      g_.updateCost(a, endWorkCost(a));

  // E. ARCS : PRICE LABELS
  for (std::vector<PriceLabelGraph> &graphs: priceLabelsGraphs_)
    for (PriceLabelGraph& plg: graphs)
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
      if (g_.nodeType(v) != PRICE_LABEL) {
        std::vector<int> ubs = g_.nodeUBs(v);
        ubs[MAX_CONS_DAYS] = maxRotationLength;
        g_.updateUBs(v, ubs);
      }
    }
  }
}

std::vector<int> SubProblem::startConsumption(int day, std::vector<int> shifts) const {
  if(pScenario_->isRestShift(shifts.back())) return {
        0,
        CDMin_,
        0,
        pContract_->minTotalShifts_,
        0
    };

  int timeDuration = 0, size = 0;
  for(int s: shifts) {
    if(pScenario_->isRestShift(s)) {
      timeDuration = 0;
      size = 0;
    } else {
      timeDuration += pScenario_->timeDurationToWork_[s];
      ++ size;
    }
  }
  std::vector<int> c = {
      size,
      CDMin_ - size,
      timeDuration,
      pContract_->minTotalShifts_ - timeDuration,
      Tools::containsWeekend(day-size+1, day)
  };
  // if need to take the historical state
  if(day == size - 1 && pLiveNurse_) {
    c[MAX_CONS_DAYS] += pLiveNurse_->pStateIni_->consDaysWorked_;
    c[MIN_CONS_DAYS] -= pLiveNurse_->pStateIni_->consDaysWorked_;
  }
  return  c;
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

bool SubProblem::canSuccStartHere(int k, const std::vector<int>& shifts) const {
  // If the starting date is forbidden, return false
  if (!startingDayStatus_[k])
    return false;
  // If the succession with the previous shift (day -1) is not allowed
  if (k == 0 and
      pScenario_->isForbiddenSuccessorShift_Shift(shifts.front(), pLiveNurse_->pStateIni_->shift_))
    return false;
  // If some day-shift is forbidden
  for (int s: shifts)
    if (!dayShiftStatus_[k++][s])
      return false;
  return true;
}

// Forbids the nodes that correspond to forbidden shifts
//
void SubProblem::forbid(const set<pair<int,int> >& forbiddenDayShifts){
	for(const pair<int,int>& p : forbiddenDayShifts){
//		std::cout << "# Forbid " << p.first << "-" << pScenario_->intToShift_[p.second] << std::endl;
		forbidDayShift(p.first,p.second);
	}
}

// Authorize the nodes that correspond to forbidden shifts
//
void SubProblem::authorize(const set<pair<int,int> >& forbiddenDayShifts){
	for(const pair<int,int> &p : forbiddenDayShifts){
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

// forbid any arc that authorizes the violation of a consecutive constraint
void SubProblem::forbidViolationConsecutiveConstraints() {
  for(PrincipalGraph& pg: principalGraphs_)
    pg.forbidViolationConsecutiveConstraints();
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

// Reset all authorizations to true--sp-type LONG --dir datasets/ --instance n005w4 --weeks 1-2-3-3 --his 0
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