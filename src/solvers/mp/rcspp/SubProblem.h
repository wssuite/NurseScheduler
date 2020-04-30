/*
 * SubProblem.h
 *
 *  Created on: 30 janv. 2015
 *      Author: samuel
 */

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include "solvers/mp/rcspp/RCGraph.h"
#include "solvers/mp/rcspp/PriceLabelsGraph.h"
#include "solvers/mp/rcspp/PrincipalGraph.h"
#include "solvers/mp/MasterProblem.h"

static int MAX_COST = 99999;

// Parameters (called in the solve function)
//
struct SubproblemParam{

	SubproblemParam(){}
	SubproblemParam(int strategy, LiveNurse* pNurse){
		initSubprobemParam(strategy, pNurse);
	}
	~SubproblemParam(){};

	void initSubprobemParam(int strategy, LiveNurse*   pNurse){
		maxRotationLength_ = pNurse->maxConsDaysWork();
		switch(strategy){

		// 0 -> [Legal only]
		//		short = day-0 and last-day,
		//		max   = CD_max
		//		sink  = one / last day
		//
		case 0: shortRotationsStrategy_=2; maxRotationLength_+=1; oneSinkNodePerLastDay_ = true; break;

		// 3 -> [Exhaustive search]
		//		short = true,
		// 		max   = LARGE
		//		sink  = one / last day
		//
		case 1:	shortRotationsStrategy_=3;	maxRotationLength_=pNurse->pStateIni_->consDaysWorked_+pNurse->nbDays_; oneSinkNodePerLastDay_ = true; break;


		// UNKNOWN STRATEGY
		default:
			std::cout << "# Unknown strategy for the subproblem (" << strategy << ")" << std::endl;
			break;
		}
	}

	// *** PARAMETERS ***
  static const int maxSubproblemStrategyLevel_ = 1;

	// 0 -> no short rotations
	// 1 -> day-0 short rotations only
	// 2 -> day-0 and last-day short rotations only
	// 3 -> all short rotations
	int shortRotationsStrategy_ = 0;

	// maximal length for a rotation
	int maxRotationLength_ = 0;

	// true  -> one sink node per day
	// false -> one single sink node for the network
	bool oneSinkNodePerLastDay_ = false;

	// Getters for the class fields
	//
	int maxRotationLength(){ return maxRotationLength_; }
	int shortRotationsStrategy(){return shortRotationsStrategy_;}
	bool oneSinkNodePerLastDay(){return oneSinkNodePerLastDay_;}

	// Setters
	//
	void maxRotationLength(int value){maxRotationLength_ = value;}
	void shortRotationsStrategy(int value){shortRotationsStrategy_ = value;}
	void oneSinkNodePerLastDay(bool value){oneSinkNodePerLastDay_ = value;}


};


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

    virtual ~SubProblem();

    // Constructor that correctly sets the resource (time + bounds), but NOT THE COST
    //
    SubProblem(Scenario* scenario, int nbDays, PConstContract contract, std::vector<State> *pInitState);

    // Initialization function for all global variables (not those of the rcspp)
    //
    virtual void init(const std::vector<State> &pInitState);

    // Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
    //
    virtual bool solve(LiveNurse*  nurse,
                       DualCosts *costs,
                       SubproblemParam param,
                       std::set<std::pair<int, int> > forbiddenDayShifts = {},
                       std::set<int> forbiddenStartingDays = {},
                       bool optimality = false,
                       double redCostBound = -1e-9);

    // Returns all rotations saved during the process of solving the SPPRC
    //
    inline const std::vector<RCSolution>& getSolutions() const { return theSolutions_; }

    virtual void build();

    // Some getters
    //
    inline Scenario* scenario() const { return pScenario_; }

    inline PConstContract contract() const { return pContract_; }

    inline const LiveNurse*  liveNurse() const { return pLiveNurse_; }

    inline int nDays() const { return nDays_; }

    inline RCGraph &g() { return g_; }

    inline int nPaths() const { return nPaths_; }

    inline int nFound() const { return nFound_; }

    // Returns true if the corresponding shift has no maximum limit of consecutive worked days
    //
    inline bool isUnlimited(int shift_type) const {
      return pScenario_->maxConsShiftsOf(shift_type) >= nDays_ + maxOngoingDaysWorked_
             or pScenario_->maxConsShiftsOf(shift_type) >= NB_SHIFT_UNLIMITED;
    }

    inline int maxCons(int shift_type) const {
      if (isUnlimited(shift_type))
        return pScenario_->minConsShiftsOf(shift_type);
      return pScenario_->maxConsShiftsOf(shift_type);
    }

    inline int addSingleNode(NodeType type, std::vector<int> lbs = {0, 0}, std::vector<int> ubs = {}) {
      if (ubs.empty())
        ubs = {maxRotationLength_, CDMin_};
      return g_.addSingleNode(type, lbs, ubs);
    }

    inline int addSingleArc(int origin, int destination, double baseCost, std::vector<int> consumptions,
                            ArcType type, int day = -1, std::vector<int> shifts = {}) {
      return g_.addSingleArc(origin, destination, baseCost, consumptions, type, day, shifts);
    }

    virtual double startWorkCost(int a) const;

    virtual double workCost(int a, bool first_day=false) const;

    virtual double endWorkCost(int a) const;

    // Print and check functions.
    //
    void printAllSolutions() const;

    void printForbiddenDayShift() const;

    void checkForbiddenDaysAndShifts(const RCSolution& sol) const;

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
    Scenario* pScenario_;

    // Number of days of the scenario (usually a multiple of 7)
    //
    int nDays_;

    // Contract type
    //
    PConstContract pContract_;

    // (Minimum) number of paths to return to the MP
    //
    int nPathsMin_;

    // Current live nurse considered
    //
    LiveNurse*  pLiveNurse_;

    // All costs from Master Problem
    //
    DualCosts *pCosts_;

    // Bound on the reduced cost: if greater than this, the rotation is not added
    //
    double maxReducedCostBound_;

    // Maximum number of consecutive days already worked by a nurse before the beginning of that period
    //
    int maxOngoingDaysWorked_;

    //----------------------------------------------------------------
    //
    // Answers: rotations, number of paths found
    //
    //----------------------------------------------------------------

    // Saved Rotations
    //
    std::vector<RCSolution> theSolutions_;

    // Number of paths found
    //
    int nPaths_;

    // Number of rotations found (that match the bound condition) at that iteration
    //
    int nFound_;

    // Best reduced cost found
    //
    double bestReducedCost_;

    //----------------------------------------------------------------
    //
    // Solving options.
    //
    //----------------------------------------------------------------

    SubproblemParam param_;

    //-----------------------
    // THE BASE COSTS
    //-----------------------
    // WARNING : for those that never change, of no use also.

    // For each day k (<= nDays_ - CDMin), contains WEIGHT_COMPLETE_WEEKEND if [it is a Saturday (resp. Sunday) AND the contract requires complete weekends]; 0 otherwise.
    std::vector<double> startWeekendCosts_, endWeekendCosts_;
    // Costs due to preferences of the nurse: for each day k (<= nDays_ - CDMin), shift s, contains WEIGHT_PREFERENCES if (k,s) is a preference of the nurse; 0 otherwise.
    vector2D<double> preferencesCosts_;

    int CDMin_;  // Minimum number of consecutive days worked for free
    int daysMin_;      // principal node network begins at this index-1;  1 if no ShortSucc, CDMin otherwis
    int nLabels_; // Number of labels to use
    int maxRotationLength_;  // MAXIMUM LENGTH OF A ROTATION (in consecutive worked days)


    //-----------------------
    // THE GRAPH
    //-----------------------
    RCGraph g_;

    // Data structures for the nodes and arcs
    std::vector<PrincipalGraph> principalGraphs_;
    std::vector<PriceLabelsGraph> priceLabelsGraphs_;
    vector4D<int> arcsFromSource_;    // Index: (shiftType, day, n, shift) of destination
    vector2D<int> arcsPrincipalToPriceLabelsIn_;  // Index: (shiftType, day) of origin

    // Creates all nodes of the rcspp (including resource window)
    void createNodes();

    // Creates all arcs of the rcspp
    void createArcs();

    // Create the specific types of arcs
    virtual void createArcsSourceToPrincipal();

    virtual void createArcsPrincipalToPrincipal();

    virtual void createArcsAllPriceLabels();

    // Updates the costs depending on the reduced costs given for the nurse
    virtual void updateArcCosts();

    double workCost(const Arc_Properties &a, bool first_day=false) const;

    // FUNCTIONS -- SOLVE
    //
    virtual bool preprocess();

    virtual bool solveRCGraph(bool optimality);

    // Function called when optimal=true in the arguments of solve
    virtual bool solveRCGraphOptimal();

    // Initializes some cost vectors that depend on the nurse
    virtual void initStructuresForSolve();

    // Resets all solutions data (rotations, number of solutions, etc.)
    void resetSolutions();

    void updatedMaxRotationLengthOnNodes(int maxRotationLentgh);

    // FORBIDDEN ARCS AND NODES
    //
    vector2D<bool> dayShiftStatus_;
    std::vector<bool> startingDayStatus_;

    // Returns true if the succession succ starting on day k does not violate any forbidden day-shift
    bool canSuccStartHere(int a) const;
    bool canSuccStartHere(const Arc_Properties& arc_prop) const;
    bool canSuccStartHere(int k, const std::vector<int>& shifts) const;

    // Forbids some days / shifts
    void forbid(const std::set<std::pair<int, int> >& forbiddenDayShifts);

    // Authorizes some days / shifts
    void authorize(const std::set<std::pair<int, int> >& forbiddenDayShifts);

    // Forbids some starting days
    void forbidStartingDays(const std::set<int>& forbiddenStartingDays);

    // Authorizes some starting days
    void authorizeStartingDays(const std::set<int>& forbiddenStartingDays);

    // Know if node
    inline bool isDayShiftForbidden(int k, int s) const { return !dayShiftStatus_[k][s]; }

    inline bool isStartingDayforbidden(int k) const { return !startingDayStatus_[k]; }

    // Forbid a node / arc
    void forbidDayShift(int k, int s);

    void forbidStartingDay(int k);

    // Authorize a node / arc
    void authorizeDayShift(int k, int s);

    void authorizeStartingDay(int k);

    void resetAuthorizations();

    // Given an arc, returns the normal travel time (i.e. travel time when authorized)
    // Test for random forbidden day-shift
    std::set<std::pair<int, int> > randomForbiddenShifts(int nbForbidden);


  public:
    // DBG
    Tools::Timer *timeInS_;
    Tools::Timer *timeInNL_;
};


#endif /* SUBPROBLEM_H_ */
