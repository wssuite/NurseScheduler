/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#ifndef SRC_SOLVERS_MP_SP_SUBPROBLEM_H_
#define SRC_SOLVERS_MP_SP_SUBPROBLEM_H_

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/BoostRCSPP.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/PrincipalGraph.h"
#include "solvers/mp/MasterProblem.h"

// Parameters (called in the solve function)
struct SubproblemParam {
  SubproblemParam() {}
  SubproblemParam(int strategy, PLiveNurse pNurse, MasterProblem *pMaster) {
    initSubproblemParam(strategy, pNurse, pMaster);
  }
  ~SubproblemParam() {}

  void initSubproblemParam(int strategy,
                           PLiveNurse pNurse,
                           MasterProblem *pMaster);

  // *** PARAMETERS ***
  int verbose_ = 0;
  double epsilon_ = 1e-5;

  SPSearchStrategy search_strategy_ = SP_BREADTH_FIRST;
  RCSPPType rcspp_type_ = BOOST_LABEL_SETTING;

  // if positive, the RC SPP stops as soon as it has generated enough
  // non-dominated path
  int nb_max_paths_ = -1;

  static const int maxSubproblemStrategyLevel_;

  // 0 -> no short rotations
  // 1 -> day-0 short rotations only
  // 2 -> day-0 and last-day short rotations only
  // 3 -> all short rotations
  int shortRotationsStrategy_ = 3;

  // maximal length for a rotation
  int maxRotationLength_ = 0;

  // true -> authorize the violation of the consecutive constraints
  // false -> forbid any arc that authorized to violate the consecutive
  // constraints
  bool violateConsecutiveConstraints_ = true;

  // true  -> one sink node per day
  // false -> one single sink node for the network
  bool oneSinkNodePerLastDay_ = true;

  // sortLabelsOption for the subProblem
  // true if labels are sorted by increasing costs before domination is checked
  bool sp_sortLabelsOption = false;
  // minimumCostFromSinksOption for the subProblem
  // true if the shortest path from each node to each sink is computed to
  // delete the labels that cannot be extended into a negative cost path
  bool sp_minimumCostFromSinksOption = false;
  // worstCaseCostOption for the subProblem
  // true if the domination rule is improved by considering the worst costs that
  // can result from the violation of soft constraints
  bool sp_worstCaseCostOption = true;
  // enumeratedSubPathOption for the subProblem
  // true if we enumerate subpaths in the RCSPP graph in order to reduce the
  // number of resources
  bool sp_enumeratedSubPathOption = false;

  // Getters for the class fields
  int maxRotationLength() { return maxRotationLength_; }
  int shortRotationsStrategy() { return shortRotationsStrategy_; }
  bool oneSinkNodePerLastDay() { return oneSinkNodePerLastDay_; }

  // Setters
  void maxRotationLength(int value) { maxRotationLength_ = value; }
  void shortRotationsStrategy(int value) { shortRotationsStrategy_ = value; }
  void oneSinkNodePerLastDay(bool value) { oneSinkNodePerLastDay_ = value; }
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

  SubProblem(PScenario scenario,
             int nDays,
             PLiveNurse pNurse,
             const SubproblemParam &param);
  virtual ~SubProblem();

  // Constructor that sets the subproblem of any nurse with a given
  // contract : this correctly sets the resource (time + bounds), but NOT THE
  // COST
  SubProblem(PScenario scenario,
             int nbDays,
             PConstContract contract,
             std::vector<State> *pInitState);

  // Initialization function for all global variables (not those of the rcspp)
  virtual void init(const std::vector<State> &pInitState);

  virtual void build() = 0;

  // Solve : Returns TRUE if negative reduced costs path were found;
  // FALSE otherwise.
  virtual bool solve(
      PLiveNurse nurse,
      DualCosts *costs,
      const SubproblemParam &param,
      const std::set<std::pair<int, int>> &forbiddenDayShifts = {},
      double redCostBound = 0);

  // Returns all rotations saved during the process of solving the SPPRC
  const std::vector<RCSolution> &getSolutions() const {
    return theSolutions_;
  }


  // Some getters
  PScenario scenario() const { return pScenario_; }

  PConstContract contract() const { return pContract_; }

  const PLiveNurse liveNurse() const { return pLiveNurse_; }

  int nDays() const { return nDays_; }

  int maxRotationLength() const { return maxRotationLength_; }

  int nPaths() const { return nPaths_; }

  int nFound() const { return nFound_; }

  double bestReducedCost() const { return bestReducedCost_; }

  double cpuInLastSolve() { return timerSolve_->dSinceStart(); }

  // Returns true if the corresponding shift has no maximum limit of
  // consecutive worked days
  bool isUnlimited(int shift_type) const {
    int maxCons = shift_type ? pScenario_->maxConsShiftsOf(shift_type)
                             : pContract_->maxConsDaysOff_;
    return maxCons
        >= std::min(nDays_ + maxOngoingDaysWorked_, NB_SHIFT_UNLIMITED);
  }

  int minCons(int shift_type) const {
    return shift_type ? pScenario_->minConsShiftsOf(shift_type)
                      : pContract_->minConsDaysOff_;
  }

  int maxCons(int shift_type) const {
    if (isUnlimited(shift_type))
      return shift_type ? pScenario_->minConsShiftsOf(shift_type)
                        : pContract_->minConsDaysOff_;
    return shift_type ? pScenario_->maxConsShiftsOf(shift_type)
                      : pContract_->maxConsDaysOff_;
  }

  Penalties initPenalties() const;

  std::vector<int> defaultLBs() const {
    return {0, 0, 0};
  }

  std::vector<int> defaultUBs() const {
    return {maxRotationLength_,
            maxTotalDuration_,
            pScenario_->nbWeeks()};
  }

  // Print and check functions.
  void printAllSolutions() const;

  void printForbiddenDayShift() const;

  void checkForbiddenDaysAndShifts(const RCSolution &sol) const;

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
  PScenario pScenario_;

  // Number of days of the scenario (usually a multiple of 7)
  int nDays_;

  // Contract type
  PConstContract pContract_;

  // (Minimum) number of paths to return to the MP
  int nPathsMin_;

  // Current live nurse considered
  PLiveNurse pLiveNurse_;

  // All costs from Master Problem
  DualCosts *pCosts_;

  // Bound on the reduced cost: if greater than this, the rotation is not added
  double maxReducedCostBound_;

  // Maximum number of consecutive days already worked by a nurse before
  // the beginning of that period
  int maxOngoingDaysWorked_;

  //----------------------------------------------------------------
  //
  // Answers: rotations, number of paths found, timers
  //
  //----------------------------------------------------------------

  // Saved solutions
  std::vector<RCSolution> theSolutions_;

  // Number of paths found
  int nPaths_;

  // Number of rotations found (that match the bound condition)
  // at that iteration
  int nFound_;

  // Best reduced cost found
  double bestReducedCost_;

  // Timer in presolve and in label setting
  Tools::Timer *timerPresolve_;
  Tools::Timer *timerSolve_;

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

  // For each day k (<= nDays_ - CDMin), contains WEIGHT_COMPLETE_WEEKEND if
  // [it is a Saturday (resp. Sunday) AND the contract requires complete
  // weekends]; 0 otherwise.
  std::vector<double> startWeekendCosts_, endWeekendCosts_;
  // Costs due to preferences of the nurse:
  // for each day k (<= nDays_ - CDMin), shift s,
  // contains WEIGHT_PREFERENCES if (k,s) is a preference of the nurse;
  // 0 otherwise.
  vector2D<double> preferencesCosts_;

  // Minimum number of consecutive days worked for free
  int CDMin_;
  // principal node network begins at this index-1;
  // 1 if no ShortSucc, CDMin otherwise
  int minConsDays_;
  // Labels to take into account when solving
  std::vector<LABEL> labels_;
  // helper  to compute the penalty associated to a label
  Penalties penalties_;
  // MAXIMUM LENGTH OF A ROTATION (in consecutive worked days)
  int maxRotationLength_;
  // Maximum time duration of a roster
  int maxTotalDuration_;

  //-----------------------
  // THE GRAPH
  //-----------------------

  // Creates all nodes of the rcspp (including resource window)
  virtual void createNodes() = 0;

  // Creates all arcs of the rcspp
  virtual void createArcs() = 0;

  // Updates the costs depending on the reduced costs given for the nurse
  virtual void updateArcDualCosts() = 0;

  // FUNCTIONS -- SOLVE
  virtual bool preprocess();

  virtual bool solveRCGraph() = 0;

  // Initializes some cost vectors that depend on the nurse
  virtual void initStructuresForSolve();

  // Resets all solutions data (rotations, number of solutions, etc.)
  void resetSolutions();

  std::vector<int> startConsumption(int day, std::vector<int> shifts) const;

  // FORBIDDEN ARCS AND NODES
  vector2D<bool> dayShiftStatus_;

  // Returns true if the succession succ starting on day k does not violate
  // any forbidden day-shift
  bool canSuccStartHere(int k, const std::vector<int> &shifts) const;

  // Forbids some days / shifts
  void forbid(const std::set<std::pair<int, int> > &forbiddenDayShifts);

  // Authorizes some days / shifts
  void authorize(const std::set<std::pair<int, int> > &forbiddenDayShifts);

  // forbid any arc that authorizes the violation of a consecutive constraint
  virtual void forbidViolationConsecutiveConstraints() = 0;

  // Know if node
  bool isDayShiftForbidden(int k, int s) const {
    return !dayShiftStatus_[k][s];
  }

  // Forbid a node / arc
  virtual void forbidDayShift(int k, int s) = 0;

  // Authorize a node / arc
  virtual void authorizeDayShift(int k, int s) = 0;

  virtual void resetAuthorizations() = 0;

  // Given an arc, returns the normal travel time (i.e. travel time when
  // authorized).
  // Test for random forbidden day-shift
  std::set<std::pair<int, int> > randomForbiddenShifts(int nbForbidden);

  virtual void updatedMaxRotationLengthOnNodes(int maxRotationLength) {}
};

class BoostSubProblem : public SubProblem {
 public:
  BoostSubProblem(): SubProblem() {}

  BoostSubProblem(PScenario scenario,
                  int nDays,
                  PLiveNurse pNurse,
                  const SubproblemParam &param):
      SubProblem(scenario, nDays, pNurse, param) {}

  BoostSubProblem(PScenario scenario, int nDays, PConstContract contract,
                  std::vector<State> *pInitState):
      SubProblem(scenario, nDays, contract, pInitState) {}

  RCGraph &g() { return g_; }

  int addSingleNode(NodeType type,
                    std::vector<int> lbs = {},
                    std::vector<int> ubs = {},
                    bool hard_lbs = false) {
    if (lbs.empty())
      lbs = defaultLBs();
    if (ubs.empty())
      ubs = defaultUBs();
    return g_.addSingleNode(type, lbs, ubs, hard_lbs);
  }

  int addSingleArc(int origin,
                   int destination,
                   double baseCost,
                   std::vector<int> consumptions,
                   ArcType type,
                   int day,
                   int shift) {
    return g_.addSingleArc(origin,
                           destination,
                           baseCost,
                           consumptions,
                           type,
                           day,
                           {shift});
  }

  void build() override;

  int addSingleArc(int origin,
                   int destination,
                   double baseCost,
                   std::vector<int> consumptions,
                   ArcType type,
                   int day = -1,
                   std::vector<int> shifts = {}) {
    return g_.addSingleArc(origin,
                           destination,
                           baseCost,
                           consumptions,
                           type,
                           day,
                           shifts);
  }

  // TODO(JO): I moved the definition here, but it seems to be related to
  //  rotations rather than to rosters and rotations. This also impacts the
  //  solution method of SubProblem, which calls this method
  // Updates the maximum arrival time at all nodes so that the maximum length
  // of a rotation is updated.
  void updatedMaxRotationLengthOnNodes(int maxRotationLength) override;

  // TODO(JO): I'm a bit lost by the way all these cost functions work, it
  //  would be nice to have more comments here
  // TODO(JO): and same comment as in MasterProblem, it seems that
  //  startWorkCost and endWorkCost belong to BoostRotationSP
  virtual double startWorkCost(int a) const;

  virtual double shiftCost(int a, bool first_day = false) const;

  virtual double endWorkCost(int a) const;

  virtual double historicalCost(int a) const;

 protected:
  //-----------------------
  // THE GRAPH
  //-----------------------
  RCGraph g_;

  // Data structures for the nodes and arcs
  std::vector<PrincipalGraph> principalGraphs_;
  // Index: (shiftType, day, n, shift) of destination
  vector4D<int> arcsFromSource_;
  // Index: (shiftType, shiftType, day)
  vector3D<int> principalToPrincipal_;
  // arcs to main sink
  std::vector<int> arcsTosink_;

  // Creates all arcs of the rcspp
  void createArcs() override;

  // Updates the costs depending on the reduced costs given for the nurse
  void updateArcDualCosts() override;

  double shiftCost(const Arc_Properties &arc_prop,
                   bool first_day = false) const;

  // Create the specific types of arcs
  virtual void createArcsSourceToPrincipal() = 0;

  virtual void createArcsPrincipalToPrincipal() = 0;

  virtual void createArcsPrincipalToSink() = 0;

  // Returns true if the succession succ starting on day k does not violate
  // any forbidden day-shift
  bool canSuccStartHere(int a) const;
  bool canSuccStartHere(const Arc_Properties &arc_prop) const;

  // forbid any arc that authorizes the violation of a consecutive constraint
  void forbidViolationConsecutiveConstraints() override;

  // Forbid a node / arc
  void forbidDayShift(int k, int s) override;

  // Authorize a node / arc
  void authorizeDayShift(int k, int s) override;

  void resetAuthorizations() override;

  // FUNCTIONS -- SOLVE
  bool solveRCGraph() override;

  // return a function that will post process any path found by the RC graph
  virtual BoostRCSPPSolver *initRCSSPSolver();
};

#endif  // SRC_SOLVERS_MP_SP_SUBPROBLEM_H_
