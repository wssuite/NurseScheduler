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
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/MasterProblem.h"

//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------
struct DualCosts {
  explicit DualCosts(const MasterProblem *pMaster) : pMaster_(pMaster) {}
  virtual ~DualCosts() = default;

  // update the dual values of every constraints based on the current solution
  void updateDuals();

  // update the dual values of every constraints randomly
  void randomUpdateDuals(
      bool useInputData = false, int nPerturbations = 10);

  // return the dual cost of a stretch based on its consumption of
  // every constraints
  double getCost(int nurseNum,
                 const Stretch &st,
                 const PAbstractShift &prevS) const;

  vector<double> getMaxDualValues() const;

  std::string toString() const;
  std::string toString(int nurseNum, const Stretch &st) const;

  const MasterProblem *pMaster_;
};

typedef std::shared_ptr<DualCosts> PDualCosts;

class SubProblem {
 public:
  SubProblem();

  SubProblem(PScenario scenario,
             int firstDayId,
             int nDays,
             PLiveNurse pNurse,
             const SubProblemParam& param);

  virtual ~SubProblem();

  // Initialization function for all global variables (not those of the rcspp)
  virtual void init(const std::vector<State> &pInitState);

  virtual void build() = 0;

  // Solve : Returns TRUE if negative reduced costs path were found;
  // FALSE otherwise.
  virtual bool solve(
      const PDualCosts &costs,
      const std::set<std::pair<int, int>> &forbiddenDayShifts = {},
      double redCostBound = 0,
      bool relaxation = false);

  // Returns all rotations saved during the process of solving the SPPRC
  const std::vector<RCSolution> &getSolutions() const {
    return theSolutions_;
  }

  // Some getters
  PScenario pScenario() const { return pScenario_; }

  PContract pContract() const { return pLiveNurse_->pContract_; }

  int nDays() const { return nDays_; }

  const PLiveNurse &pLiveNurse() const { return pLiveNurse_; }

  double bestReducedCost() const { return bestReducedCost_; }

  double cpuInLastSolve() { return timerSolve_.dSinceStart(); }

  virtual bool isLastRunOptimal() const = 0;
  virtual bool hasANextExecutionAvailable() const {
    return true;
  }

  // return a lower bound on the reduced cost. By default, worst case
  virtual double minRedCostLB() const { return -DBL_MAX; }

  // reset parameters of the subproblems. Used to give a change to the solver
  // to change their parameters. It will be used after a node has been fathomed
  virtual void updateParameters(bool useMoreTime = false) = 0;

  // Print and check functions.
  void printAllSolutions() const;

  void printForbiddenDayShift() const;

  void checkForbiddenDaysAndShifts(const RCSolution &sol) const;

  virtual void computeCost(MasterProblem *pMaster, RCSolution *rcSol) const = 0;
  void computePreferencesCost(RCSolution *rcSol) const;

  void setPreferences(const SubProblemParam &param) {
    param_ = param;
  }

  void setDefaultDomination(bool isDefault) {
    param_.rcsppImprovedDomination_ = !isDefault;
  }

  void setBidirectional(bool enable) {
    param_.rcsppBidirectional_ = enable;
  }

  void setSpMaxSolvingTimeRatioForRelaxation(double ratio) {
    param_.spMaxSolvingTimeRatioForRelaxation_ = ratio;
  }

  void attachJob(Tools::Job job) {
    job_ = job;
  }

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

  // Id of the first day of the scenario
  int firstDayId_;

  // Number of days of the scenario (usually a multiple of 7)
  int nDays_;

  // (Minimum) number of paths to return to the MP
  int nPathsMin_{};

  // Current live nurse considered
  PLiveNurse pLiveNurse_;

  // All costs from Master Problem
  PDualCosts pCosts_;

  // Bound on the reduced cost: if greater than this, the rotation is not added
  double maxReducedCostBound_{};

  // Maximum number of consecutive days already worked by a nurse before
  // the beginning of that period
  int maxOngoingDaysWorked_{};

  //----------------------------------------------------------------
  //
  // Answers: rotations, number of paths found, timers
  //
  //----------------------------------------------------------------

  // Saved solutions
  std::vector<RCSolution> theSolutions_;

  // Number of paths found
  int nPaths_{};

  // Number of rotations found (that match the bound condition)
  // at that iteration
  int nFound_{};

  // Best reduced cost found
  double bestReducedCost_{};

  // Timer in presolve and in label setting
  Tools::Timer timerPresolve_, timerSolve_, timerPostsolve_;

  //----------------------------------------------------------------
  //
  // Solving options.
  //
  //----------------------------------------------------------------

  SubProblemParam param_;

  // store the job running the solver if it's the case
  Tools::Job job_;

  //-----------------------
  // THE BASE COSTS
  //-----------------------
  // WARNING : for those that never change, of no use also.

  // For each day k (<= nDays_ - CDMin), contains completeWeekend cost if
  // [it is a Saturday (resp. Sunday) AND the contract requires complete
  // weekends]; 0 otherwise.
  std::vector<double> startWeekendCosts_, endWeekendCosts_;
  // Costs due to preferences of the nurse:
  // for each day k (<= nDays_ - CDMin), shift s,
  // contains WEIGHT_PREFERENCES if (k,s) is a preference of the nurse;
  // 0 otherwise.
  vector2D<double> preferencesCosts_;

  // Maximum time duration of a roster
  int maxTotalDuration_{};

  // random generator
  std::minstd_rand rdm_;

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
  virtual bool solve(bool initialSolve = true, bool relaxation = false);
  virtual bool presolve();
  virtual bool postprocess();

  virtual bool solveRCGraph(
      bool initialSolve = true, bool relaxation = false) = 0;

  // Initializes some cost vectors that depend on the nurse
  virtual void initStructuresForSolve();

  // Resets all solutions data (rotations, number of solutions, etc.)
  void resetSolutions();

  // FORBIDDEN ARCS AND NODES
  vector2D<bool> dayShiftStatus_;

  // Forbids some days / shifts
  void forbid(const std::set<std::pair<int, int> > &forbiddenDayShifts);

  // Authorizes some days / shifts
  void authorize(const std::set<std::pair<int, int> > &forbiddenDayShifts);

  // Know if node
  virtual bool isDayShiftForbidden(int k, int s) const {
    return !dayShiftStatus_[k][s];
  }

  // Forbid a node / arc
  virtual void forbidDayShift(int k, int s) {
    // Mark the day-shift as forbidden
    dayShiftStatus_[k][s] = false;
  }

  // Authorize a node / arc
  virtual void authorizeDayShift(int k, int s) {
    // Mark the day-shift as forbidden
    dayShiftStatus_[k][s] = true;
  }

  virtual void resetAuthorizations() {
    // set all value to true
    Tools::initVector2D(&dayShiftStatus_, nDays_, pScenario_->nShifts(), true);
  }

  // Given an arc, returns the normal travel time (i.e. travel time when
  // authorized).
  // Test for random forbidden day-shift
  std::set<std::pair<int, int> > randomForbiddenShifts(int nbForbidden);
};

#endif  // SRC_SOLVERS_MP_SP_SUBPROBLEM_H_
