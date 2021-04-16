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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_BOOST_SUBPROBLEM_H_
#define SRC_SOLVERS_MP_SP_RCSPP_BOOST_SUBPROBLEM_H_

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "RCGraph.h"
#include "RCSPP.h"
#include "PrincipalGraph.h"
#include "solvers/mp/sp/SubProblem.h"

using SP = SubProblem;

namespace boostRCSPP {

class SubProblem : public SP {
 public:
  SubProblem(PScenario scenario,
             int nDays,
             PLiveNurse pNurse,
             SubProblemParam param):
      SP(std::move(scenario), nDays, std::move(pNurse), std::move(param)),
      CDMin_(pLiveNurse_->minConsDaysWork()),
      minConsDays_(1),
      maxRotationLength_(nDays),
      defaultStrategy_(param_.strategyLevel_) {
    // set initial strategy
    initSubproblemParam(defaultStrategy_);
  }

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
                   PShift pShift) {
    return g_.addSingleArc(origin,
                           destination,
                           baseCost,
                           consumptions,
                           type,
                           day,
                           {pShift});
  }

  void build() override;

  bool solve() override;

  static const int maxSubproblemStrategyLevel_;

  void initSubproblemParam(int strategy);

  void updateParameters(bool masterFeasible) override;

  bool isLastRunOptimal() const override {
    return param_.strategyLevel_ == maxSubproblemStrategyLevel_;
  }

  int addSingleArc(int origin,
                   int destination,
                   double baseCost,
                   std::vector<int> consumptions,
                   ArcType type,
                   int day = -1,
                   std::vector<PShift> pShifts = {}) {
    return g_.addSingleArc(origin,
                           destination,
                           baseCost,
                           consumptions,
                           type,
                           day,
                           pShifts);
  }

  Penalties initPenalties() const;

  std::vector<int> startConsumption(
      int day, const std::vector<PShift> &pShifts) const;

  int maxRotationLength() const { return maxRotationLength_; }

  // Returns true if the corresponding shift has no maximum limit of
  // consecutive worked days
  bool isUnlimited(int shift_type) const {
    int maxCons = shift_type ? pScenario_->maxConsShiftsOfType(shift_type)
                             : pContract_->maxConsDaysOff_;
    return maxCons
        >= std::min(nDays_ + maxOngoingDaysWorked_, NB_SHIFT_UNLIMITED);
  }

  int minCons(int shift_type) const {
    return shift_type ? pScenario_->minConsShiftsOfType(shift_type)
                      : pContract_->minConsDaysOff_;
  }

  int maxCons(int shift_type) const {
    if (isUnlimited(shift_type))
      return shift_type ? pScenario_->minConsShiftsOfType(shift_type)
                        : pContract_->minConsDaysOff_;
    return shift_type ? pScenario_->maxConsShiftsOfType(shift_type)
                      : pContract_->maxConsDaysOff_;
  }

  std::vector<int> defaultLBs() const {
    return {0, 0, 0};
  }

  std::vector<int> defaultUBs() const {
    return {maxRotationLength_,
            maxTotalDuration_,
            pScenario_->nWeeks()};
  }

  // Updates the maximum arrival time at all nodes so that the maximum length
  // of a rotation is updated.
  void updatedMaxRotationLengthOnNodes(int maxRotationLength);

  virtual double startWorkCost(int a) const = 0;

  virtual double shiftCost(int a) const = 0;

  virtual double endWorkCost(int a) const = 0;

  // return the cost of preforming current shift on first day based
  // on the nurse historical shifts
  virtual double historicalCost(int currentShift) const;

 protected:
  //-----------------------
  // THE GRAPH
  //-----------------------
  RCGraph g_;

  // Minimum number of consecutive days worked for free
  int CDMin_;
  // principal node network begins at this index-1;
  // 1 if no ShortSucc, CDMin otherwise
  int minConsDays_;
  // MAXIMUM LENGTH OF A ROTATION (in consecutive worked days)
  int maxRotationLength_;

  // Labels to take into account when solving
  std::vector<LABEL> labels_;
  // helper  to compute the penalty associated to a label
  Penalties penalties_;

  // Data structures for the nodes and arcs
  std::vector<PrincipalGraph> principalGraphs_;
  // Index: (shiftType, day, n, shift) of destination
  vector4D<int> arcsFromSource_;
  // Index: (shiftType, shiftType, day)
  vector3D<int> principalToPrincipal_;
  // arcs to main sink
  std::vector<int> arcsTosink_;

  int defaultStrategy_;

  // Creates all arcs of the rcspp
  void createArcs() override;

  // Create the specific types of arcs
  virtual void createArcsSourceToPrincipal() = 0;

  virtual void createArcsPrincipalToPrincipal() = 0;

  virtual void createArcsPrincipalToSink() = 0;

  // Returns true if the succession succ starting on day k does not violate
  // any forbidden day-shift
  bool canSuccStartHere(int a) const;
  bool canSuccStartHere(const Arc_Properties &arc_prop) const;
  bool canSuccStartHere(int k, const std::vector<PShift> &shifts) const;

  // forbid any arc that authorizes the violation of a consecutive constraint
  vector<int> forbidViolationConsecutiveConstraints();

  // Forbid a node / arc
  void forbidDayShift(int k, int s) override;

  // Authorize a node / arc
  void authorizeDayShift(int k, int s) override;

  void resetAuthorizations() override;

  // FUNCTIONS -- SOLVE
  bool solveRCGraph() override;

  // override updateArcDualCosts to take into account end work costs when going
  // to a daily sink
  void updateArcDualCosts() override;

  // return a function that will post process any path found by the RC graph
  virtual BoostRCSPPSolver *initRCSSPSolver();
};

}  // namespace boostRCSPP

#endif  // SRC_SOLVERS_MP_SP_RCSPP_BOOST_SUBPROBLEM_H_
