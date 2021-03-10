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
  SubProblem(): SP() {}

  SubProblem(PScenario scenario,
                  int nDays,
                  PLiveNurse pNurse,
                  const SubproblemParam &param):
      SP(scenario, nDays, pNurse, param),
      CDMin_(pNurse->minConsDaysWork()),
      minConsDays_(1),
      maxRotationLength_(nDays) {}

  SubProblem(PScenario scenario, int nDays, PConstContract contract,
                  std::vector<State> *pInitState):
      SP(scenario, nDays, contract, pInitState),
      CDMin_(contract->minConsDaysWork_),
      minConsDays_(1),
      maxRotationLength_(nDays) {}

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

  Penalties initPenalties() const;

  std::vector<int> startConsumption(int day, std::vector<int> shifts) const;

  int maxRotationLength() const { return maxRotationLength_; }

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

  std::vector<int> defaultLBs() const {
    return {0, 0, 0};
  }

  std::vector<int> defaultUBs() const {
    return {maxRotationLength_,
            maxTotalDuration_,
            pScenario_->nWeeks()};
  }

  // TODO(JO): I moved the definition here, but it seems to be related to
  //  rotations rather than to rosters and rotations. This also impacts the
  //  solution method of SubProblem, which calls this method
  // Updates the maximum arrival time at all nodes so that the maximum length
  // of a rotation is updated.
  void updatedMaxRotationLengthOnNodes(int maxRotationLength) override;

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
  bool canSuccStartHere(int k, const std::vector<int> &shifts) const;

  // forbid any arc that authorizes the violation of a consecutive constraint
  void forbidViolationConsecutiveConstraints() override;

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
