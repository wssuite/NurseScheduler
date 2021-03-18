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

#ifndef SRC_SOLVERS_MP_SP_CYCLICROSTERSP_H_
#define SRC_SOLVERS_MP_SP_CYCLICROSTERSP_H_

#include <list>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "solvers/mp/sp/RosterSP.h"


/**
 * Class describing the subproblems that appear in a branch-and-price
 * approach based on the decomposition by rosters
 * The main functions are to build a RCSPP graph and call an RCSPP solver to
 * get the optimal path in the graph
 * This class differentiates of RosterSP by building cyclic rosters.
 * To do so, it creates K+1 subgraphs starting on a rest shift
 * for day 0, 1, ..., K and ending on rest shift on the start day.
 * To ensure, that the pricing of all constraints are done correctly, each
 * graph also ensures that last day is a work shift, i.e, respectively
 * day nDays_ - 1, 0, ..., K-1.
 * If K is big enough, optimality should not be lost.
 * In practice, K will be set to the Nurse maxConsWorkedShift.
 * When solving, the supbroblem will rotate on the graphs used to avoid
 * uneccessary increase in computation time. It is solving several graphs,
 * only if all of the previous ones were not able to produce rosters with a
 * negative reduced cost.
 */
class OffsetRosterSP : public RosterSP {
 public:
  OffsetRosterSP(PScenario pScenario,
                 int firstDay,
                 int nbDays,
                 PLiveNurse nurse,
                 std::vector<PResource> pResources,
                 SubproblemParam param);

  void build() override;

 protected:
  void createNodes(const PRCGraph &pRCGraph) override;
  void createArcs(const PRCGraph &pRCGraph) override;

  bool isDayShiftForbidden(int k, int s) const override {
    return !dayShiftStatus_[k % nDays_][s];
  }

  int firstDay_;  // first day of the graph
};
typedef shared_ptr<OffsetRosterSP> POffsetRosterSP;

class CyclicRosterSP : public RosterSP {
 public:
  CyclicRosterSP(PScenario pScenario,
                 int nbDays,
                 PLiveNurse nurse,
                 std::vector<PResource> pResources,
                 SubproblemParam param);

  void build() override;

  // Solve : Returns TRUE if negative reduced costs path were found;
  // FALSE otherwise.
  bool solve(
      PLiveNurse nurse,
      const PDualCosts &costs,
      const SubproblemParam &param,
      const std::set<std::pair<int, int>> &forbiddenDayShifts = {},
      double redCostBound = 0) override;

 protected:
  // offset is used to generate cyclic rosters
  // the graph is solved K+1 times with the first day from 0 to K.
  // In order to be able to link the end to the beginning,
  // the first day is enforced to be a rest and the last not a rest.
  int maxOffset_;  // K+1
  std::list<POffsetRosterSP> offsetSPs_;

  POffsetRosterSP popOffsetFront();

  void addOffsetSPs(const std::list<POffsetRosterSP> &offsetSPs) {
    offsetSPs_.insert(offsetSPs_.end(), offsetSPs.begin(), offsetSPs.end());
  }
};

#endif  // SRC_SOLVERS_MP_SP_CYCLICROSTERSP_H_
