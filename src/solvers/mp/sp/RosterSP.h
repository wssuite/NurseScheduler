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

#ifndef SRC_SOLVERS_MP_SP_ROSTERSP_H_
#define SRC_SOLVERS_MP_SP_ROSTERSP_H_

#include <vector>

#include "solvers/mp/sp/RCSPPSubProblem.h"

/**
 * Class describing the subproblems that appear in a branch-and-price
 * approach based on the decomposition by rosters
 * The main functions are to build a RCSPP graph and call an RCSPP solver to
 * get the optimal path in the graph
 */
class RosterSP : public RCSPPSubProblem {
 public:
  RosterSP(PScenario pScenario,
           int nbDays,
           PLiveNurse nurse,
           std::vector<PResource> pResources,
           SubproblemParam param);

  void computeCost(MasterProblem *pMaster, RCSolution *rcSol) const override;

 protected:
  // get the dual cost of a given stretch
  double dualCost(const PRCArc &pArc) override;

  void createNodes(const PRCGraph &pRCGraph) override;
  void createArcs(const PRCGraph &pRCGraph) override;

  // create the initial label that will be expanded from the source to the
  // sink(s)
  void createInitialLabels() override;
};

#endif  // SRC_SOLVERS_MP_SP_ROSTERSP_H_
