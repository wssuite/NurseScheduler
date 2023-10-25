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

#ifndef SRC_SOLVERS_MP_SP_RCSPPSUBPROBLEM_H_
#define SRC_SOLVERS_MP_SP_RCSPPSUBPROBLEM_H_

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "solvers/mp/sp/SubProblem.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/RCSPP.h"

#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"

/**
 * General class describing the subproblems that appear in a branch-and-price
 * approach. Contains many functions that are usefull for any decomposittions
 * that used RCSPPSolver
 */
class RCSPPSubProblem : public SubProblem {
 public:
  RCSPPSubProblem(PScenario scenario,
                  int firstDayId,
                  int nbDays,
                  PLiveNurse nurse,
                  std::vector<PResource> pResources,
                  SubProblemParam param);

  // FUNCTIONS -- SOLVE
  // run preprocessing algorithms, e.g., for pre-computing bounds on the
  // minimum cost to the sink(s) and on the minimum consumption of resources
  bool presolve() override;

  // reset parameters of the subproblems. Used to give a change to the solver
  // to change their parameters.
  void updateParameters(bool useMoreTime = false) override;

  // verify that all the parameters are compatible.
  // If not, turn off some of them.
  void fixParameters();

  bool isLastRunOptimal() const override {
    return pRcsppSolver_ && pRcsppSolver_->isLastRunOptimal();
  }

  bool hasANextExecutionAvailable() const override {
    return pRcsppSolver_ && pRcsppSolver_->hasANextExecutionAvailable();
  }

  double minRedCostLB() const override {
    if (pRcsppSolver_) return pRcsppSolver_->minRedCostLB();
    return -DBL_MAX;
  }

  // Algorithms adapted to acyclic graphs computing the costs of the shortest
  // paths from a given sink node to all the other nodes (in the reverse
  // direction of the arcs) and from the nodes of a given day (i.e., a given
  // level of the acyclic graph) to every other nodes.
  vector<double> minCostPathToSinks(
      const PRCGraph &pRCGraph, const PRCGraph &pRCGraphOriginal = nullptr);

 protected:
  PRCGraph pRCGraph_, pEnumGraph_;
  std::vector<PResource> pResources_;
  vector2D<PRCNode> pNodesPerDay_;  // list of nodes for each day
  shared_ptr<RCSPPSolver> pRcsppSolver_;  // Solver

  // Timer measuring the total time spent in enumerating sub paths in the
  // rcGraph
  Tools::Timer timerEnumerationOfSubPath_;
  // Timer measuring the total time spent in computing costs of shortest
  // paths from sink nodes to each other node in the rcGraph
  Tools::Timer timerComputeMinCostFromSink_;

  // create an arc from origin to target with the stretch {ps} starting on day
  static PRCArc addSingleArc(const PRCGraph &pRCGraph,
                      const PRCNode &pOrigin,
                      const PRCNode &pTarget,
                      const PShift &pS,
                      int day);

  // get the dual cost of a given stretch
  double dualCost(const PRCArc &pArc);

  // recover the base costs due to preferences and incomplete weekends
  void initStructuresForSolve() override {
    SubProblem::initStructuresForSolve();
  }

  // create the resource constrained graph: create nodes, arcs and resources
  void build() override;
  void build(const PRCGraph &pRCGraph);

  void createNodes() override {
    this->createNodes(pRCGraph_);
  }
  virtual void createNodes(const PRCGraph &pRCGraph) = 0;

  virtual void createArcs(const PRCGraph &pRCGraph) = 0;
  void createArcs() override {
    this->createArcs(pRCGraph_);
  }

  void createResources(const PRCGraph &pRCGraph);

  // create the RCSPP solver
  void createRCSPPSolver(int seed = -1);

  // modify the RCGraph if specified in the parameters, e.g., by enumerating
  // subpaths of consecutive identical shifts,
  void preprocessRCGraph(const PRCGraph &pRCGraph,
                         bool forceEnum);

  // update arcs costs based on the dual costs coming from the master problem
  void updateArcDualCosts() override;
  void updateArcDualCost(const PRCArc &pA);

  // create the initial label that will be expanded from the source(s) to the
  // sink(s)
  virtual void createInitialLabels() = 0;

  // Forbid a node / arc
  void forbidDayShift(int k, int s) override;

  // Authorize a node / arc
  void authorizeDayShift(int k, int s) override;

  void resetAuthorizations() override;

  // run the actual solution of the RCSPP once every preprocessing is done
  bool solveRCGraph(bool initialSolve = true, bool relaxation = false) override;

  // compute the cost of a given rcSol
  std::pair<vector<PResource>, PRCGraph> computeResourcesCosts(
      const State &initialState,
      RCSolution *rcSol) const;
  vector<PResource> computeResourcesCosts(
      const State &initialState,
      const Stretch &stretch,
      std::map<CostType, double> *costsPerType);
};

#endif  // SRC_SOLVERS_MP_SP_RCSPPSUBPROBLEM_H_
