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

#include <memory>
#include <vector>

#include "solvers/mp/sp/SubProblem.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/RCSPP.h"

/**
 * Class describing the subproblems that appear in a branch-and-price
 * approach based on the decomposition by rosters
 * The main functions are to build a RCSPP graph and call an RCSPP solver to
 * get the optimal path in the graph
 */
class RosterSP : public SubProblem {
 public:
  RosterSP(PScenario scenario,
           int nbDays,
           PLiveNurse nurse,
           std::vector<PResource> pResources,
           const SubproblemParam &param);

  // FUNCTIONS -- SOLVE
  // run preprocessing algorithms, e.g., for pre-computing bounds on the
  // minimum cost to the sink(s) and on the minimum consumption of resources
  bool preprocess() override;

  // Algorithm adapted to acyclic graphs computing the costs of the shortest
  // paths from a given sink node to all the other nodes (in the reverse
  // direction of the arcs)
  vector<double> shortestPathToSinksAcyclic(const RCGraph *pRCGraph);

  //------------------------------------------------
  // Enumeration of sub paths
  //------------------------------------------------
  // This technique consists of getting rid of some resources by adding
  // supplement arcs in the rcGraph. Typically, we can delete the resources
  // corresponding to consecutive shifts types. Arcs with stretches
  // containing several shifts are added to represent a sequence of
  // successive shifts of the same type. To set up this technique, we have to
  // add arcs representing the different possibilities of sequences of shifts
  // types and we need to set the corresponding cost on these new arcs.  To
  // set up this technique, we have to  add arcs representing the different
  // possibilities of sequences of shifts types and we need to set the
  // corresponding costs of these new arcs (base cost, dual cost and cost due
  // to penalties on soft bounds for the number consecutive shifts types). Then,
  // we need to update the cost of the existing arcs (those which existed
  // before the enumeration process). We mainly add cost due to penalties on
  // soft bounds for the number of consecutive shifts types (See Antoine’s
  // internship report for more details).


  // Principal method for the enumeration of sub paths in the rcGraph
  void enumerateSubPaths(RCGraph *pRcGraph);

  // Enumeration of sub paths from the source node
  void enumerateConsShiftTypeFromSource(RCGraph *pRCGraph);

  // Enumeration of sub paths from a PRINCIPAL_NETWORK node in the rcGraph
  // for a given day and for a given shift
  void enumerateConsShiftType(RCGraph *pRCGraph,
                              const PShift& pS,
                              int day);

  // Update of the costs of the existing arcs by adding a supplement cost
  // corresponding to penalties due to the soft bounds of the consecutive
  // shifts resources
  void updateOfExistingArcsCost(RCGraph *pRCGraph);

 private:
  RCGraph rcGraph_;
  std::vector<PResource> pResources_;
  vector2D<PRCNode> pNodesPerDayShift_;  // list of nodes for each day
  shared_ptr<MyRCSPPSolver> pRcsppSolver_;  // Solver

  // Vector containing all the consumption's upper bounds of each shift
  vector<int> consShiftsUbs_;

  // Vector containing all the consumption's lower bounds of each shift
  vector<int> consShiftsLbs_;

  // Vector containing all the cost associated with the violation of soft
  // upper bounds of each shift
  vector<double> consShiftsUbCosts_;

  // Vector containing all the cost associated with the violation of soft
  // lower bounds of each shift
  vector<double> consShiftsLbCosts_;

  // Timer measuring the total time spent in enumerating sub paths in the
  // rcGraph
  Tools::Timer timerEnumerationOfSubPath_;
  // Timer measuring the total time spent in computing costs of shortest
  // paths from sink nodes to each other node in the rcGraph
  Tools::Timer timerComputeMinCostFromSink_;

  // create an arc from origin to target with the stretch {ps} starting on day
  PRCArc addSingleArc(RCGraph* pRCGraph,
                      const PRCNode &pOrigin,
                      const PRCNode &pTarget,
                      const PShift &pS,
                      int day);

  // get the base cost of a stretch when it comes just after a given previous
  // shift
  double baseCost(const Stretch &stretch, PAbstractShift pPrevShift);

  // get the dual cost of a given stretch
  double dualCost(const Stretch &stretch);

  // recover the base costs due to preferences and incomplete weekends
  void initStructuresForSolve() override {
    SubProblem::initStructuresForSolve();
  }

  // create the resource constrained graph: create nodes, arcs and resources
  void build() override;
  void build(RCGraph *pRCGraph);
  void createNodes() override {
    this->createNodes(&rcGraph_);
  }
  void createNodes(RCGraph *pRCGraph);
  void createArcs(RCGraph *pRCGraph);
  void createArcs() override {
    this->createArcs(&rcGraph_);
  }
  void createResources(RCGraph *pRCGraph);

  // create the RCSPP solver
  void createRCSPPSolver();

  // modify the RCGraph if specified in the parameters, e.g., by enumerating
  // subpaths of consecutive identical shifts,
  void preprocessRCGraph();

  // pre-compute data for dominance and expansion of the resources
  void initializeResources();

  // update arcs costs based on the dual costs coming from the master problem
  void updateArcDualCosts() override;
  void updateArcDualCost(const PRCArc &pA);

  // create the initial label that will be expanded from the source to the
  // sink(s)
  void createInitialLabel();

  // forbid any arc that authorizes the violation of a consecutive constraint
  void forbidViolationConsecutiveConstraints() override {};

  // Forbid a node / arc
  void forbidDayShift(int k, int s) override;

  // Authorize a node / arc
  void authorizeDayShift(int k, int s) override;

  void resetAuthorizations() override;

  // run the actual solution of the RCSPP once every preprocessing is done
  bool solveRCGraph() override;
};

#endif  // SRC_SOLVERS_MP_SP_ROSTERSP_H_
