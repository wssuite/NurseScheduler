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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RCSPP_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RCSPP_H_

#include <limits>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/SubProblem.h"

/**
 * Class used to create an area of storage for labels that prepare to be
 * dominated
 */

class LabelPool {
 public:
  explicit LabelPool(PRCLabelFactory pFactory,
                     int bunchSize):
      pFactory_(pFactory),
      pLabels_() {
    pLabels_.reserve(bunchSize);
    endLabel_ = pLabels_.end();
  }

  // increase the size of the pool
  PRCLabel getNewLabel();

  // reduce by 1 the size of the pool
  void releaseLastLabel() {
    if (endLabel_ > pLabels_.begin())
      endLabel_--;
  }

  // clear the pool: allow to reuse labels
  void clear() {
    endLabel_  = pLabels_.begin();
  }

  // sort the active labels by cost in increasing or decreasing order
  void sort();

  // get a copy of all labels in the pool
  // TODO(JO): I am getting an error here on Linux
  /* vector<PRCLabel> getLabels() const {
     return std::vector<PRCLabel>(pLabels_.begin(), endLabel_);
   }
   */

  // get the number of labels in the pool
  int nLabels() const {
    return std::distance<vector<PRCLabel>::const_iterator>(
        pLabels_.begin(), endLabel_);
  }

  vector<PRCLabel>::iterator begin() {
    return pLabels_.begin();
  }

  vector<PRCLabel>::iterator end() {
    return endLabel_;
  }

 private:
  PRCLabelFactory pFactory_;
  vector<PRCLabel> pLabels_;  // all the labels in the pool
  vector<PRCLabel>::iterator endLabel_;  // end of the current used labels
};


/**
 * Class containing the label correcting algorithm to solve the RCSPP
 */

class MyRCSPPSolver {
 public:
  explicit MyRCSPPSolver(RCGraph *pRcGraph,
                         const SubproblemParam &param);

// reset the solver state without acting on the graph
  void reset();

  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm
  std::vector<RCSolution> solve(double maxReducedCostBound,
                                int verbose,
                                double epsilon,
                                SPSearchStrategy strategy,
                                int nb_max_paths);

  // initialization of  the vectors that will contain the labels of each node
  void initializeLabels() {
    int nNodes = pRcGraph_->nNodes();
    pExpandedLabelsPerNode_.resize(nNodes);
    pLabelsToExpandPerNode_.resize(nNodes);
    pLabelsNoExpandPerNode_.resize(nNodes);
  }

  // add a label to expand to a node
  void addLabelToExpand(const PRCLabel &pL) {
    pLabelsToExpandPerNode_[pL->getNode()->id].push_back(pL);
  }

  // forward label-setting algorithm on the input list of topologically
  // sorted nodes, until the input day is reached
  std::vector<PRCLabel> forwardLabelSetting(
      const std::vector<PRCNode> &sortedNodes, int finalDay = INT_MAX);

  // expand all the labels at the predecessors of the input node
  void pullLabelsFromPredecessors(const PRCNode &pN);

  // expand a label along an arc given the id of this arc
  bool expand(const PRCLabel &pLParent,
              const PRCArc &pArc,
              const PRCLabel &pLChild);

  // get the pareto-optimal set of set of labels for a given node
  void checkAllDominations(const vector<PRCLabel>::iterator &begin,
                           const vector<PRCLabel>::iterator &end,
                           int (*domFunction)(const PRCLabel &,
                                              const PRCLabel &,
                                              const vector<PResource> &),
                           vector<int> *domStatus,
                           vector<int> *domStatusNoExpand);

  // check for domination between the newly expanded labels and those that
  // are already stored at the node
  void checkDominationsWithPreviousLabels(
      const vector<PRCLabel>::iterator &begin,
      const vector<PRCLabel>::iterator &end,
      const vector<PResource> &resources,
      int (*domFunction)(const PRCLabel &,
                         const PRCLabel &,
                         const vector<PResource> &),
      vector<int> *domStatus,
      vector<int> *domStatusNoExpand);
  // check for pairwise domination between the newly expanded labels
  void checkDominationsPairwise(const vector<PRCLabel>::iterator &begin,
                                const vector<PRCLabel>::iterator &end,
                                const vector<PResource> &resources,
                                int (*domFunction)(const PRCLabel &,
                                                   const PRCLabel &,
                                                   const vector<PResource> &),
                                vector<int> *domStatus);

  // select the non-dominated labels that should be expanded
  void selectLabelsToExpand(const vector<PRCLabel>::iterator &begin,
                            const vector<PRCLabel>::iterator &end,
                            const vector<int>& domStatus,
                            const vector<int>& domStatusNoExpand,
                            const PRCNode &pN);

  // Create a solution object given a label coming from a sink node
  static RCSolution createSolution(const PRCLabel& finalLabel);

  // Algorithm adapted to acyclic graphs computing the costs of the shortest
  // paths from a given sink node to all the other nodes (in the reverse
  // direction of the arcs)
  vector<double> shortestPathToSinksAcyclic(const RCGraph *pRCGraph);

  // Computation of the costs of the shortest paths from EACH sink node to
  // all the other nodes in the graph
  void computeMinimumCostToSinks(RCGraph *pRCGraph);

  // Check if a label can produce a path to a sink node with negative cost
  bool hasPotentialImprovingPathToSinks(const PRCLabel &pl, const PRCNode& pN,
                                        double primalBound = 0) const;

  // get the number of negative cost labels in the input vector and update
  // the value of the input minimum cost pointer
  int getNbNegativeCostLabels(const vector<PRCLabel> &labels,
                              double *pMinCost);

  // if a heuristic label-setting was just executed and either optimality was
  // required or no negative-cost roster was found, we update the lsts of
  // labels and the parameters values before solving one more time
  // a heuristic may be useful when aiming at optimality, because it may
  // improve the best primal bound
  void prepareForNextExecution(const vector<PRCNode> &nodes);

  // backward label-setting from the sink nodes on the input list of
  // topologically sorted nodes and until the final day is treated
  vector<PRCLabel> backwardLabelSetting(const vector<PRCNode> &sortedNodes,
                                        int finalDay = 0);
  void pushLabelsToPredecessors(const PRCNode &pN);
  bool expandBack(const PRCLabel &pLNext,
                  const PRCArc &pArc,
                  const PRCLabel &pLPrevious);

  // display information about the graph and options if verbose >= 1
  void displaySolveInputInfo();

  // display information abot the execution of the label setting algorithm,
  // such as the number of expansions and dominations
  void displaySolveStatistics();

 private:
  // Graph where the RCSPP is solved
  const RCGraph *pRcGraph_;
  // Factory to create labels
  PRCLabelFactory pFactory_;
  // Temporary pool of labels
  LabelPool labelPool_;
  // Total number of not dominated labels during the algorithm's execution
  int nParetoLabels_;
  // Non-dominated labels at each node of the graph
  vector2D<PRCLabel> pExpandedLabelsPerNode_;
  // Non-dominated labels that will be expanded at each node of the graph
  vector2D<PRCLabel> pLabelsToExpandPerNode_;
  // Non-dominated labels that were not expanded at each node of the graph
  vector2D<PRCLabel> pLabelsNoExpandPerNode_;
  // Costs of shortest paths from sink nodes to each nodes
  vector<double> minimumCostToSinks_;
  // Cost of the best path found until now
  double bestPrimalBound_;
  // Total number of labels deleted due to the minimum cost from sinks strategy
  int number_of_infeasible_deleted_labels_;
  // Total number of generated labels during the algorithm's execution
  int total_number_of_generated_labels_;
  // Total number of domination operation during the algorithm's execution
  int total_number_of_dominations_;
  // Parameters of the solver
  SubproblemParam param_;
  // verbose
  int verbose_;
  double maxReducedCostBound_, epsilon_;
  SPSearchStrategy strategy_;
  // Maximal number of optimal path conserved after the algorithm's execution
  int nb_max_paths_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCSPP_H_
