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

  // set the initial label at the source of the graph
  void setSourceLabel(const PRCLabel& pL) {pLSource_ = pL;}

  // reset the solver state without acting on the graph
  void reset(const SubproblemParam& p) {
    setParams(p);
    resetLabels();
  }
  void resetLabels();
  void setParams(const SubproblemParam& p) {
    param_ = p;
  }

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

  // set the values of minimum costs from each node to the sinks
  void setMinimumCostToSinks(const vector<double>& minCosts) {
    minimumCostToSinks_.clear();
    minimumCostToSinks_ = minCosts;
  }

  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm
  std::vector<RCSolution> solve(double maxReducedCostBound,
                                int nb_max_paths);

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
  void pullLabelsFromSuccessors(const PRCNode &pN);
  bool expandBack(const PRCLabel &pLNext,
                  const PRCArc &pArc,
                  const PRCLabel &pLPrevious);

  // bidirectional label-setting algorithm: forward propagation from source
  // to middle-day and backward propagation from sinks to middle-day before
  // merging the labels
  vector<PRCLabel> bidirectionalLabelSetting(
      const vector<PRCNode> &sortedNodes);

  // merge the labels obtained by forward propagation with those obtained by
  // backward propagation
  void mergeLabels(vector<PRCLabel> *pForwardLabels,
                   vector<PRCLabel> *pBackwardLabels);

  // merge two labels to get a complete roster label
  bool merge(const PRCLabel &pLForward, const PRCLabel &pLBackward);

  // display information about the graph and options if verbose >= 1
  void displaySolveInputInfo();

  // display information about the execution of the label setting algorithm,
  // such as the number of expansions and dominations
  void displaySolveStatistics();

 private:
  // Graph where the RCSPP is solved
  const RCGraph *pRcGraph_;
  // Factory to create labels
  PRCLabelFactory pFactory_;
  // Temporary pool of labels
  LabelPool labelPool_;
    // Initial label at the source of the graph
  PRCLabel pLSource_;
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
  // Parameters of the solver
  SubproblemParam param_;
  // verbose
  double maxReducedCostBound_;
  // Maximal number of optimal path conserved after the algorithm's execution
  int nb_max_paths_;
  // Total number of not dominated labels during the algorithm's execution
  int total_number_of_nondominated_labels_;
  // Total number of labels deleted due to the minimum cost from sinks strategy
  int number_of_infeasible_deleted_labels_;
  // Total number of generated labels during the algorithm's execution
  int total_number_of_generated_labels_;
  // Total number of domination operation during the algorithm's execution
  int total_number_of_dominations_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCSPP_H_
