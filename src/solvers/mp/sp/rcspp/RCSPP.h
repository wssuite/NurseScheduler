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

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "Parameters.h"

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

  // get the number of labels in the pool
  int64_t nLabels() const {
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

class RCSPPSolver {
 public:
  explicit RCSPPSolver(PRCGraph pRcGraph,
                       SubProblemParam param,
                       int seed = 0,
                       RCSPPSolver* pSolver = nullptr);

  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm
  std::vector<RCSolution> solve(double maxReducedCostBound, bool relaxation);

  bool isLastRunOptimal() const {
    return isOptimal_;
  }

  bool hasANextExecutionAvailable() const {
    return !timeRanOut_ && searchLevel_ < maxSearchLevel_;
  }

  double minRedCostLB() const {
    return redCostLB_;
  }

  void setParams(const SubProblemParam &p) {
    param_ = p;
  }

  // set the initial label at the source of the graph
  void setSourceLabels(const std::vector<PRCLabel> &pLabels) {
    pLSources_ = pLabels;
  }

  // set the values of minimum costs from each node to the sinks
  void setMinimumCostToSinks(vector<double> minCosts) {
    minimumCostToSinks_ = std::move(minCosts);
  }

  // reset the solver state without acting on the graph
  void reset() {
    resetLabels();
  }

  void resetLabels();

  // initialization of  the vectors that will contain the labels of each node
  void initializeLabels() {
    int nNodes = pRcGraph_->nNodes();
    pExpandedLabelsPerNode_.resize(nNodes);
    pLabelsToExpandPerNode_.resize(nNodes);
    pLabelsNoExpandPerNode_.resize(nNodes);
  }

  // attach the job running the rcssp solver
  void attachJob(Tools::Job job) {
    job_ = job;
  }

  // reset search level and global parameters
  void resetSearchLevel() {
    maxSearchLevel_ = OPTIMAL_;
    std::lock_guard<std::recursive_mutex> l(gMutex_);
    gCanComputeLB_ = true;
    gMaxSearchLevel_ = OPTIMAL_;
  }

 protected:
  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm
  std::vector<RCSolution> solve();

  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm and ignoring some of the resources
  std::vector<RCSolution> solveRelaxation();

  // alert other RCSPPSolver that time ran out for this solver
  void timeRanOut();

  // check if there is enough time left to continue the resolution
  // if not, update timeRanOut
  bool enoughTimeLeft();

  // add a label to expand to a node
  void addLabelToExpand(const PRCLabel &pL) {
    pLabelsToExpandPerNode_[pL->getNode()->id].push_back(pL);
  }

  // compute a LB by relaxing many resources
  // if relaxedLBs, do not take into account the LBs of the resources
  std::vector<RCSolution> computeValidLB(bool relaxedLBs = false);

  // forward label-setting algorithm on the input list of topologically
  // sorted nodes, until the input day is reached
  std::vector<PRCLabel> forwardLabelSetting(
      const std::vector<PRCNode> &sortedNodes, int finalDay);

  // expand all the labels at the predecessors of the input node
  // just process arcs for which the origin day is before predecessorMaxDay
  // return true if any labels has been expanded
  void pullLabelsFromPredecessors(const PRCNode &pN, int predecessorMaxDay);

  // expand a label along an arc given the id of this arc
  bool expand(const PRCLabel &pLParent,
              const PRCArc &pArc,
              const PRCLabel &pLChild);

  // get the pareto-optimal set of set of labels for a given node
  void checkAllDominations(
      const vector<PRCLabel>::iterator &begin,
      const vector<PRCLabel>::iterator &end,
      DominationStatus (*domFunction)(RCLabel *, RCLabel *),
      vector<DominationStatus> *domStatus,
      vector<DominationStatus> *domStatusNoExpand);

  // check for domination between the newly expanded labels and those that
  // are already stored at the node
  void checkDominationsWithPreviousLabels(
      const vector<PRCLabel>::iterator &begin,
      const vector<PRCLabel>::iterator &end,
      DominationStatus (*domFunction)(RCLabel *, RCLabel *),
      vector<DominationStatus> *domStatus,
      vector<DominationStatus> *domStatusNoExpand);

  // check for pairwise domination between the newly expanded labels
  void checkDominationsPairwise(
      const vector<PRCLabel>::iterator &begin,
      const vector<PRCLabel>::iterator &end,
      DominationStatus (*domFunction)(RCLabel *, RCLabel *),
      vector<DominationStatus> *domStatus);

  // select the non-dominated labels that should be expanded
  void selectLabelsToExpand(const vector<PRCLabel>::iterator &begin,
                            const vector<PRCLabel>::iterator &end,
                            const vector<DominationStatus> &domStatus,
                            const vector<DominationStatus> &domStatusNoExpand,
                            const PRCNode &pN);

  // Create a solution object given a label coming from a sink node
  static RCSolution createSolution(
      const PRCLabel &finalLabel, const PRCGraph pRcGraph = nullptr);

  // Check if a label can produce a path to a sink node with negative cost
  bool hasPotentialImprovingPathToSinks(const PRCLabel &pl, int nodeId,
                                        double primalBound = 0) const;

  // get the number of negative cost labels in the input vector and update
  // the value of the input minimum cost pointer
  int getNbLabelsBelowMaxBound(const vector<PRCLabel> &labels,
                               double *pMinCost) const;

  // if a heuristic label-setting was just executed and either optimality was
  // required or no negative-cost roster was found, we update the lsts of
  // labels and the parameters values before solving one more time
  // a heuristic may be useful when aiming at optimality, because it may
  // improve the best primal bound
  // return true if a next execution should be done
  bool prepareForNextExecution(const vector<PRCNode> &nodes);

  // backward label-setting from the sink nodes on the input list of
  // topologically sorted nodes and until the final day is treated
  vector<PRCLabel> backwardLabelSetting(
      const vector<PRCNode> &sortedNodes,
      const std::vector<PRCNode> &mergingNodes,
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
  void displaySolveStatistics() const;

 private:
  // Graph where the RCSPP is solved
  const PRCGraph pRcGraph_;
  // Factory to create labels
  PRCLabelFactory pFactory_;
  // Temporary pool of labels
  LabelPool labelPool_;
  // Initial label at the source of the graph
  vector<PRCLabel> pLSources_;
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
  SubProblemParam param_;
  // verbose
  double maxReducedCostBound_;
  double redCostLB_;
  // store if has been solved to optimality
  bool isOptimal_;
  // level for dssr.
  // For most hard constraints, if the ub >= dssrLvl, the resouce is inactive
  // for domination.
  // WARNING: dssrLvl_ == 0 is by convention used for optimality search level
  int dssrLvl_;
  bool useDominateNoLb_;
  // Total number of not dominated labels during the algorithm's execution
  int total_number_of_nondominated_labels_;
  // Total number of labels that has stopped due to the minimum cost from sinks
  // strategy
  int number_of_infeasible_deleted_labels_;
  std::vector<int> number_of_infeasible_deleted_labels_per_resource_;
  // Total number of generated labels during the algorithm's execution
  int total_number_of_generated_labels_;
  // Total number of domination operation during the algorithm's execution
  int total_number_of_dominations_;

  std::minstd_rand rdm_;
  Tools::Timer timer_;
  bool timeRanOut_;
  double maxSolvingTime_;
  // store the job running the solver
  Tools::Job job_;

  // parameter to control the intensity of the search locally
  enum SearchLevel { NB_TO_EXPAND_, DSSR_NO_LB_, INCR_DSSR_NO_LB_,
    DSSR_LB_, INCR_DSSR_LB_, OPTIMAL_ };
  static const std::map<std::string, SearchLevel> searchLevelByName_;
  static const std::map<SearchLevel, std::string> namesBySearchLevel_;
  static std::recursive_mutex gMutex_;
  static bool gCanComputeLB_;  // true if not taking too much time
  static bool gWarningHasBeenPrinted_;
  static SearchLevel gMaxSearchLevel_;
  SearchLevel searchLevel_, maxSearchLevel_;
  bool computingRelaxation_;

  // set search level based on current parameters
  void initSearchLevel();
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCSPP_H_
