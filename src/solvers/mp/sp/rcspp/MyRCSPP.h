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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_MYRCSPP_H_
#define SRC_SOLVERS_MP_SP_RCSPP_MYRCSPP_H_

#include <limits>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/MyRCGraph.h"
#include "solvers/mp/sp/rcspp/MyRCLabel.h"

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
  explicit MyRCSPPSolver(MyRCGraph *pRcGraph);


  // reset the solver state without acting on the graph
  void reset();

  // add a label to a node
  void addLabel(const PRCLabel &pL) {
    pLabelsPerNode_[pL->getNode()->id].push_back(pL);
  }

  // Expand a label along an arc given the id of this arc
  bool expand(const PRCLabel &pLParent,
              const PRCArc &pArc,
              const PRCLabel &pLChild);

  // Get the pareto-optimal set of set of labels for a given node
  void checkAllDominations(const vector<PRCLabel>::iterator &begin,
                           const vector<PRCLabel>::iterator &end);

  // Solve the resource constrained shortest path problem in the RCGraph
  // using a label correcting algorithm
  std::vector<RCSolution> solve(double maxReducedCostBound,
                                int verbose,
                                double epsilon,
                                SPSearchStrategy strategy,
                                int nb_max_path,
                                PScenario pS,
                                vector<int> maxConsShifts,
                                vector<int> minConsShifts,
                                bool enumeratedSubPathOption = false,
                                bool sortLabelsOption = true,
                                bool minimumCostToSinksOption = false,
                                bool worstCaseCostOption = false);

  // Create a solution object given a label coming from a sink node
  RCSolution createSolution(const PRCLabel& finalLabel);

  // Display the content of a label (cost + consumptions of all the resources)
  void displayLabel(const PRCLabel &pL) const;

  // Check the domination between two labels
  // (It return true if the first label dominate the second)
  bool dominate(const PRCLabel &pL1,
                const PRCLabel &pL2,
                const vector<PResource>& resources);

  // Get the corresponding node of a label
  PRCNode getTargetNode(const PRCLabel &pl) const;

  // Display the path in the graph followed by a label coming from a sink node
  void backTracking(const PRCLabel &pl);

  // Algorithm adapted to acyclic graphs computing the costs of the shortest
  // paths from a given sink node to all the other nodes (in the reverse
  // direction of the arcs)
  vector<double> shortestPathAcyclic(const PRCNode &sink);

  // Computation of the costs of the shortest paths from EACH sink node to
  // all the other nodes in the graph
  void computeMinimumCostToSinks();

  // Check if a label can produce a path to a sink node with negative cost
  bool hasPotentialNegativePathToSinks(const PRCLabel &pl) const;


 private:
  // Graph where the RCSPP is solved
  const MyRCGraph *pRcGraph_;
  // Factory to create labels
  PRCLabelFactory pFactory_;
  // Temporary pool of labels
  LabelPool labelPool_;
  // Total number of not dominated labels during the algorithm's execution
  int nParetoLabels_;
  // Labels corresponding to each node of the graph
  vector2D<PRCLabel> pLabelsPerNode_;
  // Costs of shortest paths from sink nodes to each nodes
  vector<double> minimumCostToSinks_;
  // Total number of labels deleted due to the 'minimum cost from
  // sinks' strategy
  int number_of_infeasible_deleted_labels_;
  // Total number of generated labels during the algorithm's execution
  int total_number_of_generated_labels_;
  // Total number of domination operation during the algorithm's execution
  int total_number_of_dominations_;
  // verbose
  int verbose_;
  double maxReducedCostBound_, epsilon_;
  SPSearchStrategy strategy_;
  // Maximal number of optimal path conserved after the algorithm's execution
  int nb_max_paths_;
  // Option to sort labels before each domination operation
  bool sortLabelsOption_ = false;
  // Option to check if a label can produce a path until a sink with a
  // negative cost
  bool minimumCostToSinksOption_ = false;
  // Option to use the improved domination technique
  bool worstCaseCostOption_ = false;
  // Option to solve the RCSPP in the graph with the enumeration of sub paths
  // technique
  bool enumeratedSubPathOption_ = false;
  // Timer used to measure the execution time of the algorithm
  Tools::Timer timer_;
};


#endif  // SRC_SOLVERS_MP_SP_RCSPP_MYRCSPP_H_
