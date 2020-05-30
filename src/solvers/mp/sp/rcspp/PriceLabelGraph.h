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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_PRICELABELGRAPH_H_
#define SRC_SOLVERS_MP_SP_RCSPP_PRICELABELGRAPH_H_

#include <vector>

#include "solvers/mp/sp/rcspp/RCGraph.h"

class SubProblem;

// Subgraph to price a label based on its value
// A node is created for every label value possibility between the minimum and
// maximum levels.
// For example, for MAX_CONS_DAYS, a node is created for each value from
// maxConsDays (lb) to ub.
// For MIN_CONS_DAYS, a node is created for each value from 0 (lb) to ub.
// Then, the arcs price correctly each node
class PriceLabelGraph : public SubGraph {
 public:
  PriceLabelGraph(int day,
                  int lb,
                  int ub,
                  LABEL label,
                  SubProblem *sp = nullptr,
                  bool reset_labels_after_pricing = false);
  virtual ~PriceLabelGraph();

  void updateArcCosts();

  // link two sub graphs together:
  // 1. create an arc from the exit of inSubGraph to the current entrance
  // 2. the entrance method should now returns the entrance of the  inSubGraph
  void linkInSubGraph(SubGraph *inSubGraphs, int day = -1) override;

  // link two sub graphs together:
  // 1. create an arc from the current exit to the entrance of outSubGraph
  // 2. the exit method should now returns the exit of the outSubGraph
  void linkOutSubGraph(SubGraph *outSubGraph, int day = -1) override;

  int entrance(int day = -1) const override {
    if (!pSP_) return -1;
    if (inSubGraph_) return inSubGraph_->entrance(day);
    return entrance_;
  }

  int exit(int day = -1) const override {
    if (!pSP_) return -1;
    if (outSubGraph_) return outSubGraph_->exit(day);
    return exit_;
  }

 protected:
  SubProblem *pSP_;

  int day_, lb_, ub_;
  LABEL label_;
  bool reset_labels_after_pricing_;
// entrance node to the subnetwork
  int entrance_;
  // check nodes (measure the level of labels)
  std::vector<int> checkNodes_;
  // exit node from the subnetwork
  int exit_;

  // index of arcs going from entrance to check nodes
  std::vector<int> in_arcs_;
  // index of arcs going from check nodes to exit
  std::vector<int> out_arcs_;

  SubGraph *inSubGraph_, *outSubGraph_;

  void build();

  bool minLabel() const {
    return label_ == MIN_CONS_DAYS || label_ == MIN_DAYS;
  }

  double getLabelCost(int l) const;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_PRICELABELGRAPH_H_
