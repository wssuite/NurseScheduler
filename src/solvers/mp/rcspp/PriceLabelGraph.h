//
// Created by antoine legrain on 2020-04-11.
//

#ifndef NURSESCHEDULER_PRICELABELGRAPH_H
#define NURSESCHEDULER_PRICELABELGRAPH_H


#include <vector>
#include "RCGraph.h"

class SubProblem;

class PriceLabelGraph: public SubGraph {
  public:
    PriceLabelGraph(int lb, int ub, LABEL label,
        SubProblem* sp = nullptr, bool reset_labels_after_pricing=false);
    virtual ~PriceLabelGraph();

    void updateArcCosts();

    // link two sub graphs together:
    // 1. create an arc from the exit of inSubGraph to the current entrance
    // 2. the entrance method should now returns the entrance of the  inSubGraph
    void linkInSubGraph(SubGraph& inSubGraphs, int day=-1) override;

    // link two sub graphs together:
    // 1. create an arc from the current exit to the entrance of outSubGraph
    // 2. the exit method should now returns the exit of the outSubGraph
    void linkOutSubGraph(SubGraph& outSubGraph, int day=-1) override;

    inline int entrance(int day=-1) const override {
      if(inSubGraph_) return inSubGraph_->entrance(day);
      return entrance_;
    }

    inline int exit(int day=-1) const override {
      if(outSubGraph_) return outSubGraph_->exit(day);
      return exit_;
    }

  protected:
    SubProblem* pSP_;

    int lb_, ub_;
    LABEL label_;
    bool reset_labels_after_pricing_;

    int entrance_;				// entrance node to the subnetwork
    std::vector<int> checkNodes_;			// check nodes (measure the level of labels)
    int exit_;  // exit node from the subnetwork

    std::vector<int> in_arcs_;	// index of arcs going from entrance to check nodes
    std::vector<int> out_arcs_;	// index of arcs going from check nodes to exit

    SubGraph *inSubGraph_, *outSubGraph_;

    void build();

    bool minLabel() const {
      return label_ == MIN_CONS_DAYS || label_ == MIN_DAYS;
    }

    int getLabelCost(int l) const;
};


#endif //NURSESCHEDULER_PRICELABELGRAPH_H
