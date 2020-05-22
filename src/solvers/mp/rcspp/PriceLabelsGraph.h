//
// Created by antoine legrain on 2020-04-11.
//

#ifndef NURSESCHEDULER_PRICELABELSGRAPH_H
#define NURSESCHEDULER_PRICELABELSGRAPH_H


#include <vector>

class SubProblem;

// Subgraph to price a label based on its value
// A node is created for every label value possibility between the minimum level.
// For example, for MAX_CONS_DAYS, a node is created for each value from maxConsDays to max_ub.
// For MIN_CONS_DAYS, a node is created for each value from 0 to min_ub.
// Then, the arcs price correctly each node
class PriceLabelsGraph {
  public:
    PriceLabelsGraph(int min_ub, int max_ub, SubProblem* sp = nullptr,
        bool compute_min_cost=true, bool reset_labels_after_pricing=false);
    virtual ~PriceLabelsGraph();

    void updateArcCosts();

    inline int entrance() const { return entrance_; }

    inline int exit() const { return exit_; }

  protected:
    SubProblem* pSP_;

    int min_ub_, max_ub_;
    bool compute_min_cost_, reset_labels_after_pricing_;

    int entrance_;				// entrance node to the subnetwork
    std::vector<int> checkNodes_;			// check nodes (measure the level of labels)
    int exit_;  // exit node from the subnetwork

    std::vector<int> in_arcs_;	// index of arcs going from entrance to check nodes
    std::vector<int> out_arcs_;	// index of arcs going from check nodes to exit

    void build();
};


#endif //NURSESCHEDULER_PRICELABELSGRAPH_H
