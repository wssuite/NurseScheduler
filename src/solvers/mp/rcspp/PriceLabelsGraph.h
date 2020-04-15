//
// Created by antoine legrain on 2020-04-11.
//

#ifndef NURSESCHEDULER_PRICELABELSGRAPH_H
#define NURSESCHEDULER_PRICELABELSGRAPH_H


#include <vector>

class SubProblem;

class PriceLabelsGraph {
  public:
    PriceLabelsGraph(int min_ub, int max_ub, SubProblem* sp = nullptr, bool compute_min_cost=true);
    virtual ~PriceLabelsGraph();

    void updateArcCosts();

    inline int entrance() const { return entrance_; }

    inline int exit() const { return exit_; }

  protected:
    SubProblem* pSP_;

    int min_ub_, max_ub_;
    bool compute_min_cost_;

    int entrance_;				// entrance node to the subnetwork
    std::vector<int> checkNodes_;			// check nodes (measure the level of labels)
    int exit_;  // exit node from the subnetwork

    std::vector<int> in_arcs_;	// index of arcs going from entrance to check nodes
    std::vector<int> out_arcs_;	// index of arcs going from check nodes to exit

    void build();
};


#endif //NURSESCHEDULER_PRICELABELSGRAPH_H
