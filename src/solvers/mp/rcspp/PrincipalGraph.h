//
// Created by antoine legrain on 2020-04-11.
//

#ifndef NURSESCHEDULER_PRINCIPALGRAPH_H
#define NURSESCHEDULER_PRINCIPALGRAPH_H

#include "tools/MyTools.h"
#include "RCGraph.h"

class SubProblem;

class PrincipalGraph: public SubGraph {
  public:
    PrincipalGraph(int shift_type,  SubProblem* sp = nullptr);
    virtual ~PrincipalGraph();

    // check if feasible to link this arc at this level
    bool checkFeasibilityEntranceArc(const Arc_Properties& arc_prop, int level) const;

    void updateArcCosts();

    void forbidDayShift(int k, int s);
    void authorizeDayShift(int k, int s);
    bool checkIfShiftBelongsHere(int s, bool print_err =  false) const;

    inline int shiftType() const { return shift_type_; }

    inline int maxCons() const { return max_cons_; }

    inline int getNode(int k, int n) const { return principalNetworkNodes_[k][n]; }

    inline const std::vector<int>& getDayNodes(int k) const {
      if(pSP_) return principalNetworkNodes_[k];
      return Tools::EMPTY_INT_VECTOR;
    }

    // link two sub graphs together:
    // 1. create an arc from the exit of inSubGraph to the current entrance
    // 2. the entrance method should now returns the entrance of the  inSubGraph
    void linkInSubGraph(SubGraph& inSubGraphs, int day=-1) override;

    // link two sub graphs together:
    // 1. create an arc from the current exit to the entrance of outSubGraph
    // 2. the exit method should now returns the exit of the outSubGraph
    void linkOutSubGraph(SubGraph& outSubGraph, int day=-1) override;


    inline int entrance(int day=-1) const override {
      if(inSubGraphs_[day]) return inSubGraphs_[day]->entrance(day);
      return getNode(day, 0);
    }

    inline int exit(int day=-1) const override {
      if(outSubGraphs_[day]) return inSubGraphs_[day]->entrance(day);
      return getNode(day, max_cons_);
    }

  protected:
    SubProblem* pSP_;

    int shift_type_;
    int max_cons_;
    std::map<int,int> shifts_to_indices_;

    // Nodes of the PRINCIPAL_NETWORK subnetwork
    vector2D<int> principalNetworkNodes_;					// For each SHIFT, DAY, and # of CONSECUTIVE, the corresponding node id

    vector3D<int> arcsShiftToSameShift_;		// Index: (day, nCons, shift) of origin
    vector2D<int> arcsShiftToEndsequence_;	// Index: (day, nCons) of origin
    vector2D<int> arcsRepeatShift_;		// Index: (day, shift) of origin

    // in subgraphs
    std::vector<SubGraph*> inSubGraphs_, outSubGraphs_;

    void build();

    // return the right vector of consumption based on the day (if < 0, not performing any shift)
    std::vector<int> getConsumption(int day=-1) const;
};


#endif //NURSESCHEDULER_PRINCIPALGRAPH_H
