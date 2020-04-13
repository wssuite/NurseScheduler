//
// Created by antoine legrain on 2020-04-11.
//

#ifndef NURSESCHEDULER_PRINCIPALGRAPH_H
#define NURSESCHEDULER_PRINCIPALGRAPH_H

#include "tools/MyTools.h"
#include "RCGraph.h"

class SubProblem;

class PrincipalGraph {
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

    void build();
};


#endif //NURSESCHEDULER_PRINCIPALGRAPH_H
