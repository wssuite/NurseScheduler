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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSSHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSSHIFTRESOURCE_H_

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

/**
 * Resource corresponding to the soft min/max constraints on the number of
 * consecutive shifts: it may be used for any kind of abstract shift
 */
class ConsShift {
 public:
  ConsShift(const PAbstractShift &pShift,
            int initialConsumption = 0,
            bool endOnLastDay = false,
            bool cyclic = false,
            bool enforceLBAtStart = true):
            pAShift_(pShift),
            initialConsumption_(initialConsumption),
            lastDayEndsSequence_(endOnLastDay),
            cyclic_(cyclic),
            enforceLBAtStart_(enforceLBAtStart) {
    if (cyclic_ && lastDayEndsSequence_) {
      std::cerr << "ConsShift has disabled lastDayEndsSequence_"
                   " as cyclic_ is enable." << std::endl;
      lastDayEndsSequence_ = false;
    }
    if (cyclic_ && !enforceLBAtStart_) {
      std::cerr << "ConsShift has enabled enforceLBAtStart_"
                   " as cyclic_ is enable." << std::endl;
      enforceLBAtStart_ = true;
    }
  }

  const PAbstractShift &pAShift() const { return pAShift_; }
  bool lastDayEndsSequence() const { return lastDayEndsSequence_; }
  bool isCyclic() const { return cyclic_; }
  bool enforceLBAtStart() const { return enforceLBAtStart_; }
  int initialConsumption() const { return initialConsumption_; }

 protected:
  const PAbstractShift pAShift_;
  // consumption of the resource in the initial state
  int initialConsumption_;

  // true if the last day of horizon puts an end to a sequence of shifts (and
  // thus we need to pay LB cost if it LB is not met)
  bool lastDayEndsSequence_;

  bool cyclic_;  // true of solving the cyclic version

  // true if LB violations count on the first sequence starting on day 0
  bool enforceLBAtStart_;
};


class SoftConsShiftResource : public SoftBoundedResource, public ConsShift {
 public:
  SoftConsShiftResource(int lb,
                        int ub,
                        double lbCost,
                        double ubCost,
                        const PAbstractShift &pShift,
                        CostType costType,
                        int totalNbDays,
                        int initialConsumption,
                        bool endOnLastDay = false,
                        bool cyclic = false,
                        bool enforceLBAtStart = true,
                        std::string _name = "") :
      SoftBoundedResource(_name.empty() ?
                          "Soft Cons "+pShift->name : std::move(_name),
                          lb, ub, lbCost, ubCost),
      ConsShift(pShift, initialConsumption, endOnLastDay,
                cyclic, enforceLBAtStart) {
    totalNbDays_ = totalNbDays;
    costType_ = costType;
  }

  BaseResource* clone() const override {
    return new SoftConsShiftResource(
        lb_, ub_, lbCost_, ubCost_, pAShift_, costType_, totalNbDays_,
        initialConsumption_, lastDayEndsSequence_, cyclic_,
        enforceLBAtStart_, name);
  }

  int getConsumption(const State &initialState) const override;

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const PAbstractShift &pAShift) const override {
    return pAShift_->includes(*pAShift);
  }

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  // the worst case is when compared with a label of consumption 0
  // that will reach ub_: so we will pay consumption penalties.
  // As soon as a consumption exceeds the ub, each extra shift is immediately
  // penalized, that's why we take the min of ub_ and consumption.
  // Finally, if there are not enough days left to reach ub_ for a label of
  // consumption 0, the max that it can reach is nLeft and thus the other will
  // exceed of only std::max(ub_ - consumption - nLeft, 0) which is the default
  // worst case.
  double getWorstUbCost(int consumption) const override {
    return SoftBoundedResource::getWorstUbCost(consumption);
  }
  double getWorstUbCost(int consumption, int nLeft) const override {
    return std::min(getWorstUbCost(consumption),
                    SoftBoundedResource::getWorstUbCost(consumption, nLeft));
  }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return pAShift_->isRest(); };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  void enumerate(const PRCGraph &pRCGraph, bool forceEnum) override;
};

class HardConsShiftResource : public HardBoundedResource, public ConsShift {
 public:
  HardConsShiftResource(
      int lb, int ub, const PAbstractShift &pShift,
      int totalNbDays, int initialConsumption,
      bool endOnLastDay = false, bool cyclic = false,
      bool enforceLBAtStart = true, std::string _name = "") :
      HardBoundedResource(_name.empty() ?
                          "Hard Cons " + pShift->name : std::move(_name),
                          lb, ub),
      ConsShift(pShift, initialConsumption,
                endOnLastDay, cyclic, enforceLBAtStart) {
    totalNbDays_ = totalNbDays;
  }

  BaseResource *clone() const override {
    return new HardConsShiftResource(
        lb_, ub_, pAShift_,
        totalNbDays_, initialConsumption_,
        lastDayEndsSequence_, cyclic_, enforceLBAtStart_, name);
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  DominationStatus dominates(
      RCLabel *pL1, RCLabel *pL2, double *cost) const override;

  bool isActive(int dssrLvl) const override {
    return dssrLvl == 0 || pAShift_->isAnyWork() || pAShift_->isRest() ||
        ub_ > dssrLvl;
  }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return pAShift_->isRest(); };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/*
 * Expanders for the  resources
 */

struct ConsShiftExpander : public Expander {
  ConsShiftExpander(int rId, CostType costType,
                    bool start, bool reset,
                    int consBeforeReset, int consAfterReset,
                    bool arcToSink, int nDaysBefore = 0, int nDaysLeft = 0) :
      Expander(rId, costType),
      arcToSink(arcToSink),
      start(start),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      nDaysBefore(nDaysBefore),
      nDaysLeft(nDaysLeft) {
  }

 protected:
  bool arcToSink;  // true if the target of the arc is a sink node
  bool start;  // true if the arc starts the consumption of the resource
  bool reset;  // true if resource is reset on the arc
  int consBeforeReset;  // resource consumption before resetting the arc
  int consAfterReset;  // resource consumption after resetting the arc
  int nDaysBefore;  // number of days before the start of the stretch
  int nDaysLeft;  // total number of days left after the end of the stretch
};

struct SoftConsShiftExpander : public ConsShiftExpander {
  SoftConsShiftExpander(int indResource,
                        const SoftConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        bool arcToSink,
                        int nDaysBefore = 0,
                        int nDaysLeft = 0):
      ConsShiftExpander(indResource,
                        resource.costType(),
                        start,
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        arcToSink,
                        nDaysBefore,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const SoftConsShiftResource& resource_;
};


struct HardConsShiftExpander : public ConsShiftExpander {
  HardConsShiftExpander(int indResource,
                        const HardConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        bool arcToSink,
                        int nDaysBefore = 0,
                        int nDaysLeft = 0):
      ConsShiftExpander(indResource,
                        NO_COST,
                        start,
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        arcToSink,
                        nDaysBefore,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;

  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const HardConsShiftResource& resource_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSSHIFTRESOURCE_H_
