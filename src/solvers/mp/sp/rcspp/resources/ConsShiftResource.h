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
class SoftConsShiftResource : public SoftBoundedResource {
 public:
  SoftConsShiftResource(int lb, int ub, double lbCost, double ubCost,
                        const PAbstractShift& pShift, int totalNbDays,
                        int initialConsumption,
                        std::string _name = "") :
      SoftBoundedResource(_name.empty() ?
                          "Soft Cons "+pShift->name : std::move(_name),
                          lb, ub, lbCost, ubCost),
      pAShift_(pShift),
      initialConsumption_(initialConsumption) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const AbstractShift &aShift) const override {
    return pAShift_->includes(aShift);
  }

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  const PAbstractShift &pShift() const { return pAShift_; }

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

  const PAbstractShift &pAShift() const { return pAShift_; }

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pAShift_;
  int initialConsumption_ = 0;  // consumption of the resource in the initial
  // state
};

typedef shared_ptr<SoftConsShiftResource> PSoftConsShiftResource;

class HardConsShiftResource : public HardBoundedResource {
 public:
  HardConsShiftResource(
      int lb, int ub, const PAbstractShift& pShift, int totalNbDays) :
      HardBoundedResource("Hard Cons "+pShift->name, lb, ub),
      pAShift_(pShift) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  const PAbstractShift &pShift() const { return pAShift_; }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  const PAbstractShift &pAShift() const { return pAShift_; }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pAShift_;
};

typedef shared_ptr<HardConsShiftResource> PHardConsShiftResource;

/*
 * Expanders for the  resources
 */

struct ConsShiftExpander : public Expander {
  ConsShiftExpander(int rId, bool start, bool reset,
                    int consBeforeReset, int consAfterReset,
                    bool arcToSink, int nDaysBefore = 0, int nDaysLeft = 0) :
      Expander(rId), arcToSink(arcToSink),
      start(start),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      nDaysBefore(nDaysBefore),
      nDaysLeft(nDaysLeft) {}

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
  SoftConsShiftExpander(const SoftConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        bool arcToSink,
                        int nDaysBefore = 0,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
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
  HardConsShiftExpander(const HardConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        bool arcToSink,
                        int nDaysBefore = 0,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
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
