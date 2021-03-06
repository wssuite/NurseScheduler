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
                        const PAbstractShift pShift, int totalNbDays,
                        std::string _name = "") :
      SoftBoundedResource(_name.empty() ?
                          "Soft Cons "+pShift->name : std::move(_name),
                          lb, ub, lbCost, ubCost),
      pShift_(pShift),
      totalNbDays_(totalNbDays) {}

  int getConsumption(const State &initialState) const override;

  int getTotalNbDays() const {return totalNbDays_;}

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const AbstractShift &aShift) const override {
    return pShift_->includes(aShift);
  }

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  // the worst case is when compared with a label of consumption 0
  // that will reach ub_: so we will pay consumption penalties.
  // As soon as a consumption exceeds the ub, each extra shift is immediately
  // penalized, that's why we take the min of ub_ and consumption.
  // Finally, if there are not enough days left to reach ub_ for a label of
  // consumption 0, the max that it can reach is nLeft and thus the other will
  // exceed of only std::max(ub_ - consumption - nLeft, 0) which is the default
  // worst case.
  double getWorstUbCost(int consumption, int nLeft) const override {
    return std::min(ubCost_ * std::min(consumption, ub_),
                    SoftBoundedResource::getWorstUbCost(consumption, nLeft));
  }

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pShift_;
  int totalNbDays_;  // Total number of days in the horizon
};

class HardConsShiftResource : public HardBoundedResource {
 public:
  HardConsShiftResource(
      int lb, int ub, const PAbstractShift pShift, std::string _name = "") :
      HardBoundedResource(_name.empty() ?
                          "Hard Cons "+pShift->name : std::move(_name),
                          lb,
                          ub),
      pShift_(pShift) {}

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pShift_;
};

/*
 * Expanders for the  resources
 */

struct ConsShiftExpander : public Expander {
  ConsShiftExpander(int rId, bool start, bool reset,
                    int consBeforeReset, int consAfterReset,
                    bool arcToSink, int nDaysLeft = 0) :
      Expander(rId), arcToSink(arcToSink),
      start(start),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      nDaysLeft(nDaysLeft) {}

 protected:
  bool arcToSink;  // true if the target of the arc is a sink node
  bool start;  // true if the arc starts the consumption of the resource
  bool reset;  // true if resource is reset on the arc
  int consBeforeReset;  // resource consumption before resetting the arc
  int consAfterReset;  // resource consumption after resetting the arc
  int nDaysLeft;  // total number of days left after the end of the stretch
};

struct SoftConsShiftExpander : public ConsShiftExpander {
  SoftConsShiftExpander(const SoftConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        bool arcToSink,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
                        start,
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        arcToSink,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

 private:
  const SoftConsShiftResource& resource_;
};


struct HardConsShiftExpander : public ConsShiftExpander {
  HardConsShiftExpander(const HardConsShiftResource& resource,
                        bool start,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
                        start,
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;

  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

 private:
  const HardConsShiftResource& resource_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSSHIFTRESOURCE_H_