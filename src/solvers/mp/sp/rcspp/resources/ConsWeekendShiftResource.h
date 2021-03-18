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

#ifndef NURSESCHEDULER_SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_  // NOLINT
#define NURSESCHEDULER_SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_  // NOLINT


#include <utility>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

/**
 * Resource corresponding to the soft min/max constraints on the number of
 * consecutive shifts on weekends: it may be used for any kind of abstract shift
 */
class SoftConsWeekendShiftResource : public SoftBoundedResource {
 public:
  SoftConsWeekendShiftResource(
      int lb, int ub, double lbCost, double ubCost,
      const PAbstractShift pShift, int totalNbDays, bool cyclic = false) :
      SoftBoundedResource("Soft Weekend Cons "+pShift->name,
                          lb, ub, lbCost, ubCost),
      pShift_(pShift), cyclic_(cyclic) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const AbstractShift &aShift) const override {
    return Tools::isWeekend(dayId) && pShift_->includes(aShift);
  }

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  double getWorstUbCost(int consumption, int nLeft) const override {
    return ubCost_ * std::min(consumption, ub_);
  }

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pShift_;
  bool cyclic_;  // true of solving the cyclic version
};

class HardConsWeekendShiftResource : public HardBoundedResource {
 public:
  HardConsWeekendShiftResource(
      int lb, int ub, const PAbstractShift pShift,
      int totalNbDays, bool cyclic = false) :
      HardBoundedResource("Hard Weekend Cons "+pShift->name, lb, ub),
      pShift_(pShift), cyclic_(cyclic) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost = nullptr) override;

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

  const PAbstractShift pShift_;
  bool cyclic_;  // true of solving the cyclic version
};

/*
 * Expanders for the  resources
 */

struct ConsWeekendShiftExpander : public Expander {
  ConsWeekendShiftExpander(int rId, bool start, bool reset,
                           int consBeforeReset, int consAfterReset,
                           bool cyclic, bool arcToSink, int nWeekendsAfter) :
      Expander(rId),
      arcToSink(arcToSink),
      start(start),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      nWeekendsAfter(nWeekendsAfter),
      cyclic(cyclic) {}

 protected:
  bool arcToSink;  // true if the target of the arc is a sink node
  bool start;  // true if the arc starts the consumption of the resource
  bool reset;  // true if resource is reset on the arc
  int consBeforeReset;  // resource consumption before resetting the arc
  int consAfterReset;  // resource consumption after resetting the arc
  // number of weekends after the end of the stretch,
  // if only a portion of a weekend is left (e.g., only sunday),
  // it counts +1 if not already consumed
  int nWeekendsAfter;
  bool cyclic;  // true of solving the cyclic version
};


struct SoftConsWeekendShiftExpander : public ConsWeekendShiftExpander {
  SoftConsWeekendShiftExpander(const SoftConsWeekendShiftResource& resource,
                               bool start,
                               bool reset,
                               int consBeforeReset,
                               int consAfterReset,
                               bool cyclic,
                               bool arcToSink = false,
                               int nWeekendsAfter = 0):
      ConsWeekendShiftExpander(resource.id(),
                               start,
                               reset,
                               consBeforeReset,
                               consAfterReset,
                               cyclic,
                               arcToSink,
                               nWeekendsAfter),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const SoftConsWeekendShiftResource& resource_;
};

struct HardConsWeekendShiftExpander : public ConsWeekendShiftExpander {
  HardConsWeekendShiftExpander(const HardConsWeekendShiftResource& resource,
                               bool start,
                               bool reset,
                               int consBeforeReset,
                               int consAfterReset,
                               bool cyclic,
                               bool arcToSink = false,
                               int nWeekendsAfter = 0):
      ConsWeekendShiftExpander(resource.id(),
                               start,
                               reset,
                               consBeforeReset,
                               consAfterReset,
                               cyclic,
                               arcToSink,
                               nWeekendsAfter),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;

  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const HardConsWeekendShiftResource& resource_;
};

#endif  // NURSESCHEDULER_SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_  // NOLINT
