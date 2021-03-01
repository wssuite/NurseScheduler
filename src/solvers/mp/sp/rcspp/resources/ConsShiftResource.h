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

#include <utility>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

#include "solvers/mp/sp/rcspp/MyRCLabel.h"
#include "solvers/mp/sp/rcspp/MyRCGraph.h"

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
                        const PAbstractShift shift, int totalNbDays = 0) :
      SoftBoundedResource(lb, ub, lbCost, ubCost),
      pShift_(shift),
      totalNbDays_(totalNbDays) {}

  int getConsumption(const State &initialState) const override;

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const Shift &shift) const override {
    return pShift_->includes(shift);
  }

  double getWorstUbCost(int consumption, int nLeft = 0) const override {
    return ubCost_ * std::min(consumption, ub_);
  }

 protected:
  PExpander init(const Shift &prevShift,
                 const Stretch &stretch,
                 const RCArc &arc) override;

 private:
  const PAbstractShift pShift_;
  int totalNbDays_;  // Total number of days in the horizon
};

class HardConsShiftResource : public HardBoundedResource {
 public:
  HardConsShiftResource(int lb, int ub, const PAbstractShift pShift) :
      HardBoundedResource(lb, ub), pShift_(pShift) {}

  int getConsumption(const State &initialState) const override;

 protected:
  // initialize the expander on a given arc
  PExpander init(const Shift &prevShift,
                 const Stretch &stretch,
                 const RCArc &arc) override;

 private:
  const PAbstractShift pShift_;
};

/*
 * Expanders for the  resources
 */

struct ConsShiftExpander : public Expander {
  explicit ConsShiftExpander(int rId)
      : Expander(rId), arcToSink(false),
        reset(false), consBeforeReset(0), consAfterReset(0), cost(0),
        nDaysLeft(0) {}

  ConsShiftExpander(int rId, bool reset,
                    int consBeforeReset, int consAfterReset, double cost,
                    int nDaysLeft = 0) :
      Expander(rId), arcToSink(false),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      cost(cost),
      nDaysLeft(nDaysLeft) {}

  ConsShiftExpander(int rId, bool reset,
                    int consBeforeReset, int consAfterReset, double cost,
                    bool arcToSink, int nDaysLeft = 0) :
      Expander(rId), arcToSink(arcToSink),
      reset(reset),
      consBeforeReset(consBeforeReset),
      consAfterReset(consAfterReset),
      cost(cost),
      nDaysLeft(nDaysLeft) {}

 protected:
  bool arcToSink;  // true if the target of the arc is a sink node
  bool reset;  // true if resource is reset on the arc
  int consBeforeReset;  // resource consumption before resetting the arc
  int consAfterReset;  // resource consumption after resetting the arc
  double cost;  // cost due to soft resource costs between first and last reset
  int nDaysLeft;  // total number of days left after the end of the stretch
};

struct SoftConsShiftExpander : public ConsShiftExpander {
  explicit SoftConsShiftExpander(const SoftConsShiftResource& resource):
      ConsShiftExpander(resource.id()), resource_(resource) {}

  SoftConsShiftExpander(const SoftConsShiftResource& resource,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        double cost,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        cost,
                        nDaysLeft),
      resource_(resource) {}

  SoftConsShiftExpander(const SoftConsShiftResource& resource,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        double cost,
                        bool arcToSink,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        cost,
                        arcToSink,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const ResourceValues &vParent,
              const PRCLabel &pLChild,
              ResourceValues *vChild) override;

 private:
  const SoftConsShiftResource& resource_;
};


struct HardConsShiftExpander : public ConsShiftExpander {
  explicit HardConsShiftExpander(const HardConsShiftResource& resource):
      ConsShiftExpander(resource.id()), resource_(resource) {}

  HardConsShiftExpander(const HardConsShiftResource& resource,
                        bool reset,
                        int consBeforeReset,
                        int consAfterReset,
                        double cost,
                        int nDaysLeft = 0):
      ConsShiftExpander(resource.id(),
                        reset,
                        consBeforeReset,
                        consAfterReset,
                        cost,
                        nDaysLeft),
      resource_(resource) {}

  bool expand(const ResourceValues &vParent,
              const PRCLabel &pLChild,
              ResourceValues *vChild) override;

 private:
  const HardConsShiftResource& resource_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSSHIFTRESOURCE_H_
