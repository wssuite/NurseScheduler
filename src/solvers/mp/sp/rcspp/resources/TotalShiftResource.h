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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTRESOURCE_H_

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
 * Resource corresponding to the soft min/max constraints on the total
 * number of occurrences of a given abstract shift
 */
class SoftTotalShiftResource : public SoftBoundedResource {
 public:
  SoftTotalShiftResource(int lb, int ub, double lbCost, double ubCost,
                         const PAbstractShift shift, int totalNbDays = 0) :
      SoftBoundedResource(lb, ub, lbCost, ubCost), shift_(shift),
      totalNbDays_(totalNbDays) {}

  int getConsumption(const State &initialState) const override;

 protected:
  // initialize the expander on a given arc
  PExpander init(const Shift &prevShift,
                 const Stretch &stretch,
                 const RCArc &arc) override;

 private:
  const PAbstractShift shift_;
  int totalNbDays_;  // Total number of days in the horizon
};


/**
 * Structure storing the information that is necessary to expand the labels
 * of resources on the count of the occurrences of a given abstract shift
 */
struct SoftTotalShiftExpander : public Expander {
  explicit SoftTotalShiftExpander(const SoftTotalShiftResource& resource,
                                  int consumption) :
      Expander(resource.id()),
      resource_(resource),
      consumption(consumption),
      nDaysLeft(0),
      arcToSink(false) {}

  SoftTotalShiftExpander(
      const SoftTotalShiftResource& resource,
      int consumption,
      int nDaysLeft) :
      Expander(resource.id()),
      resource_(resource),
      consumption(consumption),
      nDaysLeft(nDaysLeft),
      arcToSink(false) {}

  SoftTotalShiftExpander(const SoftTotalShiftResource& resource,
                         int consumption,
                         int nDaysLeft,
                         bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      consumption(consumption),
      nDaysLeft(nDaysLeft),
      arcToSink(arcToSink) {}

  bool expand(const ResourceValues &vParent,
              const PRCLabel &pLChild,
              ResourceValues *vChild) override;

 protected:
  const SoftTotalShiftResource& resource_;
  int consumption;
  int nDaysLeft;  // total number of days left after the end of the stretch
  bool arcToSink;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTRESOURCE_H_
