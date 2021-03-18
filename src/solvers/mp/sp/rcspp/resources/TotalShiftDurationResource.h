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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTDURATIONRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTDURATIONRESOURCE_H_

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
 * Resource corresponding to the soft min/max constraints on the total
 * number of occurrences of a given abstract shift
 */
class SoftTotalShiftDurationResource : public SoftBoundedResource {
 public:
  SoftTotalShiftDurationResource(int lb, int ub, double lbCost, double ubCost,
                                 const PAbstractShift pShift,
                                 int totalNbDays, int defaultDuration = -1) :
      SoftBoundedResource("Soft Total "+pShift->name,
                          lb, ub, lbCost, ubCost),
      pShift_(pShift),
      defaultDuration_(defaultDuration) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  const PAbstractShift pShift() const { return pShift_; }

  bool hasDefaultDuration() const { return defaultDuration_ != -1; }
  int defaultDuration() const { return defaultDuration_; }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

 private:
  const PAbstractShift pShift_;
  int defaultDuration_;  // if defined, override duration of the PShift
};

class HardTotalShiftDurationResource : public HardBoundedResource {
 public:
  HardTotalShiftDurationResource(
      int lb, int ub, const PAbstractShift pShift,
      int totalNbDays, int defaultDuration = -1) :
      HardBoundedResource("Hard Total "+pShift->name, lb, ub),
      pShift_(pShift),
      defaultDuration_(defaultDuration) {
    totalNbDays_ = totalNbDays;
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  const PAbstractShift pShift() const { return pShift_; }

  bool hasDefaultDuration() const { return defaultDuration_ != -1; }
  int defaultDuration() const { return defaultDuration_; }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

 private:
  const PAbstractShift pShift_;
  int defaultDuration_;  // if defined, override duration of the PShift
};


/**
 * Structure storing the information that is necessary to expand the labels
 * of resources on the count of the occurrences of a given abstract shift
 */
struct SoftTotalShiftDurationExpander : public Expander {
  SoftTotalShiftDurationExpander(const SoftTotalShiftDurationResource& resource,
                                 int consumption,
                                 int nDaysBefore,
                                 int nDaysLeft,
                                 bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      consumption(consumption),
      nDaysBefore(nDaysBefore),
      nDaysLeft(nDaysLeft),
      arcToSink_(arcToSink) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const SoftTotalShiftDurationResource& resource_;
  int consumption;
  int nDaysBefore;  // number of days before the start of the stretch
  int nDaysLeft;  // total number of days left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

struct HardTotalShiftDurationExpander : public Expander {
  HardTotalShiftDurationExpander(const HardTotalShiftDurationResource& resource,
                                 int consumption,
                                 int nDaysBefore,
                                 int nDaysLeft,
                                 bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      consumption(consumption),
      nDaysBefore(nDaysBefore),
      nDaysLeft(nDaysLeft),
      arcToSink_(arcToSink) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const HardTotalShiftDurationResource& resource_;
  int consumption;
  int nDaysBefore;  // number of days before the start of the stretch
  int nDaysLeft;  // total number of days left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTDURATIONRESOURCE_H_
