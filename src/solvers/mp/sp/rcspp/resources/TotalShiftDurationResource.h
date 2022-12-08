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
class TotalShiftDuration {
 public:
  TotalShiftDuration(PAbstractShift  pShift,
                     int maxDuration,
                     int defaultDuration) :
      pAShift_(std::move(pShift)),
      maxDuration_(maxDuration),
      defaultDuration_(defaultDuration) {}

  int computeConsumption(const Stretch &stretch, bool *ready = nullptr) const;

  const PAbstractShift &pShift() const { return pAShift_; }

  int maxDuration() const { return maxDuration_; }

  int duration(const PShift &pS) const {
    return defaultDuration_ >= 0 ? defaultDuration_ : pS->duration;
  }

 protected:
  const PAbstractShift pAShift_;
  int maxDuration_;  // maximum duration of a shift: used for the dominance
  int defaultDuration_;  // if defined, override duration of the PShift
};

class SoftTotalShiftDurationResource :
    public SoftBoundedResource, public TotalShiftDuration {
 public:
  SoftTotalShiftDurationResource(int lb, int ub, double lbCost, double ubCost,
                                 const PAbstractShift& pShift,
                                 int totalNbDays,
                                 int maxDuration,
                                 int defaultDuration = -1) :
      SoftBoundedResource("Soft Total "+pShift->name,
                          lb, ub, lbCost, ubCost),
      TotalShiftDuration(pShift, maxDuration, defaultDuration) {
    totalNbDays_ = totalNbDays;
    costType_ = TOTAL_WORK_COST;
  }

  BaseResource* clone() const override {
    return new SoftTotalShiftDurationResource(
        lb_, ub_, lbCost_, ubCost_, pAShift_,
        totalNbDays_, maxDuration_, defaultDuration_);
  }

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;

  int getConsumption(const State &initialState) const override;

  int maxConsumptionPerDay() const override { return maxDuration_; }

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

class HardTotalShiftDurationResource :
    public HardBoundedResource, public TotalShiftDuration {
 public:
  HardTotalShiftDurationResource(
      int lb, int ub, const PAbstractShift& pShift,
      int totalNbDays, int maxDuration, int defaultDuration = -1) :
      HardBoundedResource("Hard Total "+pShift->name, lb, ub),
      TotalShiftDuration(pShift, maxDuration, defaultDuration) {
    totalNbDays_ = totalNbDays;
  }

  BaseResource* clone() const override {
    return new HardTotalShiftDurationResource(
        lb_, ub_, pAShift_,
        totalNbDays_, maxDuration_, defaultDuration_);
  }

  int getConsumption(const State &initialState) const override;

  int maxConsumptionPerDay() const override { return maxDuration_; }

  bool isAnyWorkShiftResource() const override { return pShift_->isAnyWork(); }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

 private:
  const PAbstractShift pShift_;
};


/**
 * Structure storing the information that is necessary to expand the labels
 * of resources on the count of the occurrences of a given abstract shift
 */
struct SoftTotalShiftDurationExpander : public Expander {
  SoftTotalShiftDurationExpander(int indResource,
                                 const SoftTotalShiftDurationResource& resource,
                                 int consumption,
                                 int maxDurationBefore,
                                 int maxDurationLeft,
                                 bool arcToSink) :
      Expander(indResource, TOTAL_WORK_COST),
      resource_(resource),
      consumption_(consumption),
      maxDurationBefore_(maxDurationBefore),
      maxDurationLeft_(maxDurationLeft),
      arcToSink_(arcToSink) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const SoftTotalShiftDurationResource& resource_;
  int consumption_;
  int maxDurationBefore_;  // total duration before the start of the stretch
  int maxDurationLeft_;  // total duration left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

struct HardTotalShiftDurationExpander : public Expander {
  HardTotalShiftDurationExpander(int indResource,
                                 const HardTotalShiftDurationResource& resource,
                                 int consumption,
                                 int maxDurationBefore,
                                 int maxDurationLeft,
                                 bool arcToSink) :
      Expander(indResource, NO_COST),
      resource_(resource),
      consumption_(consumption),
      maxDurationBefore_(maxDurationBefore),
      maxDurationLeft_(maxDurationLeft),
      arcToSink_(arcToSink) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const HardTotalShiftDurationResource& resource_;
  int consumption_;
  int maxDurationBefore_;  // total duration before the start of the stretch
  int maxDurationLeft_;  // total duration left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTDURATIONRESOURCE_H_
