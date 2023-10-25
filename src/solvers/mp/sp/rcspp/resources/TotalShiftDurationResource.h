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
  TotalShiftDuration(const PAbstractShift &pAShift, int maxDuration,
                     bool countAssignment) :
      pAShift__(pAShift), maxDuration_(maxDuration),
      countAssignment_(countAssignment) {
    if (maxDuration == 0)
      Tools::throwError("Cannot have a max duration of 0. "
                        "The constraint becomes meaningless.");
  }

  int computeConsumption(const Stretch &stretch, bool *ready = nullptr) const;

  int maxDuration() const { return maxDuration_; }

  int duration(const PShift &pS) const {
    return countAssignment_ ? 1 : pS->duration;
  }

  bool countAssignment() const {
    return countAssignment_;
  }

 protected:
  int maxDuration_;  // maximum duration of a shift: used for the dominance
  bool countAssignment_;  // if true, count assignment

 private:
  PAbstractShift pAShift__;
};

class SoftTotalShiftDurationResource :
    public SoftBoundedResource, public TotalShiftDuration {
 public:
  SoftTotalShiftDurationResource(int lb, int ub, double lbCost, double ubCost,
                                 const PAbstractShift &pShift,
                                 int totalNbDays,
                                 bool countAssignment = true,
                                 int maxDuration = 1) :
      SoftBoundedResource(pShift, "Soft Total " + pShift->name,
                          lb, ub, lbCost, ubCost),
      TotalShiftDuration(pShift, maxDuration, countAssignment) {
    totalNbDays_ = totalNbDays;
    costType_ = TOTAL_WORK_COST;
  }

  BaseResource *clone() const override {
    return new SoftTotalShiftDurationResource(
        lb_, ub_, lbCost_, ubCost_, pAShift_,
        totalNbDays_, countAssignment_, maxDuration_);
  }

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc &pA, double *cost) override;

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
      int lb, int ub,
      const PAbstractShift &pShift,
      int totalNbDays,
      bool countAssignment = true,
      int maxDuration = 1) :
      HardBoundedResource(pShift, "Hard Total " + pShift->name, lb, ub),
      TotalShiftDuration(pShift, maxDuration, countAssignment) {
    totalNbDays_ = totalNbDays;
    maxDurationForRemainingDays(totalNbDays, 1);
  }

  BaseResource *clone() const override {
    auto pR = new HardTotalShiftDurationResource(
        lb_, ub_, pAShift_,
        totalNbDays_, countAssignment_, maxDuration_);
    pR->maxDurationForRemainingDays_ = maxDurationForRemainingDays_;
    return pR;
  }

  int getConsumption(const State &initialState) const override;

  int maxConsumptionPerDay() const override { return maxDuration_; }

  bool isAnyWorkShiftResource() const override { return pAShift_->isAnyWork(); }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

  DominationStatus dominates(
      RCLabel *pL1, RCLabel *pL2, double *cost) const override;

  bool isActive(int dssrLvl) const override {
    return dssrLvl == 0 || pAShift_->isAnyWork() || pAShift_->isRest() ||
        ub_ > dssrLvl;
  }

  void maxDurationForRemainingDays(int maxConsWorkedDays, int minConsRestDays);

 protected:
  // maximum number of worked days per remaining days
  vector<int> maxDurationForRemainingDays_;
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/**
 * Structure storing the information that is necessary to expand the labels
 * of resources on the count of the occurrences of a given abstract shift
 */
struct SoftTotalShiftDurationExpander : public Expander {
  SoftTotalShiftDurationExpander(int indResource,
                                 const SoftTotalShiftDurationResource &resource,
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
  const SoftTotalShiftDurationResource &resource_;
  int consumption_;
  int maxDurationBefore_;  // total duration before the start of the stretch
  int maxDurationLeft_;  // total duration left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

struct HardTotalShiftDurationExpander : public Expander {
  HardTotalShiftDurationExpander(int indResource,
                                 const HardTotalShiftDurationResource &resource,
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
  const HardTotalShiftDurationResource &resource_;
  int consumption_;
  int maxDurationBefore_;  // total duration before the start of the stretch
  int maxDurationLeft_;  // total duration left after the end of the stretch
  bool arcToSink_;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALSHIFTDURATIONRESOURCE_H_
