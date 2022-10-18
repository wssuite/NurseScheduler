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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALWEEKENDSRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALWEEKENDSRESOURCE_H_

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
 * number of worked weekends
 * Plain expanders are used, because precomputation is more complex than
 * useful for total weekend resources
 */
class TotalWeekend {
 public:
  TotalWeekend(PAbstractShift  pShift,
               DayOfWeek first = SATURDAY, DayOfWeek last = SUNDAY) :
      pAShift_(std::move(pShift)),
      weekend_(first, last),
      totalNbWeekends_(0) {}

  virtual int computeConsumption(const Stretch &stretch,
                                 bool *ready) const;

  int computeBackConsumption(const Stretch &stretch,
                             bool *ready) const;

  const PAbstractShift& pShift() const { return pAShift_; }

  const Weekend& weekend() const { return weekend_; }

 protected:
  const PAbstractShift pAShift_;
  const Weekend weekend_;
//  BoundedResource *pR_;
  int totalNbWeekends_;
};

class SoftTotalWeekendsResource :
    public SoftBoundedResource, public TotalWeekend {
 public:
  SoftTotalWeekendsResource(
      int ub, double ubCost,
      const PAbstractShift &pShift,
      int totalNbDays,
      DayOfWeek firstWeekendDay = SATURDAY,
      DayOfWeek lastWeekendDay = SUNDAY) :
      SoftBoundedResource("Soft Weekend Total work", 0, ub, 0, ubCost),
      TotalWeekend(pShift, firstWeekendDay, lastWeekendDay) {
    totalNbDays_ = totalNbDays;
    totalNbWeekends_ = weekend_.nWeekendsInInterval(
        firstDayId(), firstDayId() + totalNbDays_ - 1);
    costType_ = TOTAL_WEEKEND_COST;
  }

  BaseResource* clone() const override {
    return new SoftTotalWeekendsResource(
        ub_, ubCost_, pAShift_, totalNbDays_,
        weekend_.firstWeekendDay().getDayOfWeek(),
        weekend_.lastWeekendDay().getDayOfWeek());
  }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  // instantiate TotalWeekEndLabel
  int getConsumption(const State &initialState) const override;


  double getWorstLbCost(int consumption) const override { return .0; }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/**
 * Expander that allows to count the consumption of the total weekend resources
 */
struct SoftTotalWeekendsExpander : public Expander {
  SoftTotalWeekendsExpander(int indResource,
                            const SoftTotalWeekendsResource& resource,
                            const Stretch& stretch,
                            int consumption,
                            int nWeekendsBefore,
                            int nWeekendsAfter,
                            bool arcToSink) :
      Expander(indResource, TOTAL_WEEKEND_COST),
      resource_(resource),
      consumption_(consumption),
      nWeekendsBefore_(nWeekendsBefore),
      nWeekendsAfter_(nWeekendsAfter),
      stretch_(stretch),
      arcToSink_(arcToSink) {}

  // expand a given resource label
  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const SoftTotalWeekendsResource& resource_;
  int consumption_;
  bool startOnWeekend_{}, endOnWeekend_{};
  // number of weekends before and including the first day of the stretch;
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsBefore_;
  // number of weekends after the end of the stretch,
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsAfter_;
  const Stretch stretch_;
  bool arcToSink_;  // true if the target of the arc is a sink node
};



class HardTotalWeekendsResource :
    public HardBoundedResource, public TotalWeekend {
 public:
  HardTotalWeekendsResource(
      int ub,
      const PAbstractShift &pShift,
      int totalNbDays,
      DayOfWeek firstWeekendDay = SATURDAY,
      DayOfWeek lastWeekendDay = SUNDAY) :
      HardBoundedResource("Hard Weekend Total work", 0, ub),
      TotalWeekend(pShift, firstWeekendDay, lastWeekendDay) {
    totalNbDays_ = totalNbDays;
  }

  BaseResource* clone() const override {
    return new HardTotalWeekendsResource(
        ub_, pAShift_, totalNbDays_,
        weekend_.firstWeekendDay().getDayOfWeek(),
        weekend_.lastWeekendDay().getDayOfWeek());
  }

  // instantiate TotalWeekEndLabel
  int getConsumption(const State &initialState) const override;

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/**
 * Expander that allows to count the consumption of the total weekend resources
 */
struct HardTotalWeekendsExpander : public Expander {
  HardTotalWeekendsExpander(int indResource,
                            const HardTotalWeekendsResource& resource,
                            const Stretch& stretch,
                            int consumption,
                            int nWeekendsBefore,
                            int nWeekendsAfter,
                            bool arcToSink) :
      Expander(indResource, NO_COST),
      consumption_(consumption),
      resource_(resource),
      nWeekendsBefore_(nWeekendsBefore),
      nWeekendsAfter_(nWeekendsAfter),
      stretch_(stretch),
      arcToSink_(arcToSink) {}

  // expand a given resource label
  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 protected:
  const HardTotalWeekendsResource& resource_;
  int consumption_;
  // number of weekends before and including the first day of the stretch;
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsBefore_;
  // number of weekends after the end of the stretch,
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsAfter_;
  const Stretch stretch_;
  bool arcToSink_;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALWEEKENDSRESOURCE_H_
