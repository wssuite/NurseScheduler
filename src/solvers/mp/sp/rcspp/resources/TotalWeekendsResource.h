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
// TODO(JO): we should generalize the concept of weekend by defining the
//  pattern of days, and adding methods that detect if a day is in the
//  pattern, if it is the first day of the pattern and if it is the last day
//  of the pattern
class SoftTotalWeekendsResource : public SoftBoundedResource {
 public:
  SoftTotalWeekendsResource(int ub, double ubCost, int totalNbDays = 0) :
      SoftBoundedResource("Soft Weekend Total work", 0, ub, 0, ubCost),
      totalNbDays_(totalNbDays) {}

  // instantiate TotalWeekEndLabel
  int getConsumption(const State &initialState) const override;

  int totalNbDays() const { return totalNbDays_; }

  double getWorstLbCost(int consumption) const override {
    return .0;
  }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

 private:
  int totalNbDays_;  // Total number of days in the horizon
};

/**
 * Expander that allows to count the consumption of the total weekend resources
 */
struct SoftTotalWeekendsExpander : public Expander {
  SoftTotalWeekendsExpander(const SoftTotalWeekendsResource& resource,
                            Stretch stretch,
                            int nWeekendsBefore,
                            int nWeekendsAfter,
                            bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      nWeekendsBefore_(nWeekendsBefore),
      nWeekendsAfter_(nWeekendsAfter),
      stretch_(std::move(stretch)),
      arcToSink_(arcToSink) {
  }

  // expand a given resource label
  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;


 protected:
  const SoftTotalWeekendsResource& resource_;
  // number of weekends before and including the first day of the stretch;
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsBefore_;
  // number of weekends after the end of the stretch,
  // if only a portion of a weekend is left (e.g., only sunday), it counts +1
  int nWeekendsAfter_;
  const Stretch stretch_;
  bool arcToSink_;  // true if the target of the arc is a sink node
};



class HardTotalWeekendsResource : public HardBoundedResource {
 public:
  HardTotalWeekendsResource(int lb, int ub, int totalNbDays = 0) :
      HardBoundedResource("Hard Weekend Total work", lb, ub),
      totalNbDays_(totalNbDays) {}

  // instantiate TotalWeekEndLabel
  int getConsumption(const State &initialState) const override;

  int totalNbDays() const { return totalNbDays_; }

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;

 private:
  int totalNbDays_;  // Total number of days in the horizon
};

/**
 * Expander that allows to count the consumption of the total weekend resources
 */
struct HardTotalWeekendsExpander : public Expander {
  HardTotalWeekendsExpander(const HardTotalWeekendsResource& resource,
                            Stretch stretch,
                            int nWeekendsBefore,
                            int nWeekendsAfter,
                            bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      nWeekendsBefore_(nWeekendsBefore),
      nWeekendsAfter_(nWeekendsAfter),
      stretch_(std::move(stretch)),
      arcToSink_(arcToSink) {
  }

  // expand a given resource label
  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;


 protected:
  const HardTotalWeekendsResource& resource_;
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
