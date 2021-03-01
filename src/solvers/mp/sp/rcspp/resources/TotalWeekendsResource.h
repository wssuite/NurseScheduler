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

#include "solvers/mp/sp/rcspp/MyRCLabel.h"
#include "solvers/mp/sp/rcspp/MyRCGraph.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

/**
 * Resource corresponding to the soft min/max constraints on the total
 * number of worked weekends
 * Plain expanders are used, because precomputation is more complex than
 * useful for total weekend resources
 * WARNING: with this resource, we implicitly assume that a week-end is
 * saturday and a sunday; otherwise there might be a lot of errors that will
 * happen with current code.
 */
class SoftTotalWeekendsResource : public SoftBoundedResource {
 public:
  SoftTotalWeekendsResource(int ub, double ubCost, int totalNbDays = 0) :
      SoftBoundedResource(0, ub, 0, ubCost), totalNbDays_(totalNbDays) {}

  // instantiate TotalWeekEndLabel
  int getConsumption(const State &initialState) const override;

  // returns true if resource label rl1 dominates rl2
  bool dominates(int conso1, int conso2) override;

  int totalNbDays() const { return totalNbDays_; }

  double getWorstLbCost(int consumption) const override {
    return .0;
  }

 protected:
  // initialize the expander on a given arc
  PExpander init(const Shift &prevShift,
                 const Stretch &stretch,
                 const RCArc &arc) override;

 private:
  int totalNbDays_;  // Total number of days in the horizon
};

/**
 * Expander that allows to count the consumption of the total weekend resources
 */
struct SoftTotalWeekendsExpander : public Expander {
  SoftTotalWeekendsExpander(const SoftTotalWeekendsResource& resource,
                            const Stretch &stretch) :
      Expander(resource.id()),
      resource_(resource),
      nSaturdaysLeft_(0),
      nSundaysLeft_(0),
      stretch_(stretch),
      arcToSink_(false) {}

  SoftTotalWeekendsExpander(const SoftTotalWeekendsResource& resource,
                            const Stretch &stretch,
                            int nSaturdaysLeft,
                            int nSundaysLeft) :
      Expander(resource.id()),
      resource_(resource),
      nSaturdaysLeft_(nSaturdaysLeft),
      nSundaysLeft_(nSundaysLeft),
      stretch_(stretch),
      arcToSink_(false) {}

  SoftTotalWeekendsExpander(const SoftTotalWeekendsResource& resource,
                            const Stretch &stretch,
                            int nSaturdaysLeft,
                            int nSundaysLeft,
                            bool arcToSink) :
      Expander(resource.id()),
      resource_(resource),
      nSaturdaysLeft_(nSaturdaysLeft),
      nSundaysLeft_(nSundaysLeft),
      stretch_(stretch),
      arcToSink_(arcToSink) {}

  // expand a given resource label
  bool expand(const ResourceValues &vParent,
              const PRCLabel &pLChild,
              ResourceValues *vChild) override;

 protected:
  const SoftTotalWeekendsResource& resource_;
  int nSaturdaysLeft_;
  int nSundaysLeft_;
  const Stretch stretch_;
  bool arcToSink_;  // true if the target of the arc is a sink node
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_TOTALWEEKENDSRESOURCE_H_
