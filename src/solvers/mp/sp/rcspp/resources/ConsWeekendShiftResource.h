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

#include "ConsShiftResource.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

/**
 * Resource corresponding to the soft min/max constraints on the number of
 * consecutive shifts on weekends: it may be used for any kind of abstract shift
 */
class SoftConsWeekendShiftResource : public SoftConsShiftResource {
 public:
  SoftConsWeekendShiftResource(
      int lb, int ub, double lbCost, double ubCost,
      const PAbstractShift pShift, int totalNbDays = 0) :
      SoftConsShiftResource(lb, ub, lbCost, ubCost, pShift, totalNbDays,
                            "Soft Weekend Cons "+pShift->name) {}

  int getConsumption(const State &initialState) const override;

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const AbstractShift &aShift) const override {
    return Tools::isWeekend(dayId) && pShift_->includes(aShift);
  }

  double getWorstUbCost(int consumption, int nLeft) const override {
    return ubCost_ * std::min(consumption, ub_);
  }

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;
};

class HardConsWeekendShiftResource : public HardConsShiftResource {
 public:
  HardConsWeekendShiftResource(int lb, int ub, const PAbstractShift pShift) :
      HardConsShiftResource(lb, ub, pShift,
                            "Hard Weekend Cons "+pShift->name) {}

  int getConsumption(const State &initialState) const override;

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const PRCArc &pArc) override;
};

#endif  // NURSESCHEDULER_SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_  // NOLINT
