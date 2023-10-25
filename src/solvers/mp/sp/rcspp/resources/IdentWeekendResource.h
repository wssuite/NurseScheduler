/*
 * Copyright (C) 2021 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_IDENTWEEKENDRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_IDENTWEEKENDRESOURCE_H_

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::vector;

/**
 * Resource corresponding to the requirement of having homogeneous
 * assignments in each weekend. May be that the complete weekend must be
 * worked or rested, or that all the assignments of the weekend must have the
 * same shift type.
 */
class IdentWeekendResource : public Resource {
 public:
  explicit IdentWeekendResource(PBaseShiftComparator pShiftCmp,
                                DayOfWeek firstWeekendDay = SATURDAY,
                                DayOfWeek lastWeekendDay = SUNDAY,
                                std::string name = "") :
      Resource(name.empty() ? "Ident Weekend" : std::move(name)),
      pShiftCmp_(std::move(pShiftCmp)),
      weekend_(firstWeekendDay, lastWeekendDay) {}

  int getConsumption(const State &initialState) const override { return  0; }

  const Weekend& weekend() const { return weekend_; }

  // this resource is taken into account in the preprocessing of the graph,
  // do we never need to check it for domination
  bool isActive(int dayId, const PAbstractShift &pAShift) const override {
    return false;
  }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return false; };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override = 0;

  const PBaseShiftComparator pShiftCmp_;
  const Weekend weekend_;
};

class SoftIdentWeekendResource : public IdentWeekendResource {
 public:
  explicit SoftIdentWeekendResource(PBaseShiftComparator pShiftCmp,
                                    double cost = 0,
                                    DayOfWeek first = SATURDAY,
                                    DayOfWeek last = SUNDAY,
                                    std::string name = "") :
      IdentWeekendResource(std::move(pShiftCmp), first, last, std::move(name)),
      cost_(cost) {
    costType_ = IDENT_WEEKEND_COST;
  }

  BaseResource* clone() const override {
    return new SoftIdentWeekendResource(
        pShiftCmp_, cost_,
        weekend_.firstWeekendDay().getDayOfWeek(),
        weekend_.lastWeekendDay().getDayOfWeek(), name);
  }

  bool isHard() const override { return false; }

  int findMaxOptimalGap() const override { return cost_; }

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;


  const double cost_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_IDENTWEEKENDRESOURCE_H_
