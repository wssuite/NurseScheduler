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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_UNWANTEDSHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_UNWANTEDSHIFTRESOURCE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data/Nurse.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

class UnwantedShiftResource : public Resource {
 public:
  // Constructor
  explicit UnwantedShiftResource(
          int nShifts,
          vector<PShift> unwantedShifts,
          double cost = HARD_COST) :
          Resource("Unwanted shift"),
          isUnwanted_(nShifts, false),
          cost_(cost) {
    costType_ = UNWANTED_SHIFT_COST;
    for (const auto &pS : unwantedShifts) isUnwanted_[pS->id] = true;
  }

  int getConsumption(const State &initialState) const override { return 0; }

  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  BaseResource *clone() const override {
    return new UnwantedShiftResource(*this);
  }

  bool isHard() const override { return isHardCost(cost_); }

  // add the cost of preference violation to all the arcs of the input graph
  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc &pA, double *cost) override;

  bool isUnwantedShift(int s) const { return isUnwanted_[s]; }
  bool isUnwantedShift(const PShift &pS) const {
    return isUnwantedShift(pS->id);
  }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return false; };


  int findMaxOptimalGap() const override {
    return isHard() ? 0 : cost_;
  }

 protected:
  std::vector<bool> isUnwanted_;
  double cost_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_UNWANTEDSHIFTRESOURCE_H_
