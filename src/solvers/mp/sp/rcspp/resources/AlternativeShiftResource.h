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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_ALTERNATIVESHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_ALTERNATIVESHIFTRESOURCE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data/Nurse.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

class AlternativeShiftResource : public Resource {
 public:
  explicit AlternativeShiftResource(const PNurse &pNurse, int nShifts) :
      Resource("Alt shift"),
      isAlternativeShift_(nShifts, false),
      cost_(pNurse->pContract_->costAlternativeShift_) {
    for (int s : pNurse->alternativeShifts_)
      isAlternativeShift_[s] = true;
    costType_ = ALTERNATIVE_COST;
  }

  // Constructor functioning for the UI input workflow
  explicit AlternativeShiftResource(
          int nShifts, vector<PShift> altShifts, double cost = LARGE_SCORE) :
      Resource("Alt shift"),
      isAlternativeShift_(nShifts, false),
      cost_(cost) {
    costType_ = ALTERNATIVE_COST;
    for (const auto &pS : altShifts) isAlternativeShift_[pS->id] = true;
  }

  int getConsumption(const State &initialState) const override { return 0; }

  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  BaseResource *clone() const override {
    return new AlternativeShiftResource(*this);
  }

  bool isHard() const override { return cost_ >= LARGE_SCORE; }

  // add the cost of preference violation to all the arcs of the input graph
  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc &pA, double *cost) override;

  bool isAlternativeShift(int s) const { return isAlternativeShift_[s]; }
  bool isAlternativeShift(const PShift &pS) const {
    return isAlternativeShift(pS->id);
  }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return false; };

 protected:
  std::vector<bool> isAlternativeShift_;
  double cost_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_ALTERNATIVESHIFTRESOURCE_H_
