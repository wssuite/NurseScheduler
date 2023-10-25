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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FREEDAYSAFTERSHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FREEDAYSAFTERSHIFTRESOURCE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

class FreeDaysAfterShiftResource : public Resource {
 public:
  explicit FreeDaysAfterShiftResource(
      PAbstractShift pAShift,
      int nbDays = 0,
      std::string _name = "") :
      Resource(_name.empty() ?
               "Rest after work " : std::move(_name)),
      pAShift_(std::move(pAShift)),
      nbFreeDays_(nbDays) {}

  int getConsumption(const State &initialState) const override { return 0;}

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const PAbstractShift &pAShift) const override {
    return true;
  }

  int getUb() const {
    return nbFreeDays_ + 1;
  }

  bool isInRosterMaster() const override { return false; };
  // TODO(AL): I guess that this must be handled by a modifications of the
  //  rotation graph: not sure as it depends on the shift worked.
  //  Modifications will need to be done on both master and subproblem.
  bool isInRotationMaster() const override {
    Tools::throwError("FreeDaysAfterShiftResource cannot be handled "
                      "in the pattern if using a rotation decomposition.");
    return false;
  };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  const PAbstractShift pAShift_;
  const int nbFreeDays_;
};

class SoftFreeDaysAfterShiftResource : public FreeDaysAfterShiftResource {
 public:
  explicit SoftFreeDaysAfterShiftResource(
      PAbstractShift pAShift, int nbDays, double cost, std::string _name = ""):
      FreeDaysAfterShiftResource(std::move(pAShift),
                                 nbDays, std::move(_name)),
      cost_(cost) {
    costType_ = REST_AFTER_SHIFT_COST;
  }

  BaseResource* clone() const override {
    return new SoftFreeDaysAfterShiftResource(
        pAShift_, nbFreeDays_, cost_, name);
  }

  bool isHard() const override {return false;}
  double getCost() const {return cost_;}

  int findMaxOptimalGap() const override {
    return cost_;
  }

  void enumerate(const PRCGraph &pRCGraph, bool forceEnum) override;

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override {
    std::cerr << "SoftFreeDaysAfterShiftResource cannot be handled "
                 "if using a rotation decomposition."
              << std::endl;
    return true;
  };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  const double cost_;
};


#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FREEDAYSAFTERSHIFTRESOURCE_H_
