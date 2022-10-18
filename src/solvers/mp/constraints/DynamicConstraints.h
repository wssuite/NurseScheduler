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

#ifndef SRC_SOLVERS_MP_CONSTRAINTS_DYNAMICCONSTRAINTS_H_
#define SRC_SOLVERS_MP_CONSTRAINTS_DYNAMICCONSTRAINTS_H_

#include <string>
#include <utility>
#include <vector>

#include "ConstraintsMP.h"
#include "ResourceConstraints.h"


class DynamicConstraints : public ConstraintMP {
 public:
  explicit DynamicConstraints(
      MasterProblem *pMaster,
      TotalShiftDurationConstraint *totalShiftDurationConstraint,
      TotalWeekendConstraint *totalWeekendConstraint);

  // update the right hand side of the constraints based on the dynamic weights
  void update() override;

  // update the dual values of the constraints based on the current solution
  void updateDuals() override;

  // return the dual cost of a stretch based on its consumption of
  // the constraints
  double getDualCost(int nurseNum,
                     const Stretch &st,
                     const PAbstractShift &prevS) const override;

  // add a given constraint to the column
  void addConsToCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Column &col) const override;

  std::string toString(int nurseNum, const Stretch &st) const override;

  double getTotalCost() const override {
    return pModel()->getTotalCost(minWorkedDaysAvgVars_) +
        pModel()->getTotalCost(maxWorkedDaysAvgVars_) +
        pModel()->getTotalCost(maxWorkedWeekendAvgVars_) +
        pModel()->getTotalCost(minWorkedDaysContractAvgVars_) +
        pModel()->getTotalCost(maxWorkedDaysContractAvgVars_) +
        pModel()->getTotalCost(maxWorkedWeekendContractAvgVars_);
  }

 protected:
  int dynamicWeightsVersion_;
  TotalShiftDurationConstraint *totalShiftDurationConstraint_;
  std::vector<SoftTotalShiftDurationResource*> totalShiftDurationResources_;
  TotalWeekendConstraint *totalWeekendConstraint_;
  std::vector<SoftTotalWeekendsResource*> totalWeekendResources_;
  vector<double> dualValues_, weekendDualValues_;

  /*
   * Variables
   */
  // count the number of missing worked days from average per nurse
  std::vector<MyVar *> minWorkedDaysAvgVars_;
  // count the number of exceeding worked days from average per nurse
  std::vector<MyVar *> maxWorkedDaysAvgVars_;
  // count the number of exceeding worked weekends from average per nurse
  std::vector<MyVar *> maxWorkedWeekendAvgVars_;

  // count the number of missing worked days from average per contract
  std::vector<MyVar *> minWorkedDaysContractAvgVars_;
  // count the number of exceeding worked days from average per contract
  std::vector<MyVar *> maxWorkedDaysContractAvgVars_;
  // count the number of exceeding worked weekends from average per contract
  std::vector<MyVar *> maxWorkedWeekendContractAvgVars_;

  /*
   * Constraints
   */
  // count the number of missing worked days from average per nurse
  std::vector<MyCons *> minWorkedDaysAvgCons_;
  // count the number of exceeding worked days from average per nurse
  std::vector<MyCons *> maxWorkedDaysAvgCons_;
  // count the number of exceeding worked weekends from average per nurse
  std::vector<MyCons *> maxWorkedWeekendAvgCons_;

  // count the number of missing worked days from average per contract
  std::vector<MyCons *> minWorkedDaysContractAvgCons_;
  // count the number of exceeding worked days from average per contract
  std::vector<MyCons *> maxWorkedDaysContractAvgCons_;
  //  the number of exceeding worked weekends from average per contract
  std::vector<MyCons *> maxWorkedWeekendContractAvgCons_;

  // build the constraints
  void build();

  // return the consumption for the duration and the weekend resources
  std::pair<int, int> getConsumptions(
      int nurseNum, const Stretch &st, const PAbstractShift &prevS) const;
};

#endif  // SRC_SOLVERS_MP_CONSTRAINTS_DYNAMICCONSTRAINTS_H_
