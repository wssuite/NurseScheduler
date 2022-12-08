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

#ifndef SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_
#define SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_

#include "ConstraintsMP.h"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

template <class H, class S>
class BoundedResourceConstraint : public ConstraintMP {
 public:
  explicit BoundedResourceConstraint(MasterProblem *pMaster, std::string name);

  // compute the consumption of the resource H and S
  // (must be the same if both defined)
  int computeConsumption(const Stretch &st,
                         const std::pair<H*, S*> &p,
                         const PAbstractShift &prevS) const;

  // update the dual values of the constraints based on the current solution
  void updateDuals() override;

  // update the bounds and costs of slack variables as well as
  // bounds of the constraints
  void update() override;

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

  const std::vector<std::map<std::pair<H*, S*>, std::pair<MyVar*, MyVar*>>>&
  getVariables() const {
    return resourceVars_;
  }

  const std::vector<std::map<std::pair<H*, S*>, std::pair<MyCons*, MyCons*>>>&
  getConstraints() const {
    return resourceCons_;
  }

  double getTotalCost() const override {
    return pModel()->getTotalCost(resourceVars_);
  }

  std::string individualCostToString() const;

 protected:
  // variables associated to the resources
  std::vector<std::map<std::pair<H*, S*>, std::pair<MyVar*, MyVar*>>>
      resourceVars_;
  // constraints associated to resources
  std::vector<std::map<std::pair<H*, S*>, std::pair<MyCons*, MyCons*>>>
      resourceCons_;
  // dual values per nurse and resources
  std::vector<std::map<std::pair<H*, S*>, double>> dualValues_;

  // update the dual values of the constraints randomly
  void createRandomUpdateDuals(double weight);

  // build a constraint for the LB and the UB of these bounded resources:
  // pHR = pointer to a Hard Resource, pSR = pointer to a Soft Resource
  // WARNING: The two resources MUST represent the same resource if both defined
  // The goal of having both is to be able to define both soft and hard bounds
  // on the same constraint
  // If there is only a soft constraint (SR), set pHR to nullptr
  // If there is only a hard constraint (HR), set pSR to nullptr
  void build(H *pHR, S *pSR, const LiveNurse &pN);

  // compute bounds for constraint and slack
  std::pair<std::pair<int, int>, std::pair<int, int>>
  computeBounds(H *pHR, S *pSR) const;
};

class TotalShiftDurationConstraint : public BoundedResourceConstraint<
    HardTotalShiftDurationResource, SoftTotalShiftDurationResource> {
 public:
  explicit TotalShiftDurationConstraint(MasterProblem *pMaster);

  void addConstraintFor(const shared_ptr<HardTotalShiftDurationResource> &pHR,
                        const shared_ptr<SoftTotalShiftDurationResource> &pSR,
                        const LiveNurse &pN);

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;

  std::string writeIndividualCost() const;
};

class TotalWeekendConstraint : public BoundedResourceConstraint<
    HardTotalWeekendsResource, SoftTotalWeekendsResource> {
 public:
  explicit TotalWeekendConstraint(MasterProblem *pMaster);

  void addConstraintFor(const shared_ptr<HardTotalWeekendsResource> &pHR,
                        const shared_ptr<SoftTotalWeekendsResource> &pSR,
                        const LiveNurse &pN);

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;

  std::string writeIndividualCost() const;
};

// class ConsWeekendConstraint : public BoundedResourceConstraint<
//    HardConsWeekendShiftResource, SoftConsWeekendShiftResource> {
// public:
//  explicit ConsWeekendConstraint(MasterProblem *pMaster);
//
//  void addConstraintFor(const shared_ptr<HardConsWeekendShiftResource> &pHR,
//                        const shared_ptr<SoftConsWeekendShiftResource> &pSR,
//                        const LiveNurse &pN);
//
//  // update the dual values of the constraints randomly
//  void randomUpdateDuals(bool useInputData, int nPerturbations) override;
//
//  // compute consumption
//  int computeConsumption(
//      const Stretch &st,
//      const std::pair<HardConsWeekendShiftResource*,
//      SoftConsWeekendShiftResource*> &p,
//      const PAbstractShift &prevS) const;
//};
#endif  // SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_
