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

#ifndef SRC_SOLVERS_MP_ROTATIONMP_H_
#define SRC_SOLVERS_MP_ROTATIONMP_H_

#include "MasterProblem.h"

#include <algorithm>
#include <set>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/constraints/AssignmentConstraints.h"
#include "solvers/mp/constraints/ResourceConstraints.h"
#include "solvers/mp/constraints/DynamicConstraints.h"


//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct RotationColumn : public Column {
  // Specific constructors and destructors
  //
  RotationColumn(int firstDay,
                  std::vector<PShift> pShifts,
                  int nurseNum,
                  double cost = DBL_MAX,
                  double dualCost = DBL_MAX) :
      RotationColumn(RCSolution(firstDay, std::move(pShifts), cost, dualCost),
                      nurseNum) {}

  RotationColumn(RCSolution sol, int nurseNum) :
      Column(std::move(sol), nurseNum) {}

  RotationColumn(MyVar *var, const PScenario &pScenario) :
      Column(var, pScenario) {}

  ~RotationColumn() override = default;

  // when branching on this column, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid all the shifts that would be worked on a day
  // that is already covered by this column.
  // Moreover, there needs to be a resting day before and after each rotation,
  // so the shifts can also be forbidden on these two days
  // (if the rotation is not at an extremity of the horizon).
  void addForbiddenShifts(
          std::set<std::pair<int, int> > *forbiddenShifts,
          int firstDayId, int nDays, int nbShifts) const override;

  //  Compute the dual cost of a column
  void checkReducedCost(const DualCosts &dualCosts,
                        bool printBadPricing) const override;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------
class RotationMP : public MasterProblem {
 public:
  RotationMP(const PScenario& pScenario,
             SolverType solver);
  RotationMP(const PScenario& pScenario,
             SolverType solver,
             const SolverParam &param);
  ~RotationMP() override;

  PColumn getPColumn(MyVar *var) const override;
  PColumn getPColumn(const RCSolution &st, int nurseNum) const override;

  MyVar *addColumn(int nurseNum, const RCSolution &solution) override;

  // get a reference to the restsPerDay_ for a Nurse
  std::vector<MyVar *> getRestVarsPerDay(
      PLiveNurse pNurse, int day) const override {
    return rotationGraphConstraint_->getVariables(pNurse, day);
  }

  // build the, possibly fractional, roster corresponding to the solution
  // currently stored in the model
  vector3D<double> fractionalRoster() const override;

 protected:
  // Main method to build the rostering problem for a given input
  void build(const SolverParam &parameters) override;

  // update the demand with a new one of the same size
  // change the rhs of the constraints minDemandCons_ and optDemandCons_
  void update(vector<PDemand> pDemands) override;

  // Provide an initial solution to the solver. If empty, add artificial columns
  void initializeSolution(const std::vector<Roster> &solution) override;

  void buildResourceCons(const SolverParam &param);

  // return the costs of all active columns associated to the type
  std::map<CostType, double> getColumnsCosts() const override;

  // separate the resources between the ones that will be managed by
  // the master, the master rotation graph, and the sub problems
  // must initialize spResources_
  void splitPResources() override;

  // return the value V used to choose the number of columns on which to branch.
  // Choose as many columns as possible such that: sum (1 - value(column)) < V
  double getBranchColumnValueMax() const override {
    return std::max(1, nDays() / 14) - pModel_->epsilon();
  }

  /*
   * RPresources
   */
  vector2D<PBoundedResource> masterConstraintResources_,
      masterRotationGraphResources_;

  // Rotation graph
  std::unique_ptr<RotationGraphConstraint> rotationGraphConstraint_;

  // constraints associated to resources
  std::unique_ptr<TotalShiftDurationConstraint> totalShiftDurationConstraint_;
  std::unique_ptr<TotalWeekendConstraint> totalWeekendConstraint_;
//  std::unique_ptr<ConsWeekendConstraint> consWeekendConstraints_;
  std::unique_ptr<DynamicConstraints> dynamicConstraints_;
};

#endif  // SRC_SOLVERS_MP_ROTATIONMP_H_
