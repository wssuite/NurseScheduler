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
#include <climits>  // for INT_MAX
#include <set>
#include <map>
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

struct RotationPattern : public Pattern {
  // Specific constructors and destructors
  //

  RotationPattern(int firstDay,
                  std::vector<PShift> pShifts,
                  int nurseNum,
                  double cost = DBL_MAX,
                  double dualCost = DBL_MAX) :
      RotationPattern(RCSolution(firstDay, pShifts, cost, dualCost),
                      nurseNum) {}

  RotationPattern(RCSolution sol, int nurseNum) :
      Pattern(std::move(sol), nurseNum) {}

  RotationPattern(MyVar *var, const PScenario &pScenario) :
      Pattern(var, pScenario) {}

  ~RotationPattern() {}

  // when branching on this pattern, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid all the shifts that would be worked on a day
  // that is already covered by this pattern.
  // Moreover, there needs to be a resting day before and after each rotation,
  // so the shifts can also be forbidden on these two days
  // (if the rotation is not at an extremity of the horizon).
  void addForbiddenShifts(std::set<std::pair<int, int> > *forbiddenShifts,
                          int nbShifts,
                          PDemand pDemand) const override;

  // Level of the branch and bound tree where the rotation has been generated
  int treeLevel_ = 0;

  //  Compute the dual cost of a column
  void checkReducedCost(const DualCosts &dualCosts,
                        bool printBadPricing = true) const override;
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
  virtual ~RotationMP();

  PPattern getPattern(MyVar *var) const override;

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

  // Provide an initial solution to the solver. If empty, add artificial columns
  void initializeSolution(const std::vector<Roster> &solution) override;

  void buildResourceCons(const SolverParam &param);

  // return the costs of all active columns associated to the type
  double getColumnsCost(CostType costType) const override;
  double getDaysCost() const override;
  double getWeekendCost() const override;

  // separate the resources between the ones that will be managed by
  // the master, the master rotation graph, and the sub problems
  // must initialize spResources_
  void splitPResources() override;

  // return the value V used to choose the number of columns on which to branch.
  // Choose as many columns as possible such that: sum (1 - value(column)) < V
  double getBranchColumnValueMax() const override {
    return std::max(1.0, pScenario_->nWeeks() / 2.0);
  }

  /*
   * RPresources
   */
  vector2D<PBoundedResource> masterConstraintResources_,
      masterRotationGraphResources_;

  // Rotation graph
  RotationGraphConstraint *rotationGraphConstraint_;

  // constraints associated to resources
  TotalShiftDurationConstraint *totalShiftDurationConstraint_;
  TotalWeekendConstraint *totalWeekendConstraint_;
  DynamicConstraints *dynamicConstraints_;
};

#endif  // SRC_SOLVERS_MP_ROTATIONMP_H_
