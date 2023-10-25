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

#ifndef SRC_SOLVERS_MP_ROSTERMP_H_
#define SRC_SOLVERS_MP_ROSTERMP_H_

#include "MasterProblem.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/constraints/AssignmentConstraints.h"
#include "solvers/mp/constraints/DynamicConstraints.h"


//-----------------------------------------------------------------------------
//
//  S t r u c t   R o s t e r C o l u m n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct RosterColumn : public Column {
  // Specific constructors and destructors
  RosterColumn(std::vector<PShift> pShifts,
                int nurseNum,
                double cost = DBL_MAX,
                double dualCost = DBL_MAX) :
      Column(RCSolution(0, std::move(pShifts), cost, dualCost), nurseNum) {}

  RosterColumn(RCSolution rcSol, int nurseNum) :
      Column(std::move(rcSol), nurseNum) {}

  RosterColumn(MyVar *var, const PScenario &pScenario) :
      Column(var, pScenario) {}

  // When branching on this column, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid any shifts on any days as the nurse already has a roster.
  void addForbiddenShifts(
          std::set<std::pair<int, int> > *forbiddenShifts,
          int firstDayId, int nDays, int nbShifts) const override;

  // Compute the reduced cost of a roster and compare it to the one found
  // by the subproblem
  void checkReducedCost(const DualCosts &dualCosts,
                        bool printBadPricing) const override;

  // Returns true if both columns are disjoint (needRest not used)
  bool isDisjointWith(const PColumn &pat, bool needRest) const override {
    if (pat->nurseNum() == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < nDays(); k++)
      // work both
      if (pShift(k)->isWork() && pat->pShift(k)->isWork())
        return false;
    return true;
  };

  // Returns true if both columns are disjoint for shifts
  bool isShiftDisjointWith(const PColumn &pat, bool needRest) const override {
    if (pat->nurseNum() == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < nDays(); k++) {
      const PShift &pS1 = pShift(k), &pS2 = pat->pShift(k);
      if (pS1->isWork() && pS2->isWork() && pS1->id == pS2->id)
        return false;  // work on same shift
    }
    return true;
  }
};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------
class RosterMP : public MasterProblem {
 public:
  RosterMP(const PScenario& pScenario,
           SolverType solver);
  virtual ~RosterMP();

  PColumn getPColumn(MyVar *var) const override;
  PColumn getPColumn(const RCSolution &st, int nurseNum) const override;

  MyVar *addColumn(int nurseNum, const RCSolution &solution) override;

  // get a reference to the restsPerDay_ for a Nurse
  std::vector<MyVar *> getRestVarsPerDay(PLiveNurse pNurse,
                                         int day) const override;

  // STAB: compute the lagrangian bound
  double computeLagrangianBound(double objVal) const override;

  // return the value V used to choose the number of columns on which to branch.
  // Choose as many columns as possible such that: sum (1 - value(column)) < V
  double getBranchColumnValueMax() const override {
    return 1 - pModel_->epsilon();
  }

 protected:
  // Main method to build the rostering problem for a given input
  void build(const SolverParam &parameters) override;

  // Provide an initial solution to the solver.
  // If empty, add artificial columns
  void initializeSolution(const std::vector<Roster> &solution) override;

  // split the resources between the master and the subproblem
  // must initialize spResources_
  void splitPResources() override;

  /*
  * Constraints
  */
  // Ensure that is nurse has a roster assigned to her
  std::unique_ptr<RosterAssignmentConstraint> assignmentConstraint_;
  std::unique_ptr<DynamicConstraints> dynamicConstraints_;
};

#endif  // SRC_SOLVERS_MP_ROSTERMP_H_
