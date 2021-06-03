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
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/constraints/AssignmentConstraints.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o s t e r C o l u m n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct RosterPattern : public Pattern {
  // Specific constructors and destructors
  RosterPattern(std::vector<PShift> pShifts,
                int nurseNum,
                double cost = DBL_MAX,
                double dualCost = DBL_MAX) :
      Pattern(RCSolution(0, pShifts, cost, dualCost), nurseNum) {}

  RosterPattern(RCSolution rcSol, int nurseNum) :
      Pattern(std::move(rcSol), nurseNum) {}

  RosterPattern(MyVar *var,
                const PScenario &pScenario) :
      Pattern(var, pScenario) {}

  ~RosterPattern() = default;

  // When branching on this pattern, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid any shifts on any days as the nurse already has a roster.
  void addForbiddenShifts(std::set<std::pair<int, int> > *forbidenShifts,
                          int nbShifts,
                          PDemand pDemand) const override;

  // Level of the branch and bound tree where the rotation has been generated
  int treeLevel_ = 0;

  // compact the rotation in a vector
  std::vector<double> getCompactPattern() const override {
    std::vector<double> pattern = Pattern::getCompactPattern();
    return pattern;
  }

  // Compute the reduced cost of a roster and compare it to the one found
  // by the subproblem
  void checkReducedCost(const DualCosts &dualCosts,
                        bool printBadPricing = true) const override;

  // Returns true if both columns are disjoint (needRest not used)
  bool isDisjointWith(PPattern pat, bool needRest = true) const override {
    if (pat->nurseNum() == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < nDays(); k++)
      // work both
      if (pShift(k)->isWork() && pat->pShift(k)->isWork())
        return false;
    return true;
  };

  // Returns true if both columns are disjoint for shifts
  bool isShiftDisjointWith(PPattern pat, bool needRest = true) const override {
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

  PPattern getPattern(MyVar *var) const override;

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

  double getDaysCost() const override;
  double getWeekendCost() const override;

  // split the resources between the master and the subproblem
  // must initialize spResources_
  void splitPResources() override {
    // put all the resources in the sub problem
    spResources_.clear();
    for (const auto &m : pResources_) {
      vector<PResource> pResources;
      for (const auto &p : m) {
        p.first->setId(pResources.size());
        pResources.push_back(p.first);
      }
      spResources_.push_back(pResources);
    }
  }

  /*
  * Constraints
  */
  // Ensure that is nurse has a roster assigned to her
  RosterAssignmentConstraint *assignmentConstraint_;
};

#endif  // SRC_SOLVERS_MP_ROSTERMP_H_
