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

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o s t e r C o l u m n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct RosterPattern : Pattern {
  // Specific constructors and destructors
  RosterPattern(const std::vector<int> &shifts,
                const PScenario &pScenario,
                int nurseNum = -1,
                double cost = DBL_MAX,
                double dualCost = DBL_MAX) :
      Pattern(0, shifts, pScenario, nurseNum, cost, dualCost),
      nbWeekends_(-1) {}

  RosterPattern(const std::vector<double> &compactPattern,
                const PScenario &pScenario) :
      Pattern(compactPattern, pScenario),
      nbWeekends_(static_cast<int>(compactPattern.back())) {}

  RosterPattern(const RosterPattern &roster, int nurseNum) :
      Pattern(roster, nurseNum),
      nbWeekends_(roster.nbWeekends_) {}

  ~RosterPattern() = default;

  // When branching on this pattern, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid any shifts on any days as the nurse already has a roster.
  void addForbiddenShifts(std::set<std::pair<int, int> > *forbidenShifts,
                          int nbShifts,
                          PDemand pDemand) const override;

  // Level of the branch and bound tree where the rotation has been generated
  int treeLevel_ = 0;

  int nbWeekends_;

  // compact the rotation in a vector
  std::vector<double> getCompactPattern() const override {
    std::vector<double> pattern = Pattern::getCompactPattern();
    pattern.push_back(nbWeekends_);
    return pattern;
  }

  // Compute the cost of a roster
  void computeCost(const MasterProblem *pMaster,
      const PLiveNurse &pNurse) override;

  // Compute the reduced cost of a roster and compare it to the one found
  // by the subproblem
  void checkReducedCost(
      const PDualCosts &costs, PScenario Scenario, bool printBadPricing = true);

  // Returns true if both columns are disjoint (needRest not used)
  bool isDisjointWith(PPattern pat, bool needRest = true) const override {
    if (pat->nurseNum_ == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < stretch_.nDays(); k++)
      // work both
      if (stretch_.pShift(k)->isWork() && pat->stretch_.pShift(k)->isWork())
        return false;
    return true;
  };

  // Returns true if both columns are disjoint for shifts
  bool isShiftDisjointWith(PPattern pat, bool needRest = true) const override {
    if (pat->nurseNum_ == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < stretch_.nDays(); k++) {
      const PShift &pS1 = stretch_.pShift(k), &pS2 = pat->stretch_.pShift(k);
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
  RosterMP(PScenario pScenario,
           PDemand pDemand,
           PPreferences pPreferences,
           std::vector<State> *pInitState,
           SolverType solver);
  virtual ~RosterMP();

  PPattern getPattern(MyVar *var) const override;

  MyVar *addColumn(int nurseNum, const RCSolution &solution) override;

  // define the resources used for the sub problem
  std::vector<PResource> createResources(
      const PLiveNurse &pN,
      std::map<int, CostType> *resourceCostType) const override;

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

  // Create a new rotation variable
  // add the correct constraints and coefficients for the nurse i working
  // on a rotation
  // if s=-1, the nurse works on all shifts
  MyVar *addRoster(const RosterPattern &rotation,
                   const char *baseName,
                   bool coreVar = false);

  /* Build each set of constraints
   * Add also the coefficient of a column for each set
   */
  void buildAssignmentCons(const SolverParam &parameters);
  int addRosterConsToCol(std::vector<MyCons *> *cons,
                         std::vector<double> *coeffs,
                         int i);

  // return the costs of all active columns associated to the type
  double getColumnsCost(CostType costType) const override;
  double getColumnsCost(CostType costType,
                        const std::vector<MyVar *> &vars) const;

  double getDaysCost() const override;
  double getWeekendCost() const override;

  /* retrieve the dual values */
  double getConstantDualvalue(PLiveNurse pNurse) const override;

  /*
  * Constraints
  */
  // Ensure that is nurse has a roster assigned to her
  std::vector<MyCons *> assignmentCons_;
};

#endif  // SRC_SOLVERS_MP_ROSTERMP_H_
