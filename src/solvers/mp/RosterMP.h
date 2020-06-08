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

#include <utility>
#include <set>
#include <string>
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
  RosterPattern(std::vector<int> shifts,
                int nurseNum = -1,
                double cost = DBL_MAX,
                double dualCost = DBL_MAX) :
      Pattern(0, shifts.size(), nurseNum, cost, dualCost),
      shifts_(shifts),
      consShiftsCost_(0),
      consDaysOffCost_(0),
      consDaysWorkedCost_(0),
      completeWeekendCost_(0),
      preferenceCost_(0),
      initRestCost_(0),
      initWorkCost_(0),
      minDaysCost_(0),
      maxDaysCost_(0),
      maxWeekendCost_(0),
      timeDuration_(-1),
      nbWeekends_(-1) {}

  explicit RosterPattern(const std::vector<double> &compactPattern) :
      Pattern(compactPattern),
      shifts_(length_),
      consShiftsCost_(0),
      consDaysOffCost_(0),
      consDaysWorkedCost_(0),
      completeWeekendCost_(0),
      preferenceCost_(0),
      initRestCost_(0),
      initWorkCost_(0),
      minDaysCost_(0),
      maxDaysCost_(0),
      maxWeekendCost_(0),
      timeDuration_(static_cast<int>(
          compactPattern[compactPattern.size() - 2])),
      nbWeekends_(static_cast<int>(compactPattern.back())) {
    for (int k = 0; k < length_; k++)
      shifts_[firstDay_ + k] = static_cast<int>(compactPattern[k + 3]);
  }

  RosterPattern(const RosterPattern &roster, int nurseNum) :
      Pattern(roster, nurseNum),
      shifts_(roster.shifts_),
      consShiftsCost_(roster.consShiftsCost_),
      consDaysOffCost_(roster.consDaysOffCost_),
      consDaysWorkedCost_(roster.consDaysWorkedCost_),
      completeWeekendCost_(roster.completeWeekendCost_),
      preferenceCost_(roster.preferenceCost_),
      initRestCost_(roster.initRestCost_),
      initWorkCost_(roster.initWorkCost_),
      minDaysCost_(roster.minDaysCost_),
      maxDaysCost_(roster.maxDaysCost_),
      maxWeekendCost_(roster.maxWeekendCost_),
      timeDuration_(roster.timeDuration_),
      nbWeekends_(roster.nbWeekends_) {}

  ~RosterPattern() {}

  int getShift(int day) const override {
    return shifts_[day];
  }

  // When branching on this pattern, this method add the corresponding
  // forbidden shifts to the set.
  // It will forbid any shifts on any days as the nurse already has a roster.
  void addForbiddenShifts(std::set<std::pair<int, int> > *forbidenShifts,
                          int nbShifts,
                          PDemand pDemand) const override;

  // Shifts to be performed
  std::vector<int> shifts_;

  // Cost
  double consShiftsCost_, consDaysOffCost_, consDaysWorkedCost_,
      completeWeekendCost_, preferenceCost_,
      initRestCost_, initWorkCost_, minDaysCost_, maxDaysCost_, maxWeekendCost_;

  // Level of the branch and bound tree where the rotation has been generated
  int treeLevel_ = 0;

  // Time duration (in a certain unit: day, hours, half-hours, ...)
  int timeDuration_;
  int nbWeekends_;

  // compact the rotation in a vector
  std::vector<double> getCompactPattern() const override {
    std::vector<double> pattern = Pattern::getCompactPattern();
    pattern.insert(pattern.end(), shifts_.begin(), shifts_.end());
    pattern.push_back(timeDuration_);
    pattern.push_back(nbWeekends_);
    return pattern;
  }

  // Compute the cost of a roster
  void computeCost(PScenario pScenario,
                   const std::vector<PLiveNurse> &liveNurses);

  // Compute the reduced cost of a roster and compare it to the one found
  // by the subproblem
  void checkReducedCost(const DualCosts &costs, PScenario Scenario);

  // Returns true if both columns are disjoint (needRest not used)
  bool isDisjointWith(PPattern pat, bool needRest = true) const override {
    if (pat->nurseNum_ == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < shifts_.size(); k++)
      if (getShift(k) > 0 && pat->getShift(k) > 0)  // work both
        return false;
    return true;
  };

  // Returns true if both columns are disjoint for shifts !!
  bool isShiftDisjointWith(PPattern pat, bool needRest = true) const override {
    if (pat->nurseNum_ == nurseNum_)
      return false;  // cannot be disjoint if same nurse

    for (int k = 0; k < shifts_.size(); k++) {
      int s1 = getShift(k), s2 = pat->getShift(k);
      if (s1 > 0 && s2 > 0 && s1 == s2)
        return false;  // work on same shift
    }
    return true;
  }

  std::string toString(
      int nbDays = -1,
      std::vector<int> shiftIDToShiftTypeID = {}) const override;

  std::string costsToString(bool initialCosts = true) const;
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
           MySolverType solver);
  virtual ~RosterMP();

  PPattern getPattern(const std::vector<double> &pattern) const override;

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
  double getColumnsCost(CostType costType, bool historicalCosts) const override;
  double getColumnsCost(CostType costType,
                        const std::vector<MyVar *> &vars,
                        bool historicalCosts) const;

  double getMinDaysCost() const override;
  double getMaxDaysCost() const override;
  double getMaxWeekendCost() const override;

  /* retrieve the dual values */
  double getConstantDualvalue(PLiveNurse pNurse) const override;

  /*
  * Constraints
  */
  // Ensure that is nurse has a roster assigned to her
  std::vector<MyCons *> assignmentCons_;
};

#endif  // SRC_SOLVERS_MP_ROSTERMP_H_
