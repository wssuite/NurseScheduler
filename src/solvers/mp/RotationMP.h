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

#include <climits>  // for INT_MAX

#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct RotationDualCosts : public DualCosts {
 public:
  RotationDualCosts(vector2D<double> workedShiftsCosts,
                    std::vector<double> startWorkCosts,
                    std::vector<double> endWorkCosts,
                    double workedWeekendCost,
                    double constant = 0) :
      DualCosts(std::move(workedShiftsCosts), constant),
      startWorkCosts_(std::move(startWorkCosts)),
      endWorkCosts_(std::move(endWorkCosts)),
      workedWeekendCost_(workedWeekendCost) {}

  // GETTERS
  //
  double startWorkCost(int day) const override {
    return (startWorkCosts_[day]);
  }
  double endWorkCost(int day) const override {
    return (endWorkCosts_[day]);
  }
  double workedWeekendCost() const override {
    return workedWeekendCost_;
  }

 protected:
  // Indexed by : day
  std::vector<double> startWorkCosts_;

  // Indexed by : day
  std::vector<double> endWorkCosts_;

  // Reduced cost of the weekends
  double workedWeekendCost_;
};

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct RotationPattern : Pattern {
  // Specific constructors and destructors
  //
  RotationPattern(std::map<int, int> shifts,
                  const PScenario &pScenario,
                  int nurseNum = -1,
                  double cost = DBL_MAX,
                  double dualCost = DBL_MAX) :
      Pattern(shifts, pScenario, nurseNum, cost, dualCost),
      consShiftsCost_(0),
      consDaysWorkedCost_(0),
      completeWeekendCost_(0),
      preferenceCost_(0),
      initRestCost_(0) {}

  RotationPattern(int firstDay,
                  std::vector<int> shiftSuccession,
                  const PScenario &pScenario,
                  int nurseNum = -1,
                  double cost = DBL_MAX,
                  double dualCost = DBL_MAX) :
      Pattern(firstDay, shiftSuccession, pScenario, nurseNum, cost, dualCost),
      consShiftsCost_(0),
      consDaysWorkedCost_(0),
      completeWeekendCost_(0),
      preferenceCost_(0),
      initRestCost_(0) {}

  explicit RotationPattern(const std::vector<double> &compactPattern,
                           const PScenario &pScenario) :
      Pattern(compactPattern, pScenario),
      consShiftsCost_(0),
      consDaysWorkedCost_(0),
      completeWeekendCost_(0),
      preferenceCost_(0),
      initRestCost_(0) {}

  RotationPattern(const RotationPattern &rotation, int nurseNum) :
      Pattern(rotation, nurseNum),
      consShiftsCost_(rotation.consShiftsCost_),
      consDaysWorkedCost_(rotation.consDaysWorkedCost_),
      completeWeekendCost_(rotation.completeWeekendCost_),
      preferenceCost_(rotation.preferenceCost_),
      initRestCost_(rotation.initRestCost_) {}

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

  // Cost
  //
  double consShiftsCost_, consDaysWorkedCost_, completeWeekendCost_,
      preferenceCost_, initRestCost_;

  // Level of the branch and bound tree where the rotation has been generated
  int treeLevel_ = 0;

  // Compute the cost of a rotation
  void computeCost(const MasterProblem *pMaster,
                   const PLiveNurse &pNurse) override;

  //  Compute the dual cost of a rotation
  void checkReducedCost(const PDualCosts &costs, bool printBadPricing = true);
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
             PDemand pDemand,
             PPreferences pPreferences,
             std::vector<State> *pInitState,
             SolverType solver);
  virtual ~RotationMP();

  PPattern getPattern(MyVar *var) const override;

  MyVar *addColumn(int nurseNum, const RCSolution &solution) override;

  // get a reference to the restsPerDay_ for a Nurse
  std::vector<MyVar *> getRestVarsPerDay(PLiveNurse pNurse,
                                         int day) const override {
    return restsPerDay_[pNurse->num_][day];
  }

  // build the, possibly fractional, roster corresponding to the solution
  // currently stored in the model
  vector3D<double> fractionalRoster() const override;

  const std::vector<MyVar *> &getMinWorkedDaysVars() const {
    return minWorkedDaysVars_;
  }
  const std::vector<MyVar *> &getMaxWorkedDaysVars() const {
    return maxWorkedDaysVars_;
  }
  const std::vector<MyVar *> &getMaxWorkedWeekendVars() const {
    return maxWorkedWeekendVars_;
  }

 protected:
  // Main method to build the rostering problem for a given input
  void build(const SolverParam &parameters) override;

  // Provide an initial solution to the solver. If empty, add artificial columns
  void initializeSolution(const std::vector<Roster> &solution) override;

  // Create a new rotation variable
  // add the correct constraints and coefficients for the nurse i working
  // on a rotation
  // if s=-1, the nurse works on all shifts
  MyVar *addRotation(const RotationPattern &rotation,
                     const char *baseName,
                     bool coreVar = false);

  // compute and add the last rotation finishing on the day just before
  // the first one
  RotationPattern computeInitStateRotation(const PLiveNurse& pNurse);

  /* Build each set of constraints
   * Add also the coefficient of a column for each set */
  void buildRotationCons(const SolverParam &parameters);
  int addRotationConsToCol(std::vector<MyCons *> *cons,
                           std::vector<double> *coeffs,
                           int i,
                           int k,
                           bool firstDay,
                           bool lastDay);

  void buildMinMaxCons(const SolverParam &parameters);
  int addMinMaxConsToCol(std::vector<MyCons *> *cons,
                         std::vector<double> *coeffs,
                         int i,
                         int nbDays,
                         int nbWeekends);

  // return the costs of all active columns associated to the type
  double getColumnsCost(CostType costType) const override;
  double getColumnsCost(CostType costType,
                        const std::vector<MyVar *> &vars) const;

  double getDaysCost() const override;
  double getWeekendCost() const override;

  /* retrieve the dual values */
  PDualCosts buildDualCosts(PLiveNurse pNurse) const override;
  vector2D<double> getShiftsDualValues(PLiveNurse pNurse) const override;
  std::vector<double> getStartWorkDualValues(const PLiveNurse& pNurse) const;
  std::vector<double> getEndWorkDualValues(const PLiveNurse& pNurse) const;
  double getWorkedWeekendDualValue(const PLiveNurse& pNurse) const;

  PDualCosts buildRandomDualCosts(bool optimalDemandConsidered,
                                  int NDaysShifts) const override;

  // Functions to generate the resources for a given nurse
  std::map<PResource, CostType>
  defaultgeneratePResources(const PLiveNurse &pN) const override {
    Tools::throwError("defaultgeneratePResources is not implemented.");
    return {};
  }

  /*
  * Variables
  */
  // stores all the arcs that are resting on a day for each nurse
  vector3D<MyVar *> restsPerDay_;
  // binary variables for the resting arcs in the rotation network
  vector2D<MyVar *> restingVars_;
  // binary variables for the resting arcs in the rotation network
  vector3D<MyVar *> longRestingVars_;
  // stores all the initial rotations finishing on the first day
  std::vector<MyVar *> initialStateVars_;

  // count the number of missing worked days per nurse
  std::vector<MyVar *> minWorkedDaysVars_;
  // count the number of exceeding worked days per nurse
  std::vector<MyVar *> maxWorkedDaysVars_;
  // count the number of exceeding worked weekends per nurse
  std::vector<MyVar *> maxWorkedWeekendVars_;

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
  // transmission of the flow on the resting nodes
  // initialization of the flow constraint at the first position of
  // each restFlowCons_[i] (i=nurse)
  vector2D<MyCons *> restFlowCons_;
  // transmission of the flow on the working nodes
  // end of the flow constraint at the last position of each
  // workFlowCons_[i] (i=nurse)
  vector2D<MyCons *> workFlowCons_;

  // count the number of missing worked days per nurse
  std::vector<MyCons *> minWorkedDaysCons_;
  // count the number of exceeding worked days per nurse
  std::vector<MyCons *> maxWorkedDaysCons_;
  // count the number of exceeding worked weekends per nurse
  std::vector<MyCons *> maxWorkedWeekendCons_;
  // count the total number of weekends that will be penalized
  MyCons *sumMaxWorkedWeekendCons_;

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

  // vectors of booleans indicating whether some above constraints are present
  // in the model
  std::vector<bool> isMinWorkedDaysAvgCons_,
      isMaxWorkedDaysAvgCons_,
      isMaxWorkedWeekendAvgCons_,
      isMinWorkedDaysContractAvgCons_,
      isMaxWorkedDaysContractAvgCons_,
      isMaxWorkedWeekendContractAvgCons_;
};

#endif  // SRC_SOLVERS_MP_ROTATIONMP_H_
