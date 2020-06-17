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

#ifndef SRC_SOLVERS_MP_MASTERPROBLEM_H_
#define SRC_SOLVERS_MP_MASTERPROBLEM_H_

#include <memory>
#include <algorithm>
#include <utility>
#include <set>
#include <vector>
#include <string>

#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/Solver.h"
#include "tools/MyTools.h"
#include "data/Nurse.h"
#include "solvers/mp/modeler/Modeler.h"
#include "OsiSolverInterface.hpp"

enum CostType {
  ROTATION_COST, CONS_SHIFTS_COST, CONS_WORKED_DAYS_COST,
  COMPLETE_WEEKEND_COST, PREFERENCE_COST, REST_COST,
  MIN_DAYS_COST, MAX_DAYS_COST, MAX_WEEKEND_COST
};

struct DualCosts {
 public:
  DualCosts(const vector2D<double> &workedShiftsCosts,
            const std::vector<double> &startWorkCosts,
            const std::vector<double> &endWorkCosts,
            double workedWeekendCost, double constant = 0) :
      workedShiftsCosts_(workedShiftsCosts),
      startWorkCosts_(startWorkCosts),
      endWorkCosts_(endWorkCosts),
      workedWeekendCost_(workedWeekendCost),
      constant_(constant) {}

  // GETTERS
  //
  int nDays() const { return startWorkCosts_.size(); }
  double workedDayShiftCost(int day, int shift) const {
    if (!shift) return 0;
    return (workedShiftsCosts_[day][shift - 1]);
  }
  double startWorkCost(int day) const { return (startWorkCosts_[day]); }
  double endWorkCost(int day) const { return (endWorkCosts_[day]); }
  double workedWeekendCost() const { return workedWeekendCost_; }
  double constant() const { return constant_; }

 protected:
  // Indexed by : (day, shift) !! 0 = shift 1 !!
  vector2D<double> workedShiftsCosts_;

  // Indexed by : day
  std::vector<double> startWorkCosts_;

  // Indexed by : day
  std::vector<double> endWorkCosts_;

  // Reduced cost of the weekends
  double workedWeekendCost_;

  // constant part of the reduced cost
  double constant_;
};

struct Pattern;
typedef std::shared_ptr<Pattern> PPattern;

struct Pattern {
  Pattern(int firstDay,
          int length,
          int nurseNum = -1,
          double cost = DBL_MAX,
          double dualCost = DBL_MAX) :
      nurseNum_(nurseNum), firstDay_(firstDay), length_(length),
      cost_(cost), reducedCost_(dualCost) {}

  Pattern(const Pattern &pat, int nurseNum) :
      nurseNum_(nurseNum),
      firstDay_(pat.firstDay_),
      length_(pat.length_),
      cost_(pat.cost_),
      reducedCost_(pat.reducedCost_) {
    if (pat.nurseNum_ != nurseNum_) {
      cost_ = DBL_MAX;
      reducedCost_ = DBL_MAX;
    }
  }

  explicit Pattern(const std::vector<double> &pattern) :
      nurseNum_(static_cast<int>(pattern[0])),
      firstDay_(static_cast<int>(pattern[1])),
      length_(static_cast<int>(pattern[2])),
      cost_(DBL_MAX),
      reducedCost_(DBL_MAX) {}

  virtual ~Pattern() {}

  virtual bool equals(PPattern pat) const {
    if (nurseNum_ != pat->nurseNum_) return false;
    if (firstDay_ != pat->firstDay_) return false;
    if (length_ != pat->length_) return false;
    for (int k = firstDay_; k < firstDay_ + length_; ++k) {
      if (getShift(k) != pat->getShift(k)) return false;
    }
    return true;
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isDisjointWith(PPattern pat, bool needRest = true) const {
    return ((firstDay_ + length_ < pat->firstDay_ - needRest)
        || (pat->firstDay_ + pat->length_ < firstDay_ - needRest));
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isShiftDisjointWith(PPattern pat, bool needRest = true) const {
    if (isDisjointWith(pat, needRest))
      return true;

    int commomFirstDay = std::max(firstDay_, pat->firstDay_),
        commomLastDay =
        std::min(firstDay_ + length_, pat->firstDay_ + pat->length_) - 1;
    for (int k = commomFirstDay; k <= commomLastDay; ++k)
      if (getShift(k) == pat->getShift(k)) return false;

    return true;
  }

  virtual std::string toString(
      int nbDays = -1,
      std::vector<int> shiftIDToShiftTypeID = {}) const {
    return "";
  }

  virtual int getShift(int day) const = 0;

  // when branching on this pattern,
  // this method add the corresponding forbidden shifts to the set
  virtual void addForbiddenShifts(
      std::set<std::pair<int, int> > *forbidenShifts,
      int nbShifts,
      PDemand pDemand) const = 0;

  // need to be able to write the pattern as a vector and
  // to create a new one from it
  virtual std::vector<double> getCompactPattern() const {
    std::vector<double> compact;
    compact.push_back(nurseNum_);
    compact.push_back(firstDay_);
    compact.push_back(length_);
    return compact;
  }

  int endDay() const {
    return firstDay_ + length_ - 1;
  }

  int nurseNum_;
  int firstDay_;
  int length_;

  // Cost
  double cost_;

  // Dual cost as found in the subproblem
  double reducedCost_;

  // Compare rotations on cost
  static bool compareCost(PPattern pat1, PPattern pat2);

  // Compare rotations on dual cost
  static bool compareDualCost(PPattern pat1, PPattern pat2);
};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

class MasterProblem : public Solver, public PrintSolution {
 public:
  // Specific constructor and destructor
  MasterProblem(PScenario pScenario,
                PDemand pDemand,
                PPreferences pPreferences,
                std::vector<State> *pInitState,
                MySolverType solver);
  ~MasterProblem();

  // solve the rostering problem
  double solve(const std::vector<Roster> &solution = {}) override;

  // solve the rostering problem or just the relaxation(root node)
  double solve(const std::vector<Roster> &solution, bool rebuild);

  // Solve with parameters
  double solve(const SolverParam &parameters,
               const std::vector<Roster> &solution = {}) override;

  // Resolve the problem with another demand and keep the same preferences
  double resolve(PDemand pDemand,
                 const SolverParam &parameters,
                 const std::vector<Roster> &solution = {}) override;

  // needs to be specialized: add a column  to the master from a solution of
  // the subproblem
  virtual MyVar *addColumn(int nurseNum, const RCSolution &solution) = 0;

  // retrieve the object represented ny the  vector pattern
  virtual PPattern getPattern(MyVar *var) const = 0;

  // throw an error if pattern is already present as an active column
  void checkIfPatternAlreadyPresent(const std::vector<double> &pattern) const;

  virtual std::vector<MyVar *> getRestVarsPerDay(PLiveNurse pNurse,
                                                 int day) const = 0;

  // get the pointer to the model
  Modeler *getModel() {
    return pModel_;
  }

  MyPricer *getPricer() {
    return pPricer_;
  }

  // override save and currentSolToString virtual method fom printFunction
  void save(const std::vector<int> &weekIndices, std::string outdir) override;
  std::string currentSolToString() const override;

  bool isSolutionInteger() const override {
    return pModel_->isSolutionInteger();
  }

  double getLB() const override {
    return pModel_->getBestLB();
  }

  // build the, possibly fractional, roster corresponding to the solution
  // currently stored in the model
  vector3D<double> getFractionalRoster() const override;

  // build a DualCosts structure
  DualCosts buildDualCosts(PLiveNurse pNurse) const;

  // build a random DualCosts structure
  DualCosts buildRandomDualCosts() const;

  // return the value V used to choose the number of columns on which to branch.
  // Choose as many columns as possible such that: sum (1 - value(column)) < V
  virtual double getBranchColumnValueMax() const {
    return std::max(1.0, pScenario_->nbWeeks_ / 2.0);
  }

  //------------------------------------------------
  // Solution with rolling horizon process
  //------------------------------------------------

  // relax/unrelax the integrality constraints of the variables corresponding
  // to input days
  void relaxDays(const std::vector<bool> &isRelax) override;
  void unrelaxDays(const std::vector<bool> &isUnrelax) override;

  // set availability for the days in fixDays based on
  // the current solution
  void fixAvailabilityBasedOnSolution(
      const std::vector<bool> &fixDays) override;

  // set the available days per nurse
  void nursesAvailabilities(
      const vector3D<bool> &availableNursesDaysShifts) override;

  // remove any column that does not respect the availabilities
  void filterColumnsBasedOnAvailability();

  // fix/unfix all the variables corresponding to the input vector of nurses
  void fixNurses(const std::vector<bool> &isFixNurse) override;
  void unfixNurses(const std::vector<bool> &isUnfixNurse) override;

  // Solve the problem with a method that allows for a warm start
  double rollingSolve(const SolverParam &parameters,
                      int firstDay,
                      const std::vector<Roster> &solution = {}) override;

  // Special solve function for LNS
  // It is a priori the same as a regular, but it might be modified if needed
  double LNSSolve(const SolverParam &parameters,
                  const std::vector<Roster> &solution = {}) override;

  // Initialization of the master problem with/without solution
  void initialize(const SolverParam &parameters,
                  const std::vector<Roster> &solution = {}) override;

  // STAB: compute the lagrangian bound
  bool lagrangianBoundAvailable() const { return lagrangianBoundAvail_; }
  virtual double computeLagrangianBound(double objVal) const;
  virtual double computeApproximateDualUB(double objVal) const;

  /*
  * Solving parameter doubles
  */
  const char *PB_NAME = "GenCol";
  int solvingTime;

  const vector3D<MyVar *> &getOptDemandVars() { return optDemandVars_; }

  const vector4D<MyVar *> &getSkillsAllocVars() { return skillsAllocVars_; }

  const std::vector<int> &getPositionsForSkill(int sk) const {
    return positionsPerSkill_[sk];
  }

  /* Display functions */
  std::string costsConstrainstsToString() const override;
  std::string allocationToString(bool printInteger = true) const;
  std::string coverageToString(bool printInteger = true) const;

 protected:
  Modeler *pModel_;
  vector2D<int> positionsPerSkill_;  // link positions to skills
  vector2D<int> skillsPerPosition_;  // link skills to positions
  MyPricer *pPricer_;  // prices the rotations
  MyTree *pTree_;  // store the tree information
  MyBranchingRule *pRule_;  // choose the variables on which we should branch
  MySolverType solverType_;  // which solver is used
  bool lagrangianBoundAvail_ = false;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // IMPORTANT:  WORKED DAYS CAN ALSO BE USED WORKED HOURS
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
  * Variables
  */

  // count the number of missing nurse to reach the optimal
  vector3D<MyVar *> optDemandVars_;
  // count the number of nurses by position on each day, shift
  vector3D<MyVar *> numberOfNursesByPositionVars_;
  // makes the allocation of the skills
  vector4D<MyVar *> skillsAllocVars_;

  /*
  * Constraints
  */

// ensure a minimal coverage per day, per shift, per skill
  vector3D<MyCons *> minDemandCons_;
  // count the number of missing nurse to reach the optimal
  vector3D<MyCons *> optDemandCons_;
  // ensure there are enough nurses for numberOfNursesByPositionVars_
  vector3D<MyCons *> numberOfNursesByPositionCons_;
  // ensures that each nurse works with the good skill
  vector3D<MyCons *> feasibleSkillsAllocCons_;

  /*
  * Methods
  */

  // Initialize the solver at construction
  void initializeSolver(MySolverType solverType);

  // Main method to build the rostering problem for a given input
  virtual void build(const SolverParam &parameters);

  // Provide an initial solution to the solver
  virtual void initializeSolution(const std::vector<Roster> &solution) = 0;

  // solve method to catch execption
  void solveWithCatch();

  // solve a solution in the output
  void storeSolution() override;

  // set parameters and update printFuntion pointer with this
  void setParameters(const SolverParam &param) {
    pModel_->setParameters(param, this);
    param_ = pModel_->getParameters();
  }

  // return the costs of all active columns associated to the type
  virtual double getColumnsCost(CostType costType,
                                bool justHistoricalCosts) const = 0;

  virtual double getMinDaysCost() const = 0;
  virtual double getMaxDaysCost() const = 0;
  virtual double getMaxWeekendCost() const = 0;

  // update the demand with a new one of the same size
  // change the rhs of the constraints minDemandCons_ and optDemandCons_
  void updateDemand(PDemand pDemand);

  /* Build each set of constraints
   * Add also the coefficient of a column for each set
   */
  void buildSkillsCoverageCons(const SolverParam &parameters);
  int addSkillsCoverageConsToCol(std::vector<MyCons *> *cons,
                                 std::vector<double> *coeffs,
                                 const Pattern &pat) const;

  /* retrieve the dual values */
  virtual vector2D<double> getShiftsDualValues(PLiveNurse pNurse) const;
  virtual std::vector<double> getStartWorkDualValues(PLiveNurse pNurse) const;
  virtual std::vector<double> getEndWorkDualValues(PLiveNurse pNurse) const;
  virtual double getWorkedWeekendDualValue(PLiveNurse pNurse) const;
  virtual double getConstantDualvalue(PLiveNurse pNurse) const;

  //---------------------------------------------------------------------------
  //
  // Methods required to implement stabilization in the column generation
  //
  //---------------------------------------------------------------------------

 public:
  // STAB
  // Update the stabilization variables based on the dual solution
  // 1- When the column generation has found a new columns:
  //     - increase the cost of the duals (the bound for the primal)
  //     - increase the bound of the duals (the cost for the primal) if the
  //        the dual lays outside otherwise reduce it.
  // 2- When the column generation has not found new columns:
  //     - recenter the bounds on the current duals
  //     - keep the same range for the bounds of the duals
  //     - reduce the cost of the duals
  void stabUpdate(OsiSolverInterface *solver, bool recenter = true);

  // STAB
  // initialize the stabilization variables and center them on the current duals
  void stabInitializeBoundAndCost(OsiSolverInterface *solver);

  // STAB
  // deactivate the stabilization variables
  void stabDeactivateBoundAndCost(OsiSolverInterface *solver);

  // STAB
  // Check the stopping criterion of the relaxation solution specific to the
  // the stabilization
  // The point is that current solution can be infeasible if  stabilization
  // variables are non zero
  bool stabCheckStoppingCriterion() const;

  // STAB
  // return the current cost of the stabilization variables
  double getStabCost() const;

 private:
  std::vector<MyVar *> stabVariables_;
  // constraints associated with each two stab variables
  std::vector<MyCons *> stabConstraints_;
  std::vector<double> stabBoxCenters_;

 protected:
  // add stabilization variables:
  // dual: obj += - c_+ z_+ - c_- z_-, s.t.: b_- - z_-<= Pi <= b_+ + z_+
  // primal: obj += -b_- y_- + b_+ y_+, s.t.: y_- <= c_-, y_+ <= c_+
  // if primal constraint is <= -> create just plus var
  // if primal constraint is >= -> create just minus var
  // WARNING: they are inactive at the beginning
  void addStabVariables(const SolverParam &param,
                        const char *name,
                        MyCons *cons,
                        bool LECons,
                        bool GECons);

  // STAB
  // Multiply the upper bound of the input variable by the input factor
  void multiplyUbInSolver(MyVar *pVar,
                          OsiSolverInterface *solver,
                          double factor);
  // Set the bound of the input variable to the input value
  void updateVarUbInSolver(MyVar *pVar,
                           OsiSolverInterface *solver,
                           double value);

  // STAB
  // Set the cost of the input variable to the input value
  void updateVarCostInSolver(MyVar *pVar,
                             OsiSolverInterface *solver,
                             double value);
};

#endif  // SRC_SOLVERS_MP_MASTERPROBLEM_H_
