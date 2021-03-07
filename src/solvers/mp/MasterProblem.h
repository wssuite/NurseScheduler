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

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>


#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/Solver.h"
#include "tools/Tools.h"
#include "data/Nurse.h"
#include "solvers/mp/modeler/Modeler.h"
#include "OsiSolverInterface.hpp"

enum CostType {
  ROTATION_COST, CONS_SHIFTS_COST, CONS_WORK_COST,
  COMPLETE_WEEKEND_COST, PREFERENCE_COST, CONS_REST_COST,
  DAYS_COST, WEEKEND_COST
};

//---------------------------------------------------------------------------
//
// C l a s s   R C S o l u t i o n
//
// Contains a solution of the RCSPP
//
//---------------------------------------------------------------------------
// TODO(AL): if shifts and shift types do not match, the vector of shifts of
//  a solution should not contain shift types
struct RCSolution {
  RCSolution(int firstDay, const std::vector<int> &shifts, double c) :
      firstDay(firstDay), shifts(shifts), cost(c) {}
  RCSolution() : firstDay(-1), cost(0) {}

  int firstDay;
  std::vector<int> shifts;
  double cost;

  std::string toString(std::vector<int> shiftIDToShiftTypeID = {}) const;
};

struct DualCosts {
 public:
  DualCosts(vector2D<double> workedShiftsCosts,
            double constant) :
      workedShiftsCosts_(std::move(workedShiftsCosts)),
      constant_(constant) {}
  virtual ~DualCosts() = default;

  // GETTERS
  //
  virtual int nDays() const {
    return workedShiftsCosts_.size();
  }
  virtual double workedDayShiftCost(int day, int shift) const {
    return shift == 0 ? 0 : workedShiftsCosts_[day][shift - 1];
  }
  virtual double startWorkCost(int day) const { return 0; }
  virtual double endWorkCost(int day) const { return 0; }
  virtual double workedWeekendCost() const { return 0; }
  virtual double constant() const { return constant_; }

 protected:
  // TODO(JO): this indexation is VERY confusing and prone to error. Since we
  //  have decided to treat rest shifts as any shift, we should never have to
  //  assume that the rest shift is the one with index 0 (except in the
  //  isWork()/isRest() functions. And we could still include a cost for the
  //  rest shift in this vector, with a zero value.
  //  AL: indeed, but at the same time, it's protected and
  //  the workedDayShiftCost method takes care of that. The rest shift not 0
  //  would be fixed another time
  // Indexed by : (day, shift) !! 0 = shift 1 !!
  vector2D<double> workedShiftsCosts_;

  // constant part of the reduced cost
  double constant_;
};
typedef std::shared_ptr<DualCosts> PDualCosts;

struct Pattern;
typedef std::shared_ptr<Pattern> PPattern;

class MasterProblem;

struct Pattern {
  /* Static helpers */
  static int nurseNum(const std::vector<double> &pattern) {
    return static_cast<int>(pattern[0]);
  }
  static int firstDay(const std::vector<double> &pattern) {
    return static_cast<int>(pattern[1]);
  }
  static int lastDay(const std::vector<double> &pattern) {
    return Pattern::firstDay(pattern) + Pattern::nDays(pattern) - 1;
  }
  static int nDays(const std::vector<double> &pattern) {
    return static_cast<int>(pattern[2]);
  }
  static int shift(const std::vector<double> &pattern, int d) {
    return static_cast<int>(pattern[3+d-Pattern::firstDay(pattern)]);
  }
  static int duration(const std::vector<double> &pattern) {
    return static_cast<int>(pattern[3+Pattern::nDays(pattern)]);
  }

  static int nurseNum(MyVar *var) {
    return Pattern::nurseNum(var->getPattern());
  }
  static int firstDay(MyVar *var) {
    return Pattern::firstDay(var->getPattern());
  }
  static int lastDay(MyVar *var) {
    return Pattern::lastDay(var->getPattern());
  }
  static int nDays(MyVar *var) {
    return Pattern::nDays(var->getPattern());
  }
  static int shift(MyVar *var, int d) {
    return Pattern::shift(var->getPattern(), d);
  }
  static int duration(MyVar *var) {
    return Pattern::duration(var->getPattern());
  }

  /* class functions */
  Pattern(int firstDay,
          const std::vector<int> &shifts,
          const PScenario &pScenario,
          int nurseNum = -1,
          double cost = DBL_MAX,
          double dualCost = DBL_MAX);

  Pattern(const std::map<int, int> &shifts,
          const PScenario &pScenario,
          int nurseNum = -1,
          double cost = DBL_MAX,
          double dualCost = DBL_MAX);

  Pattern(const Pattern &pat, int nurseNum) :
      nurseNum_(nurseNum),
      stretch_(pat.stretch_),
      cost_(pat.cost_),
      costs_(pat.costs_),
      reducedCost_(pat.reducedCost_) {
    if (pat.nurseNum_ != nurseNum_) {
      cost_ = DBL_MAX;
      reducedCost_ = DBL_MAX;
    }
  }

  Pattern(const std::vector<double> &pattern, const PScenario &pScenario);

  virtual ~Pattern() = default;

  virtual bool equals(PPattern pat) const {
    if (nurseNum_ != pat->nurseNum_) return false;
    return !(stretch_ != pat->stretch_);
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isDisjointWith(PPattern pat, bool needRest = true) const {
    return ((lastDay() < pat->firstDay() - needRest)
        || (pat->lastDay() < firstDay() - needRest));
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isShiftDisjointWith(PPattern pat, bool needRest = true) const {
    if (isDisjointWith(pat, needRest))
      return true;

    int commomFirstDay =
        std::max(stretch_.firstDay(), pat->stretch_.firstDay()),
        commomLastDay =
        std::min(stretch_.lastDay(), pat->stretch_.lastDay());
    for (int k = commomFirstDay; k <= commomLastDay; ++k)
      if (shift(k) == pat->shift(k)) return false;

    return true;
  }

  virtual std::string toString(int nbDays = -1) const;

  virtual int shift(int day) const {
    return pShift(day)->id;
  }

  virtual const PShift& pShift(int day) const {
    return stretch_.pShift(day-firstDay());
  }

  // when branching on this pattern,
  // this method add the corresponding forbidden shifts to the set
  virtual void addForbiddenShifts(
      std::set<std::pair<int, int> > *forbidenShifts,
      int nbShifts,
      PDemand pDemand) const = 0;

  // compute the cost of a pattern based on the master
  virtual void computeCost(const MasterProblem *pMaster,
                           const PLiveNurse &pNurse) = 0;

  void computeResourcesCosts(const MasterProblem *pMaster,
                             const State &initialState);

  void addCost(double c, CostType t)  {
    cost_ += c;
    costs_[t] += c;
  }

  double cost(CostType t) const {
    if (t == ROTATION_COST)
      return costs_.at(CONS_SHIFTS_COST) + costs_.at(CONS_WORK_COST) +
          costs_.at(PREFERENCE_COST) + costs_.at(COMPLETE_WEEKEND_COST);
    else
      return costs_.at(t);
  }

  std::string costsToString() const;

  // need to be able to write the pattern as a vector and
  // to create a new one from it
  virtual std::vector<double> getCompactPattern() const {
    std::vector<double> pattern;
    pattern.reserve(stretch_.nDays()+10);
    pattern.push_back(nurseNum_);
    pattern.push_back(stretch_.firstDay());
    pattern.push_back(stretch_.nDays());
    for (const PShift& pS : stretch_.pShifts())
      pattern.push_back(pS->id);
    pattern.push_back(stretch_.duration());
    return pattern;
  }

  int firstDay() const { return stretch_.firstDay(); }

  int lastDay() const { return stretch_.lastDay(); }

  int nDays() const { return stretch_.nDays(); }

  int duration() const { return stretch_.duration(); }

  const int nurseNum_;
  const Stretch stretch_;

  // Cost
  double cost_;
  std::map<CostType, double> costs_;

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
                SolverType solver);
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

  // define the resources used for the sub problem
  std::vector<PResource> createResources(const PLiveNurse &pN) const {
    return createResources(pN, nullptr);
  }
  virtual std::vector<PResource>
  createResources(const PLiveNurse &pN,
                  std::map<int, CostType> *resourceCostType) const = 0;

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

  double LB() const override {
    return pModel_->getBestLB();
  }

  // build the, possibly fractional, roster corresponding to the solution
  // currently stored in the model
  vector3D<double> fractionalRoster() const override;

  // build a DualCosts structure
  virtual PDualCosts buildDualCosts(PLiveNurse pNurse) const;

  // build a random DualCosts structure
  // If the option 'optimalDemandConsidered' is selected, the dual costs will
  // randomly represent the optimal demand over all the horizon. Otherwise,
  // the dual costs will be created totally randomly.
  virtual PDualCosts buildRandomDualCosts(
      bool optimalDemandConsidered = false, int NDaysShifts = 10) const;

  vector2D<double> getRandomWorkedDualCosts(
      bool optimalDemandConsidered, int NDaysShifts) const;

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

  // Compute an approximation of the dual UB based on the lagrangian bound
  // It could be useful to measure the quality of a dual solution (used when
  // stabilizing).
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
  SolverType solverType_;  // which solver is used
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

  // STAB
  // Stabilization variables for each constraint
  // Two variables are needed for equality constraints and one for inequalities
  // The constraints on average values are not stabilized yet
  // The position and allocation constraints do not require stabilization

  // ensure a minimal coverage per day, per shift, per skill
  vector3D<MyVar *> stabMinDemandPlus_;
  // count the number of missing nurse to reach the optimal
  vector3D<MyVar *> stabOptDemandPlus_;

  /*
  * Methods
  */

  // Initialize the solver at construction
  void initializeSolver(SolverType solverType);

  // Main method to build the rostering problem for a given input
  virtual void build(const SolverParam &parameters);

  // Provide an initial solution to the solver
  virtual void initializeSolution(const std::vector<Roster> &solution) = 0;

  // solve method to catch execption
  void solveWithCatch();

  // solve a solution in the output
  void storeSolution() override;

  // set parameters and update printFunction pointer with this
  void setParameters(const SolverParam &param) {
    pModel_->setParameters(param, this);
    param_ = pModel_->getParameters();
  }

  // return the costs of all active columns associated to the type
  virtual double getColumnsCost(CostType costType) const = 0;

  virtual double getDaysCost() const = 0;
  virtual double getWeekendCost() const = 0;

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
  virtual double getConstantDualvalue(PLiveNurse pNurse) const;

  //---------------------------------------------------------------------------
  //
  // STAB: Methods required to implement stabilization in the column generation
  //
  // Ref: Lübbecke, Marco E., and Jacques Desrosiers.
  // "Selected topics in column generation."
  // Operations research 53.6 (2005): 1007-1023.
  //
  //---------------------------------------------------------------------------

 public:
// STAB
// Update the stabilization variables based on the dual solution
// 1- When the dual lays inside the box:
//     - increase the penalty of the duals (the bound for the primal)
//     - decrease the radius of the duals (the cost for the primal).
// 2- When the dual lays outside the box:
//     - decrease the penalty of the duals (the bound for the primal)
//     - increase the radius of the duals (the cost for the primal).
// When a dual solution (of the original problem) of better quality
// is obtained, recenter the box.
// The issue here is that the  dual solution is not available as the lagrangian
// bound needs to be computed (and available) and all sub problems need to
// have been solved to optimality.
// Instead, the solution is recenter when asked (recenter=true).
// Currently, the box is recentered when no more columns are generated.
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
  std::vector<MyVar *> stabVariablesPlus_, stabVariablesMinus_;
  // constraints associated with each two stab variables
  std::vector<MyCons *> stabConstraints_;
  std::vector<double> stabBoxCenters_;

 protected:
  // Add stabilization variables z for the box [b_, b+] with the penalties c
  // if getting outside of the box:
  // dual = obj += - c_+ z_+ - c_- z_-, s.t.: b_- - z_-<= Pi <= b_+ + z_+
  // primal = obj += -b_- y_- + b_+ y_+, s.t.: y_- <= c_-, y_+ <= c_+
  // if primal constraint is <= -> create just minus var
  // if primal constraint is >= -> create just plus var
  // WARNING: they are inactive at the beginning
  //
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
