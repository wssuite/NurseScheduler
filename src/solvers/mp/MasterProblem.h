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
#include "solvers/mp/modeler/Stabilization.h"
#include "solvers/mp/constraints/ConstraintsMP.h"


// forward declaration
class MasterProblem;

//---------------------------------------------------------------------------
//
// C l a s s   R C S o l u t i o n
//
// Contains a solution of the RCSPP
//
//---------------------------------------------------------------------------
struct RCSolution : public Stretch {
  explicit RCSolution(int firstDay = -1,
                      std::vector<PShift> pShifts = {},
                      double cost = DBL_MAX,
                      double reducedCost = DBL_MAX) :
      Stretch(firstDay, std::move(pShifts)),
      cost_(cost),
      reducedCost_(reducedCost) {}

  RCSolution(Stretch stretch, double cost, double reducedCost) :
      Stretch(std::move(stretch)),
      cost_(cost),
      reducedCost_(reducedCost) {}

  std::string toString() const override;

  double cost() const { return cost_; }

  double reducedCost() const { return reducedCost_; }

  void addCost(double c, CostType t)  {
    cost_ += c;
    costs_[t] += c;
  }

  double costByType(CostType t) const {
    if (t == ROTATION_COST)
      return costs_.at(CONS_SHIFTS_COST) + costs_.at(CONS_WORK_COST) +
          costs_.at(PREFERENCE_COST) + costs_.at(COMPLETE_WEEKEND_COST);
    else
      return costs_.at(t);
  }

  std::string costsToString() const;

  void resetCosts() {
    cost_ = 0;
    costs_.clear();
    for (int t=CONS_SHIFTS_COST; t <= TOTAL_WEEKEND_COST; t++)
      costs_[(CostType)t] = 0;
  }

  // Compare rotations on cost
  static bool compareCost(const RCSolution &sol1, const RCSolution &sol2);

  // Compare rotations on dual cost
  static bool compareReducedCost(
      const RCSolution &sol1, const RCSolution &sol2);

  static void sort(std::vector<RCSolution> *solutions) {
    std::stable_sort(solutions->begin(), solutions->end(),
                     [](const RCSolution &sol1, const RCSolution &sol2) {
                       return sol1.reducedCost() < sol2.reducedCost();
                     });
  }

 protected:
  double cost_, reducedCost_;
  std::map<CostType, double> costs_;
};


struct Pattern;
typedef std::shared_ptr<Pattern> PPattern;

struct DualCosts;

struct Pattern : public RCSolution {
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
  Pattern(): nurseNum_(-1) {}

  Pattern(RCSolution rcSol,
          int nurseNum) :
          RCSolution(std::move(rcSol)),
          nurseNum_(nurseNum) {}

  Pattern(MyVar *var, const PScenario &pScenario);

  virtual ~Pattern() = default;

  virtual bool equals(const PPattern &pat) const {
    if (nurseNum_ != pat->nurseNum_) return false;
    return !(*this != *pat);
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

    int commomFirstDay = std::max(firstDay(), pat->firstDay()),
        commomLastDay = std::min(lastDay(), pat->lastDay());
    for (int k = commomFirstDay; k <= commomLastDay; ++k)
      if (shift(k) == pat->shift(k)) return false;

    return true;
  }

  // when branching on this pattern,
  // this method add the corresponding forbidden shifts to the set
  virtual void addForbiddenShifts(
      std::set<std::pair<int, int> > *forbidenShifts,
      int nbShifts,
      PDemand pDemand) const = 0;

  std::string toString() const override;

  // need to be able to write the pattern as a vector and
  // to create a new one from it
  virtual std::vector<double> getCompactPattern() const {
    std::vector<double> pattern;
    pattern.reserve(nDays()+10);
    pattern.push_back(nurseNum_);
    pattern.push_back(firstDay());
    pattern.push_back(nDays());
    for (const PShift& pS : pShifts())
      pattern.push_back(pS->id);
    pattern.push_back(duration());
    return pattern;
  }

  //  Compute the dual cost of a column
  virtual void checkReducedCost(const DualCosts &dualCosts,
                                bool printBadPricing = true) const = 0;

  int nurseNum() const { return nurseNum_; }

  // need to redefine this function because of the static functions
  // having the same name
  int firstDay() const override { return Stretch::firstDay(); }

  int lastDay() const override { return Stretch::lastDay(); }

  int nDays() const override { return Stretch::nDays(); }

  int shift(int day) const override { return Stretch::shift(day); }

  int duration() const override { return Stretch::duration(); }

 protected:
  const int nurseNum_;
};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------
class RCPricer;

class MasterProblem : public Solver, public PrintSolution {
 public:
  // Specific constructor and destructor
  MasterProblem(PScenario pScenario,
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

  void computePatternCost(Pattern *pat) const;

  CostType resourceCostType(const PResource &pR, const PLiveNurse &pN) const {
    return pResourceCostTypes_[pN->num_].at(pR);
  }

  // throw an error if pattern is already present as an active column
  void checkIfPatternAlreadyPresent(const Pattern &pat) const;

  virtual std::vector<MyVar *> getRestVarsPerDay(
      PLiveNurse pNurse, int day) const = 0;

  // get the pointer to the model
  Modeler * pModel() const {
    return pModel_;
  }

  MyPricer * pPricer() const {
    return pModel_->pPricer();
  }

  MyTree * pTree() const {
    return pModel_->pTree();
  }

  MyBranchingRule *pBranchingRule() const {
    return pModel_->pBranchingRule();
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

  // return the value V used to choose the number of columns on which to branch.
  // Choose as many columns as possible such that: sum (1 - value(column)) < V
  virtual double getBranchColumnValueMax() const = 0;

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

  // fix/unfix all the variables corresponding to the input vector of nurses
  void fixNurses(const std::vector<bool> &isFixNurse) override;
  void unfixNurses(const std::vector<bool> &isUnfixNurse) override;

  // remove any column that does not respect the availabilities
  void filterInitialColumnsBasedOnAvailability();

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

  const vector3D<MyVar *> &getOptDemandVars() const {
    return optDemandConstraint_->getVariables();
  }

  const vector3D<MyVar *> &getNursePositionCountVars() const {
    return nursePositionCountConstraint_->getVariables();
  }

  const vector4D<MyVar *> &getSkillsAllocVars() const {
    return allocationConstraint_->getVariables();
  }

  const vector3D<MyCons *> &getOptDemandCons() const {
    return optDemandConstraint_->getConstraints();
  }

  /* Display functions */
  std::string costsConstrainstsToString() const override;
  std::string allocationToString(bool printInteger = true) const;
  std::string coverageToString(bool printInteger = true) const;

  vector<PResource> getSPResources(const PLiveNurse &pN) const {
    return spResources_[pN->num_];
  }

 protected:
  Modeler *pModel_;
  RCPricer *pRCPricer_;
  SolverType solverType_;  // which solver is used
  bool lagrangianBoundAvail_ = false;

  // PResources for the problem and subproblems
  vector<std::map<PResource, CostType>> pResources_;  // per  nurse
  vector2D<PResource> spResources_;  // per nurse

  // create the resources used for the sub problem
  void createPResources();

  // split the resources between the master and the subproblem
  // must initialize spResources_
  virtual void splitPResources() = 0;

  // Functions to generate the default resources for a given nurse
  virtual std::map<PResource, CostType>
  defaultGeneratePResources(const PLiveNurse &pN) const;

  // create resources and the map that associates a PRessource to its cost type
  // create resources and the map that associates a PRessource to its cost type
  virtual std::map<PResource, CostType>
  generatePResources(const PLiveNurse &pN) {
    if (pResourceCostTypes_[pN->num_].empty()) {
      // if defined, do not use the default function
      if (generatePResourcesFunc_)
        pResourceCostTypes_[pN->num_] = generatePResourcesFunc_(pN);
      else
        pResourceCostTypes_[pN->num_] = defaultGeneratePResources(pN);
    }
    return pResourceCostTypes_[pN->num_];
  }

  vector<std::map<PResource, CostType>> pResourceCostTypes_;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // IMPORTANT:  WORKED DAYS CAN ALSO BE USED WORKED HOURS
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*
  * Master Problem Constraints
  */
  NursePositionCountConstraint *nursePositionCountConstraint_;
  AllocationConstraint *allocationConstraint_;
  DemandConstraint *minDemandConstraint_, *optDemandConstraint_;

  // vector containing the constraints that are involved
  // in the column generation
  vector<ConstraintMP*> columnConstraints_;

 public:
  const vector<ConstraintMP*> & columnConstraints() const {
    return columnConstraints_;
  }

  void addColumnConstraint(ConstraintMP* pC) {
    columnConstraints_.push_back(pC);
  }

  // add the column to the problem
  MyVar *createColumn(const Pattern &col, const char *baseName);

  // add a given constraint to the column
  void addConstoCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Pattern &col) const;

 protected:
  /*
  * Methods
  */

  // Main method to build the rostering problem for a given input
  virtual void build(const SolverParam &parameters);

  // Provide an initial solution to the solver
  virtual void initializeSolution(const std::vector<Roster> &solution) = 0;

  // solve method to catch exception
  void solveWithCatch();

  // solve a solution in the output
  void storeSolution() override;

  // set parameters and update printFunction pointer with this
  void setParameters(const SolverParam &param) {
    pModel_->setParameters(param, this);
    param_ = pModel_->getParameters();
  }

  // return the costs of all active columns associated to the type
  // return the costs of all active columns associated to the type
  virtual double getColumnsCost(CostType costType) const;
  double getColumnsCost(CostType costType,
                        const std::vector<MyVar *> &vars) const;

  virtual double getDaysCost() const = 0;
  virtual double getWeekendCost() const = 0;

  // update the demand with a new one of the same size
  // change the rhs of the constraints minDemandCons_ and optDemandCons_
  void updateDemand(PDemand pDemand);
};

#endif  // SRC_SOLVERS_MP_MASTERPROBLEM_H_
