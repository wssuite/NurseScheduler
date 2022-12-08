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
// C l a s s   C o l u m n
//
// Data of the variables of the column generation algorithm
//
//---------------------------------------------------------------------------
struct DualCosts;

struct Column : public RCSolution {
  /* Static helpers */
  static int nurseNum(const std::vector<double> &column) {
    return static_cast<int>(column[0]);
  }
  static int firstDay(const std::vector<double> &column) {
    return static_cast<int>(column[1]);
  }
  static int lastDay(const std::vector<double> &column) {
    return Column::firstDay(column) + Column::nDays(column) - 1;
  }
  static int nDays(const std::vector<double> &column) {
    return static_cast<int>(column[2]);
  }
  static int shift(const std::vector<double> &column, int d) {
    return static_cast<int>(column[3+d-Column::firstDay(column)]);
  }
  static int duration(const std::vector<double> &column) {
    return static_cast<int>(column[3+Column::nDays(column)]);
  }

  static int nurseNum(MyVar *var) {
    return Column::nurseNum(var->getCompactColumn());
  }
  static int firstDay(MyVar *var) {
    return Column::firstDay(var->getCompactColumn());
  }
  static int lastDay(MyVar *var) {
    return Column::lastDay(var->getCompactColumn());
  }
  static int nDays(MyVar *var) {
    return Column::nDays(var->getCompactColumn());
  }
  static int shift(MyVar *var, int d) {
    return Column::shift(var->getCompactColumn(), d);
  }
  static int duration(MyVar *var) {
    return Column::duration(var->getCompactColumn());
  }

  /* class functions */
  Column(): nurseNum_(-1) {}

  Column(RCSolution rcSol, int nurseNum) :
      RCSolution(std::move(rcSol)),
      nurseNum_(nurseNum) {}

  Column(MyVar *var, const PScenario &pScenario);

  virtual bool equals(const std::shared_ptr<Column> &pat) const {
    if (nurseNum_ != pat->nurseNum_) return false;
    return !(*this != *pat);
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isDisjointWith(
      const std::shared_ptr<Column> &pat, bool needRest) const {
    return ((lastDayId() < pat->firstDayId() - needRest)
        || (pat->lastDayId() < firstDayId() - needRest));
  }

  // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest)
  virtual bool isShiftDisjointWith(
      const std::shared_ptr<Column> &pat, bool needRest) const {
    if (isDisjointWith(pat, needRest))
      return true;

    int commonFirstDay = std::max(firstDayId(), pat->firstDayId()),
        commonLastDay = std::min(lastDayId(), pat->lastDayId());
    for (int k = commonFirstDay; k <= commonLastDay; ++k)
      if (shift(k) == pat->shift(k)) return false;

    return true;
  }

  // when branching on this column,
  // this method add the corresponding forbidden shifts to the set
  virtual void addForbiddenShifts(
      std::set<std::pair<int, int> > *forbidenShifts,
      int nbShifts,
      PDemand pDemand) const = 0;

  std::string toString() const override;

  // need to be able to write the column as a vector and
  // to create a new one from it
  virtual std::vector<double> getCompactColumn() const {
    std::vector<double> column;
    column.reserve(nDays()+10);
    column.push_back(nurseNum_);
    column.push_back(firstDayId());
    column.push_back(nDays());
    for (const PShift& pS : pShifts())
      column.push_back(pS->id);
    column.push_back(duration());
    return column;
  }

  //  Compute the dual cost of a column
  virtual void checkReducedCost(const DualCosts &dualCosts,
                                bool printBadPricing) const = 0;

  int nurseNum() const { return nurseNum_; }

  // need to redefine this function because of the static functions
  // having the same name
  int firstDayId() const override { return Stretch::firstDayId(); }

  int lastDayId() const override { return Stretch::lastDayId(); }

  int nDays() const override { return Stretch::nDays(); }

  int shift(int day) const override { return Stretch::shift(day); }

  int duration() const override { return Stretch::duration(); }

 protected:
  const int nurseNum_;
};

typedef std::shared_ptr<Column> PColumn;


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
  MasterProblem(const PScenario& pScenario,
                SolverType solver);
  ~MasterProblem() override;

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

  // retrieve the object represented ny the  vector column
  virtual PColumn getPColumn(MyVar *var) const = 0;

  void computeColumnCost(Column *col) const;

  // throw an error if column is already present as an active column
  bool checkIfColumnAlreadyPresent(
      const Column &column, bool printErr = false) const;

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
  void saveSolution() override;
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
      const std::vector<bool> &fixDays,
      const std::vector<Roster> &solution = {}) override;

  // fix/unfix all the variables corresponding to the input vector of nurses
  void fixNurses(const std::vector<bool> &isFixNurse) override;
  void unfixNurses(const std::vector<bool> &isUnfixNurse) override;

  // remove any column that does not respect the availabilities
  void filterInitialColumnsBasedOnAvailability();

  // Solve the problem with a method that allows for a warm start
  double rollingSolve(const SolverParam &parameters,
                      int firstDay,
                      const std::vector<Roster> &solution) override;

  // Special solve function for LNS
  double LNSSolve(const SolverParam &parameters) override;

  // Initialization of the master problem
  void initialize(const SolverParam &parameters);

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

  // solve a solution in the output
  // return false if fails, true otherwise
  bool storeSolution() override;

  /* Display functions */
  std::string costsConstraintsToString() const override;
  std::string costsColumnsToString() const;
  string allocationToString() const;
  string coverageToString() const;

  vector<PResource> getSPResources(const PLiveNurse &pN) const {
    return spResources_[pN->num_];
  }

  void addNewSPResources(const PLiveNurse &pN, const PResource &pR) {
    pR->setId(static_cast<int>(spResources_[pN->num_].size()));
    spResources_[pN->num_].push_back(pR);
    pResources_[pN->num_].push_back(pR);
  }

 protected:
  Modeler *pModel_;
  RCPricer *pRCPricer_;
  bool lagrangianBoundAvail_ = false;

  // PResources for the problem and subproblems
  vector2D<PResource> pResources_, spResources_;  // per nurse

  // create the resources used for the sub problem
  void createPResources();

  void addPResources();

  // split the resources between the master and the subproblem
  // must initialize spResources_
  virtual void splitPResources() = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // IMPORTANT:  WORKED DAYS CAN ALSO BE USED WORKED HOURS
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*
  * Master Problem Constraints
  */
  NursePositionCountConstraint *nursePositionCountConstraint_{};
  AllocationConstraint *allocationConstraint_{};
  DemandConstraint *minDemandConstraint_{}, *optDemandConstraint_{};

  // vector containing the constraints that are involved
  // in the column generation
  vector<ConstraintMP*> constraints_, columnConstraints_;

 public:
  const vector<ConstraintMP*> & columnConstraints() const {
    return columnConstraints_;
  }

  void addConstraint(ConstraintMP* pC, bool affectColumns) {
    constraints_.push_back(pC);
    if (affectColumns) columnConstraints_.push_back(pC);
  }

  // add the column to the problem
  MyVar *createColumn(const Column &col, const char *baseName);

  // add a given constraint to the column
  void addConsToCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Column &col) const;

 protected:
  /*
  * Methods
  */

  // Main method to build the rostering problem for a given input
  virtual void build(const SolverParam &parameters);

  // Provide an initial solution to the solver
  virtual void initializeSolution(const std::vector<Roster> &solution) = 0;

  // solve method to catch exception
  void solveWithCatch(const vector<Roster> &solution);

  // set parameters and update printFunction pointer with this
  void setParameters(const SolverParam &param) override {
    pModel_->setParameters(param, this);
    param_ = pModel_->getParameters();
  }

  // return the costs of all active columns associated to the type
  virtual std::map<CostType, double> getColumnsCosts() const;
  std::map<CostType, double> getColumnsCosts(
      const std::vector<MyVar *> &vars) const;

  // update the demand with a new one of the same size
  // change the rhs of the constraints minDemandCons_ and optDemandCons_
  virtual void update(const PDemand& pDemand);
};

#endif  // SRC_SOLVERS_MP_MASTERPROBLEM_H_
