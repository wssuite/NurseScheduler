/*
*Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
*All Rights Reserved.
 *
*You may use, distribute and modify this code under the terms of the MIT
*license.
 *
*Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
*full license detail.
 */

#include "MasterProblem.h"

#include <list>

#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

/* namespace usage */
using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;


// P a t t e r n   s t a t i c   m e t h o d s
std::vector<PShift> getPShifts(const std::vector<double> &column,
                               const PScenario &pScenario) {
  int length = Column::nDays(column);
  std::vector<PShift> pShifts(length);
  for (int k = 0; k < length; k++)
    pShifts[k] = pScenario->pShift(static_cast<int>(column[k + 3]));
  return pShifts;
}

Column::Column(MyVar *var, const PScenario &pScenario) :
    RCSolution(Column::firstDay(var->getCompactColumn()),
               getPShifts(var->getCompactColumn(), pScenario),
               var->getCost()),
    nurseNum_(Column::nurseNum(var->getCompactColumn())) {}

std::string Column::toString() const {
  std::stringstream rep;
  rep << "Nurse " << nurseNum_ << " - " << RCSolution::toString();
  return rep.str();
}

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

// Default constructor
MasterProblem::MasterProblem(const PScenario& pScenario,
                             SolverType solverType)  :
    Solver(pScenario),
    PrintSolution(),
    pModel_(nullptr),
    pRCPricer_(nullptr),
    pResources_(pScenario->nNurses()) {
  // build the model
  pModel_ = new BcpModeler(this, PB_NAME, solverType);
  this->preprocessData();
}

MasterProblem::~MasterProblem() {
  delete nursePositionCountConstraint_;
  delete allocationConstraint_;
  delete minDemandConstraint_;
  delete optDemandConstraint_;
  delete pModel_;
}

// solve the rostering problem
double MasterProblem::solve(const vector<Roster> &solution) {
  return solve(solution, true);
}

// Solve the rostering problem with parameters

double MasterProblem::solve(
    const SolverParam &param, const vector<Roster> &solution) {
  setParameters(param);
  return solve(solution, true);
}

// solve the rostering problem
double MasterProblem::solve(const vector<Roster> &solution, bool rebuild) {
  // build the model first
  if (rebuild) {
    pModel_->clear();
    this->build(param_);
    if (param_.isStabilization_) pModel_->stab().initAllStabVariables();
  } else {
    // reset model
    pModel_->reset();
  }

  // input an initial solution
  this->initializeSolution(solution);

  this->solveWithCatch();

  if (pModel_->getParameters().printBranchStats_) {
    pModel_->printStats();
  }

  // if was solving just the relaxation: stopAfterXSolution_ = 0
  if (param_.stopAfterXSolution_ == 0 && pModel()->nbSolutions() == 0) {
    // if column generation has been solved to optimality
    if (status_ != INFEASIBLE)
      status_ = OPTIMAL;

    // print constraints' costs
    std::cout << costsConstraintsToString() << std::endl;

    return pModel_->getRelaxedObjective();
  }

  // print best solution
  pModel_->printBestSol();

  // store solution
  storeSolution();

  // print constraints' costs
  std::cout << costsConstraintsToString() << std::endl;

  return pModel_->getObjective();
}

void MasterProblem::solveWithCatch() {
  try {
    // initialize pricer (useful to authorize/forbid shifts for nurses)
    pPricer()->initNursesAvailabilities();
    // filter out initial columns that do not respect nurses' availabilities
    filterInitialColumnsBasedOnAvailability();
    // then solve
    pModel_->solve();
  } catch (const std::exception &e) {
    std::stringstream buff;
    buff << "MasterProblem::solveWithCatch() caught an exception: "
         << e.what() << std::endl;

    if (strcmp(e.what(), "std::bad_alloc") == 0) {
      double memGB = Tools::getResidentMemoryGB();
      buff << "There is a high probability that the program is OUT OF MEMORY;"
              " it has consumed " << std::setprecision(3) << memGB << " GB."
              << std::endl;
      buff << "You may increase maxRelativeLPGapToKeepChild and "
              "decrease maxLevelDifference to reduce the number of in-memory "
              "branching nodes, and thus the memory used." << std::endl;
      buff << "You may also modify the global search strategy by enabling/"
              "disabling diving, strong branching or/and the heuristic."
              << std::endl;
    }
    if (!pModel_->logfile().empty()) {
      Tools::LogOutput log(pModel_->logfile(), true);
      log << buff.str();
    }
    std::cerr << buff.str();
  }
}

// Resolve the problem with another demand and keep the same preferences
double MasterProblem::resolve(PDemand pDemand,
                              const SolverParam &param,
                              const std::vector<Roster> &solution) {
  setParameters(param);
  update(pDemand);
  return solve(solution, false);
}

// Initialization of the master problem with/without solution
void MasterProblem::initialize(const SolverParam &param,
                               const std::vector<Roster> &solution) {
  build(param);
  initializeSolution(solution);
  setParameters(param);
}

// build the rostering problem
void MasterProblem::build(const SolverParam &param) {
  /* initialize resources */
  createPResources();

  /* Skills coverage constraints */
  nursePositionCountConstraint_ = new NursePositionCountConstraint(this);
  allocationConstraint_ = new AllocationConstraint(this);
  minDemandConstraint_ = new DemandConstraint(this, true);
  if (pDemand_->isOptDemand_)
    optDemandConstraint_ = new DemandConstraint(
        this, false, true, pScenario_->weights().optimalDemand);


  /* Initialize the objects used in the branch and price
   * unless the CBC is used to solve the problem
   */
//  if (solverType_ != CBC) {
    /* Rotation pricer */
    pRCPricer_ = new RCPricer(this, "pricer", param);
    pModel_->addPricer(pRCPricer_);  // transfer ownership

    /* Tree */
    RestTree *pTree =
        new RestTree(pScenario_, pDemand_, param.epsilon_, param.verbose_ > 0);
    pModel_->addTree(pTree);

    /* Branching rule */
    auto *pRule = new DiveBranchingRule(this, pTree, "branching rule");
    pModel_->addBranchingRule(pRule);
//  }
}

void MasterProblem::createPResources() {
  pResources_.clear();
  // create the resources for all nurses
  for (const PLiveNurse &pN : theLiveNurses_)
    pResources_.push_back(pN->pResources());
  // split the resources between the master and the subproblem
  // must initialize spResources_
  splitPResources();
}

//------------------------------------------------
// Solution with rolling horizon process
//------------------------------------------------

// relax/unrelax the integrality constraints of the variables corresponding to
// input days
void MasterProblem::relaxDays(const std::vector<bool> &isRelax) {
  if (isRelaxDay_.empty()) {
    isRelaxDay_.insert(isRelaxDay_.begin(), pDemand_->nDays_, false);
    isPartialRelaxDays_ = true;
  }

  for (int day = 0; day < pDemand_->nDays_; day++) {
    isRelaxDay_[day] = isRelaxDay_[day] || isRelax[day];
  }
}

void MasterProblem::unrelaxDays(const std::vector<bool> &isUnrelax) {
  if (isRelaxDay_.empty()) {
    isPartialRelaxDays_ = false;
  } else {
    for (int day = 0; day < pDemand_->nDays_; day++)
      isRelaxDay_[day] = isRelaxDay_[day] && !isUnrelax[day];
  }
}

// set availability for the days before fixBefore (<=) based on
// the current solution
void MasterProblem::fixAvailabilityBasedOnSolution(
    const std::vector<bool> &fixDays) {
  // check if a solution has been loaded
  if (pModel_->getActiveColumns().empty())
    Tools::throwError("Cannot call fixAvailabilityBasedOnSolution "
                      "if no solutions in the solver");

  // by default set all availabilities to true
  vector3D<bool> availableNursesDaysShifts;
  Tools::initVector3D(&availableNursesDaysShifts,
                      nNurses(), nDays(), nShifts(), true);
  // then get the solution and set to true only the shift belonging to
  // a solution
  vector3D<double> roster = fractionalRoster();
  for (int k = 0; k < nDays(); ++k)
    if (fixDays[k]) {
      for (int n = 0; n < nNurses(); ++n)
        for (int s = 0; s < nShifts(); ++s)
          availableNursesDaysShifts[n][k][s] = roster[n][k][s] > epsilon();
    }

  nursesAvailabilities(availableNursesDaysShifts);
}

// fix/unfix all the variables corresponding to the input vector of nurse ids
void MasterProblem::fixNurses(const std::vector<bool> &isFix) {
  // initialize the list of fixed nurses if empty
  if (isFixNurse_.empty()) {
    Tools::initVector(&isFixNurse_, nNurses(), false);
    isPartialFixNurses_ = true;
  }
  // set the list of fixed day
  // + forbid the generation of rotations of the input nurses
  for (const PLiveNurse& pNurse : theLiveNurses_) {
    int n = pNurse->num_;
    isFixNurse_[n] = isFixNurse_[n] || isFix[n];
    if (isFixNurse_[n]) pPricer()->forbidNurse(n);
  }
}

void MasterProblem::unfixNurses(const std::vector<bool> &isUnfix) {
  // easy to treat the case where no nurse is fixed yet
  if (isFixNurse_.empty()) {
    isPartialFixNurses_ = false;
    pPricer()->clearForbiddenNurses();
  } else {
    // set the list of unfixed nurses
    // + authorize the generation of rotations for the input nurses
    for (const PLiveNurse& pNurse : theLiveNurses_) {
      int n = pNurse->num_;
      isFixNurse_[n] = isFixNurse_[n] && !isUnfix[n];
      if (!isFixNurse_[n]) pPricer()->authorizeNurse(n);
    }
  }
}

// remove any column that does not respect the availabilities
// TODO(legraina): try to keep the column in a pool
//  to perhaps add them back later on
void MasterProblem::filterInitialColumnsBasedOnAvailability() {
  // if all available, all columns can remain in the pool
  if (!isPartialAvailable_) return;

  // iterate through the pool of columns
  std::vector<MyVar *> initialColumns;
  for (MyVar *var : pModel_->getInitialColumns()) {
    // check if column is feasible
    PColumn pat = getPColumn(var);
    for (int k = pat->firstDayId(); k <= pat->lastDayId(); k++)
      if (!isNurseAvailableOnDayShift(pat->nurseNum(), k, pat->shift(k))) {
        // delete infeasible column
        delete var;
        var = nullptr;
        break;
      }
    // if not deleted, keep it
    if (var)
      initialColumns.push_back(var);
  }

  // reset the pool
  pModel_->setInitialColumns(initialColumns);
}


//------------------------------------------------------------------------------
// Solve the problem with a method that allows for a warm start
//------------------------------------------------------------------------------
double MasterProblem::rollingSolve(const SolverParam &param,
                                   int firstDay,
                                   const std::vector<Roster> &solution) {
  // build the model and initialize with artificial columns
  // at the first iteration
  if (firstDay == 0) {
    initialize(param, solution);
  } else {
    // copy active columns (the solution) from previous run to initial and
    // clear solution
    pModel_->copyActiveToInitialColumns();
    // reset the model
    pModel_->reset();
    setParameters(param);
    // add the solution in the model
    initializeSolution(solution);
  }

  // solve the problem
  solveWithCatch();
  pModel_->loadBestSol(false);

  // output information and save the solution
  if (pModel_->getParameters().printBranchStats_)
    pModel_->printStats();

  return pModel_->getObjective();
}

//------------------------------------------------------------------------------
// Solve the problem with a method that can be specific to our implementation
// of LNS
//------------------------------------------------------------------------------
double MasterProblem::LNSSolve(const SolverParam &param,
                               const std::vector<Roster> &solution) {
  // copy active columns (the solution) to initial from previous run and
  // clear solution
  pModel_->copyActiveToInitialColumns();

  // reset the model
  pModel_->reset();
  setParameters(param);
  // add the solution
  initializeSolution(solution);

  // solve the problem
  solveWithCatch();

  // load and store best solution
  pModel_->loadBestSol();
  storeSolution();

  if (pModel_->getVerbosity() >= 1)
    std::cout << costsConstraintsToString() << std::endl;

  // output information and save the solution
  if (pModel_->getParameters().printBranchStats_)
    pModel_->printStats();

  return pModel_->getObjective();
}


//------------------------------------------------------------------------------
// Store the current solution of the master problem in each LiveNurse
//------------------------------------------------------------------------------
bool MasterProblem::storeSolution() {
  // check if a solution if loaded
  if (pModel_->getActiveColumns().empty()) {
    if (!pModel_->loadBestSol()) {
      std::cerr << "No solution available to store." << std::endl;
      return false;
    }
  }

  // retrieve a feasible allocation of skills as the current value of the
  // allocation variables in the master problem
  vector4D<double> skillsAllocation;
  Tools::initVector4D(&skillsAllocation,
                      pDemand_->nDays_,
                      pScenario_->nShifts(),
                      pScenario_->nSkills(),
                      0,
                      .0);


  for (int k = 0; k < pDemand_->nDays_; ++k)
    for (int s = 0; s < pScenario_->nShifts(); ++s)
      for (int sk = 0; sk < pScenario_->nSkills(); ++sk)
        skillsAllocation[k][s][sk] =
            pModel_->getVarValues(getSkillsAllocVars()[k][s][sk]);

  // build the rosters: while doing so, substract the value of the active
  // columns of each nurse to each corresponding allocation in vector
  // skillsAllocation; if there is no error the skills allocation vector
  // should contain only zeros at the end
  solution_.clear();
  for (const PLiveNurse& pNurse : theLiveNurses_) {
    pNurse->roster_.reset(pScenario_->pRestShift());
    pNurse->columns_.clear();
    pNurse->colVals_.clear();
  }

  std::list<MyVar*> activeColumns(
      pModel_->getActiveColumns().begin(), pModel_->getActiveColumns().end());
  int size = static_cast<int>(activeColumns.size());
  int n = 0;
  while (!activeColumns.empty()) {
    MyVar *var = activeColumns.front();
    activeColumns.pop_front();
    ++n;
    double v = pModel_->getVarValue(var);
    if (v > epsilon()) {
      PColumn pCol = getPColumn(var);
      PLiveNurse pNurse = theLiveNurses_[pCol->nurseNum()];
      bool fractional = v < 1 - epsilon();
      // if first round, do not process fractional columns
      // store them for next round
      if (fractional && n < size) {
        activeColumns.push_back(var);
        continue;
      }
      // check that each assignment is valid and initialize the assignment
      // table in the roster_ variable of the LiveNurse
      for (int k = pCol->firstDayId(); k <= pCol->lastDayId(); ++k) {
        if (fractional && !isRelaxDay(k)) {
          std::cerr << "Column has a fractional value (" << v
                    << ") while it should be integer: " << std::endl
                    << pCol->toString();
          return false;
        }
        int s = pCol->shift(k);
        if (s == 0) continue;  // nothing to do
        // assign a skill to the nurse for the shift
        bool assigned = false;
        double vday = v;
        for (int sk = 0; sk < pScenario_->nSkills(); ++sk)
          if (skillsAllocation[k][s][sk][pNurse->pPosition_->id_]
              > epsilon()) {
#ifdef DBG
            if (!pNurse->hasSkill(sk)) {
              std::cout << currentSolToString() << std::endl;
              std::cout << allocationToString() << std::endl;
              std::cout << coverageToString() << std::endl;
              Tools::throwError(
                  "Trying to assign skill %s that Nurse %d (%s) hasn't.",
                  pScenario_->skillName(sk).c_str(),
                  pNurse->num_, pNurse->name_.c_str());
            }
#endif
            assigned = true;
            pNurse->roster_.assignTask(k, pScenario_->pShift(s), sk);
            if (skillsAllocation[k][s][sk][pNurse->pPosition_->id_] >
                vday - epsilon()) {
              skillsAllocation[k][s][sk][pNurse->pPosition_->id_] -= vday;
              break;
            }
            if (!fractional)
              printf("Nurse %d (%s) must be allocated to different skills as "
                     "the coverage is fractional (%.9f) on day %d "
                     "on shift %d (%s) for the skill %s.\n",
                     pNurse->num_, pNurse->name_.c_str(),
                     skillsAllocation[k][s][sk][pNurse->pPosition_->id_],
                     k, s, pScenario_->shiftName(s).c_str(),
                     pScenario_->skillName(sk).c_str());
            vday -= skillsAllocation[k][s][sk][pNurse->pPosition_->id_];
            skillsAllocation[k][s][sk][pNurse->pPosition_->id_] = 0;
          }
        if (!assigned) {
          std::cout << currentSolToString() << std::endl;
          std::cout << allocationToString() << std::endl;
          std::cout << coverageToString() << std::endl;
          std::cerr << "No skill found for Nurse " << pNurse->num_ << " ("
                    << pNurse->name_ << ") on day " << k << " on shift " << s
                    << " (" << pScenario_->shiftName(s) << ")" << std::endl;
          return false;
        }
      }

      // record the column (and its value) as one of those active for the
      // LiveNurse
      computeColumnCost(pCol.get());
      pNurse->columns_.push_back(*pCol);
      pNurse->colVals_.push_back(v);
    }
  }

#ifdef DBG
  //
  for (int k = 0; k < pDemand_->nDays_; ++k)
    for (int s = 0; s < pScenario_->nShifts(); ++s)
      for (int sk = 0; sk < pScenario_->nSkills(); ++sk)
        for (double v : skillsAllocation[k][s][sk])
          if (v > epsilon()) {
            std::cout << currentSolToString() << std::endl;
            std::cout << allocationToString() << std::endl;
            std::cout << coverageToString() << std::endl;
            Tools::throwError(
                "Allocation on day %d on shift %d (%s) for skill %s is not "
                "covered (missing %.9f).", k, s,
                pScenario_->shiftName(s).c_str(),
                pScenario_->skillName(sk).c_str(), v);
          }
#endif

  // build the states of each nurse
  for (const auto& pNurse : theLiveNurses_) {
    pNurse->buildStates();
    solution_.push_back(pNurse->roster_);
  }

  // set the corresponding cost
  solutionCost_ = pModel_->getBestUB();

  return true;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void MasterProblem::computeColumnCost(Column *col) const {
  pRCPricer_->computeCost(col);
}


void MasterProblem::saveSolution() {
  if (!storeSolution()) return;
  if (pScenario_->isINRC2_)
    // add the number of nurses only when using the default values
    // for the week indices (i.e. weekIndices_ is empty)
    displaySolutionMultipleWeeks(param_.weekIndices_.empty());
  else
    solutionToXmlINRC();
}

std::string MasterProblem::currentSolToString() const {
//  rep << allocationToString();
//  rep << coverageToString();
  // fetch the active columns and store them by nurses
  vector2D<MyVar *> colsByNurses(pScenario_->nNurses());
  for (MyVar *var : pModel_->getActiveColumns())
    if (pModel_->getVarValue(var) > epsilon())
      colsByNurses[var->getNurseNum()].push_back(var);

  // print the solutions
  std::stringstream rep;
  for (const auto &vec : colsByNurses)
    for (MyVar *var : vec) {
      double v = pModel_->getVarValue(var);
      rep << var->getNurseNum() << ": " << v << std::endl;
      PColumn pCol = getPColumn(var);
      computeColumnCost(pCol.get());
      rep << pCol->toString() << std::endl;
    }
  return rep.str();
}

//------------------------------------------------------------------------------
// Build the, possibly fractional, roster corresponding to the solution
// currently stored in the model
//------------------------------------------------------------------------------
vector3D<double> MasterProblem::fractionalRoster() const {
  vector3D<double> fractionalRoster;
  Tools::initVector3D(&fractionalRoster,
                      nNurses(),
                      nDays(),
                      pDemand_->nShifts_,
                      .0);

  // Retrieve current fractional roster for each nurse
  for (MyVar *var : pModel_->getActiveColumns()) {
    if (var->getCompactColumn().empty()) continue;
    double value = pModel_->getVarValue(var);
    if (value < epsilon()) continue;
    PColumn pat = getPColumn(var);
    vector2D<double> &fractionalRoster2 = fractionalRoster[pat->nurseNum()];
    for (int k = pat->firstDayId(); k <= pat->lastDayId(); ++k)
      fractionalRoster2[k][pat->shift(k)] += value;
  }

  return fractionalRoster;
}

bool MasterProblem::checkIfColumnAlreadyPresent(
    const Column &col, bool printErr) const {
  std::vector<double> column = col.getCompactColumn();
  for (MyVar *var : pModel_->getActiveColumns()) {
    bool equal = true;
    for (int j = 0; j < column.size(); ++j)
      if (std::abs(column[j] - var->getCompactColumn()[j]) > epsilon()) {
        equal = false;
        break;
      }
    if (equal) {
      // check if current pattern is forbidden
      if (printErr) {
        if (var->getUB() < epsilon()) {
          std::cerr << "Column already present as a forbidden column (ub=0): "
                    << var->name_ << std::endl;
          // print branched nodes
          std::cerr << pTree()->getCurrentNode()->writeInheritance()
                    << std::endl;
        } else {
          std::cerr << "Column already present as an active column: "
                    << var->name_ << std::endl;
        }
        std::cerr << col.toString() << std::endl;

        DualCosts duals(this);
        std::cerr << duals.toString(col.nurseNum(), col) << std::endl;
        break;
      }
      return true;
    }
  }
  return false;
}

// add the column to the problem
MyVar * MasterProblem::createColumn(const Column &col, const char *baseName) {
  // add all constraints to the column
  vector<MyCons *> cons;
  vector<double> coeffs;
  addConsToCol(&cons, &coeffs, col);

  // create the column variable
  MyVar *var;
  char name[255];
  snprintf(name, sizeof(name), "%s_N%d", baseName, col.nurseNum());
  pModel_->createIntColumn(&var,
                           name,
                           col.cost(),
                           col.getCompactColumn(),
                           col.reducedCost(),
                           cons,
                           coeffs);
  return var;
}


// add a given constraint to the column
void MasterProblem::addConsToCol(std::vector<MyCons *> *cons,
                                 std::vector<double> *coeffs,
                                 const Column &col) const {
  for (ConstraintMP *pC : columnConstraints_)
    pC->addConsToCol(cons, coeffs, col);
}

/******************************************************
*Get the duals values per day for a nurse
 ******************************************************/
void MasterProblem::update(const PDemand& pDemand) {
  if (pDemand->nDays_ != pDemand_->nDays_)
    Tools::throwError("The new demand must have the same "
                      "size than the old one, so that's ");

  // set the pointer
  pDemand_ = pDemand;

  // update the other constraints that may have changed
  for (ConstraintMP* pC : constraints_)
    pC->update();
}

// return the costs of all active columns associated to the type
std::map<CostType, double> MasterProblem::getColumnsCosts() const {
  return getColumnsCosts(pModel_->getActiveColumns());
}

std::map<CostType, double> MasterProblem::getColumnsCosts(
    const std::vector<MyVar *> &vars) const {
  std::map<CostType, double> costs;
  for (MyVar *var : vars) {
    double value = pModel_->getVarValue(var);
    if (value > epsilon()) {
      PColumn pCol = getPColumn(var);
      computeColumnCost(pCol.get());
      for (const auto &p : pCol->costs())
        costs[p.first] += p.second * value;
    }
  }
  return costs;
}

string MasterProblem::costsConstraintsToString() const {
  std::stringstream rep;
  char buffer[100];

  rep << costsColumnsToString();
  rep << "-----------------------------------------\n";

  for (auto pC : constraints_) {
    if (!pC->printInSolutionCosts()) continue;
    double total = pC->getTotalCost();
    if (abs(total) > epsilon()) {
      std::string n = pC->name + " costs";
      snprintf(buffer, sizeof(buffer), "%-40s %10.0f \n", n.c_str(), total);
      rep << buffer;
    }
  }
  rep << "-----------------------------------------\n";
  rep << "\n";

  return rep.str();
}

string MasterProblem::costsColumnsToString() const {
  std::map<CostType, double> costs = getColumnsCosts();

  std::stringstream rep;
  char buffer[100];

  double total = 0;
  for (const auto &p : costs) total += p.second;
  snprintf(buffer, sizeof(buffer), "%-40s %10.0f \n", "Columns costs", total);
  rep << buffer;
  rep << "-----------------------------------------\n";

  for (const auto &p : costs) {
    if (abs(p.second) < epsilon()) continue;
    snprintf(buffer, sizeof(buffer),
             "%5s%-35s %10.0f \n", "",
             prettyNamesByCostType.at(p.first).c_str(), p.second);
    rep << buffer;
  }

  return rep.str();
}

string MasterProblem::allocationToString() const {
  std::stringstream rep;

  int nbNurses = pScenario_->nNurses();
  int nbShifts = pScenario_->nShifts();
  int firstDay = pDemand_->firstDayId_, nbDays = pDemand_->nDays_;

  rep << std::endl;
  rep << "Allocations of the (potentially fractional) current solution:"
      << std::endl;
  char buff[100];
  snprintf(buff, sizeof(buff), "%20s", "");
  rep << buff;
  for (int day = firstDay; day < firstDay + nbDays; day++) {
    if (day != firstDay && Day::isFirstDayOfWeek(day)) rep << "| ";
    rep << "|" << Day::toDayOfWeekShortName(day) << " ";
  }
  rep << "|" << std::endl;
  rep << "-------------------------------------" << std::endl;

  vector3D<double> fractionalRoster = this->fractionalRoster();
  for (int n = 0; n < nbNurses; n++) {
    PLiveNurse pNurse = theLiveNurses_[n];
    const vector2D<double> &fnurseFractionalRoster = fractionalRoster[n];
    snprintf(buff, sizeof(buff), "%-12s", pNurse->name_.c_str());
    rep << buff;
    for (int s = 1; s < nbShifts; ++s) {
      if (s > 1) {
        snprintf(buff, sizeof(buff), "%12s", "");
        rep << buff;
      }
      snprintf(buff, sizeof(buff), "%-8s", pScenario_->shiftName(s).c_str());
      rep << buff;
      for (int day = firstDay; day < firstDay + nbDays; day++) {
        if (day != firstDay && Day::isFirstDayOfWeek(day))
          rep << "| ";
        double shiftValue = fnurseFractionalRoster[day][s];
        if (shiftValue > 1 - epsilon()) {
          rep << "| 1  ";
        } else if (shiftValue > epsilon()) {
          char buffer[100];
          snprintf(buffer, sizeof(buffer), "|%4.2f", shiftValue);
          rep << buffer;
        } else {
          rep << "| -  ";
        }
      }
      rep << "|" << std::endl;
    }
    rep << std::endl;
  }
  rep << std::endl;

  return rep.str();
}

string MasterProblem::coverageToString() const {
  std::stringstream rep;

  int nbShifts = pScenario_->nShifts();
  int nbSkills = pScenario_->nSkills();
  int firstDay = pDemand_->firstDayId_, nbDays = pDemand_->nDays_;

  rep << std::endl;
  rep << "Coverage of the (potentially fractional) current solution:"
      << std::endl;
  char buff[100];
  snprintf(buff, sizeof(buff), "%21s", "");
  rep << buff;
  for (int day = firstDay; day < firstDay + nbDays; day++) {
    if (day != firstDay && Day::isFirstDayOfWeek(day)) rep << "| ";
    rep << "|" << Day::toDayOfWeekShortName(day) << " ";
  }
  rep << "|" << std::endl;
  rep << "-------------------------------------" << std::endl;

  for (int s = 1; s < nbShifts; ++s) {
    snprintf(buff, sizeof(buff), "%-8s", pScenario_->shiftName(s).c_str());
    rep << buff;
    for (int sk = 0; sk < nbSkills; sk++) {
      if (sk > 0) {
        snprintf(buff, sizeof(buff), "%8s", "");
        rep << buff;
      }
      snprintf(
          buff, sizeof(buff), "%-12s", pScenario_->skillName(sk).c_str());
      rep << buff;

      for (int day = firstDay; day < firstDay + nbDays; day++) {
        if (day != firstDay && Day::isFirstDayOfWeek(day))
          rep << "| ";
        double shiftValue =
            pModel_->getVarValue(getSkillsAllocVars()[day][s][sk]);
        char buffer[100];
        if (std::fabs(shiftValue - round(shiftValue)) < epsilon())
          snprintf(buffer, sizeof(buffer), "|%4.0f", round(shiftValue));
        else
          snprintf(buffer, sizeof(buffer), "|%4.2f", shiftValue);
        rep << buffer;
      }

      rep << "|" << std::endl;
    }
  }
  rep << std::endl;

  return rep.str();
}

// Compute the lagrangian bound
double MasterProblem::computeLagrangianBound(double objVal) const {
  Tools::throwError("Lagrangian bound not implemented "
                    "for this master problem.");
  return -LARGE_SCORE;
// return objVal+sumRedCost-getStabCost();
}

// Compute an approximation of the dual UB based on the lagrangian bound
// It could be useful to measure the quality of a dual solution (used when
// stabilizing).
double MasterProblem::computeApproximateDualUB(double objVal) const {
  double sumRedCost = 0, minRedCost = pPricer()->getLastMinReducedCost();
  for (double v : pPricer()->getLastMinReducedCosts()) {
    if (v < minRedCost) v = minRedCost;
    sumRedCost += v;
  }
  return objVal + sumRedCost;
}
