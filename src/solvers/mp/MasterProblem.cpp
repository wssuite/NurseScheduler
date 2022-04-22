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
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"


/* namespace usage */
using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;


std::string RCSolution::toString() const {
  std::stringstream rep;
  rep.setf(std::ios_base::fixed, std::ios_base::floatfield);
  rep << "RCSolution: ";
  if (cost_ < DBL_MAX - 1) rep << "cost=" << std::setprecision(0) << cost_;
  if (reducedCost_ < DBL_MAX - 1)
    rep << "  dualCost=" << std::setprecision(2) << reducedCost_;
  rep << "; " << Stretch::toString();
  return rep.str();
}

std::string RCSolution::costsToString() const {
  if (costs_.empty()) {
    std::cout << "WARNING: costs are not computed by cost types." << std::endl;
    return "";
  }
  std::stringstream rep;
  rep.setf(std::ios_base::fixed, std::ios_base::floatfield);
  rep.precision(0);
  rep << "#       | Cons. shifts     : "
      << costs_.at(CONS_SHIFTS_COST) << std::endl;
  rep << "#       | Cons. weekend    : "
      << costs_.at(CONS_WEEKEND_SHIFTS_COST) << std::endl;
  rep << "#       | Cons. days off   : "
      << costs_.at(CONS_REST_COST) << std::endl;
  rep << "#       | Cons. days       : "
      << costs_.at(CONS_WORK_COST) << std::endl;
  rep << "#       | Complete weekends: "
      << costs_.at(COMPLETE_WEEKEND_COST) << std::endl;
  rep << "#       | Preferences      : "
      << costs_.at(PREFERENCE_COST) << std::endl;
  rep << "#       | Worked days      : "
      << costs_.at(TOTAL_WORK_COST) << std::endl;
  rep << "#       | Worked weekend   : "
      << costs_.at(TOTAL_WEEKEND_COST) << std::endl;
  return rep.str();
}

// Compare rotations on cost
bool RCSolution::compareCost(
    const RCSolution &sol1, const RCSolution &sol2) {
  if (sol1.cost_ == DBL_MAX || sol2.cost_ == DBL_MAX)
    Tools::throwError("Pattern cost not computed.");
  return (sol1.cost_ < sol2.cost_);
}

// Compare rotations on dual cost
bool RCSolution::compareReducedCost(
    const RCSolution &sol1, const RCSolution &sol2) {
  if (sol1.reducedCost_ == DBL_MAX || sol2.reducedCost_ == DBL_MAX)
    Tools::throwError("Pattern dual cost not computed.");
  return (sol1.reducedCost_ < sol2.reducedCost_);
}

// P a t t e r n   s t a t i c   m e t h o d s
std::vector<PShift> getPShifts(const std::vector<double> &pattern,
                               const PScenario &pScenario) {
  int length = Pattern::nDays(pattern);
  std::vector<PShift> pShifts(length);
  for (int k = 0; k < length; k++)
    pShifts[k] = pScenario->pShift(static_cast<int>(pattern[k + 3]));
  return pShifts;
}

Pattern::Pattern(MyVar *var, const PScenario &pScenario) :
    RCSolution(Pattern::firstDay(var->getPattern()),
               getPShifts(var->getPattern(), pScenario),
               var->getCost()),
    nurseNum_(Pattern::nurseNum(var->getPattern())) {}

std::string Pattern::toString() const {
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
MasterProblem::MasterProblem(PScenario pScenario, SolverType solverType) :
    Solver(pScenario),
    PrintSolution(),
    pModel_(nullptr),
    pRCPricer_(nullptr),
    solverType_(solverType),
    pResourceCostTypes_(pScenario->nNurses()) {
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
    return pModel_->getRelaxedObjective();
  }

  // print best solution
  pModel_->printBestSol();

  // store solution
  storeSolution();

  // print constraints' costs
  std::cout << costsConstrainstsToString() << std::endl;

  return pModel_->getObjective();
}

void MasterProblem::solveWithCatch() {
  // initialize pricer (useful to authorize/forbid shifts for nurses)
  pPricer()->initNursesAvailabilities();
  // filter out initial columns that does not respect nurses' availabilities
  filterInitialColumnsBasedOnAvailability();
  // then solve
  pModel_->solve();
}

// Resolve the problem with another demand and keep the same preferences
//
double MasterProblem::resolve(PDemand pDemand,
                              const SolverParam &param,
                              const std::vector<Roster> &solution) {
  updateDemand(pDemand);
  setParameters(param);
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
  columnConstraints_.push_back(nursePositionCountConstraint_);
  allocationConstraint_ = new AllocationConstraint(this);
  minDemandConstraint_ = new DemandConstraint(this, true);
  optDemandConstraint_ = new DemandConstraint(
      this, false, true, pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);


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
    pResources_.push_back(generatePResources(pN));
  // split the resources between the master and the subproblem
  // must initialize spResources_
  splitPResources();
}

// Functions to generate the resources for a given nurse
std::map<PResource, CostType>
MasterProblem::defaultGeneratePResources(const PLiveNurse &pN) const {
  const Weights &w = pScenario_->weights();
  PAbstractShift pWork = std::make_shared<AnyWorkShift>();
  std::map<PResource, CostType> mResources = {
      // initialize resource on the total number of working days
      {std::make_shared<SoftTotalShiftDurationResource>(
          minTotalShifts_[pN->num_],
          maxTotalShifts_[pN->num_],
          weightTotalShiftsMin_[pN->num_],
          weightTotalShiftsMax_[pN->num_],
          pWork,
          nDays(),
          pScenario_->maxDuration()), TOTAL_WORK_COST},
      // initialize resource on the total number of week-ends
      {std::make_shared<SoftTotalWeekendsResource>(
          maxTotalWeekends_[pN->num_],
          weightTotalWeekendsMax_[pN->num_],
          pWork,
          nDays()), TOTAL_WEEKEND_COST},
      // initialize resource on the number of consecutive worked days
      {std::make_shared<SoftConsShiftResource>(
          pN->minConsDaysWork(),
          pN->maxConsDaysWork(),
          w.WEIGHT_CONS_DAYS_WORK,
          w.WEIGHT_CONS_DAYS_WORK,
          pWork,
          nDays(),
          pN->pStateIni_->consDaysWorked_),
       CONS_WORK_COST}
  };

  // initialize resources on the number of consecutive shifts of each type
  for (int st = 0; st < pScenario_->nShiftTypes(); st++) {
    shared_ptr<AbstractShift> absShift =
        std::make_shared<AnyOfTypeShift>(st, pScenario_->shiftType(st));
    if (absShift->isWork()) {
      int consShiftsInitial =
          absShift->includes(*pN->pStateIni_->pShift_) ?
          pN->pStateIni_->consShifts_ : 0;
      mResources[std::make_shared<SoftConsShiftResource>(
          pScenario_->minConsShiftsOfType(st),
          pScenario_->maxConsShiftsOfType(st),
          w.WEIGHT_CONS_SHIFTS,
          w.WEIGHT_CONS_SHIFTS,
          absShift,
          nDays(),
          consShiftsInitial)] = CONS_SHIFTS_COST;
    } else if (absShift->isRest()) {
      // needed to evaluate thee historical cost
      mResources[std::make_shared<SoftConsShiftResource>(
          pN->minConsDaysOff(),
          pN->maxConsDaysOff(),
          w.WEIGHT_CONS_DAYS_WORK,
          w.WEIGHT_CONS_DAYS_WORK,
          absShift,
          nDays(),
          pN->pStateIni_->consDaysOff_)] = CONS_REST_COST;
    }
  }

  return mResources;
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
    isRelaxDay_[day] = isRelaxDay_[day] ? true : isRelax[day];
  }
}

void MasterProblem::unrelaxDays(const std::vector<bool> &isUnrelax) {
  if (isRelaxDay_.empty()) {
    isPartialRelaxDays_ = false;
  } else {
    for (int day = 0; day < pDemand_->nDays_; day++)
      isRelaxDay_[day] = isRelaxDay_[day] ? !isUnrelax[day] : false;
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
  for (PLiveNurse pNurse : theLiveNurses_) {
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
    for (PLiveNurse pNurse : theLiveNurses_) {
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
    PPattern pat = getPattern(var);
    for (int k = pat->firstDay(); k <= pat->lastDay(); k++)
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
    std::cout << costsConstrainstsToString() << std::endl;

  // output information and save the solution
  if (pModel_->getParameters().printBranchStats_)
    pModel_->printStats();

  return pModel_->getObjective();
}


//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void MasterProblem::storeSolution() {
  // check if a solution if loaded
  if (pModel_->getActiveColumns().empty()) {
    if (!pModel_->loadBestSol()) {
      std::cerr << "No solution available to store." << std::endl;
      return;
    }
  }

  // retrieve a feasible allocation of skills
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

  // build the rosters
  for (PLiveNurse pNurse : theLiveNurses_)
    pNurse->roster_.reset(pScenario_->pRestShift());

  std::list<MyVar*> activeColumns(
      pModel_->getActiveColumns().begin(), pModel_->getActiveColumns().end());
  int size = activeColumns.size(), n = 0;
  while (!activeColumns.empty()) {
    MyVar *var = activeColumns.front();
    activeColumns.pop_front();
    ++n;
    double v = pModel_->getVarValue(var);
    if (v > epsilon()) {
      PPattern pat = getPattern(var);
      PLiveNurse pNurse = theLiveNurses_[pat->nurseNum()];
      bool fractional = v < 1 - epsilon();
      // if first round, do not process fractional columns
      // store them for next round
      if (fractional && n < size) {
        activeColumns.push_back(var);
        continue;
      }
      for (int k = pat->firstDay(); k <= pat->lastDay(); ++k) {
        if (fractional && !isRelaxDay(k))
          Tools::throwError("Column has a fractional value (%.9f) while "
                            "it should be integer: %s.",
                            v, pat->toString().c_str());
        int s = pat->shift(k);
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
                  pScenario_->skill(sk).c_str(),
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
                     k, s, pScenario_->shift(s).c_str(),
                     pScenario_->skill(sk).c_str());
            vday -= skillsAllocation[k][s][sk][pNurse->pPosition_->id_];
            skillsAllocation[k][s][sk][pNurse->pPosition_->id_] = 0;
          }
        if (!assigned) {
          std::cout << currentSolToString() << std::endl;
          std::cout << allocationToString() << std::endl;
          std::cout << coverageToString() << std::endl;
          Tools::throwError(
              "No skill found for Nurse %d (%s) on day %d on shift %d (%s)",
              pNurse->num_, pNurse->name_.c_str(),
              k, s, pScenario_->shift(s).c_str());
        }
      }
    }
  }

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
                "covered (missing %.9f).", k, s, pScenario_->shift(s).c_str(),
                pScenario_->skill(sk).c_str(), v);
          }

  // build the states of each nurse
  solution_.clear();
  for (PLiveNurse pNurse : theLiveNurses_) {
    pNurse->buildStates();
    solution_.push_back(pNurse->roster_);
  }

  // set the corresponding cost
  solutionCost_ = pModel_->getBestUB();
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void MasterProblem::computePatternCost(Pattern *pat) const {
  pRCPricer_->computeCost(pat);
}


void MasterProblem::save(const vector<int> &weekIndices, string outfile) {
  storeSolution();
  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = weekIndices.size();

  // write separately the solutions of each week in the required output format
  int firstDay = pDemand_->firstDay_;
  for (int w = 0; w < nbWeeks; ++w) {
    // string solutionFile = outdir+std::to_string(weekIndices[w])+".txt";
    Tools::LogOutput solutionStream(outfile);
    solutionStream << solutionToString(firstDay, 7, pScenario_->thisWeek() + w);
    firstDay += 7;
  }
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
      rep << "N" << var->getNurseNum() << ": " << v << std::endl;
      PPattern pat = getPattern(var);
      computePatternCost(pat.get());
      rep << pat->toString();
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
    if (var->getPattern().empty()) continue;
    double value = pModel_->getVarValue(var);
    if (value < epsilon()) continue;
    PPattern pat = getPattern(var);
    vector2D<double> &fractionalRoster2 = fractionalRoster[pat->nurseNum()];
    for (int k = pat->firstDay(); k <= pat->lastDay(); ++k)
      fractionalRoster2[k][pat->shift(k)] += value;
  }

  return fractionalRoster;
}

void MasterProblem::checkIfPatternAlreadyPresent(
    const Pattern &pat) const {
  std::vector<double> pattern = pat.getCompactPattern();
  for (MyVar *var : pModel_->getActiveColumns()) {
    bool equal = true;
    for (int j = 0; j < pattern.size(); ++j)
      if (std::fabs(pattern[j] - var->getPattern()[j]) > epsilon()) {
        equal = false;
        break;
      }
    if (equal) {
      // check if current pattern is forbidden
      if (var->getUB() < epsilon()) {
        std::cerr << "Pattern already present as a forbidden column (ub=0): "
                  << var->name_ << std::endl;
        // print branched nodes
        std::cerr << pTree()->getCurrentNode()->writeInheritance() << std::endl;
      } else {
        std::cerr << "Pattern already present as an active column: "
                  << var->name_ << std::endl;
      }
      std::cerr << pat.toString() << std::endl;
    }
  }
}

// add the column to the problem
MyVar * MasterProblem::createColumn(const Pattern &col, const char *baseName) {
  // add all constraints to the column
  vector<MyCons *> cons;
  vector<double> coeffs;
  addConstoCol(&cons, &coeffs, col);

  // create the column variable
  MyVar *var;
  char name[255];
  snprintf(name, sizeof(name), "%s_N%d", baseName, col.nurseNum());
  pModel_->createIntColumn(&var,
                           name,
                           col.cost(),
                           col.getCompactPattern(),
                           col.reducedCost(),
                           cons,
                           coeffs);
  return var;
}


// add a given constraint to the column
void MasterProblem::addConstoCol(std::vector<MyCons *> *cons,
                  std::vector<double> *coeffs,
                  const Pattern &col) const {
  for (ConstraintMP *pC : columnConstraints_)
    pC->addConsToCol(cons, coeffs, col);
}

/******************************************************
*Get the duals values per day for a nurse
 ******************************************************/
void MasterProblem::updateDemand(PDemand pDemand) {
  if (pDemand->nDays_ != pDemand_->nDays_)
    Tools::throwError("The new demand must have the same "
                      "size than the old one, so that's ");

  // set the pointer
  pDemand_ = pDemand;

  // modify the associated constraints
  minDemandConstraint_->updateDemand();
  optDemandConstraint_->updateDemand();
}

// return the costs of all active columns associated to the type
double MasterProblem::getColumnsCost(CostType costType) const {
  return getColumnsCost(costType, pModel_->getActiveColumns());
}

double MasterProblem::getColumnsCost(CostType costType,
                                const std::vector<MyVar *> &vars) const {
  double cost = 0;
  for (MyVar *var : vars) {
    double value = pModel_->getVarValue(var);
    if (value > epsilon()) {
      PPattern pat = getPattern(var);
      computePatternCost(pat.get());
      cost += pat->costByType(costType) * value;
    }
  }
  return cost;
}

string MasterProblem::costsConstrainstsToString() const {
  std::stringstream rep;

  char buffer[100];
  snprintf(buffer, sizeof(buffer),
           "%-40s %10.0f \n",
           "Rotation costs",
           getColumnsCost(ROTATION_COST));
  rep << buffer;
  rep << "-----------------------------------------\n";
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Cons. shifts costs",
           getColumnsCost(CONS_SHIFTS_COST));
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Cons. weekend shifts costs",
           getColumnsCost(CONS_WEEKEND_SHIFTS_COST));
  rep << buffer;
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Cons. worked days costs",
           getColumnsCost(CONS_WORK_COST));
  rep << buffer;
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Complete weekend costs",
           getColumnsCost(COMPLETE_WEEKEND_COST));
  rep << buffer;
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Preferences costs",
           getColumnsCost(PREFERENCE_COST));
  rep << buffer;
  rep << "-----------------------------------------\n";
  snprintf(buffer, sizeof(buffer),
           "%5s%-35s %10.0f \n",
           "",
           "Cons. rest days costs",
           getColumnsCost(CONS_REST_COST));
  rep << buffer;
  snprintf(buffer,
           sizeof(buffer),
           "%-40s %10.0f \n",
           "Worked days costs",
           getDaysCost());
  rep << buffer;
  snprintf(buffer, sizeof(buffer),
           "%-40s %10.0f \n",
           "Worked weekend costs",
           getWeekendCost());
  rep << buffer;
  snprintf(buffer, sizeof(buffer),
           "%-40s %10.0f \n",
           "Coverage costs",
           pModel_->getTotalCost(optDemandConstraint_->getVariables()));
  rep << buffer;
  double c = pModel_->getTotalCost(allocationConstraint_->getVariables());
  if (c > epsilon()) {
    snprintf(buffer, sizeof(buffer), "%-40s %10.0f \n",
        "Alternative skills costs", c);
    rep << buffer;
  }
  rep << "-----------------------------------------\n";
  rep << "\n";

  return rep.str();
}

string MasterProblem::allocationToString(bool printInteger) const {
  std::stringstream rep;

  int nbNurses = pScenario_->nNurses();
  int nbShifts = pScenario_->nShifts();
  int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nDays_;

  rep << std::endl;
  rep << "Allocations of the (potentially fractional) current solution:"
      << std::endl;
  char buff[100];
  snprintf(buff, sizeof(buff), "%20s", "");
  rep << buff;
  for (int day = firstDay; day < firstDay + nbDays; day++) {
    rep << "|" << Tools::intToDay(day) << " ";
    if (Tools::isSunday(day)) rep << "| ";
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
      snprintf(buff, sizeof(buff), "%-8s", pScenario_->shift(s).c_str());
      rep << buff;
      for (int day = firstDay; day < firstDay + nbDays; day++) {
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
        if (Tools::isSunday(day)) rep << "| ";
      }
      rep << "|" << std::endl;
    }
    rep << std::endl;
  }
  rep << std::endl;

  return rep.str();
}

string MasterProblem::coverageToString(bool printInteger) const {
  std::stringstream rep;

  int nbShifts = pScenario_->nShifts();
  int nbSkills = pScenario_->nSkills();
  int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nDays_;

  rep << std::endl;
  rep << "Coverage of the (potentially fractional) current solution:"
      << std::endl;
  char buff[100];
  snprintf(buff, sizeof(buff), "%21s", "");
  rep << buff;
  for (int day = firstDay; day < firstDay + nbDays; day++) {
    rep << "|" << Tools::intToDay(day) << " ";
    if (Tools::isSunday(day)) rep << "| ";
  }
  rep << "|" << std::endl;
  rep << "-------------------------------------" << std::endl;

  for (int s = 1; s < nbShifts; ++s) {
    snprintf(buff, sizeof(buff), "%-8s", pScenario_->shift(s).c_str());
    rep << buff;
    for (int sk = 0; sk < nbSkills; sk++) {
      if (sk > 0) {
        snprintf(buff, sizeof(buff), "%8s", "");
        rep << buff;
      }
      snprintf(
          buff, sizeof(buff), "%-12s", pScenario_->skill(sk).c_str());
      rep << buff;

      for (int day = firstDay; day < firstDay + nbDays; day++) {
        double shiftValue =
            pModel_->getVarValue(getSkillsAllocVars()[day][s][sk]);
        char buffer[100];
        if (std::fabs(shiftValue - round(shiftValue)) < epsilon())
          snprintf(buffer, sizeof(buffer), "|%4.0f", round(shiftValue));
        else
          snprintf(buffer, sizeof(buffer), "|%4.2f", shiftValue);
        rep << buffer;
        if (Tools::isSunday(day)) rep << "| ";
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
