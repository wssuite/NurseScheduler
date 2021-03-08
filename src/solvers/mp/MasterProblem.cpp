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

#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"


#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"
#endif

#ifdef USE_CBC
#include "CbcModeler.h"
#endif

#ifdef USE_SCIP
#include "ScipModeler.h"
#endif

/* namespace usage */
using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;


// print a rc solution
std::string RCSolution::toString(std::vector<int> shiftIDToShiftTypeID) const {
  std::stringstream buff;
  buff << "RC solution of cost " << cost
       << " starting on day " << firstDay << ":" << std::endl;
  for (int k = 0; k < firstDay; k++) buff << "|     ";

  char buffer[500];
  if (shiftIDToShiftTypeID.empty()) {
    for (int s : shifts) {
      snprintf(buffer, sizeof(buffer), "| -:%2d", s);
      buff << buffer;
    }
  } else {
    for (int s : shifts) {
      snprintf(buffer, sizeof(buffer), "|%2d:%2d", shiftIDToShiftTypeID[s], s);
      buff << buffer;
    }
  }
  buff << "|" << std::endl;
  return buff.str();
}

// P a t t e r n   s t a t i c   m e t h o d s

Stretch getStretch(int firstDay,
                   const std::vector<int> &shifts,
                   const PScenario &pScenario) {
  std::vector<PShift> pShifts;
  pShifts.reserve(shifts.size());
  for (int s : shifts) pShifts.emplace_back(pScenario->pShifts_[s]);
  return Stretch(pShifts, firstDay);
}

Stretch getStretch(const std::vector<double> &pattern,
                   const PScenario &pScenario) {
  int firstDay = Pattern::firstDay(pattern),
      length = Pattern::nDays(pattern);
  std::vector<PShift> pShifts(length);
  for (int k = 0; k < length; k++)
    pShifts[k] = pScenario->pShifts_[static_cast<int>(pattern[k + 3])];
  return Stretch(pShifts, firstDay);
}

Stretch getStretch(const std::map<int, int> &shifts,
                   const PScenario &pScenario) {
  if (shifts.empty()) return Stretch(vector<PShift>(), -1);
  // find first day
  int firstDay = shifts.begin()->first;
  for (const auto &p : shifts)
    if (p.first < firstDay) firstDay = p.first;
  // build the vector of PShifts
  std::vector<PShift> pShifts(shifts.size());
  for (const auto &p : shifts)
    pShifts[p.first - firstDay] = pScenario->pShifts_[p.second];
  return Stretch(pShifts, firstDay);
}

Pattern::Pattern(int firstDay,
                 const std::vector<int> &shifts,
                 const PScenario &pScenario,
                 int nurseNum,
                 double cost,
                 double dualCost) :
    nurseNum_(nurseNum), stretch_(getStretch(firstDay, shifts, pScenario)),
    cost_(cost), reducedCost_(dualCost) {}

Pattern::Pattern(const std::map<int, int> &shifts,
                 const PScenario &pScenario,
                 int nurseNum,
                 double cost,
                 double dualCost) :
    nurseNum_(nurseNum), stretch_(getStretch(shifts, pScenario)),
    cost_(cost), reducedCost_(dualCost) {}

Pattern::Pattern(const std::vector<double> &pattern,
                 const PScenario &pScenario) :
    nurseNum_(Pattern::nurseNum(pattern)),
    stretch_(getStretch(pattern, pScenario)),
    cost_(DBL_MAX),
    reducedCost_(DBL_MAX) {}

void Pattern::computeResourcesCosts(const MasterProblem *pMaster,
                                    const State &initialState) {
  // reset costs
  costs_.clear();
  for (int t=CONS_SHIFTS_COST; t <= WEEKEND_COST; t++)
    costs_[(CostType)t] = 0;
  std::map<int, CostType> resourceCostType;
  // create origin, destination and arc
  PRCNode pSource = std::make_shared<RCNode>(
      0, SOURCE_NODE, firstDay(), initialState.pShift_),
      pSink = std::make_shared<RCNode>(
      1, SINK_NODE, nDays(), pShift(lastDay()));
  PRCArc pArc = std::make_shared<RCArc>(
      0, pSource, pSink, stretch_, 0, TO_SINK);
  // create resources
  auto resources = pMaster->createResources(pMaster->liveNurses()[nurseNum_],
                                            &resourceCostType);
  // create expander for stretch
  double c = pArc->cost;
  for (const auto &pR : resources) {
    pR->initialize(*pSource->pAShift, stretch_, pArc);
    addCost(pArc->cost-c, resourceCostType[pR->id()]);
    c = pArc->cost;
  }
  // expand
  PRCLabel pL = std::make_shared<RCLabel>(resources, initialState);
  c = pL->cost();
  for (const auto &pE : pArc->expanders) {
    ResourceValues &v = pL->getResourceValues(pE->resourceId);
    pE->expand(pL, &v);
    addCost(pL->cost()-c, resourceCostType[pE->resourceId]);
    c = pL->cost();
  }
}

std::string Pattern::costsToString() const {
  std::stringstream rep;
  rep << "#       | Consecutive shifts  : "
      << costs_.at(CONS_SHIFTS_COST) << std::endl;
  rep << "#       | Consecutive days off: "
      << costs_.at(CONS_REST_COST) << std::endl;
  rep << "#       | Consecutive days    : "
      << costs_.at(CONS_WORK_COST) << std::endl;
  rep << "#       | Complete weekends   : "
      << costs_.at(COMPLETE_WEEKEND_COST) << std::endl;
  rep << "#       | Preferences         : "
      << costs_.at(PREFERENCE_COST) << std::endl;
  rep << "#       | Worked days         : "
      << costs_.at(DAYS_COST) << std::endl;
  rep << "#       | Worked weekend  : "
      << costs_.at(WEEKEND_COST) << std::endl;
  return rep.str();
}

std::string Pattern::toString(int nbDays) const {
  if (nbDays == -1) nbDays = lastDay()+1;
  std::stringstream rep;
  rep << "#   | PATTERN: N=" << nurseNum_ << "  cost=" << cost_
      << "  dualCost=" << reducedCost_ << "  firstDay=" << firstDay()
      << "  length=" << nDays() << "  duration=" << duration() << std::endl;
  rep << "#   |";
  std::vector<PShift> allTasks(nbDays);
  for (int k = 0; k < nbDays; ++k) {
    if (k < firstDay() || k > lastDay() || pShift(k)->isRest())
      rep << "   |";
    else
      rep << pShift(k)->type << ":" << pShift(k)->id << "|";
  }
  rep << std::endl;
  return rep.str();
}

// Compare rotations on cost
bool Pattern::compareCost(PPattern pat1, PPattern pat2) {
  if (pat1->cost_ == DBL_MAX || pat2->cost_ == DBL_MAX)
    Tools::throwError("Pattern cost not computed.");
  return (pat1->cost_ < pat2->cost_);
}

// Compare rotations on dual cost
bool Pattern::compareDualCost(PPattern pat1, PPattern pat2) {
  if (pat1->reducedCost_ == DBL_MAX || pat2->reducedCost_ == DBL_MAX)
    Tools::throwError("Pattern dual cost not computed.");
  return (pat1->reducedCost_ < pat2->reducedCost_);
}


//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

// Default constructor
MasterProblem::MasterProblem(PScenario pScenario,
                             PDemand pDemand,
                             PPreferences pPreferences,
                             vector<State> *pInitState,
                             SolverType solverType) :

    Solver(pScenario, pDemand, pPreferences, pInitState),
    PrintSolution(),
    pModel_(nullptr),
    positionsPerSkill_(pScenario->nbSkills_),
    skillsPerPosition_(pScenario->nbPositions()),
    pPricer_(nullptr),
    pTree_(nullptr),
    pRule_(nullptr),
    solverType_(solverType),
    optDemandVars_(pDemand_->nDays_),
    numberOfNursesByPositionVars_(pDemand_->nDays_),
    skillsAllocVars_(pDemand_->nDays_),
    minDemandCons_(pDemand_->nDays_),
    optDemandCons_(pDemand_->nDays_),
    numberOfNursesByPositionCons_(pScenario->nbPositions()),
    feasibleSkillsAllocCons_(pDemand_->nDays_) {
  // build the model
  this->initializeSolver(solverType);
}

MasterProblem::~MasterProblem() {
  delete pPricer_;
  delete pRule_;
  delete pTree_;
  delete pModel_;
}

// Initialize the solver at construction
//
void MasterProblem::initializeSolver(SolverType solverType) {
  // This OSI interface is created only to retrieve the proper value of
  // infinity for the solver
  double inf = -1;
  switch (solverType) {
    case S_CLP:pModel_ = new BcpModeler(this, PB_NAME, CLP);
      inf = OsiClpSolverInterface().getInfinity();
      break;
    case S_Gurobi:
#ifdef USE_GUROBI
      pModel_ = new BcpModeler(this, PB_NAME, Gurobi);
      inf = OsiGrbSolverInterface().getInfinity();
#else
      Tools::throwError("BCP has not been built with Gurobi.");
#endif
      break;
    case S_Cplex:
#ifdef USE_CPLEX
      pModel_ = new BcpModeler(this, PB_NAME, Cplex);
      inf = OsiCpxSolverInterface().getInfinity();
#else
      Tools::throwError("BCP has not been built with Cplex.");
#endif
      break;
    case S_CBC:pModel_ = new BcpModeler(this, PB_NAME);
      inf = OsiClpSolverInterface().getInfinity();
      break;
    default:
      Tools::throwError("MasterProblem::initializeSolver: "
                        "the requested solver is not supported presently.");
      break;
  }

  pModel_->setInfinity(inf);

  this->preprocessData();

  /*
  *Build the two vectors linking positions and skills
   */
  for (unsigned int p = 0; p < skillsPerPosition_.size(); p++) {
    vector<int> skills(pScenario_->pPositions()[p]->skills_.size());
    for (unsigned int sk = 0; sk < pScenario_->pPositions()[p]->skills_.size();
         ++sk)
      skills[sk] = pScenario_->pPositions()[p]->skills_[sk];
    skillsPerPosition_[p] = skills;
  }
  for (unsigned int sk = 0; sk < positionsPerSkill_.size(); sk++) {
    vector<int> positions(pScenario_->nbPositions());
    int i(0);
    for (unsigned int p = 0; p < positions.size(); p++)
      if (find(skillsPerPosition_[p].begin(), skillsPerPosition_[p].end(), sk)
          != skillsPerPosition_[p].end()) {
        positions[i] = p;
        ++i;
      }
    positions.resize(i);
    positionsPerSkill_[sk] = positions;
  }
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
  } else {
    pModel_->reset();
  }

  // input an initial solution
  this->initializeSolution(solution);

  // DBG
#ifdef  DBG
  pModel_->writeProblem("master_model");
#endif

  this->solveWithCatch();

  if (pModel_->getParameters().printBranchStats_) {
    pModel_->printStats();
  }

  // if was solving just the relaxation: stopAfterXSolution_ = 0
  if (param_.stopAfterXSolution_ == 0 && getModel()->nbSolutions() == 0) {
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
  // ROLLING/LNS use it
  pPricer_->initNursesAvailabilities();
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
  /* Skills coverage constraints */
  buildSkillsCoverageCons(param);

  /* Initialize the objects used in the branch and price unless the CBC is used
      to solve the problem
   */
  if (solverType_ != S_CBC) {
    /* Rotation pricer */
    pPricer_ = new RCPricer(this, "pricer", param);
    pModel_->addPricer(pPricer_);

    /* Tree */
    RestTree *pTree = new RestTree(pScenario_, pDemand_, param.epsilon_);
    pTree_ = pTree;
    pModel_->addTree(pTree_);

    /* Branching rule */
    pRule_ = new DiveBranchingRule(this, pTree, "branching rule");
    pModel_->addBranchingRule(pRule_);
  }
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
  for (int n = 0; n < nNurses(); ++n)
    for (int k = 0; k <= nDays(); ++k)
      if (fixDays[k]) {
        for (int s = 0; s < nShifts(); ++s)
          availableNursesDaysShifts[n][k][s] = roster[n][k][s] > epsilon();
      }

  nursesAvailabilities(availableNursesDaysShifts);
}

void MasterProblem::nursesAvailabilities(
    const vector3D<bool> &nursesAvailabilities) {
  Solver::nursesAvailabilities(nursesAvailabilities);

  // remove from the pool any column that does not respect the availability
  filterColumnsBasedOnAvailability();
}

// remove any column that does not respect the availabilities
// TODO(legraina): try to keep the column in a pool
//  to perhaps add them back later on
void MasterProblem::filterColumnsBasedOnAvailability() {
  // if all available, all columns can remain in the pool
  if (!isPartialAvailable_) return;

  // iterate through the pool of columns
  std::vector<MyVar *> initialColumns;
  for (MyVar *var : pModel_->getInitialColumns()) {
    // check if column is feasible
    PPattern pat = getPattern(var);
    bool feasible = true;
    for (int k = pat->firstDay(); k <= pat->lastDay(); k++)
      if (!isNurseAvailableOnDayShift(pat->nurseNum_, k, pat->shift(k))) {
        feasible = false;
        break;
      }
    // if feasible, keep it, otherwise delete it
    if (feasible)
      initialColumns.push_back(var);
    else
      delete var;
  }

  // reset the pool
  pModel_->setInitialColumns(initialColumns);
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
    isFixNurse_[n] = isFixNurse_[n] ? true : isFix[n];
    if (isFixNurse_[n]) pPricer_->forbidNurse(n);
  }

  // actually fix the active columns of the input nurses
  for (MyVar *var : pModel_->getInitialColumns())
    if (isFixNurse(var->getNurseNum()))
      var->setLB(1);
}

void MasterProblem::unfixNurses(const std::vector<bool> &isUnfix) {
  // easy to treat the case where no nurse is fixed yet
  if (isFixNurse_.empty()) {
    isPartialFixNurses_ = false;
    pPricer_->clearForbiddenNurses();
  } else {
    // set the list of unfixed nurses
    // + authorize the generation of rotations for the input nurses
    for (PLiveNurse pNurse : theLiveNurses_) {
      int n = pNurse->num_;
      isFixNurse_[n] = isFixNurse_[n] ? !isUnfix[n] : false;
      if (!isFixNurse_[n]) pPricer_->authorizeNurse(n);
    }

    // actually unfix the active columns of the input nurses
    for (MyVar *var : pModel_->getInitialColumns())
      if (!isFixNurse(var->getNurseNum()))
        var->setLB(0);
  }
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
    // reset the model
    pModel_->reset();
    setParameters(param);
    // add the solution in the model
    initializeSolution(solution);
  }

  // solve the problem
  solveWithCatch();
  pModel_->loadBestSol(false);
  // copy active columns (the solution) to initial for next run and
  // clear solution
  pModel_->copyActiveToInitialColumns();
  // reload solution
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
  // reset the model
  pModel_->reset();
  setParameters(param);
  // add the solution
  initializeSolution(solution);

  // solve the problem
  solveWithCatch();
  if (pModel_->loadBestSol())
    storeSolution();
  costsConstrainstsToString();

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
                      pScenario_->nbShifts_ - 1,
                      pScenario_->nbSkills_,
                      0,
                      .0);

  for (int k = 0; k < pDemand_->nDays_; ++k)
    for (int s = 0; s < pScenario_->nbShifts_ - 1; ++s)
      for (int sk = 0; sk < pScenario_->nbSkills_; ++sk)
        skillsAllocation[k][s][sk] =
            pModel_->getVarValues(skillsAllocVars_[k][s][sk]);

  // build the rosters
  for (PLiveNurse pNurse : theLiveNurses_)
    pNurse->roster_.reset();

  for (MyVar *var : pModel_->getActiveColumns()) {
    if (pModel_->getVarValue(var) > epsilon()) {
      PPattern pat = getPattern(var);
      PLiveNurse pNurse = theLiveNurses_[pat->nurseNum_];
      for (int k = pat->firstDay(); k <= pat->lastDay(); ++k) {
        int s = pat->shift(k);
        if (s == 0) continue;  // nothing to do
        // assign a skill to the nurse for the shift
        bool assigned = false;
        for (int sk = 0; sk < pScenario_->nbSkills_; ++sk)
          if (skillsAllocation[k][s - 1][sk][pNurse->pPosition_->id_]
              > epsilon()) {
            pNurse->roster_.assignTask(k, s, sk);
            skillsAllocation[k][s - 1][sk][pNurse->pPosition_->id_]--;
            assigned = true;
            break;
          }
        if (!assigned)
          Tools::throwError("No skill found for Nurse %d on day %d on shift %d",
                            pNurse->num_,
                            k,
                            s);
      }
    }
  }

  for (int k = 0; k < pDemand_->nDays_; ++k)
    for (int s = 0; s < pScenario_->nbShifts_ - 1; ++s)
      for (int sk = 0; sk < pScenario_->nbSkills_; ++sk)
        for (double v : skillsAllocation[k][s][sk])
          if (v > epsilon())
            Tools::throwError("Some allocation are not covered.");

  // build the states of each nurse
  solution_.clear();
  for (PLiveNurse pNurse : theLiveNurses_) {
    pNurse->buildStates();
    solution_.push_back(pNurse->roster_);
  }
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

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
  vector2D<MyVar *> colsByNurses(pScenario_->nbNurses());
  for (MyVar *var : pModel_->getActiveColumns()) {
    double v = pModel_->getVarValue(var);
    if (v < epsilon()) continue;
    colsByNurses[var->getNurseNum()].push_back(var);
    std::stringstream rep;
    rep << var->getNurseNum() << ": " << v << std::endl;
    PPattern pat = getPattern(var);
    rep << pat->toString(nDays())
        << std::endl;
  }

  // print the solutions
  std::stringstream rep;
  for (const auto &vec : colsByNurses)
    for (MyVar *var : vec) {
      double v = pModel_->getVarValue(var);
      rep << var->getNurseNum() << ": " << v << std::endl;
      PPattern pat = getPattern(var);
      rep << pat->toString(nDays())
          << std::endl;
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
    vector2D<double> &fractionalRoster2 = fractionalRoster[pat->nurseNum_];
    for (int k = pat->firstDay(); k <= pat->lastDay(); ++k)
      fractionalRoster2[k][pat->shift(k)] += value;
  }

  return fractionalRoster;
}

void MasterProblem::checkIfPatternAlreadyPresent(
    const std::vector<double> &pattern) const {
  for (MyVar *var : pModel_->getActiveColumns()) {
    bool equal = true;
    for (int j = 0; j < pattern.size(); ++j)
      if (std::fabs(pattern[j] - var->getPattern()[j]) > epsilon()) {
        equal = false;
        break;
      }
    if (equal) {
      string name = var->name_;
      Tools::throwError("Pattern already present as column: " + name);
    }
  }
}

/******************************************************
*Get the duals values per day for a nurse
 ******************************************************/
// build a DualCosts structure
PDualCosts MasterProblem::buildDualCosts(PLiveNurse pNurse) const {
  return std::make_shared<DualCosts>(getShiftsDualValues(pNurse),
                                     getConstantDualvalue(pNurse));
}

// return the dual values associated to the demand
vector2D<double> MasterProblem::getShiftsDualValues(PLiveNurse pNurse) const {
  vector2D<double> dualValues(pDemand_->nDays_);
  for (int k = 0; k < pDemand_->nDays_; ++k)
    dualValues[k] = pModel_->getDuals(
        numberOfNursesByPositionCons_[pNurse->pPosition_->id_][k]);
  return dualValues;
}

double MasterProblem::getConstantDualvalue(PLiveNurse pNurse) const {
  return 0;
}

PDualCosts MasterProblem::buildRandomDualCosts(bool optimalDemandConsidered,
                                               int NDaysShifts) const {
  return std::make_shared<DualCosts>(
      getRandomWorkedDualCosts(optimalDemandConsidered, NDaysShifts),
      Tools::randomDouble(-10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND,
                          10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND));
}

vector2D<double> MasterProblem::getRandomWorkedDualCosts(
    bool optimalDemandConsidered, int NDaysShifts) const {
  // TODO(JO): below, every dual cost is initialized, so I fear that they
  //  will all be used to modify arcs costs, even in the roster-based
  //  decomposition
  vector2D<double> workDualCosts;
  if (optimalDemandConsidered) {
    // This following 2D vector will contain all the dual costs corresponding
    // to each day and to each shift associated with this day.
    Tools::initVector2D(
        &workDualCosts, pDemand_->nDays_, pScenario_->nbShifts_ - 1, 0.0);
    // All the values in the vector are initialized to 0 and only a few of
    // these will be randomly replaced by a value corresponding to the
    // optimal demand weight multiplied by a random coefficient between 1 and 3.
    // The number of days-shifts whose dual cost will be changed is given by
    // the parameter 'NDaysShifts'.
    std::set<int> daysToBeModified;
    while (daysToBeModified.size() < NDaysShifts)
      daysToBeModified.insert(Tools::randomInt(0, pDemand_->nDays_-1));
    for (auto day : daysToBeModified) {
      int idShift = Tools::randomInt(0, pScenario_->nbShifts_-2);
      int coeff = Tools::randomInt(1, 3);
      workDualCosts[day][idShift] =
          coeff*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND;
    }
  } else {
    workDualCosts = Tools::randomDoubleVector2D(
        pDemand_->nDays_, pScenario_->nbShifts_ - 1,
        0, 3*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
  }
  return workDualCosts;
}

/*
*Skills coverage constraints
 */
void MasterProblem::buildSkillsCoverageCons(const SolverParam &param) {
  char name[255];
  // initialize vectors
  Tools::initVector3D<MyVar *>(&optDemandVars_,
                               pDemand_->nDays_,
                               pScenario_->nbShifts_ - 1,
                               pScenario_->nbSkills_,
                               nullptr);
  Tools::initVector3D<MyVar *>(&numberOfNursesByPositionVars_,
                               pDemand_->nDays_,
                               pScenario_->nbShifts_ - 1,
                               pScenario_->nbPositions(),
                               nullptr);
  Tools::initVector4D<MyVar *>(&skillsAllocVars_,
                               pDemand_->nDays_,
                               pScenario_->nbShifts_ - 1,
                               pScenario_->nbSkills_,
                               pScenario_->nbPositions(),
                               nullptr);
  Tools::initVector3D<MyCons *>(&minDemandCons_,
                                pDemand_->nDays_,
                                pScenario_->nbShifts_ - 1,
                                pScenario_->nbSkills_,
                                nullptr);
  Tools::initVector3D<MyCons *>(&optDemandCons_,
                                pDemand_->nDays_,
                                pScenario_->nbShifts_ - 1,
                                pScenario_->nbSkills_,
                                nullptr);
  Tools::initVector3D<MyCons *>(&numberOfNursesByPositionCons_,
                                pScenario_->nbPositions(),
                                pDemand_->nDays_,
                                pScenario_->nbShifts_ - 1,
                                nullptr);
  Tools::initVector3D<MyCons *>(&feasibleSkillsAllocCons_,
                                pDemand_->nDays_,
                                pScenario_->nbShifts_ - 1,
                                pScenario_->nbPositions(),
                                nullptr);

  for (int k = 0; k < pDemand_->nDays_; k++) {
    // forget s=0, it's a resting shift
    for (int s = 1; s < pScenario_->nbShifts_; s++) {
      for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
        // create variables
        snprintf(name, sizeof(name), "optDemandVar_%d_%d_%d", k, s, sk);
        pModel_->createPositiveVar(&optDemandVars_[k][s - 1][sk],
                                   name,
                                   pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
        for (int p = 0; p < pScenario_->nbPositions(); p++) {
          snprintf(name, sizeof(name),
                   "skillsAllocVar_%d_%d_%d_%d", k, s, sk, p);
          pModel_->createPositiveVar(&skillsAllocVars_[k][s - 1][sk][p],
                                     name,
                                     0);
        }
        // adding variables and building minimum demand constraints
        vector<MyVar *> vars1(positionsPerSkill_[sk].size());
        vector<double> coeffs1(positionsPerSkill_[sk].size());
        for (unsigned int p = 0; p < positionsPerSkill_[sk].size(); ++p) {
          vars1[p] = skillsAllocVars_[k][s - 1][sk][positionsPerSkill_[sk][p]];
          coeffs1[p] = 1;
        }

        MyVar *vFeasibility;
        snprintf(name, sizeof(name),
                 "minDemandFeasibilityVar_%d_%d_%d", k, s, sk);
        pModel_->createPositiveFeasibilityVar(&vFeasibility, name);
        vars1.push_back(vFeasibility);
        coeffs1.push_back(1);

        snprintf(name, sizeof(name), "minDemandCons_%d_%d_%d", k, s, sk);
        pModel_->createGEConsLinear(&minDemandCons_[k][s - 1][sk],
                                    name,
                                    pDemand_->minDemand_[k][s][sk],
                                    vars1,
                                    coeffs1);

        // STAB:Add stabilization variable
        if (param.isStabilization_) {
          snprintf(name, sizeof(name), "stabMinDemand_%d_%d_%d", k, s, sk);
          addStabVariables(
              param, name, minDemandCons_[k][s - 1][sk], false, true);
        }

        // adding variables and building optimal demand constraints
        vars1.push_back(optDemandVars_[k][s - 1][sk]);
        coeffs1.push_back(1);
        snprintf(name, sizeof(name), "optDemandCons_%d_%d_%d", k, s, sk);
        pModel_->createGEConsLinear(&optDemandCons_[k][s - 1][sk],
                                    name,
                                    pDemand_->optDemand_[k][s][sk],
                                    vars1,
                                    coeffs1);

        // STAB:Add stabilization variable
        if (param.isStabilization_) {
          snprintf(name, sizeof(name), "stabOptDemand_%d_%d_%d", k, s, sk);
          addStabVariables(
              param, name, optDemandCons_[k][s - 1][sk], false, true);
        }
      }

      for (int p = 0; p < pScenario_->nbPositions(); p++) {
        // creating variables
        snprintf(name, sizeof(name), "nursesNumber_%d_%d_%d", k, s, p);
        // DBG
        // pModel_->createIntVar(&numberOfNursesByPositionVars2[p], name, 0);
        pModel_->createPositiveVar(
            &numberOfNursesByPositionVars_[k][s - 1][p], name, 0);
        // adding variables and building number of nurses constraints
        vector<MyVar *> vars3;
        vector<double> coeff3;
        vars3.push_back(numberOfNursesByPositionVars_[k][s - 1][p]);
        coeff3.push_back(-1);
        snprintf(name, sizeof(name), "nursesNumberCons_%d_%d_%d", k, s, p);
        pModel_->createEQConsLinear(&numberOfNursesByPositionCons_[p][k][s - 1],
                                    name,
                                    0,
                                    vars3,
                                    coeff3);

        // adding variables and building skills allocation constraints
        int const nonZeroVars4(1 + skillsPerPosition_[p].size());
        vector<MyVar *> vars4(nonZeroVars4);
        vector<double> coeff4(nonZeroVars4);
        vars4[0] = numberOfNursesByPositionVars_[k][s - 1][p];
        coeff4[0] = 1;
        for (int sk = 1; sk < nonZeroVars4; ++sk) {
          vars4[sk] =
              skillsAllocVars_[k][s - 1][skillsPerPosition_[p][sk - 1]][p];
          coeff4[sk] = -1;
        }
        snprintf(name, sizeof(name),
                 "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
        pModel_->createEQConsLinear(&feasibleSkillsAllocCons_[k][s - 1][p],
                                    name,
                                    0,
                                    vars4,
                                    coeff4);
      }
    }
  }
}

int MasterProblem::addSkillsCoverageConsToCol(vector<MyCons *> *cons,
                                              vector<double> *coeffs,
                                              const Pattern &pat) const {
  int nbCons(0);

  int p = theLiveNurses_[pat.nurseNum_]->pPosition()->id_;
  for (int k = pat.firstDay(); k <= pat.lastDay(); ++k) {
    int s = pat.shift(k);
    if (pScenario_->isAnyShift(s)) {
      for (int s0 = 1; s0 < pScenario_->nbShifts_; ++s0) {
        ++nbCons;
        cons->push_back(numberOfNursesByPositionCons_[p][k][s0 - 1]);
        coeffs->push_back(1.0);
      }
    } else if (pScenario_->isWorkShift(s)) {  // if work
      ++nbCons;
      cons->push_back(numberOfNursesByPositionCons_[p][k][s - 1]);
      coeffs->push_back(1.0);
    }
  }

  return nbCons;
}

void MasterProblem::updateDemand(PDemand pDemand) {
  if (pDemand->nDays_ != pDemand_->nDays_)
    Tools::throwError("The new demand must have the same "
                      "size than the old one, so that's ");

  // set the pointer
  pDemand_ = pDemand;

  // modify the associated constraints
  for (int k = 0; k < pDemand_->nDays_; k++)
    for (int s = 1; s < pScenario_->nbShifts_; s++)
      for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
        minDemandCons_[k][s - 1][sk]->setLhs(pDemand_->minDemand_[k][s][sk]);
        optDemandCons_[k][s - 1][sk]->setLhs(pDemand_->optDemand_[k][s][sk]);
      }
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
           pModel_->getTotalCost(optDemandVars_));  // , true));
  rep << buffer;
  rep << "-----------------------------------------\n";
  rep << "\n";

  return rep.str();
}

string MasterProblem::allocationToString(bool printInteger) const {
  std::stringstream rep;

  int nbNurses = pScenario_->nbNurses_;
  int nbShifts = pScenario_->nbShifts_;
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
      snprintf(buff, sizeof(buff), "%-8s", pScenario_->intToShift_[s].c_str());
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

  int nbShifts = pScenario_->nbShifts_;
  int nbSkills = pScenario_->nbSkills_;
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
    snprintf(buff, sizeof(buff), "%-8s", pScenario_->intToShift_[s].c_str());
    rep << buff;
    for (int sk = 0; sk < nbSkills; sk++) {
      if (sk > 0) {
        snprintf(buff, sizeof(buff), "%8s", "");
        rep << buff;
      }
      snprintf(
          buff, sizeof(buff), "%-12s", pScenario_->intToSkill_[sk].c_str());
      rep << buff;

      for (int day = firstDay; day < firstDay + nbDays; day++) {
        double shiftValue =
            pModel_->getVarValue(skillsAllocVars_[day][s - 1][sk]);
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
  double sumRedCost = 0, minRedCost = pPricer_->getLastMinReducedCost();
  for (double v : pPricer_->getLastMinReducedCosts()) {
    if (v < minRedCost) v = minRedCost;
    sumRedCost += v;
  }
  return objVal + sumRedCost;
}

//---------------------------------------------------------------------------
//
// STAB: Methods required to implement stabilization in the column generation
//
// Ref: Lübbecke, Marco E., and Jacques Desrosiers.
// "Selected topics in column generation."
// Operations research 53.6 (2005): 1007-1023.
//
//---------------------------------------------------------------------------

// Add stabilization variables z for the box [b_, b+] with the penalties c
// if getting outside of the box:
// dual = obj += - c_+ z_+ - c_- z_-, s.t.: b_- - z_-<= Pi <= b_+ + z_+
// primal = obj += - b_- y_- + b_+ y_+, s.t.: y_- <= c_-, y_+ <= c_+
// if primal constraint is <= -> dual <= 0 -> just need the LB of the box
//                            -> create just minus var
// if primal constraint is >= -> dual >= 0 -> just need the UB of the box
//                            -> create just plus var
// WARNING: they are inactive at the beginning
//
void MasterProblem::addStabVariables(
    const SolverParam &param,
    const char *name,
    MyCons *cons,
    bool LECons,
    bool GECons) {
  MyVar *var;
  char n[255];

  // The lower side of the box
  if (LECons) {
    snprintf(n, sizeof(n), "%s_minus", name);
    pModel_->createPositiveVar(&var, n, LARGE_SCORE, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, -1);
    stabVariablesMinus_.push_back(var);
  } else {
    stabVariablesMinus_.push_back(nullptr);
  }

  // The  upper side of the box
  if (GECons) {
    snprintf(n, sizeof(n), "%s_plus", name);
    pModel_->createPositiveVar(&var, n, LARGE_SCORE, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, 1);
    stabVariablesPlus_.push_back(var);
  } else {
    stabVariablesPlus_.push_back(nullptr);
  }

  stabConstraints_.push_back(cons);
  stabBoxCenters_.push_back(0);
}

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
void MasterProblem::stabUpdate(OsiSolverInterface *solver, bool recenter) {
  // stabilization variables corresponding to the cover constraints
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    double center = stabBoxCenters_[i],
        radius = varPlus ? varPlus->getCost() - center :
                 center + varMinus->getCost();
    double dual = pModel_->getDual(stabConstraints_[i], true);
    double penaltyFactor = 1;

    // if dual within the box, decrease radius and increase cost
    if (dual > center + radius + epsilon() &&
        dual < center - radius - epsilon()) {
      // decrease radius
      if (param_.isStabUpdateBoxRadius_)
        radius /= param_.stabBoxRadiusFactor_;
      // increase penalty
      if (param_.isStabUpdatePenalty_)
        penaltyFactor *= param_.stabBoxRadiusFactor_;
    } else {
      // increase radius
      if (param_.isStabUpdateBoxRadius_)
        radius *= param_.stabPenaltyFactor_;
      // decrease penalty
      if (param_.isStabUpdatePenalty_)
        penaltyFactor /= param_.stabPenaltyFactor_;
    }

    if (radius > param_.stabBoxRadiusMax_)
      radius = param_.stabBoxRadiusMax_;

    // recenter box
    if (recenter) {
      center = dual;
      stabBoxCenters_[i] = dual;
    }

    // update box
    if (varPlus) updateVarCostInSolver(varPlus, solver, center + radius);
    if (varMinus) updateVarCostInSolver(varMinus, solver, -center + radius);
    // update penalty
    if (param_.isStabUpdatePenalty_) {
      if (varPlus) multiplyUbInSolver(varPlus, solver, penaltyFactor);
      if (varMinus) multiplyUbInSolver(varMinus, solver, penaltyFactor);
    }
  }
}

// STAB
// activate the stabilization variables and center them on the current duals
void MasterProblem::stabInitializeBoundAndCost(OsiSolverInterface *solver) {
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    double center = pModel_->getDual(stabConstraints_[i], true);
    stabBoxCenters_[i] = center;
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver,
                            center + param_.stabBoxRadiusIni_);
      updateVarUbInSolver(varPlus, solver, param_.stabPenaltyIni_);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver,
                            -center + param_.stabBoxRadiusIni_);
      updateVarUbInSolver(varMinus, solver, param_.stabPenaltyIni_);
    }
  }
}

// STAB
// deactivate the stabilization variables
void MasterProblem::stabDeactivateBoundAndCost(OsiSolverInterface *solver) {
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    stabBoxCenters_[i] = 0;
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver, LARGE_SCORE);
      updateVarUbInSolver(varPlus, solver, 0);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver, LARGE_SCORE);
      updateVarUbInSolver(varMinus, solver, 0);
    }
  }
}

// STAB
// Stop when the stabilization variables are all null
bool MasterProblem::stabCheckStoppingCriterion() const {
  if (!pModel_->getParameters().isStabilization_)
    return true;

  for (MyVar *var : stabVariablesPlus_)
    if (var && pModel_->getVarValue(var) > epsilon())
      return false;
  for (MyVar *var : stabVariablesMinus_)
    if (var && pModel_->getVarValue(var) > epsilon())
      return false;

  return true;
}

// STAB
// return the current cost of the stabilization variables
double MasterProblem::getStabCost() const {
  if (!pModel_->getParameters().isStabilization_)
    return 0;
  return pModel_->getTotalCost(stabVariablesPlus_) +
      pModel_->getTotalCost(stabVariablesMinus_);
}

// STAB
// Multiply the upper bound of the input variable by the input factor
void MasterProblem::multiplyUbInSolver(MyVar *pVar,
                                       OsiSolverInterface *solver,
                                       double factor) {
  int varind = pVar->getIndex();
  double ub = pVar->getUB();

  if (ub != solver->getColUpper()[varind]) {
    Tools::throwError("multiplyUbInSolver: the upper bound stored in the "
                      "variable is not the same as that in the solver!");
  }

  ub *= factor;
  if (ub > param_.stabBoxBoundMax_) ub = param_.stabBoxBoundMax_;

  solver->setColUpper(varind, ub);
  pVar->setUB(ub);
}

// STAB
// Set the bound of the input variable to the input value
void MasterProblem::updateVarUbInSolver(MyVar *pVar,
                                        OsiSolverInterface *solver,
                                        double value) {
  int varind = pVar->getIndex();
  double ub = pVar->getUB();

  if (ub != solver->getColUpper()[varind]) {
    Tools::throwError("updateVarUbInSolver: the upper bound stored in the "
                      "variable is not the same as that in the solver!");
  }

  if (value > param_.stabBoxBoundMax_) value = param_.stabBoxBoundMax_;

  solver->setColUpper(varind, value);
  pVar->setUB(value);
}

// STAB
// Set the cost of the input variable to the input value
void MasterProblem::updateVarCostInSolver(MyVar *pVar,
                                          OsiSolverInterface *solver,
                                          double value) {
  int varind = pVar->getIndex();
  double cost = pVar->getCost();

  if (cost != solver->getObjCoefficients()[varind]) {
    Tools::throwError("updateVarCostInSolver: the cost stored in the variable "
                      "is not the same as that in the solver!");
  }

  if (value > param_.stabPenaltyMax_) value = param_.stabPenaltyMax_;
  else if (-value > param_.stabPenaltyMax_) value = -param_.stabPenaltyMax_;

  solver->setObjCoeff(varind, value);
  pVar->setCost(value);
}
