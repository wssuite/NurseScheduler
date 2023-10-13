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

#include "DeterministicSolver.h"

#include <algorithm>
#include <list>
#include <map>
#include <memory>

#include "InitializeInstance.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "ReadWrite.h"


using std::string;
using std::vector;
using std::map;
using std::pair;


//-----------------------------------------------------------------------------
//
//  C l a s s   D e t e r m i n i s t i c S o l v e r
//
//  Solves the problem with deterministic demand
//  From a given problem (number of weeks, nurses, etc.), can compute a
//  solution.
//
//-----------------------------------------------------------------------------
DeterministicSolver::DeterministicSolver(
    const PScenario& pScenario, const DeterministicSolver &solver) :
    Solver(pScenario),
    options_(solver.options_),
    rollingParameters_(solver.getRollingParameters()),
    lnsParameters_(solver.getLnsParameters()),
    completeParameters_(solver.getCompleteParameters()),
    pCompleteSolver_(nullptr),
    pRollingSolver_(nullptr),
    pLNSSolver_(nullptr) {
  setGeneratePResourcesFunction(solver.generatePResourcesFunc());
  // The nurses must be preprocessed to use their positions
  if (!isPreprocessedNurses_) this->preprocessTheNurses();
}

DeterministicSolver::DeterministicSolver(const PScenario& pScenario,
                                         const InputPaths &inputPaths) :
    Solver(pScenario),
    pCompleteSolver_(nullptr), pRollingSolver_(nullptr), pLNSSolver_(nullptr) {
#ifdef NS_DEBUG
  std::cout << "# New deterministic solver created!" << std::endl;
#endif

  // The nurses must be preprocessed to use their positions
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  // set the options of the deterministic solver
  // (the corresponding method needs to be changed manually for tests)
  //
#ifdef NS_DEBUG
  std::cout << "# Set the options" << std::endl;
#endif
  this->initializeOptions(inputPaths);
  std::cout << std::endl;
}

DeterministicSolver::~DeterministicSolver() {
  stop();
}

void DeterministicSolver::stop(bool wait) {
  // force to stop if was paused
  if (pCompleteSolver_) {
    pCompleteSolver_->getJob().askStop();
    if (wait) pCompleteSolver_->getJob().wait();
  }
}
//----------------------------------------------------------------------------
//
// PARAMETERS FUNCTIONS
// Set the parameters of every solver
//
//----------------------------------------------------------------------------

#define PARAM(field, value) { completeParameters_.field = value; \
            lnsParameters_.field = value; \
            rollingParameters_.field = value; }

// Initialize deterministic options with default values
//
void DeterministicSolver::initializeOptions(const InputPaths &inputPaths) {
  // check if any path in input
  param_.outdir_ = inputPaths.solutionPath();
  if (param_.outdir_.empty()) param_.outdir_ = inputPaths.logPath();
  if (param_.outdir_.empty()) param_.outdir_ = "outfiles/";
  if (inputPaths.logPath().empty())
    param_.logfile_ = param_.outdir_+"ns_log.txt";
  Tools::mkdirs(param_.outdir_);  // create directories

  // initialize the log files when no path is specified
  string logDeterministic =
      inputPaths.logPath().empty() ? "" : inputPaths.logPath()
          + "LogDeterministic.txt";
  string logSolver = inputPaths.logPath().empty() ? "" : inputPaths.logPath()
      + "LogSolver.txt";

  // global options are set to default values in the .h

  // random seed is initialized
  options_.randomSeed_ = inputPaths.randSeed();

  // default parameters for the param
  completeParameters_ =
      SolverParam(options_.verbose_, TWO_DIVES, logDeterministic, logSolver);
  completeParameters_.outdir_ = param_.outdir_;
  lnsParameters_ = completeParameters_;
  rollingParameters_ = completeParameters_;

  if (pScenario_->isINRC2_) {
    PARAM(optimalAbsoluteGap_, 5)
    PARAM(absoluteGap_, 5)
  }

  completeParameters_.optimalityLevel(OPTIMALITY);

  // read any options defined in the param file
  if (!inputPaths.paramFile().empty())
    readOptionsFromFile(inputPaths);

  // override default values with argument values
  if (inputPaths.verbose() >= 0) {
    options_.verbose_ = inputPaths.verbose();
    completeParameters_.verbose(options_.verbose_);
    rollingParameters_.verbose(options_.verbose_);
    lnsParameters_.verbose(options_.verbose_);
  }

  if (inputPaths.timeOut() >= 0)
    options_.totalTimeLimitSeconds_ = static_cast<int>(inputPaths.timeOut());
  PARAM(maxSolvingTimeSeconds_, options_.totalTimeLimitSeconds_)

  if (!inputPaths.SPType().empty()) {
    std::string type = Tools::toUpperCase(inputPaths.SPType());
    SPType t = spTypesByName.at(type);
    PARAM(spType_, t)
  }

  if (inputPaths.SPStrategy() >= 0) {
    PARAM(spParam_.strategyLevel_, inputPaths.SPStrategy())
  }

  if (!inputPaths.RCSPPType().empty()) {
    std::string type = Tools::toUpperCase(inputPaths.RCSPPType());
    RCSPPType t = rcsppTypesByName.at(type);
    PARAM(rcsppType_, t)
  }

  if (inputPaths.nCandidates() >= 1) {
    PARAM(nCandidates_, inputPaths.nCandidates())
  }

  // set the number of threads
  if (inputPaths.nThreads() >= 0)
    options_.nThreads_ = inputPaths.nThreads();
  Tools::ThreadsPool::setMaxGlobalThreads(options_.nThreads_);

  // override default values if not solving relaxation to optimality
  if (!completeParameters_.solveRelaxationToOptimality_) {
    if (completeParameters_.spParam_.spMaxSolvingTimeSeconds_
        > LARGE_TIME - 1) {
      double availTime = completeParameters_.maxSolvingTimeSeconds_
          * completeParameters_.maxSolvingTimeRatioForSubproblems_
          / nNurses() * options_.nThreads_;
      PARAM(spParam_.spMaxSolvingTimeSeconds_, availTime)
    }
  }

  if (completeParameters_.spType_ != ROSTER &&
      completeParameters_.spParam_.rcsppBidirectional_) {
    std::cout << "Disable bidirectional search when solving "
                 "a rotation sub problem." << std::endl;
    PARAM(spParam_.rcsppBidirectional_, false)
  }

  if (completeParameters_.spType_ == ROSTER &&
      completeParameters_.spParam_.rcsppRandomStartDay_) {
    std::cout << "Disable RandomStartDay when solving "
                 "a roster sub problem." << std::endl;
    PARAM(spParam_.rcsppRandomStartDay_, false)
  }
}

// Read deterministic options from a file
void DeterministicSolver::readOptionsFromFile(const InputPaths &inputPaths) {
  // load all the options
  string content = Tools::loadOptions(
      inputPaths.paramFile(),
      [this](const std::string &field, std::fstream *file) {
        // Here and below, read the parameters that refer to branch and price
        // solution.
        // They are not options of the deterministic solver, but they need to be
        // set at this stage
        if (Tools::strEndsWith(field, "verbose")) {
          *file >> options_.verbose_;
          completeParameters_.verbose(options_.verbose_);
          rollingParameters_.verbose(options_.verbose_);
          lnsParameters_.verbose(options_.verbose_);
        } else if (options_.setParameter(field, file)) {
          // nothing
        } else if (completeParameters_.setParameter(
            field, file, {&rollingParameters_, &lnsParameters_})) {
          // nothing
        } else if (completeParameters_.spParam_.setParameter(
            field, file, {&rollingParameters_.spParam_,
                          &lnsParameters_.spParam_})) {
          // nothing
        } else if (Tools::strEndsWith(field, "rollingOptimalityLevel")) {
          std::string strOpt;
          *file >> strOpt;
          rollingParameters_.optimalityLevel(strToOptimalityLevel(strOpt));
        } else if (Tools::strEndsWith(field, "lnsOptimalityLevel")) {
          std::string strOpt;
          *file >> strOpt;
          lnsParameters_.optimalityLevel(strToOptimalityLevel(strOpt));
        } else if (Tools::strEndsWith(field, "completeOptimalityLevel")) {
          std::string strOpt;
          *file >> strOpt;
          completeParameters_.optimalityLevel(strToOptimalityLevel(
              strOpt));
        } else {
          return false;
        }
        return true;
      });

  Tools::LogOutput log(options_.logfile_, true, true);
  log << "===================================================" << std::endl;
  log.addCurrentTime() << std::endl << "Deterministic options : "
                       << inputPaths.paramFile() << std::endl;
  log << content << std::endl;
  log << "===================================================" << std::endl;
}

//----------------------------------------------------------------------------
//
// STATISTICS OF THE OVERALL SOLUTION PROCESS
//
//----------------------------------------------------------------------------

void DeterministicSolver::updateInitialStats(Solver *pSolver) {
  auto *pMaster = dynamic_cast<MasterProblem *>(pSolver);
  if (!pMaster) return;
  auto *pModel = dynamic_cast<BcpModeler *>(pMaster->pModel());

  stats_.bestUBInitial_ = pMaster->computeSolutionCost();
  stats_.bestUB_ = stats_.bestUBInitial_;
  stats_.rootLB_ = pModel->getRootLB();
  double f = param_.absoluteGap_;
  stats_.bestLB_ = std::min(
      stats_.bestUB_, f * ceil((pModel->getBestLB() - pModel->epsilon()) / f));
  stats_.timeInitialSol_ = pMaster->timerTotal().dSinceInit();
  if (pMaster->getJob().finish())
    pMaster->resetTimerTotal();

  // Details on Branch and price related runtimes
  //
  stats_.timeGenColRoot_ = pModel->getTimeFirstRoot();
  stats_.timeGenColMaster_ = pModel->getTimeStats().time_lp_solving;
  stats_.timeGenSubProblems_ = pModel->getTimeStats().time_var_generation;

  // details on Branch and price iterations
  //
  stats_.itGenColInitial_ = pModel->getNbLpIterations();
  stats_.nodesBBInitial_ = pModel->getNbNodes();
  if (!stats_.nodesBBInitial_) stats_.timeGenColRoot_ = stats_.timeInitialSol_;
}

void DeterministicSolver::updateImproveStats(Solver *pSolver) {
  auto *pMaster = dynamic_cast<MasterProblem *>(pSolver);
  if (!pMaster) return;
  auto *pModel = dynamic_cast<BcpModeler *>(pMaster->pModel());

  stats_.bestUB_ = pMaster->computeSolutionCost();
  stats_.timeImproveSol_ = pMaster->timerTotal().dSinceInit();

  // Details on Branch and price related runtimes
  stats_.timeGenColMaster_ += pModel->getTimeStats().time_lp_solving;
  stats_.timeGenSubProblems_ += pModel->getTimeStats().time_var_generation;

  // details on Branch and price iterations
  stats_.itGenColImprove_ = pModel->getNbLpIterations();
  stats_.nodesBBImprove_ = pModel->getNbNodes();

  // number of problems solved in each phase
  //
  if (options_.withLNS_) {
    stats_.itImproveSol_ = 0;
  }
}


//----------------------------------------------------------------------------
//
// SOLVE FUNCTIONS
// The one is general for the whole process
//
//----------------------------------------------------------------------------

// Main function
double DeterministicSolver::solve(const std::vector<Roster> &solution) {
  objValue_ = XLARGE_SCORE;

  // set the random seed to the same fixed value for every solution of a new
  // scenario
  // this is absolutely necessary for reproducibility of the results
  Tools::initializeRandomGenerator(this->options_.randomSeed_);

  // First divide the scenario into connected positions if requested
  // and if instance is big enough
  if (options_.divideIntoConnectedPositions_ && !useCompleteSolverAnyway()) {
    options_.divideIntoConnectedPositions_ = false;
    objValue_ = this->solveByConnectedPositions();
  } else {
    // If the scenario is divided into connected positions, the solution
    // of the subproblems goes in the "else" below.
    // Find a good feasible solution using a rolling horizon planning or
    // solving directly the complete horizon.
    // When providing an initial solution with the LNS, skip this step
    solveOneComponent(solution);
  }
  return objValue_;
}

// After dividing into connected components, solve one component
double DeterministicSolver::solveOneComponent(
    const std::vector<Roster> &solution) {
  objValue_ = XLARGE_SCORE;

  // Always solve small problems to optimality; this can actually save time
  if (useCompleteSolverAnyway()) {
    // if not solving the linear relaxation
    if (completeParameters_.stopAfterXSolution_ > 0)
      completeParameters_.optimalityLevel(OPTIMALITY);
    objValue_ = solveCompleteHorizon(solution);
    return objValue_;
  }
  // Find a good feasible solution using a rolling horizon planning or
  // solving directly the complete horizon.
  if (options_.withRollingHorizon_)
    objValue_ = solveWithRollingHorizon(solution);
  else
    objValue_ = this->solveCompleteHorizon(solution);
  // do not bother improving the solution if it is already optimal
  if (status_ == OPTIMAL)
    return objValue_;

  // Improve the solution with an LNS
  if (options_.withLNS_)
    objValue_ = this->solveWithLNS(solution);

  return objValue_;
}

//------------------------------------------------------------------------
//
// Solve the complete horizon using the input algorithm
//
//------------------------------------------------------------------------

double DeterministicSolver::solveCompleteHorizon(
    const std::vector<Roster> &solution) {
  // Initialize solver and solve
  bool firstSolve = !pCompleteSolver_;
  if (firstSolve || pCompleteSolver_->getJob().finish()) {
    pCompleteSolver_ = std::unique_ptr<Solver>(newSolverWithInputAlgorithm(
            options_.solutionAlgorithm_, completeParameters_));
    Tools::Job job([this, solution](Tools::Job job) {
      pCompleteSolver_->solve(completeParameters_, solution);
    });
    pCompleteSolver_->attachJob(job);
  } else {
    firstSolve = true;  // if resuming, it's still the first solve
    pCompleteSolver_->setParameters(completeParameters_);
  }
  // add LNS heuristic if scenario is big enough
  // WARNING: doesn't work as 2 BCP cannot run in parallel
//  if (!useCompleteSolverAnyway() && completeParameters_.performLNSHeuristic_){
//    auto pLNS = new DeterministicSolver(pScenario_, *this);
//    pLNS->options_.withLNS_ = true;  // activate LNS
//    pLNS->options_.nThreads_ = 1;
//    // deactivate heuristic
//    pLNS->lnsParameters_.performHeuristicAfterXNode_ = -1;
//    pLNS->options_.divideIntoConnectedPositions_ = false;
//    pCompleteSolver_->attachHeuristic(pLNS);
//  }
  // start or resume the resolution
  Tools::Job job = pCompleteSolver_->getJob();
  if (options_.pauseSolveOnTimeout_ || job.active()) {
    pThreadsPool_->run(job, true);
    pThreadsPool_->wait();
  } else {
    job.run();
  }

  if (options_.verbose_ > 1 && pCompleteSolver_->isSolutionInteger()) {
    std::cout << pCompleteSolver_->currentSolToString() << std::endl;
    std::cout << pCompleteSolver_->solutionToString() << std::endl;
  }

  // store stat
  if (firstSolve)
    updateInitialStats(pCompleteSolver_.get());
  else
    updateImproveStats(pCompleteSolver_.get());

  return treatResults(pCompleteSolver_.get());
}


//----------------------------------------------------------------------------
// After the end of a solution process: retrieve status, solution, etc.
//----------------------------------------------------------------------------

double DeterministicSolver::treatResults(Solver *pSolver) {
  // Change status back to feasible a feasible solution was found
  status_ = pSolver->status();
  double timeSinceStart = timerTotal_.dSinceStart();
  std::cout << "Time spent until then: " << timeSinceStart << " s ";
  std::cout << "(time limit is " << options_.totalTimeLimitSeconds_ << " s)"
            << std::endl;

  // Print the solution if required
  if (completeParameters_.printIntermediarySol_) {
    std::cout << pSolver->currentSolToString() << std::endl;
  }

  // Update nurses' states when a feasible solution has been found
  if (pSolver->isSolutionInteger())
    copySolution(pSolver);

  if (pSolver->parameters().stopAfterXSolution_ > 0)
    objValue_ = this->computeSolutionCost();
  else
    objValue_ = pSolver->LB();

  costsConstraints_ = pSolver->costsConstraintsByName();

  return objValue_;
}


//------------------------------------------------------------------------
//
// Solve the problem using a decomposition of the set nurses by connected
// components of the rcspp of positions
//
//------------------------------------------------------------------------

double DeterministicSolver::solveByConnectedPositions() {
  // DIVIDE THE SCENARIO INTO CONNECTED COMPONENTS AND PRINT THE RESULT
  vector<PScenario> scenariosPerComponent =
      divideScenarioIntoConnectedPositions(pScenario_);

  // SOLVE THE PROBLEM COMPONENT-WISE
  int nbNursesToSolve = pScenario_->nNurses();
  std::list<PScenario> scenariosToSolve(scenariosPerComponent.begin(),
                                        scenariosPerComponent.end());
  std::map<PScenario, std::unique_ptr<DeterministicSolver>> solverForScenario;
  while (!scenariosToSolve.empty()) {
    // check if need to recompute the number of nurses to solve
    // for another phase.
    // The idea is to wait that all scenarios has been solved before restarting
    // a new phase and divide the time left proportionally to the number of
    // nurses in each component not solved at optimality.
    // change also the settings to try to consume all the solving time provided
    if (nbNursesToSolve == 0) {
      for (const auto& pSc : scenariosToSolve)
        nbNursesToSolve += pSc->nNurses();
    }

    // retrieve and remove first scenario
    PScenario pScenario = scenariosToSolve.front();
    scenariosToSolve.pop_front();

    // SET THE SOLVER AND SOLVE THE SUBPROBLEM
    // set allowed time proportionally to the number of nurses in the
    // component compared to the total number of nurses left to solve
    double timeLeft =
        options_.totalTimeLimitSeconds_ - timerTotal_.dSinceStart();
    // if no more time left, stop
    if (timeLeft <= 10)
      break;
    double allowedTime = pScenario->nNurses() * timeLeft / nbNursesToSolve;

    if (options_.verbose_ > 0) {
      std::cout << "COMPONENT-WISE SCENARIO" << std::endl;
      std::cout << pScenario->toStringINRC2() << std::endl;
    }

    // fetch the pSolver if already existing or create anew one
    DeterministicSolver *pSolver;
    auto it = solverForScenario.find(pScenario);
    if (it == solverForScenario.end()) {
      pSolver = new DeterministicSolver(pScenario, *this);
      solverForScenario[pScenario] =
              std::unique_ptr<DeterministicSolver>(pSolver);
      // do not pause last component on timeout. There is no need as either:
      // it's a timeout and thus the last solve or the solver finishes properly
      if (solverForScenario.size() != scenariosPerComponent.size())
        pSolver->options_.pauseSolveOnTimeout_ = true;
    } else {
      pSolver = it->second.get();
    }

    // solve the component
    pSolver->setTotalTimeLimit(allowedTime);
    pSolver->solveOneComponent(pSolver->solution());

    // Print the stat of the pSolver
    std::cout << pSolver->getGlobalStat().toString() << std::endl;

    // break if job finished and if the status is not that of
    // a normally finished solution process
    Status newStatus = pSolver->status();
    if (newStatus == UNSOLVED || newStatus == INFEASIBLE) {
      status_ = newStatus;
      std::cout << "Solution process by connected component did not "
                   "terminate normally" << std::endl;
      return XLARGE_SCORE;
    }

    // STORE THE SOLUTION
    // Be particularly cautious that the nurse indices are not the same in the
    // initial scenario and in the solvers per component
    for (const PLiveNurse &pN : pSolver->theLiveNurses_) {
      PLiveNurse pNurse = theLiveNurses_[pN->nurseNum_];
      pNurse->roster(pN->roster_);
    }

    // update the number of nurses to solve
    nbNursesToSolve -= pScenario->nNurses();

    // If not optimal or a timeout has been specified or not trying to
    // solve at optimality,
    // push back the scenario in the list to be solved again in case
    // there is some time left at the end
    if (newStatus == TIME_LIMIT)
      scenariosToSolve.push_back(pScenario);
  }

  // 1- update status only if it's the first component (UNSOLVED) or
  // if all previous ones were optimal
  // 2- update the global stat
  // 3- delete the solver
  for (auto &p : solverForScenario) {
    if (status_ == UNSOLVED || status_ == OPTIMAL)
      status_ = p.second->status(true);
    stats_.add(p.second->getGlobalStat());
    updateCostsConstraints(p.second.get());
  }

  // update nurses' states
  solution_.clear();
  for (int n = 0; n < pScenario_->nNurses(); ++n) {
    solution_.push_back(theLiveNurses_[n]->roster_);
    theLiveNurses_[n]->buildStates();
  }

  return computeSolutionCost();
}

//------------------------------------------------------------------------
// Solve the problem with a rolling horizon algorithm
// The sample period is the number of days for which column variables
// must be integer
//------------------------------------------------------------------------

double DeterministicSolver::solveWithRollingHorizon(
    const std::vector<Roster> &solution) {
  std::cout << "SOLVE WITH ROLLING HORIZON" << std::endl;

  int samplePeriod = options_.rollingSamplePeriod_;
  int controlPeriod = options_.rollingControlHorizon_;

  // Initialize the solver that will handle the iterative solution of the
  // rolling horizon
  //
  bool firstSolve = !pRollingSolver_;
  pRollingSolver_ = std::unique_ptr<Solver>(newSolverWithInputAlgorithm(
          options_.solutionAlgorithm_, rollingParameters_));
  pRollingSolver_->solution(solution);

  // Solve the instance iteratively with a rolling horizon
  //
  int firstDay = 0;  // first day of the current horizon
  double LB = 0;  // LB obtained on the first solve (this is only real LB)
  while (firstDay < pDemand_->nDays_) {
    std::cout << "FIRST DAY = " << firstDay << std::endl << std::endl;

    // last days of the sample horizon and of the control horizon
    //
    int lastDaySample =
        std::min(firstDay + samplePeriod - 1, pDemand_->nDays_ - 1);
    int lastDayControl =
        std::min(firstDay + controlPeriod - 1, pDemand_->nDays_ - 1);

    // Relax the integrality constraints on the variables outside the horizon
    // control (and inside the prediction horizon, but at this stage the
    // prediction horizon includes the whole horizon)
    //
    vector<bool> isRelaxDay(pDemand_->nDays_, false);
    for (int day = lastDayControl + 1; day < pDemand_->nDays_; day++)
      isRelaxDay[day] = true;
    pRollingSolver_->relaxDays(isRelaxDay);

    // Solve the problem with a method that allows for a warm start
    pRollingSolver_->rollingSolve(
        rollingParameters_, firstDay, pRollingSolver_->solution());
    if (firstDay == 0)
      LB = pRollingSolver_->LB();

    // check if solution is optimal and integer
    if (pRollingSolver_->status() == OPTIMAL &&
        pRollingSolver_->isSolutionInteger())
      break;
    if (pRollingSolver_->status(true) == INFEASIBLE)
      break;

    if (rollingParameters_.printIntermediarySol_)
      std::cout << pRollingSolver_->currentSolToString() << std::endl;

    // fix the days of the sample horizon to the values of the solution
    pRollingSolver_->fixAvailabilityBasedOnSolution(lastDaySample);

    // update the first and last day of the sample period
    firstDay = firstDay + samplePeriod;
    lastDayControl = std::min(firstDay + controlPeriod - 1,
                              pDemand_->nDays_ - 1);

    // Set back the integrality constraints for next sampling period
    vector<bool> isUnrelaxDay(pDemand_->nDays_, false);
    for (int day = 0; day <= lastDayControl; day++) isUnrelaxDay[day] = true;
    pRollingSolver_->unrelaxDays(isUnrelaxDay);

    // stop rolling horizon if runtime is exceeded
    //
    double timeSinceStart = timerTotal_.dSinceStart();
    std::cout << "Time spent until then: " << timeSinceStart << " s"
              << std::endl;
    if (timeSinceStart > options_.totalTimeLimitSeconds_) {
      std::cout << "Stop the rolling horizon: time limit is reached!"
                << std::endl;
      break;
    }
  }
  // put the status to optimal if really the case:
  // 1- solution must be integer;
  // 2- solution must  be closed enough of LB (i.e. the LB obtained
  //    on the first solve).
  if (pRollingSolver_->status() == OPTIMAL &&
      pRollingSolver_->isSolutionInteger()) {
    double UB = computeSolutionCost();
    // If UB is not optimal over the whole horizon
    if (UB + pRollingSolver_->epsilon() > LB + rollingParameters_.absoluteGap_)
      pRollingSolver_->status(FEASIBLE);
  }

  // clean the solver and store the solution
  pRollingSolver_->resetNursesAvailabilities();
  if (pRollingSolver_->isSolutionInteger())
    pRollingSolver_->storeSolution();
  if (pRollingSolver_->status() != INFEASIBLE)
    pRollingSolver_->costsConstraintsToString();

  std::cout << "END OF ROLLING HORIZON" << std::endl << std::endl;

  // store stat
  if (firstSolve)
    updateInitialStats(pRollingSolver_.get());
  else
    updateImproveStats(pRollingSolver_.get());

  return treatResults(pRollingSolver_.get());
}

// Set the optimality level of the rolling horizon solver
// This function needs to be called before each new solution, and the behavior
// depends on the first day of the horizon
//
// void DeterministicSolver::rollingSetOptimalityLevel(int firstDay) {
//  rollingParameters_.setOptimalityLevel(TWO_DIVES);
//
//  if (pDemand_->nbDays_ - firstDay <= 21) {
//    rollingParameters_.setOptimalityLevel(options_.);
//  }
// }


//------------------------------------------------------------------------
//
// Solve the problem with a large neighborhood search
// Iteratively fix the schedule of complete weeks or of a set of nurses and
// solve the rest to improve a initial feasible solution
//
//------------------------------------------------------------------------


// Perform the LNS
//
double DeterministicSolver::solveWithLNS(const std::vector<Roster> &solution) {
  std::cout << "SOLVE WITH LNS" << std::endl << std::endl;

  // Initialize data structures for LNS
  //
  this->initializeLNS();
  std::vector<double>
      nursesSelectionWeights(nursesSelectionOperators_.size(), 3.0);
  std::vector<double> daysSelectionWeights(daysSelectionOperators_.size(), 3.0);
  std::vector<double> repairWeights(repairOperators_.size(), 0.1);

  // Initialize the solver that will handle the repair problems
  if (options_.withRollingHorizon_) {
    pLNSSolver_ = pRollingSolver_.get();
  } else {
    pLNSSolver_ = pCompleteSolver_.get();
  }

  // load solution
  if (!solution.empty()) pLNSSolver_->solution(solution);

  // Perform destroy/repair iterations until a given number of iterations
  // without improvement is reached
  //
  int nbItWithoutImprovement = 0;
  double bestObjVal = this->computeSolutionCost();
  std::vector<Roster> bestSolution = solution;
  while (!getJob().shouldStop()) {
    // nbItWithoutImprovement < options_.lnsMaxItWithoutImprovement_) {
    // draw the next destroy operators randomly according to the weights
    int nurseIndex = Tools::drawRandomWithWeights(nursesSelectionWeights);
    int dayIndex = Tools::drawRandomWithWeights(daysSelectionWeights);
    int repairIndex = Tools::drawRandomWithWeights(repairWeights);
    NursesSelectionOperator
        nurseOperator = nursesSelectionOperators_[nurseIndex];
    DaysSelectionOperator dayOperator = daysSelectionOperators_[dayIndex];

    // apply the destroy operator
    adaptiveDestroy(nurseOperator, dayOperator);

    // set the computational time left for the lns after the initialization
    double timeSinceStart = timerTotal_.dSinceStart();
    lnsParameters_.maxSolvingTimeSeconds_ =
        options_.totalTimeLimitSeconds_ - static_cast<int>(timeSinceStart);

    // run the repair operator
    double currentObjVal = pLNSSolver_->LNSSolve(lnsParameters_);

    // stop lns if runtime is exceeded
    timeSinceStart = timerTotal_.dSinceStart();
    std::cout << "Time spent until then: " << timeSinceStart << " s "
              << "(time limit is " << options_.totalTimeLimitSeconds_ << " s)"
              << std::endl;
    if (pLNSSolver_->status() == TIME_LIMIT
        && timeSinceStart <= options_.totalTimeLimitSeconds_ - 5.0) {
      Tools::throwError("Error with the timers in LNS!");
    }

    // store the solution
    if (pLNSSolver_->isSolutionInteger()) {
      status_ = pLNSSolver_->status();
      if (lnsParameters_.printIntermediarySol_) {
        std::cout << pLNSSolver_->currentSolToString() << std::endl;
      }

      // update the weight of adaptive lns
      //
      if (currentObjVal < bestObjVal - lnsParameters_.epsilon_) {
        // update stats
        stats_.lnsImprovementValueTotal_ += bestObjVal - currentObjVal;
        stats_.lnsNbIterationsWithImprovement_++;

        // update best solution
        bestObjVal = currentObjVal;
        bestSolution = pLNSSolver_->solution();

        // update weights of the destroy operators
        nbItWithoutImprovement = 0;
        lnsParameters_.optimalityLevel(TWO_DIVES);
        nursesSelectionWeights[nurseIndex] =
            nursesSelectionWeights[nurseIndex] + 1.0;
        daysSelectionWeights[nurseIndex] =
            daysSelectionWeights[nurseIndex] + 1.0;

        // update weights of the repair operators
        double timeIteration = timerTotal_.dSinceInit() - timeSinceStart;
        repairWeights[repairIndex] += 10.0 / timeIteration;


        // update the counters of iterations with improvement
        stats_.nbImprovementsWithNursesSelection_[nurseIndex]++;
        stats_.nbImprovementsWithDaysSelection_[dayIndex]++;
        stats_.nbImprovementsWithRepair_[repairIndex]++;
      } else {
        nbItWithoutImprovement++;

        if (nbItWithoutImprovement
            > std::min(30, options_.lnsMaxItWithoutImprovement_ / 2)) {
          lnsParameters_.optimalityLevel(OPTIMALITY);
        } else if (nbItWithoutImprovement
            > std::min(10, options_.lnsMaxItWithoutImprovement_ / 4)) {
          lnsParameters_.optimalityLevel(REPEATED_DIVES);
        }
      }
    }

    std::cout << "**********************************************" << std::endl
              << "LNS iteration: " << stats_.lnsNbIterations_
              << "\t" << "Best solution: " << bestObjVal << std::endl
              << "**********************************************" << std::endl;

    if (pLNSSolver_->status() == TIME_LIMIT
        || timeSinceStart > options_.totalTimeLimitSeconds_) {
      std::cout << "Stop the lns: time limit is reached" << std::endl;
      break;
    }

    // unfix every nurse and/or days for next iteration
    std::vector<bool> isUnfixNurse(pScenario_->nNurses(), true);
    pLNSSolver_->unfixNurses(isUnfixNurse);
    pLNSSolver_->resetNursesAvailabilities();

    stats_.lnsNbIterations_++;
  }

  std::cout << "END OF LNS" << std::endl << std::endl;

  // reload best solution
  pLNSSolver_->loadSolution(bestSolution);

  // treat results
  this->updateImproveStats(pLNSSolver_);

  return treatResults(pLNSSolver_);
}

// Prepare data structures for LNS
//
void DeterministicSolver::initializeLNS() {
  // initialize the set of destroy operators of the solver
  if (options_.lnsNursesRandomDestroy_) {
    nursesSelectionOperators_.push_back(NURSES_RANDOM);
  }
  if (options_.lnsNursesPositionDestroy_) {
    nursesSelectionOperators_.push_back(NURSES_POSITION);
    this->organizeTheLiveNursesByPosition();
  }
  if (options_.lnsNursesContractDestroy_) {
    nursesSelectionOperators_.push_back(NURSES_CONTRACT);
    this->organizeTheLiveNursesByContract();
  }
  if (options_.lnsDestroyOverTwoWeeks_) {
    daysSelectionOperators_.push_back(TWO_WEEKS);
  }
  if ((options_.lnsDestroyOverFourWeeks_ && nDays() > 28) ||
      (options_.lnsDestroyOverAllWeeks_ && nDays() <= 28)) {
    daysSelectionOperators_.push_back(FOUR_WEEKS);
  }
  if (options_.lnsDestroyOverAllWeeks_ && nDays() > 28) {
    daysSelectionOperators_.push_back(ALL_WEEKS);
  }

  // initialize the set of repair operators of the lns
  repairOperators_.push_back(REPAIR_TWO_DIVES);
  repairOperators_.push_back(REPAIR_REPEATED_DIVES);
  repairOperators_.push_back(REPAIR_OPTIMALITY);

  // initialize the counters of improvements
  stats_.nbImprovementsWithRepair_.insert(
      stats_.nbImprovementsWithRepair_.begin(),
      repairOperators_.size(),
      0);
  stats_.nbImprovementsWithNursesSelection_.insert(
      stats_.nbImprovementsWithNursesSelection_.begin(),
      nursesSelectionOperators_.size(),
      0);
  stats_.nbImprovementsWithDaysSelection_.insert(
      stats_.nbImprovementsWithDaysSelection_.begin(),
      daysSelectionOperators_.size(),
      0);
}

// Application of the destroy operator
//
void DeterministicSolver::adaptiveDestroy(NursesSelectionOperator nurseOp,
                                          DaysSelectionOperator dayOp) {
  // apply the destroy operator
  std::vector<bool> isFixNurse(pScenario_->nNurses(), true);
  std::vector<bool> isFixDay(pScenario_->nDays(), true);
  std::vector<int> randIndVector;

  // FIRST SET THE NUMBER OF NURSES AND DAYS THAT MUST BE FIXED
  int nbNursesDestroy = 0;
  int nbDaysDestroy = 0;
  switch (dayOp) {
    case TWO_WEEKS: nbNursesDestroy = options_.lnsNbNursesDestroyOverTwoWeeks_;
      nbDaysDestroy = 14;
      break;
    case FOUR_WEEKS:nbNursesDestroy = options_.lnsNbNursesDestroyOverFourWeeks_;
      nbDaysDestroy = 28;
      break;
    case ALL_WEEKS: nbNursesDestroy = options_.lnsNbNursesDestroyOverAllWeeks_;
      nbDaysDestroy = nDays();
      break;
  }

  // SECOND GENERATE THE NURSES WHOSE PLANNING WILL BE FIXED
  //
  switch (nurseOp) {
    case NURSES_RANDOM: {
      randIndVector = Tools::drawRandomIndices(
          nbNursesDestroy, 0, pScenario_->nNurses() - 1);
      for (int ind : randIndVector)
        isFixNurse[ind] = false;
      break;
    }
    case NURSES_POSITION: {
      int randPos = Tools::drawRandomWithWeights(positionWeights_);
      int nbNursesWithPos =
          static_cast<int>(theLiveNursesByPosition_[randPos].size());
      randIndVector =
          Tools::drawRandomIndices(nbNursesDestroy, 0, nbNursesWithPos - 1);
      for (int ind : randIndVector)
        isFixNurse[theLiveNursesByPosition_[randPos][ind]->num_] = false;
      break;
    }
    case NURSES_CONTRACT: {
      int randContract = Tools::drawRandomWithWeights(contractWeights_);
      int nbNursesWithContract =
          static_cast<int>(theLiveNursesByContract_[randContract].size());
      randIndVector = Tools::drawRandomIndices(
          nbNursesDestroy, 0, nbNursesWithContract - 1);
      for (int ind : randIndVector)
        isFixNurse[theLiveNursesByContract_[randContract][ind]->num_] = false;
      break;
    }
  }
  // Fix the nurses that are not destroyed
  pLNSSolver_->fixNurses(isFixNurse);

  // GENERATE THE DAYS THAT WILL BE DESTROYED AND FIX THE OTHERS
  // fix no day if the number of days in the scenario is small

  if (nbDaysDestroy < this->nDays()) {
    // draw the first day of the relaxed interval

    std::vector<double>
        weightDays(pScenario_->nDays() - nbDaysDestroy  - 1, 1.0);
    weightDays[0] = 7;
    for (int i = 1; i < std::max(nDays() - nbDaysDestroy - 1, 6); i++) {
      weightDays[i] = 0.1;
    }
    // Tools::randomInt(0, getNbDays()-nbDaysDestroy-1);
    int firstDay = Tools::drawRandomWithWeights(weightDays);
    for (int day = 0; day < nbDaysDestroy; day++) {
      isFixDay[firstDay + day] = false;
    }
    pLNSSolver_->fixAvailabilityBasedOnSolution(
        isFixDay, pLNSSolver_->solution());
  }

  // DBG
  std::cout << "REPAIR NURSES: ";
  for (unsigned int i = 0; i < isFixNurse.size(); i++) {
    if (!isFixNurse[i]) std::cout << i << "\t";
  }
  std::cout << std::endl;
  std::cout << "REPAIR DAYS: ";
  for (unsigned int i = 0; i < isFixDay.size(); i++) {
    if (!isFixDay[i]) std::cout << i << "\t";
  }
  std::cout << std::endl;
}

// Initialize the organized vectors of live nurses
//
void DeterministicSolver::organizeTheLiveNursesByPosition() {
  nbPositions_ = 0;
  std::vector<PLiveNurse> copyTheLiveNurses = theLiveNurses_;

  // Transfer the nurses with same positions from the copy vector of live nurse
  // to the organized vector of live nurses
  while (!copyTheLiveNurses.empty()) {
    nbPositions_++;
    std::vector<PLiveNurse> theLiveNursesAtPosition;

    // initialize with the first nurse
    theLiveNursesAtPosition.push_back(copyTheLiveNurses[0]);
    const PPosition thisPosition = copyTheLiveNurses[0]->pPosition_;
    copyTheLiveNurses.erase(copyTheLiveNurses.begin());

    // search for the nurses with same position in the remaining nurses
    int nbNursesLeft = static_cast<int>(copyTheLiveNurses.size());
    for (int n = nbNursesLeft - 1; n >= 0; n--) {
      if (copyTheLiveNurses[n]->pPosition_->id_ == thisPosition->id_) {
        theLiveNursesAtPosition.push_back(copyTheLiveNurses[n]);
        copyTheLiveNurses.erase(copyTheLiveNurses.begin() + n);
      }
    }
    theLiveNursesByPosition_.push_back(theLiveNursesAtPosition);
  }

  // Initialize the position weights according to the number of nurses with
  // each position
  //
  for (int p = 0; p < nbPositions_; p++) {
    positionWeights_.push_back(
        static_cast<int>(theLiveNursesByPosition_[p].size()));
  }
}

void DeterministicSolver::organizeTheLiveNursesByContract() {
  nbContracts_ = 0;
  std::vector<PLiveNurse> copyTheLiveNurses = theLiveNurses_;

  // Transfer the nurses with same contract from the copy vector of live nurse
  // to the organized vector of live nurses
  while (!copyTheLiveNurses.empty()) {
    nbContracts_++;
    std::vector<PLiveNurse> theLiveNursesWithContract;

    // initialize with the first nurse
    theLiveNursesWithContract.push_back(copyTheLiveNurses[0]);
    PContract thisContract = copyTheLiveNurses[0]->pContract();
    copyTheLiveNurses.erase(copyTheLiveNurses.begin());

    // search for the nurses with same contract in the remaining nurses
    int nbNursesLeft = copyTheLiveNurses.size();
    for (int n = nbNursesLeft - 1; n >= 0; n--) {
      if (copyTheLiveNurses[n]->pContract()->id_ == thisContract->id_) {
        theLiveNursesWithContract.push_back(copyTheLiveNurses[n]);
        copyTheLiveNurses.erase(copyTheLiveNurses.begin() + n);
      }
    }
    theLiveNursesByContract_.push_back(theLiveNursesWithContract);
  }

  // Initialize the contract weights according to the number of nurses with
  // each contract
  //
  for (int c = 0; c < nbContracts_; c++) {
    contractWeights_.push_back(theLiveNursesByContract_[c].size());
  }
}


//----------------------------------------------------------------------------
//
// Construct solvers that will be used to really make the job
//
//----------------------------------------------------------------------------

// Return a solver with the input algorithm
Solver *DeterministicSolver::newSolverWithInputAlgorithm(
    Algorithm algorithm, const SolverParam &param) {
  Solver *pSolver = newSolver(pScenario_, algorithm,
                              param.spType_,
                              options_.mySolverType_);
  // override the default function if generatePResourcesFunc_
  if (pSolver && generatePResourcesFunc())
    pSolver->setGeneratePResourcesFunction(generatePResourcesFunc());
  // attach job if any
  pSolver->attachJob(getJob());
  return pSolver;
}
