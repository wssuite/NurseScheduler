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

#include "solvers/mp/RotationMP.h"
#include "solvers/mp/RosterMP.h"
#include "InitializeSolver.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "tools/ReadWrite.h"


// #define COMPARE_EVALUATIONS


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

DeterministicSolver::DeterministicSolver(PScenario pScenario) :
    Solver(pScenario),
    pCompleteSolver_(nullptr),
    pRollingSolver_(nullptr),
    pLNSSolver_(nullptr) {
#ifdef DBG
  std::cout << "# New deterministic solver created!" << std::endl;
#endif

  // The nurses must be preprocessed to use their positions
  if (!isPreprocessedNurses_) this->preprocessTheNurses();
}

DeterministicSolver::DeterministicSolver(PScenario pScenario,
                                         const InputPaths &inputPaths) :
    Solver(pScenario),
    pCompleteSolver_(nullptr), pRollingSolver_(nullptr), pLNSSolver_(nullptr) {
#ifdef DBG
  std::cout << "# New deterministic solver created!" << std::endl;
#endif

  // The nurses must be preprocessed to use their positions
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  // set the options of the deterministic solver
  // (the corresponding method needs to be changed manually for tests)
  //
#ifdef DBG
  std::cout << "# Set the options" << std::endl;
#endif
  this->initializeOptions(inputPaths);
  std::cout << std::endl;

  if (!options_.logfile_.empty()) {
    FILE *pFile;
    pFile = fopen(options_.logfile_.c_str(), "w");
    fprintf(pFile, "HERE COME THE LOGS OF THE MASTER PROBLEM\n\n");
    fclose(pFile);
  }
}

DeterministicSolver::~DeterministicSolver() {
  // delete also the reusable solvers
  delete pCompleteSolver_;
  delete pRollingSolver_;
  // DBG if (pLNSSolver_) delete pLNSSolver_;
}


//----------------------------------------------------------------------------
//
// PARAMETERS FUNCTIONS
// Set the parameters of every solver
//
//----------------------------------------------------------------------------

// Initialize deterministic options with default values
//
void DeterministicSolver::initializeOptions(const InputPaths &inputPaths) {
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
  lnsParameters_ =
      SolverParam(options_.verbose_, TWO_DIVES, logDeterministic, logSolver);
  rollingParameters_ =
      SolverParam(options_.verbose_, TWO_DIVES, logDeterministic, logSolver);

  // read any options defined in the param file
  if (!inputPaths.paramFile().empty())
    readOptionsFromFile(inputPaths);

  // override default values with argument values
  if (inputPaths.verbose() >= 0)
    options_.verbose_ = inputPaths.verbose();
  completeParameters_.verbose(options_.verbose_);
  rollingParameters_.verbose(options_.verbose_);
  lnsParameters_.verbose(options_.verbose_);

  if (inputPaths.timeOut() >= 0)
    options_.totalTimeLimitSeconds_ = inputPaths.timeOut();
  completeParameters_.maxSolvingTimeSeconds_ = options_.totalTimeLimitSeconds_;
  rollingParameters_.maxSolvingTimeSeconds_ = options_.totalTimeLimitSeconds_;
  lnsParameters_.maxSolvingTimeSeconds_ = options_.totalTimeLimitSeconds_;

  if (!inputPaths.SPType().empty()) {
    std::string type = Tools::toUpperCase(inputPaths.SPType());
    SPType t = SPTypesByName.at(type);
    completeParameters_.sp_type_ = t;
    rollingParameters_.sp_type_ = t;
    lnsParameters_.sp_type_ = t;
  }

  if (inputPaths.SPStrategy() >= 0) {
    completeParameters_.spParam_.strategyLevel_ = inputPaths.SPStrategy();
    rollingParameters_.spParam_.strategyLevel_ = inputPaths.SPStrategy();
    lnsParameters_.spParam_.strategyLevel_ = inputPaths.SPStrategy();
  }

  if (!inputPaths.RCSPPType().empty()) {
    std::string type = Tools::toUpperCase(inputPaths.RCSPPType());
    RCSPPType t = RCSPPTypesByName.at(type);
    completeParameters_.rcspp_type_ = t;
    rollingParameters_.rcspp_type_ = t;
    lnsParameters_.rcspp_type_ = t;
  }

  if (inputPaths.nCandidates() >= 1) {
    completeParameters_.nCandidates_ = inputPaths.nCandidates();
    rollingParameters_.nCandidates_ = inputPaths.nCandidates();
    lnsParameters_.nCandidates_ = inputPaths.nCandidates();
  }

  // set the number of threads
  if (inputPaths.nThreads() >= 0)
    options_.nThreads_ = inputPaths.nThreads();
  Tools::ThreadsPool::setMaxGlobalThreads(options_.nThreads_);
}

// Read deterministic options from a file
//
void DeterministicSolver::readOptionsFromFile(const InputPaths &inputPaths) {
  // open the file
  //
  std::fstream file;
  ReadWrite::openFile(inputPaths.paramFile(), &file);


  // go through all the lines of the parameter file
  std::string title;
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '#', '=', &title);
    // ignore line
    if (title.empty() || Tools::strEndsWith(title, "\n")) {
      // read line
      Tools::readUntilChar(&file, '\n', &title);
      continue;
    }


    // Read the name of the scenario
    //
    if (Tools::strEndsWith(title, "divideIntoConnectedPositions")) {
      file >> options_.divideIntoConnectedPositions_;
    } else if (Tools::strEndsWith(title, "withRollingHorizon")) {
      file >> options_.withRollingHorizon_;
    } else if (Tools::strEndsWith(title, "rollingSamplePeriod")) {
      file >> options_.rollingSamplePeriod_;
    } else if (Tools::strEndsWith(title, "rollingControlHorizon")) {
      file >> options_.rollingControlHorizon_;
    } else if (Tools::strEndsWith(title, "rollingPredictionHorizon")) {
      file >> options_.rollingPredictionHorizon_;
    } else if (Tools::strEndsWith(title, "withLNS")) {
      file >> options_.withLNS_;
    } else if (Tools::strEndsWith(title, "lnsMaxItWithoutImprovement")) {
      file >> options_.lnsMaxItWithoutImprovement_;
    } else if (Tools::strEndsWith(title, "lnsNursesRandomDestroy")) {
      file >> options_.lnsNursesRandomDestroy_;
    } else if (Tools::strEndsWith(title, "lnsNursesPositionDestroy")) {
      file >> options_.lnsNursesPositionDestroy_;
    } else if (Tools::strEndsWith(title, "lnsNursesContractDestroy")) {
      file >> options_.lnsNursesContractDestroy_;
    } else if (Tools::strEndsWith(title, "lnsDestroyOverTwoWeeks")) {
      file >> options_.lnsDestroyOverTwoWeeks_;
    } else if (Tools::strEndsWith(title, "lnsDestroyOverFourWeeks")) {
      file >> options_.lnsDestroyOverFourWeeks_;
    } else if (Tools::strEndsWith(title, "lnsDestroyOverAllWeeks")) {
      file >> options_.lnsDestroyOverAllWeeks_;
    } else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverTwoWeeks")) {
      file >> options_.lnsNbNursesDestroyOverTwoWeeks_;
    } else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverFourWeeks")) {
      file >> options_.lnsNbNursesDestroyOverFourWeeks_;
    } else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverAllWeeks")) {
      file >> options_.lnsNbNursesDestroyOverAllWeeks_;
    } else if (Tools::strEndsWith(title, "solutionAlgorithm")) {
      std::string algoName;
      file >> algoName;
      algoName = Tools::toUpperCase(algoName);
      options_.solutionAlgorithm_ = AlgorithmsByName.at(algoName);
    } else if (Tools::strEndsWith(title, "solverType")) {
      std::string solverName;
      file >> solverName;
      solverName = Tools::toUpperCase(solverName);
      options_.MySolverType_ = SolverTypesByName.at(solverName);
#ifdef DBG
      std::cout << "LP solver :" << solverName << std::endl;
#endif
    } else if (Tools::strEndsWith(title, "verbose")) {
      file >> options_.verbose_;
    } else if (Tools::strEndsWith(title, "isStabilization")) {
      // Here and below, read the parameters that refer to branch and price
      // solution.
      // They are not options of the deterministic solver, but they need to be
      // set at this stage
      file >> completeParameters_.isStabilization_;
      lnsParameters_.isStabilization_ = completeParameters_.isStabilization_;
      rollingParameters_.isStabilization_ =
          completeParameters_.isStabilization_;
    } else if (Tools::strEndsWith(title, "isStabUpdateBoxRadius")) {
      file >> completeParameters_.isStabUpdateBoxRadius_;
      lnsParameters_.isStabUpdateBoxRadius_ =
          completeParameters_.isStabUpdateBoxRadius_;
      rollingParameters_.isStabUpdateBoxRadius_ =
          completeParameters_.isStabUpdateBoxRadius_;
    } else if (Tools::strEndsWith(title, "isStabUpdatePenalty")) {
      file >> completeParameters_.isStabUpdatePenalty_;
      lnsParameters_.isStabUpdatePenalty_ =
          completeParameters_.isStabUpdatePenalty_;
      rollingParameters_.isStabUpdatePenalty_ =
          completeParameters_.isStabUpdatePenalty_;
    } else if (Tools::strEndsWith(title, "nbDiveIfBranchOnColumns")) {
      file >> completeParameters_.nbDiveIfBranchOnColumns_;
      lnsParameters_.nbDiveIfBranchOnColumns_ =
          completeParameters_.nbDiveIfBranchOnColumns_;
      rollingParameters_.nbDiveIfBranchOnColumns_ =
          completeParameters_.nbDiveIfBranchOnColumns_;
    } else if (Tools::strEndsWith(title, "branchColumnDisjoint")) {
      file >> completeParameters_.branchColumnDisjoint_;
      lnsParameters_.branchColumnDisjoint_ =
          completeParameters_.branchColumnDisjoint_;
      rollingParameters_.branchColumnDisjoint_ =
          completeParameters_.branchColumnDisjoint_;
    } else if (Tools::strEndsWith(title, "branchColumnUntilValue")) {
      file >> completeParameters_.branchColumnUntilValue_;
      lnsParameters_.branchColumnUntilValue_ =
          completeParameters_.branchColumnUntilValue_;
      rollingParameters_.branchColumnUntilValue_ =
          completeParameters_.branchColumnUntilValue_;
    } else if (Tools::strEndsWith(title, "stopAfterXDegenerateIt")) {
      file >> completeParameters_.stopAfterXDegenerateIt_;
      lnsParameters_.stopAfterXDegenerateIt_ =
          completeParameters_.stopAfterXDegenerateIt_;
      rollingParameters_.stopAfterXDegenerateIt_ =
          completeParameters_.stopAfterXDegenerateIt_;
    } else if (Tools::strEndsWith(title, "heuristicMinIntegerPercent")) {
      file >> completeParameters_.heuristicMinIntegerPercent_;
      lnsParameters_.heuristicMinIntegerPercent_ =
          completeParameters_.heuristicMinIntegerPercent_;
      rollingParameters_.heuristicMinIntegerPercent_ =
          completeParameters_.heuristicMinIntegerPercent_;
    } else if (Tools::strEndsWith(title, "performHeuristicAfterXNode")) {
      file >> completeParameters_.performHeuristicAfterXNode_;
      lnsParameters_.performHeuristicAfterXNode_ =
          completeParameters_.performHeuristicAfterXNode_;
      rollingParameters_.performHeuristicAfterXNode_ =
          completeParameters_.performHeuristicAfterXNode_;
    } else if (Tools::strEndsWith(title, "nCandidates")) {
      file >> completeParameters_.nCandidates_;
      lnsParameters_.nCandidates_ = completeParameters_.nCandidates_;
      rollingParameters_.nCandidates_ = completeParameters_.nCandidates_;
    } else if (Tools::strEndsWith(title, "rollingOptimalityLevel")) {
      std::string strOpt;
      file >> strOpt;
      strOpt = Tools::toUpperCase(strOpt);
      rollingParameters_.optimalityLevel(stringToOptimalityLevel[strOpt]);
    } else if (Tools::strEndsWith(title, "lnsOptimalityLevel")) {
      std::string strOpt;
      file >> strOpt;
      strOpt = Tools::toUpperCase(strOpt);
      lnsParameters_.optimalityLevel(stringToOptimalityLevel[strOpt]);
    } else if (Tools::strEndsWith(title, "completeOptimalityLevel")) {
      std::string strOpt;
      file >> strOpt;
      strOpt = Tools::toUpperCase(strOpt);
      completeParameters_.optimalityLevel(stringToOptimalityLevel[strOpt]);
    } else if (Tools::strEndsWith(title, "spDefaultStrategy")) {
      file >> completeParameters_.spParam_.strategyLevel_;
    } else if (Tools::strEndsWith(title, "spNbRotationsPerNurse")) {
      file >> completeParameters_.spParam_.sp_nbrotationspernurse_;
    } else if (Tools::strEndsWith(title, "spNbNursesToPrice")) {
      file >> completeParameters_.spParam_.sp_nbnursestoprice_;
    } else if (Tools::strEndsWith(title, "spMaxReducedCostBound")) {
      file >> completeParameters_.spParam_.sp_max_reduced_cost_bound_;
    } else if (Tools::strEndsWith(title, "rcsppToOptimality")) {
      file >> completeParameters_.spParam_.rcsppToOptimality_;
    } else if (Tools::strEndsWith(title, "rcsppSortLabels")) {
      file >> completeParameters_.spParam_.rcsppSortLabels_;
    } else if (Tools::strEndsWith(title, "rcsppMinCostToSinks")) {
      file >> completeParameters_.spParam_.rcsppMinCostToSinks_;
    } else if (Tools::strEndsWith(title, "rcsppImprovedDomination")) {
      file >> completeParameters_.spParam_.rcsppImprovedDomination_;
    } else if (Tools::strEndsWith(title, "rcsppEnumerateSubpaths")) {
      file >> completeParameters_.spParam_.rcsppEnumSubpaths_;
    } else if (Tools::strEndsWith(title,
                                  "rcsppEnumerateSubpathsForMinCostToSinks")) {
      file >> completeParameters_.spParam_
           .rcsppEnumSubpathsForMinCostToSinks_;
    } else if (Tools::strEndsWith(title, "rcsppDssr")) {
      file >> completeParameters_.spParam_.rcsppDssr_;
    } else if (Tools::strEndsWith(title, "rcsppNbToExpand")) {
      file >> completeParameters_.spParam_.rcsppNbToExpand_;
    } else if (Tools::strEndsWith(title, "rcsppBidirectional")) {
      file >> completeParameters_.spParam_.rcsppBidirectional_;
    } else if (Tools::strEndsWith(title, "rcsppPulse")) {
      file >> completeParameters_.spParam_.rcsppPulse_;
    } else if (Tools::strEndsWith(title, "rcsppWithRotationGraph")) {
      file >> completeParameters_.spParam_.rcsppWithRotationGraph_;
    }
  }
  lnsParameters_.spParam_ = completeParameters_.spParam_;
  rollingParameters_.spParam_ = completeParameters_.spParam_;
}

//----------------------------------------------------------------------------
//
// STATISTICS OF THE OVERALL SOLUTION PROCESS
//
//----------------------------------------------------------------------------

void DeterministicSolver::updateInitialStats(Solver *pSolver) {
  MasterProblem *pMaster = dynamic_cast<MasterProblem *>(pSolver);
  if (!pMaster) return;
  BcpModeler *pModel = dynamic_cast<BcpModeler *>(pMaster->pModel());

  stats_.bestUBInitial_ = pMaster->computeSolutionCost();
  stats_.bestUB_ = stats_.bestUBInitial_;
  stats_.rootLB_ = pModel->getRootLB();
  stats_.bestLB_ = std::min(stats_.bestUB_,
                            5.0 * ceil((pModel->getBestLB() - pModel->epsilon())
                                           / 5.0));
  stats_.timeInitialSol_ = pMaster->timerTotal()->dSinceStart();
  pMaster->timerTotal()->reset();

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
  MasterProblem *pMaster = dynamic_cast<MasterProblem *>(pSolver);
  if (!pMaster) return;
  BcpModeler *pModel = dynamic_cast<BcpModeler *>(pMaster->pModel());

  stats_.bestUB_ = pMaster->computeSolutionCost();
  stats_.timeImproveSol_ = pMaster->timerTotal()->dSinceStart();

  // Details on Branch and price related runtimes
  //
  stats_.timeGenColMaster_ += pModel->getTimeStats().time_lp_solving;
  stats_.timeGenSubProblems_ += pModel->getTimeStats().time_var_generation;

  // details on Branch and price iterations
  //
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
  objValue_ = 0.0;

  // set the random seed to the same fixed value for every solution of a new
  // scenario
  // this is absolutely necessary for reproducibility of the results
  //
  Tools::initializeRandomGenerator(this->options_.randomSeed_);
  srand(this->options_.randomSeed_);

  // Always solve small problems to optimality
  // This can actually save time
  //
  if ((pScenario_->horizon() <= 28 && pScenario_->nNurses() <= 8)
      || (pScenario_->horizon() <= 56 && pScenario_->nNurses() <= 5)) {
    completeParameters_.optimalityLevel(OPTIMALITY);
    objValue_ = solveCompleteHorizon(solution);
    return objValue_;
  }

  // First divide the scenario into connected positions if requested
  //
  if (options_.divideIntoConnectedPositions_) {
    options_.divideIntoConnectedPositions_ = false;
    objValue_ = this->solveByConnectedPositions();
  } else {
    // If the the scenario is divided into connected positions, the solution
    // of the subproblems goes in the "else" below.
    // Find a good feasible solution using a rolling horizon planning or
    // solving directly the complete horizon.
    //
    if (options_.withRollingHorizon_)
      objValue_ = solveWithRollingHorizon(solution);
    else
      objValue_ = this->solveCompleteHorizon(solution);
    // do not bother improving the solution if it is already optimal
    if (status_ == OPTIMAL)
      return objValue_;

    // Improve the solution with an LNS
    //
    if (options_.withLNS_)
      objValue_ = this->solveWithLNS(solution);
  }
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
  //
  bool firstSolve = !pCompleteSolver_;
  delete pCompleteSolver_;
  pCompleteSolver_ = setSolverWithInputAlgorithm(options_.solutionAlgorithm_,
                                                 completeParameters_);
  pCompleteSolver_->solve(completeParameters_, solution);
  if (options_.verbose_ > 0 && pCompleteSolver_->isSolutionInteger()) {
    std::cout << pCompleteSolver_->currentSolToString() << std::endl;
    std::cout << pCompleteSolver_->solutionToString() << std::endl;
  }

  // store stat
  if (firstSolve)
    updateInitialStats(pCompleteSolver_);
  else
    updateImproveStats(pCompleteSolver_);

  return treatResults(pCompleteSolver_);
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
  if (pSolver->isSolutionInteger()) {
    solution_ = pSolver->solution();
    for (int n = 0; n < pScenario_->nNurses(); ++n) {
      theLiveNurses_[n]->roster_ = solution_[n];
      theLiveNurses_[n]->buildStates();
    }
  }

  objValue_ = this->computeSolutionCost();

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
  std::map<PScenario, DeterministicSolver *> solverForScenario;
  while (!scenariosToSolve.empty()) {
    // check if need to recompute the number of nurses to solve
    // for another phase.
    // The idea is to wait that all scenarios has been solved before restarting
    // a new phase and divide the time left proportionally to the number of
    // nurses in each component not solved at optimality.
    // change also the settings to try to consume all the solving time provided
    if (nbNursesToSolve == 0) {
      for (auto pSc : scenariosToSolve)
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
      std::cout << pScenario->toString() << std::endl;
    }

    // fetch the solver if already existing or create anew one
    DeterministicSolver *solver = nullptr;
    auto it = solverForScenario.find(pScenario);
    if (it == solverForScenario.end()) {
      solver = new DeterministicSolver(pScenario);
      solver->copyParameters(this);
      solverForScenario[pScenario] = solver;
    } else {
      solver = it->second;
    }

    // solve the component
    solver->setTotalTimeLimit(allowedTime);
    solver->solve(solver->solution());

    // Print the stat of the solver
    std::cout << solver->getGlobalStat().toString() << std::endl;

    // break if the status is not that of a normally finished solution process
    Status newStatus = solver->status();
    if (newStatus == UNSOLVED || newStatus == INFEASIBLE) {
      status_ = newStatus;
      std::cout << "Solution process by connected component did not "
                   "terminate normally" << std::endl;
      return LARGE_SCORE;
    }

    // STORE THE SOLUTION
    // Be particularly cautious that the nurse indices are not the same in the
    // initial scenario and in the solvers per component
    for (PLiveNurse pNurse : solver->theLiveNurses_)
      theLiveNurses_[pNurse->nurseNum_]->roster_ =
          solver->solution()[pNurse->num_];

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
    delete p.second;
  }

  // update nurses' states
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
  delete pRollingSolver_;
  pRollingSolver_ = setSolverWithInputAlgorithm(options_.solutionAlgorithm_,
                                                rollingParameters_);
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
    //
//    this->rollingSetOptimalityLevel(firstDay);
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
    //
    pRollingSolver_->fixAvailabilityBasedOnSolution(lastDaySample);

    // update the first and last day of the sample period
    //
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
    if (UB + pRollingSolver_->epsilon()
        > LB + rollingParameters_.optimalAbsoluteGap_)
      pRollingSolver_->status(FEASIBLE);
  }

  // clean the solver and store the solution
  pRollingSolver_->resetNursesAvailabilities();
  if (pRollingSolver_->isSolutionInteger())
    pRollingSolver_->storeSolution();
  if (pRollingSolver_->status() != INFEASIBLE)
    pRollingSolver_->costsConstrainstsToString();

  std::cout << "END OF ROLLING HORIZON" << std::endl << std::endl;

  // store stat
  if (firstSolve)
    updateInitialStats(pRollingSolver_);
  else
    updateImproveStats(pRollingSolver_);

  return treatResults(pRollingSolver_);
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
  //
  if (options_.withRollingHorizon_) {
    pLNSSolver_ = pRollingSolver_;
  } else {
    pLNSSolver_ = pCompleteSolver_;
  }
  if (!solution.empty()) pLNSSolver_->solution(solution);

  // set the computational time left for the lns after the initialization
  double timeSinceStart = timerTotal_.dSinceStart();
  lnsParameters_.maxSolvingTimeSeconds_ =
      options_.totalTimeLimitSeconds_ - timeSinceStart;

  // Perform destroy/repair iterations until a given number of iterations
  // without improvement is reached
  //
  int nbItWithoutImprovement = 0;
  double bestObjVal = this->computeSolutionCost();
  while (true) {
    // nbItWithoutImprovement < options_.lnsMaxItWithoutImprovement_) {
    // draw the next destroy operators randomly according to the weights
    int nurseIndex = Tools::drawRandomWithWeights(nursesSelectionWeights);
    int dayIndex = Tools::drawRandomWithWeights(daysSelectionWeights);
    int repairIndex = Tools::drawRandomWithWeights(repairWeights);
    NursesSelectionOperator
        nurseOperator = nursesSelectionOperators_[nurseIndex];
    DaysSelectionOperator dayOperator = daysSelectionOperators_[dayIndex];

    // apply the destroy operator
    this->adaptiveDestroy(nurseOperator, dayOperator);

    // run the repair operator
    //
    double currentObjVal = pLNSSolver_->LNSSolve(lnsParameters_);

    // stop lns if runtime is exceeded
    //
    timeSinceStart = timerTotal_.dSinceStart();
    std::cout << "Time spent until then: " << timeSinceStart << " s";
    std::cout << "(time limit is " << options_.totalTimeLimitSeconds_ << " s)"
              << std::endl;
    if (pLNSSolver_->status() == TIME_LIMIT
        && timeSinceStart <= options_.totalTimeLimitSeconds_ - 5.0) {
      Tools::throwError("Error with the timers in LNS!");
    }
    if (pLNSSolver_->status() == TIME_LIMIT
        || timeSinceStart > options_.totalTimeLimitSeconds_) {
      std::cout << "Stop the lns: time limit is reached" << std::endl;
      break;
    }

    // store the solution
    //
    solution_ = pLNSSolver_->solution();
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

      // update weights of the destroy operators
      bestObjVal = currentObjVal;
      nbItWithoutImprovement = 0;
      lnsParameters_.optimalityLevel(TWO_DIVES);
      nursesSelectionWeights[nurseIndex] =
          nursesSelectionWeights[nurseIndex] + 1.0;
      daysSelectionWeights[nurseIndex] = daysSelectionWeights[nurseIndex] + 1.0;

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

    std::cout << "**********************************************" << std::endl
              << "LNS iteration: " << stats_.lnsNbIterations_
              << "\t" << "Best solution: " << bestObjVal << std::endl
              << "**********************************************" << std::endl;

    // unfix every nurse and/or days for next iteration
    //
    std::vector<bool> isUnfixNurse(pScenario_->nNurses(), true);
    pLNSSolver_->unfixNurses(isUnfixNurse);
    pLNSSolver_->resetNursesAvailabilities();

    stats_.lnsNbIterations_++;
  }

  std::cout << "END OF LNS" << std::endl << std::endl;

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
  std::vector<bool> isFixDay(pScenario_->horizon(), true);
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
      int nbNursesWithPos = theLiveNursesByPosition_[randPos].size();
      randIndVector =
          Tools::drawRandomIndices(nbNursesDestroy, 0, nbNursesWithPos - 1);
      for (int ind : randIndVector)
        isFixNurse[theLiveNursesByPosition_[randPos][ind]->num_] = false;
      break;
    }
    case NURSES_CONTRACT: {
      int randContract = Tools::drawRandomWithWeights(contractWeights_);
      int nbNursesWithContract = theLiveNursesByContract_[randContract].size();
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
        weightDays(pScenario_->horizon() - nbDaysDestroy - 1, 1.0);
    weightDays[0] = 7;
    for (int i = 1; i < std::max(nDays() - nbDaysDestroy - 1, 6); i++) {
      weightDays[i] = 0.1;
    }
    // Tools::randomInt(0, getNbDays()-nbDaysDestroy-1);
    int firstDay = Tools::drawRandomWithWeights(weightDays);
    for (int day = 0; day < nbDaysDestroy; day++) {
      isFixDay[firstDay + day] = false;
    }
    pLNSSolver_->fixAvailabilityBasedOnSolution(isFixDay);
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
    int nbNursesLeft = copyTheLiveNurses.size();
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
    positionWeights_.push_back(theLiveNursesByPosition_[p].size());
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
    PConstContract thisContract = copyTheLiveNurses[0]->pContract_;
    copyTheLiveNurses.erase(copyTheLiveNurses.begin());

    // search for the nurses with same contract in the remaining nurses
    int nbNursesLeft = copyTheLiveNurses.size();
    for (int n = nbNursesLeft - 1; n >= 0; n--) {
      if (copyTheLiveNurses[n]->pContract_->id_ == thisContract->id_) {
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
Solver *DeterministicSolver::setSolverWithInputAlgorithm(
    Algorithm algorithm, const SolverParam &param) {
  Solver *pSolver = nullptr;
  switch (algorithm) {
    case GENCOL:
      switch (param.sp_type_) {
        case LONG_ROTATION:
        case ALL_ROTATION:
          pSolver = new RotationMP(pScenario_, options_.MySolverType_);
          break;
        case ROSTER:
          pSolver = new RosterMP(pScenario_, options_.MySolverType_);
          break;
        default:Tools::throwError("The subproblem type is not handled yet");
      }
      break;
    default:Tools::throwError("The algorithm is not handled yet");
      break;
  }
  // override the default function if generatePResourcesFunc_
  if (generatePResourcesFunc_)
    pSolver->setGeneratePResourcesFunction(generatePResourcesFunc_);
  return pSolver;
}
