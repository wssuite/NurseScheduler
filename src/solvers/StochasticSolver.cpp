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

#include "StochasticSolver.h"

#include <algorithm>
#include <memory>
#include <utility>

#include "tools/DemandGenerator.h"
#include "solvers/mp/RotationMP.h"
#include "tools/ReadWrite.h"

// #define COMPARE_EVALUATIONS

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;


/******************************************************************************
 * Set the options of the stochastic solver
 * This is not automated, so the options need to be changed inside the code
 * during the tests
 * The solution time depends on the number of nurses
 ******************************************************************************/

void setStochasticSolverOptions(StochasticSolverOptions *options,
                                PScenario pScenario,
                                string solPath,
                                string logPathIni,
                                double timeout) {
  string logStochastic =
      logPathIni.empty() ? "" : logPathIni + "LogStochastic.txt";
  string logSolver = logPathIni.empty() ? "" : logPathIni + "LogSolver.txt";

  options->withIterativeDemandIncrease_ = false;
  options->withEvaluation_ = true;
  options->generationCostPerturbation_ = true;
  options->evaluationCostPerturbation_ = true;
  options->withResolveForGeneration_ = false;
  options->generationAlgorithm_ = GENCOL;
  options->withResolveForEvaluation_ = true;
  options->evaluationAlgorithm_ = GENCOL;
  options->rankingStrategy_ = RK_SCORE;
  options->totalTimeLimitSeconds_ = timeout;
  options->nExtraDaysGenerationDemands_ = 7;
  options->nEvaluationDemands_ = 2;
  options->nDaysEvaluation_ = 14;
  options->nGenerationDemandsMax_ = 100;
  options->logfile_ = logStochastic;
  options->rankingStrategy_ = RK_MEAN;
  options->demandingEvaluation_ = true;
  options->verbose_ = 0;
#ifdef DBG
  options->verbose_ = 1;
#endif

  SolverParam generationParameters;
  generationParameters.maxSolvingTimeSeconds_ =
      options->totalTimeLimitSeconds_ - 1.0;
  generationParameters.printEverySolution_ = false;
  generationParameters.outfile_ = solPath;
  generationParameters.logfile_ = logSolver;
  generationParameters.absoluteGap_ = 5;
  generationParameters.minRelativeGap_ = .05;
  generationParameters.relativeGap_ = .1;
  generationParameters.nbDiveIfMinGap_ = 1;
  generationParameters.nbDiveIfRelGap_ = 2;
  generationParameters.solveToOptimality_ = false;
  generationParameters.weightStrategy_ = RANDOMMEANMAX;

  options->generationParameters_ = generationParameters;

  SolverParam evaluationParameters;
  evaluationParameters.maxSolvingTimeSeconds_ =
      options->totalTimeLimitSeconds_ - 1.0;
  evaluationParameters.printEverySolution_ = false;
  // evaluationParameters.outfile_ = "";
  evaluationParameters.logfile_ = logSolver;
  evaluationParameters.absoluteGap_ = 5;
  evaluationParameters.minRelativeGap_ = .05;
  evaluationParameters.relativeGap_ = .1;
  evaluationParameters.nbDiveIfMinGap_ = 1;
  evaluationParameters.nbDiveIfRelGap_ = 2;
  evaluationParameters.solveToOptimality_ = false;
  evaluationParameters.weightStrategy_ = BOUNDRATIO;
  evaluationParameters.stopAfterXSolution_ = 0;

  options->evaluationParameters_ = evaluationParameters;
}

void setStochasticSolverOptions(
    StochasticSolverOptions *stochasticSolverOptions,
    string instanceName,
    string solPath,
    string logPathIni,
    string stochasticOptionsFile,
    string generationOptionsFile,
    string evaluationOptionsFile) {
  string logStochastic =
      logPathIni.empty() ? "" : logPathIni + "LogStochastic.txt";
  string logSolver = logPathIni.empty() ? "" : logPathIni + "LogSolver.txt";

  ReadWrite::readStochasticSolverOptions(stochasticOptionsFile,
                                         stochasticSolverOptions);
  stochasticSolverOptions->logfile_ = logStochastic;

  SolverParam generationParameters;
  ReadWrite::readSolverOptions(generationOptionsFile, &generationParameters);
  generationParameters.outfile_ = solPath;
  stochasticSolverOptions->generationParameters_ = generationParameters;
  stochasticSolverOptions->generationParameters_.verbose_ =
      stochasticSolverOptions->verbose_;

  SolverParam evaluationParameters;
  ReadWrite::readSolverOptions(evaluationOptionsFile, &evaluationParameters);
  evaluationParameters.logfile_ = logSolver;
  stochasticSolverOptions->evaluationParameters_ = evaluationParameters;
  stochasticSolverOptions->evaluationParameters_.verbose_ =
      stochasticSolverOptions->verbose_;
}

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.),
//  can compute a solution.
//
//-----------------------------------------------------------------------------

StochasticSolver::StochasticSolver(PScenario pScenario,
                                   StochasticSolverOptions options,
                                   vector<PDemand> demandHistory,
                                   double costPreviousWeeks) :
    Solver(pScenario,
           pScenario->pWeekDemand(),
           pScenario->pWeekPreferences(),
           pScenario->pInitialState()),
    options_(options),
    demandHistory_(demandHistory),
    pReusableGenerationSolver_(nullptr),
    costPreviousWeeks_(costPreviousWeeks) {
  std::cout << "# New stochastic solver created!" << std::endl;

  int remainingDays = (pScenario_->nbWeeks_ - pScenario_->thisWeek() - 1) * 7;
  options_.nDaysEvaluation_ =
      std::min(options_.nDaysEvaluation_, remainingDays);

  options_.generationParameters_.maxSolvingTimeSeconds_ =
      options_.totalTimeLimitSeconds_;
  options_.generationParameters_.weekIndices_ = {pScenario_->thisWeek()};

  options_.generationParameters_.verbose_ = options_.verbose_;
  options_.evaluationParameters_.verbose_ = options_.verbose_;

  bestScore_ = LARGE_SCORE;
  bestSchedule_ = -1;
  nGenerationDemands_ = 0;
  nSchedules_ = 0;

  // initialize the log output
  pLogStream_ = new Tools::LogOutput(options_.logfile_);

  // initialize random of tools
  Tools::initializeRandomGenerator();

  if (!options_.generationParameters_.logfile_.empty()) {
    FILE *pFile;
    pFile = fopen(options_.generationParameters_.logfile_.c_str(), "w");
    fprintf(pFile, "HERE COME THE LOGS OF THE MASTER PROBLEM\n\n");
    fclose(pFile);
  }
}

StochasticSolver::~StochasticSolver() {
  // Delete the output
  delete pLogStream_;

  // delete the solvers used for generation
  for (Solver *pSolver : pGenerationSolvers_)
    delete pSolver;

  // delete the solvers used for evaluation
  for (unsigned int n = 0; n < pEvaluationSolvers_.size(); n++)
    for (Solver *pSolver : pEvaluationSolvers_[n])
      delete pSolver;

  // delete also the reusable solvers
  delete pReusableGenerationSolver_;
  for (Solver *pSolver : pReusableEvaluationSolvers_)
    delete pSolver;
}

//----------------------------------------------------------------------------
//
// SOLVE FUNCTIONS
// The one is general for the whole process
//
//----------------------------------------------------------------------------

// Main function
double StochasticSolver::solve(const vector<Roster> &initialSolution) {
  options_.nExtraDaysGenerationDemands_ =
      std::min(options_.nExtraDaysGenerationDemands_,
               7 * (pScenario_->nbWeeks() - (pScenario_->thisWeek() + 1)));
  options_.nDaysEvaluation_ = std::min(options_.nDaysEvaluation_,
                                       7 * (pScenario_->nbWeeks()
                                           - (pScenario_->thisWeek() + 1)));
  // Special case of the last week -> always to optimality with no time limit
  //
  if (pScenario_->nbWeeks() - 1 == pScenario_->thisWeek()) {
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " as the LAST WEEK (hence, to optimality !)" << std::endl;
    // General options
    options_.nExtraDaysGenerationDemands_ = 0;
    options_.withEvaluation_ = false;
    options_.generationCostPerturbation_ = false;
    // Options for the generation algo
    // (-> optimality, no time limit, write every solution)
    options_.generationParameters_.solveToOptimality_ = true;
    options_.generationParameters_.printEverySolution_ = true;
  }

  if (options_.withEvaluation_) {
    // A. Generation-evaluation
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with GENERATION-EVALUATION." << std::endl;
    solveOneWeekGenerationEvaluation();
  } else if (options_.withIterativeDemandIncrease_) {
    // B. Iterative increase in the demand
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with ITERATIVE DEMAND INCREASE." << std::endl;
    solveIterativelyWithIncreasingDemand();
  } else {
    // C. No generation-evaluation
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with PERTURBATIONS." << std::endl;
    solveOneWeekNoGenerationEvaluation();
    while (status_ == INFEASIBLE || status_ == UNSOLVED) {
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] Status is INFEASIBLE || UNSOLVED..." << std::endl;
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] Solving week no. " << pScenario_->thisWeek()
                     << " with PERTURBATIONS -> trying again." << std::endl;
      solveOneWeekNoGenerationEvaluation();
    }
  }

  /* update nurse States */
  if (solution_.empty())
    Tools::throwError("No feasible schedule has been found "
                      "in the available computational time.");
  for (int n = 0; n < pScenario_->nbNurses_; ++n) {
    theLiveNurses_[n]->roster_ = solution_[n];
    theLiveNurses_[n]->buildStates();
  }

  return computeSolutionCost();
}

// Does everything for the one week and only keeps the best schedule for it
void StochasticSolver::solveOneWeekNoGenerationEvaluation() {
  // Need to extend the current demand?
  //
  if (options_.nExtraDaysGenerationDemands_ > 0) {
    generateSingleGenerationDemand();
    pReusableGenerationSolver_ = setSubSolverWithInputAlgorithm(
        pGenerationDemands_[0],
        options_.generationAlgorithm_);
  } else {
    pReusableGenerationSolver_ =
        setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(),
                                       options_.generationAlgorithm_);
  }

  // Need to perturb the costs?
  //
  if (options_.generationCostPerturbation_) {
    pReusableGenerationSolver_->setBoundsAndWeights(
        options_.generationParameters_.weightStrategy_);
  }

  // Solve
  //
  (*pLogStream_) << "# Solve without evaluation\n";
  pReusableGenerationSolver_->solve(options_.generationParameters_);
  if (!options_.withRealDemand_) {
    solution_ = pReusableGenerationSolver_->getSolutionAtDay(6);
  } else {
    solution_ = pReusableGenerationSolver_->getSolutionAtDay(13);
  }
  status_ = pReusableGenerationSolver_->getStatus(true);
}

// Special case of the last week
void StochasticSolver::solveOneWeekWithoutPenalties() {
  Solver *pSolver =
      setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(), GENCOL);
  pSolver->solve(options_.generationParameters_);
  solution_ = pSolver->getSolution();
  status_ = pSolver->getStatus(true);
  delete pSolver;
}

// Solves the problem by generation + evaluation of scenarios
void StochasticSolver::solveOneWeekGenerationEvaluation() {
  while (nSchedules_ < options_.nGenerationDemandsMax_) {
    // get the time left to solve another schedule
    double
        timeLeft = options_.totalTimeLimitSeconds_ - pTimerTotal_->dSinceInit();
    if (nSchedules_ > 0) {
      if (timeLeft < 1.0) break;
      // options_.generationParameters_.maxSolvingTimeSeconds_  =
      // (timeLeft-1.0)/2.0;
      // options_.evaluationParameters_.maxSolvingTimeSeconds_  =
      // (timeLeft-1.0)/(2.0*options_.nEvaluationDemands_);
    }
    (*pLogStream_) << "# Time left: " << timeLeft << std::endl;

    // This the main function that finds a new schedule and evaluates it
    bool printOption = options_.generationParameters_.printEverySolution_;
    if (nSchedules_ == 0) {
      options_.generationParameters_.printEverySolution_ = true;
      options_.generationParameters_.weekIndices_ = {pScenario_->thisWeek()};
    }

    if (addAndSolveNewSchedule()) {
      // Get the new best schedule
      //
      int newBestSchedule = -1;
      double newBestScore = LARGE_SCORE;
      double bestBaseCost = 0;
      for (int i = 0; i < nSchedules_; i++) {
        if ((options_.demandingEvaluation_ && theScores_[i] + 30 < newBestScore)
            || (!options_.demandingEvaluation_
                && theScores_[i] < newBestScore)) {
          newBestScore = theScores_[i];
          newBestSchedule = i;
          bestBaseCost = theBaseCosts_[i];
        } else if (theScores_[i] == newBestScore
            && theBaseCosts_[i] < bestBaseCost) {
          newBestScore = theScores_[i];
          newBestSchedule = i;
          bestBaseCost = theBaseCosts_[i];
        }
      }

      // write the output NOW so that it is not lost
      //
      bestScore_ = newBestScore;
      if (newBestSchedule != bestSchedule_) {
        bestSchedule_ = newBestSchedule;

        // solution_ = pGenerationSolvers_[bestSchedule_]->getSolutionAtDay(6);
        solution_ = schedules_[bestSchedule_];
        loadSolution(solution_);
        Tools::LogOutput outStream(options_.generationParameters_.outfile_);
        outStream << solutionToString();

        (*pLogStream_) << "# New best is schedule n°" << bestSchedule_
                       << " (score: " << bestScore_ << ")" << std::endl;
        (*pLogStream_) << "# The new best solution was written in "
                       << options_.generationParameters_.outfile_ << std::endl;

      } else {
        (*pLogStream_) << "# Best schedule did not change and is no. "
                       << bestSchedule_ << " (score: " << bestScore_ << ")"
                       << std::endl;
      }
    } else {
      (*pLogStream_) << "# Time has run out." << std::endl;
    }

    options_.generationParameters_.printEverySolution_ = printOption;
  }
#ifdef COMPARE_EVALUATIONS
  for (int i = 0; i < nSchedules_; i++) {
    (*pLogStream_) << " The score of schedule " << i << ". GENCOL : "
                   << theScores_[i] << " ; GREEDY : " << theScoresGreedy_[i]
                   << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------
//
// Iterative solution process in which the week is first solved by itself,
// before adding one perturbated week demand and solving the new extended
// demand demand until no time is left
//
//----------------------------------------------------------------------------

void StochasticSolver::solveIterativelyWithIncreasingDemand() {
  // Set the options corresponding to this algorithm
  options_.withEvaluation_ = false;
  options_.generationCostPerturbation_ = true;

  // Initialize the values that intervene in the stopping criterion
  double timeLeft =
      options_.totalTimeLimitSeconds_ - pTimerTotal_->dSinceInit();
  Tools::Timer *timerSolve = new Tools::Timer();
  double timeLastSolve = 0.0;
  int maxNbAddedWeeks = pScenario_->nbWeeks() - (pScenario_->thisWeek() + 1);
  int nbAddedWeeks = 0;

  // Launch the iterative process
  vector<Roster> previousSolution;
  while (timeLeft > timeLastSolve && nbAddedWeeks <= maxNbAddedWeeks) {
    (*pLogStream_) << "# Solve with " << nbAddedWeeks
                   << " additional weeks to the demand" << std::endl;
    (*pLogStream_) << "# Time left: " << timeLeft << std::endl;

    // Update the properties of the solver
    options_.generationParameters_.maxSolvingTimeSeconds_ = timeLeft - 1.0;
    options_.nExtraDaysGenerationDemands_ = 7 * nbAddedWeeks;

    // Solve the week with no evaluation and
    // nbAddedWeek extra weeks in the demand
    timerSolve->start();
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with PERTURBATIONS." << std::endl;
    solveOneWeekNoGenerationEvaluation();
    if (nbAddedWeeks > 0) {
      while (status_ == INFEASIBLE || status_ == UNSOLVED) {
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Status is INFEASIBLE || UNSOLVED..." << std::endl;
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Solving week no. " << pScenario_->thisWeek()
                       << " with PERTURBATIONS -> trying again." << std::endl;
        solveOneWeekNoGenerationEvaluation();
      }
      // Go back to the last solution if the solver was interrupted
      timeLeft = options_.totalTimeLimitSeconds_ - pTimerTotal_->dSinceInit();
      if (timeLeft <= 1.0) {
        (*pLogStream_)
            << "# The execution had to be interrupted, "
               "so the solution is not kept"
            << std::endl;
        loadSolution(previousSolution);
      } else {
        (*pLogStream_) << "# New schedule based on extended demand with "
                       << nbAddedWeeks << " extra weeks" << std::endl;
        loadSolution(solution_);
        previousSolution = solution_;
        Tools::LogOutput outStream(options_.generationParameters_.outfile_);
        outStream << solutionToString();
      }

      // Delete and popback the last generation demand to be consistent
      // with the implementation of solveOneWeekNoGenerationEvaluation
      pGenerationDemands_.pop_back();
      nGenerationDemands_--;
    } else if (status_ == INFEASIBLE || status_ == UNSOLVED) {
      Tools::throwError("# solveIterativelyWithIncreasingDemand: "
                        "no solution was found for this instance!");
    } else {
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] The demand is feasible, "
                        "write the schedule based on one week"
                     << std::endl;
      previousSolution = solution_;
      Tools::LogOutput outStream(options_.generationParameters_.outfile_);
      outStream << solutionToString();
    }

    timerSolve->stop();
    timeLastSolve = timerSolve->dSinceStart();
    nbAddedWeeks++;
  }
}


//----------------------------------------------------------------------------
//
// GENERIC FUNCTION TO DO EVERYTHING FOR ONE SCHEDULE
// Includes generation, evaluation of the score, and update of the rankings
// and data.
//
//----------------------------------------------------------------------------

// Do everything for the new schedule (incl. generation, score, ranking)
bool StochasticSolver::addAndSolveNewSchedule() {
  generateNewSchedule();

//  std::cout << pReusableGenerationSolver_->solutionToLogString() << std::endl;

  if (nSchedules_ == 1)
    generateAllEvaluationDemands();
  return evaluateSchedule(nSchedules_ - 1);
}


//----------------------------------------------------------------------------
//
// GENERATION OF DEMANDS FOR THE CURRENT WEEK (=FOR SCHEDULE GENERATION)
// Note that these demands share a common first week which is the week
// we currently try to solve.
//
//----------------------------------------------------------------------------

// Generate a new demand for generation
void StochasticSolver::generateSingleGenerationDemand() {
  int nDaysInDemand = options_.nExtraDaysGenerationDemands_;
  bool isFeasible = false;
  PDemand pCompleteDemand;

  (*pLogStream_) << "# Generating new generation demand..." << std::endl;

  // use the real demand instead of generating a random one
  if (options_.withRealDemand_) {
    pCompleteDemand = std::make_shared<Demand>(*(pScenario_->pWeekDemand()));
    PDemand pFutureDemand = std::make_shared<Demand>(*(demandHistory_[0]));
    pFutureDemand->keepFirstNDays(nDaysInDemand);
    pCompleteDemand = pCompleteDemand->append(pFutureDemand);
  } else {
    PDemand pSingleDemand;
    while (!isFeasible) {
      DemandGenerator dg(1, nDaysInDemand, demandHistory_, pScenario_);
      // no feasibility check here
      pSingleDemand = dg.generateSinglePerturbatedDemand(false);
      pCompleteDemand = pScenario_->pWeekDemand()->append(pSingleDemand);
      isFeasible = true;
//         isFeasible = dg.checkDemandFeasibility(pCompleteDemand);
//         if(!isFeasible){
//            (*pLogStream_) << "# Demand has been deleted because "
//                              "it was infeasible." << std::endl;
//            delete pCompleteDemand;
//            delete pSingleDemand;
//         }
    }
  }

  pGenerationDemands_.push_back(pCompleteDemand);
  nGenerationDemands_++;
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Generation demand no. " << (nGenerationDemands_ - 1)
                 << " created (over "
                 << pGenerationDemands_[nGenerationDemands_ - 1]->nbDays_
                 << " days)." << std::endl;
}



//----------------------------------------------------------------------------
//
// GENERATION OF SCENARIOS FOR THE FUTURE (=FOR SCHEDULE EVALUATION)
//
//----------------------------------------------------------------------------

// Generate the schedules that are used for evaluation
void StochasticSolver::generateAllEvaluationDemands() {
  DemandGenerator dg(options_.nEvaluationDemands_,
                     options_.nDaysEvaluation_,
                     demandHistory_,
                     pScenario_);
  pEvaluationDemands_ = dg.generatePerturbedDemands();
  // Initialize structures for scores
  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
    map<double, set<int> > m;
    schedulesFromObjectiveByEvaluationDemand_.push_back(m);

#ifdef COMPARE_EVALUATIONS
    schedulesFromObjectiveByEvaluationDemandGreedy_.push_back(m);
#endif

    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Evaluation demand no. " << j << " created (over "
                   << options_.nDaysEvaluation_ << " days)." << std::endl;
  }
}



//----------------------------------------------------------------------------
//
// GENERATION OF SCHEDULES
// A solution is a potential candidate to be the chosen schedule
// for the week we are solving.
// A result is, given a solution and a potential future, the value obtained
// for that couple (solution,demand) [i.e. LP bound for instance]
// A score is, given a solution, the average score it obtains, compared to the
// other solutions (the precise meaning of "score" should be better defined)
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule GENERATION
Solver *StochasticSolver::setGenerationSolverWithInputAlgorithm(
    PDemand pDemand) {
  switch (options_.generationAlgorithm_) {
    case GENCOL:
      return new RotationMP(pScenario_,
                            pDemand,
                            pScenario_->pWeekPreferences(),
                            pScenario_->pInitialState(),
                            S_CLP);
    default:Tools::throwError("The algorithm is not handled yet");
      break;
  }
  return nullptr;
}

// Generate a new schedule
void StochasticSolver::generateNewSchedule() {
  bool hasFoundFeasible = false;

  while (!hasFoundFeasible) {
    // A. Generate a demand that will be the origin of the scenario generation
    //
    generateSingleGenerationDemand();
    PDemand newDemand = pGenerationDemands_[nGenerationDemands_ - 1];

    // B. Solve this schedule (in a way that should be defined)
    // so as to have a schedule
    //
    // Create a new solver if first schedule || if RE-solve is forbidden
    if (nSchedules_ == 0 || !(options_.withResolveForGeneration_)) {
      if (nSchedules_ > 0) delete pReusableGenerationSolver_;
      pReusableGenerationSolver_ =
          setGenerationSolverWithInputAlgorithm(newDemand);
    }

    if (options_.generationCostPerturbation_) {
      pReusableGenerationSolver_->setBoundsAndWeights(
          options_.generationParameters_.weightStrategy_);
    }

    // If first || no RE-solve, solve normally.
    // Otherwise, re-solve with a new demand
    if (nSchedules_ == 0 || !(options_.withResolveForGeneration_))
      pReusableGenerationSolver_->solve(options_.generationParameters_);
    else
      pReusableGenerationSolver_->resolve(newDemand,
                                          options_.generationParameters_);

    if (pReusableGenerationSolver_->getStatus(true) == FEASIBLE
        || pReusableGenerationSolver_->getStatus() == OPTIMAL) {
      hasFoundFeasible = true;
    } else {
      nGenerationDemands_--;
      pGenerationDemands_.pop_back();
    }
  }

  // C. Store the solution
  //
  schedules_.push_back(pReusableGenerationSolver_->getSolutionAtDay(6));
  finalStates_.push_back(pReusableGenerationSolver_->getStatesOfDay(6));

  // D. Update the data
  //
  nSchedules_++;

  // E. Display
  //
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Candidate schedule no. " << (nSchedules_ - 1)
                 << " generated: (length: "
                 << pReusableGenerationSolver_->getNbDays() << " days)"
                 << std::endl;
}



//----------------------------------------------------------------------------
//
// EVALUATION OF SCHEDULES
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule EVALUATION
Solver *StochasticSolver::setEvaluationWithInputAlgorithm(
    PDemand pDemand, vector<State> *stateEndOfSchedule) {
  Solver *pSolver(nullptr);
  PScenario pScen = std::make_shared<Scenario>(*pScenario_);

  // update the scenario to treat next week
  PPreferences pEmptyPref = std::make_shared<Preferences>(
      pScenario_->nbNurses(),
      options_.nDaysEvaluation_,
      pScenario_->nbShifts());
  pScen->updateNewWeek(pDemand, pEmptyPref, *stateEndOfSchedule);

  switch (options_.evaluationAlgorithm_) {
    case GENCOL:
      pSolver = new RotationMP(pScen,
                               pDemand,
                               pEmptyPref,
                               stateEndOfSchedule,
                               S_CLP);
      break;
    default:Tools::throwError("The algorithm is not handled yet");
      break;
  }
  return pSolver;
}

// Initialization
void StochasticSolver::initScheduleEvaluation(int sched) {
  // Extend pEvaluationSolvers_
  vector<Solver *> v;
//  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
//    Solver *s;
//    v.push_back(s);
//  }
//  pEvaluationSolvers_.push_back(v);
  Solver *so(0);
  pReusableEvaluationSolvers_.push_back(so);
}

// Evaluate 1 schedule on all evaluation instances
bool StochasticSolver::evaluateSchedule(int sched) {
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Evaluation of the schedule no. " << sched << std::endl;

#ifdef COMPARE_EVALUATIONS
  vector<Solver *> pGreedyEvaluators;
  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
    Solver *s;
    pGreedyEvaluators.push_back(s);
  }
#endif

  initScheduleEvaluation(sched);
  vector<State> initialStates = finalStates_[nSchedules_ - 1];
  for (int i = 0; i < pScenario_->nbNurses_; i++)
    initialStates[i].dayId_ = 0;

  int baseCost = pReusableGenerationSolver_->computeSolutionCost(7);
  theBaseCosts_.push_back(baseCost);

  // set the time per evaluation to the ratio of the time left
  // over the number of evaluations
  // double timeLeft =
  // options_.totalTimeLimitSeconds_-pTimerTotal_->dSinceInit();
  // options_.evaluationParameters_.maxSolvingTimeSeconds_ =
  // (timeLeft-1.0)/(double)options_.nEvaluationDemands_;
  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
    double timeLeft =
        options_.totalTimeLimitSeconds_ - pTimerTotal_->dSinceInit();
    if (nSchedules_ > 0) {
      if (timeLeft < 1.0) {
        std::cout << "# Time has run out when evaluating schedule no."
                  << (nSchedules_ - 1) << std::endl;
        std::cout << options_.totalTimeLimitSeconds_ << "; "
                  << pTimerTotal_->dSinceInit() << std::endl;
        return false;
      }
    }

    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Starting evaluation of schedule no. " << sched
                   << " over evaluation demand no. " << j << std::endl;

    if (j == 0) {
      pReusableEvaluationSolvers_[sched] = setEvaluationWithInputAlgorithm(
          pEvaluationDemands_[j], &initialStates);
    }

#ifdef COMPARE_EVALUATIONS
    options_.evaluationAlgorithm_ = GREEDY;
    pGreedyEvaluators[j] = setEvaluationWithInputAlgorithm(
        pEvaluationDemands_[j], & initialStates);
    options_.evaluationAlgorithm_ = GENCOL;
#endif

    if (options_.evaluationCostPerturbation_) {
      if (pReusableEvaluationSolvers_[sched]->getNbDays()
          + (7 * pScenario_->thisWeek() + 1) < 7 * pScenario_->nbWeeks_) {
        pReusableEvaluationSolvers_[sched]->setBoundsAndWeights(
            options_.evaluationParameters_.weightStrategy_);
#ifdef COMPARE_EVALUATIONS
        pGreedyEvaluators[j]->computeWeightsTotalShiftsForStochastic();
#endif
      }
    }

    // Only perform the evaluation if the schedule is feasible and
    // there is time for more than one schedule
    double currentCost = costPreviousWeeks_ + baseCost;

#ifdef COMPARE_EVALUATIONS
    double currentCostGreedy = costPreviousWeeks_ + baseCost;
#endif

    if (pReusableGenerationSolver_->getStatus(true) == INFEASIBLE) {
      currentCost = 1.0e6;
#ifdef COMPARE_EVALUATIONS
      currentCostGreedy = 1.0e6;
#endif
    } else {
      // Perform the actual evaluation on demand j
      // by running the chosen algorithm
      if (j == 0) {
        currentCost += pReusableEvaluationSolvers_[sched]->solve(
            options_.evaluationParameters_);
      } else {
        currentCost += pReusableEvaluationSolvers_[sched]->resolve(
            pEvaluationDemands_[j],
            options_.evaluationParameters_);
      }

#ifdef COMPARE_EVALUATIONS
      pGreedyEvaluators[j]->solve();
      currentCostGreedy += pGreedyEvaluators[j]->computeSolutionCost();
#endif
    }

    // Display
    //
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Schedule no. "
                   << sched << " evaluated over evaluation demand no. " << j
                   << " (solution cost: " << currentCost << ")." << std::endl;

    // Insert the solution cost and solution
    //
    // If already in the costs -> add it to the set of schedules
    // that found that cost
    if (schedulesFromObjectiveByEvaluationDemand_[j].find(currentCost)
        != schedulesFromObjectiveByEvaluationDemand_[j].end()) {
      schedulesFromObjectiveByEvaluationDemand_[j].at(
          currentCost).insert(sched);
    } else {
      // Otherwise, add a new pair
      set<int> s;
      s.insert(sched);
      schedulesFromObjectiveByEvaluationDemand_[j].insert(
          pair<double, set<int> >(currentCost, s));
    }
#ifdef COMPARE_EVALUATIONS
    if (schedulesFromObjectiveByEvaluationDemandGreedy_[j].find(
        currentCostGreedy) !=
        schedulesFromObjectiveByEvaluationDemandGreedy_[j].end()) {
      schedulesFromObjectiveByEvaluationDemandGreedy_[j].at(
          currentCostGreedy).insert(sched);
    } else {
      set<int> s;
      s.insert(sched);
      schedulesFromObjectiveByEvaluationDemandGreedy_[j].insert(
          pair<double, set<int> >(currentCostGreedy, s));
    }
#endif
  }

  (*pLogStream_) << "# Evaluation of schedule no. " << sched << " done!"
                 << std::endl;

  updateRankingsAndScores(options_.rankingStrategy_);

  return true;
}

// Recompute all scores after one schedule evaluation
void StochasticSolver::updateRankingsAndScores(RankingStrategy strategy) {
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Starting the update of the scores and ranking."
                 << std::endl;
  vector<double> theNewScores(nSchedules_, 0);

  if (options_.nEvaluationDemands_ == 0) {
    for (int sched = 0; sched < nSchedules_; sched++) {
      theNewScores[sched] = theBaseCosts_[sched];
    }
  }

  switch (strategy) {
    case RK_SCORE:
      for (int j = 0; j < options_.nEvaluationDemands_; j++) {
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Solution costs for demand no. " << j << std::endl;
        int localRank = 1;
        map<double, set<int> > localCosts =
            schedulesFromObjectiveByEvaluationDemand_[j];
        for (auto it = localCosts.begin(); it != localCosts.end(); ++it) {
          for (int sched : it->second) {
            theNewScores[sched] += localRank + (it->second.size() - 1)
                * 1.0 / it->second.size();
            (*pLogStream_) << "#     | sched " << sched << " -> " << it->first
                           << " (score += " << localRank
                               + (it->second.size() - 1)
                                   * 1.0 / it->second.size() << ")"
                           << std::endl;
          }
          localRank += it->second.size();
        }
      }
      break;
    case RK_MEAN:
      for (int j = 0; j < options_.nEvaluationDemands_; j++) {
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Solution costs for demand no. " << j << std::endl;
        map<double, set<int> >
            localCosts = schedulesFromObjectiveByEvaluationDemand_[j];
        for (pair<double, set<int> > p : localCosts)
          for (int sched : p.second) {
            theNewScores[sched] +=
                static_cast<int>(p.first / options_.nEvaluationDemands_);
            (*pLogStream_) << "#     | sched " << sched << " -> " << p.first
                           << std::endl;
          }
      }
      break;
    case RK_NONE:Tools::throwError("Ranking strategy set to NONE.");
      break;
    default:Tools::throwError("Ranking strategy not defined.");
  }

#ifdef COMPARE_EVALUATIONS
  vector<double> theNewScoresGreedy;
  Tools::initVector(&theNewScoresGreedy, nSchedules_, .0);
  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solution costs for demand no. " << j << std::endl;
    int localRank = 1;
    map<double, set<int> > localCosts =
        schedulesFromObjectiveByEvaluationDemandGreedy_[j];
    for (auto it = localCosts.begin(); it != localCosts.end(); ++it) {
      for (int sched : it->second) {
        theNewScoresGreedy[sched] +=
            localRank + (it->second.size() - 1) * 1.0 / it->second.size();
        // If is infeasible -> double that amount to get more robust
        if ((pEvaluationSolvers_[sched][j])->getStatus() == INFEASIBLE) {
          theNewScoresGreedy[sched] +=
              localRank + (it->second.size() - 1) * 1.0 / it->second.size();
        }
        (*pLogStream_) << "#     | sched " << sched << " -> " << it->first
                       << " (score += " << localRank
                           + (it->second.size() - 1) * 1.0 / it->second.size()
                       << ")" << std::endl;
      }
      localRank += it->second.size();
    }
  }
  theScoresGreedy_ = theNewScoresGreedy;
#endif

  theScores_ = theNewScores;
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Update of the scores and ranking done!" << std::endl;
}

Solver *StochasticSolver::setSubSolverWithInputAlgorithm(
    PDemand pDemand, Algorithm algorithm) {
  Solver *pSolver(nullptr);
  switch (algorithm) {
    case GENCOL:
      pSolver = new RotationMP(pScenario_,
                               pDemand,
                               pScenario_->pWeekPreferences(),
                               pScenario_->pInitialState(),
                               S_CLP);
      break;
    default:Tools::throwError("The algorithm is not handled yet");
      break;
  }
  return pSolver;
}
