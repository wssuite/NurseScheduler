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
    Solver(pScenario),
    options_(options),
    demandHistory_(std::move(demandHistory)),
    pGenerationSolver_(nullptr),
    nGeneratedSolutions_(0),
    costPreviousWeeks_(costPreviousWeeks) {
  if (!pScenario_->isINRC2_)
    Tools::throwError("StochasticSolver works only with INRC2 scenarios.");

  std::cout << "# New stochastic solver created!" << std::endl;

  int remainEvalDays = (pScenario_->nWeeks() - pScenario_->thisWeek() - 1) * 7;
  options_.nDaysEvaluation_ =
      std::min(options_.nDaysEvaluation_, remainEvalDays);

  options_.generationParameters_.maxSolvingTimeSeconds_ =
      options_.totalTimeLimitSeconds_;
  options_.generationParameters_.weekIndices_ = {pScenario_->thisWeek()};

  options_.generationParameters_.verbose_ = options_.verbose_;
  options_.evaluationParameters_.verbose_ = options_.verbose_;

  options_.generationParameters_.minRelativeGap_ = 1e-6;
  options_.generationParameters_.relativeGap_ = 1e-6;

  if (options_.rankingStrategy_ == RK_SCORE && options_.demandingEvaluation_) {
    options_.demandingEvaluation_ = false;
    std::cout <<"Deactivate parameter demandingEvaluation_ as it is not "
                "compatible with the ranking strategy RK_SCORE." << std::endl;
  }

  bestScore_ = XLARGE_SCORE;
  bestSchedule_ = -1;

  // initialize the log output
  pLogStream_ =
          std::make_unique<Tools::LogOutput>(options_.logfile_, true, true);

  if (!options_.generationParameters_.logfile_.empty()) {
    Tools::LogOutput log(options_.generationParameters_.logfile_);
    log << "HERE COME THE LOGS OF THE MASTER PROBLEM" << std::endl << std::endl;
  }
}

StochasticSolver::~StochasticSolver() {}

//----------------------------------------------------------------------------
//
// SOLVE FUNCTIONS
// The one is general for the whole process
//
//----------------------------------------------------------------------------

// Main function
double StochasticSolver::solve(const vector<Roster> &initialSolution) {
  // set the number of days in the horizon for generation and evaluation
  options_.nExtraDaysGenerationDemands_ =
      std::min(options_.nExtraDaysGenerationDemands_,
               7 * (pScenario_->nWeeks() - (pScenario_->thisWeek() + 1)));
  options_.nDaysEvaluation_ = std::min(options_.nDaysEvaluation_,
                                       7 * (pScenario_->nWeeks()
                                           - (pScenario_->thisWeek() + 1)));

  // Special case of the last week -> always to optimality with no time limit
  //
  if (pScenario_->nWeeks() - 1 == pScenario_->thisWeek()) {
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " as the LAST WEEK (hence, to optimality !)" << std::endl;
    // General options
    options_.nExtraDaysGenerationDemands_ = 0;
    options_.withEvaluation_ = false;
    options_.generationCostPerturbation_ = false;
    options_.withIterativeDemandIncrease_ = false;
    // Options for the generation algo
    // (-> optimality, no time limit, write every solution)
    options_.generationParameters_.optimalityLevel(OPTIMALITY);
  }

  if (options_.withEvaluation_) {
    // A. Generation-evaluation
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with GENERATION-EVALUATION." << std::endl;
    solveOneWeekGenerationEvaluation();
    (*pLogStream_) << "# Best is schedule n°" << bestSchedule_
                   << " (score: " << bestScore_ << ")" << std::endl;
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
                   << " with"
                   << (!options_.generationCostPerturbation_ ? "out" : "")
                   << " PERTURBATIONS." << std::endl;
    solveOneWeekNoGenerationEvaluation();
    while (status_ == INFEASIBLE || status_ == UNSOLVED) {
      // get the time left to solve another schedule
      double timeLeft =
          options_.totalTimeLimitSeconds_ - timerTotal_.dSinceInit();
      if (timeLeft < 1.0) {
        (*pLogStream_) << "# Time has run out." << std::endl;
        break;
      }
      (*pLogStream_) << "# Time left: " << timeLeft << std::endl;
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] Status is INFEASIBLE || UNSOLVED..." << std::endl;
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] Solving week no. " << pScenario_->thisWeek()
                     << " with PERTURBATIONS -> trying again." << std::endl;
      options_.generationParameters_.maxSolvingTimeSeconds_ = timeLeft;
      solveOneWeekNoGenerationEvaluation();
    }
    loadSolution(solution_);
    solutionToTxt(options_.generationParameters_.outdir_);
  }

  /* update nurse States */
  if (solution_.empty())
    Tools::throwError("No feasible schedule has been found "
                      "in the available computational time.");
  for (int n = 0; n < pScenario_->nNurses(); ++n) {
    theLiveNurses_[n]->roster_ = solution_[n];
    theLiveNurses_[n]->buildStates();
  }

  return costPreviousWeeks_ + computeSolutionCost();
}

// Does everything for the one week and only keeps the best schedule for it
void StochasticSolver::solveOneWeekNoGenerationEvaluation() {
  // Need to extend the current demand?
  auto pDemand = options_.nExtraDaysGenerationDemands_ > 0 ?
      generateSingleGenerationDemand() : pScenario_->pDemand();
  pGenerationSolver_ = std::unique_ptr<Solver>(
          setGenerationSolverWithInputAlgorithm(pDemand));

  // Need to perturb the costs?
  int endSchedule = pGenerationSolver_->pDemand()->nDays_
      + 7 * pScenario_->thisWeek();
  if (options_.generationCostPerturbation_ &&
      endSchedule < 7 * pScenario_->nWeeks()) {
    pGenerationSolver_->boundsAndWeights(
        options_.generationParameters_.weightStrategy_);
  }

  // Solve
  (*pLogStream_) << "# Solve without evaluation\n";
  pGenerationSolver_->solve(options_.generationParameters_);
  copySolution(pGenerationSolver_.get());
  solution_ = solutionAtDay(options_.withRealDemand_ ? 13 : 6);
  status_ = pGenerationSolver_->status(true);
}

// Special case of the last week
void StochasticSolver::solveOneWeekWithoutPenalties() {
  auto pSolver = std::unique_ptr<Solver>(
      setGenerationSolverWithInputAlgorithm(pScenario_->pDemand()));
  pSolver->solve(options_.generationParameters_);
  copySolution(pSolver.get());
  status_ = pSolver->status(true);
}

// Solves the problem by generation + evaluation of scenarios
void StochasticSolver::solveOneWeekGenerationEvaluation() {
  // ensure to keep enough time for evaluation for the first schedule
  options_.generationParameters_.maxSolvingTimeSeconds_ =
      options_.totalTimeLimitSeconds_ / 2;
  // generate schedules while enough time
  double timeLeft = 2;
  while (timeLeft >= 1.0) {
    // get the time left to solve another schedule
    timeLeft = options_.totalTimeLimitSeconds_ - timerTotal_.dSinceInit();
    (*pLogStream_) << "# Time left: " << timeLeft << std::endl;

    // print the first solution found ASAP
    bool printOption = options_.generationParameters_.printEverySolution_;
    if (nGeneratedSolutions_ == 0) {
      options_.generationParameters_.printEverySolution_ = true;
      options_.generationParameters_.weekIndices_ = {pScenario_->thisWeek()};
    }

    // This the main function that finds a new schedule and evaluates it
    if (addAndSolveNewSchedule()) {
      // Get the new best schedule
      int newBestSchedule = -1;
      double newBestScore = XLARGE_SCORE;
      double bestBaseCost = 0;
      double minGap = options_.demandingEvaluation_ ?
                      pScenario_->weights().underCoverage : epsilon();
      for (int i = 0; i < nGeneratedSolutions_; i++) {
        if (theScores_[i] + minGap < newBestScore) {
          newBestScore = theScores_[i];
          newBestSchedule = i;
          bestBaseCost = theBaseCosts_[i];
        } else if (!options_.demandingEvaluation_ &&
            theScores_[i] < newBestScore + epsilon()
            && theBaseCosts_[i] + epsilon() < bestBaseCost) {
          newBestScore = theScores_[i];
          newBestSchedule = i;
          bestBaseCost = theBaseCosts_[i];
        }
      }

      // write the output NOW so that it is not lost
      bestScore_ = newBestScore;
      if (newBestSchedule != bestSchedule_) {
        bestSchedule_ = newBestSchedule;
        copySolution(pGenerationSolver_.get());
        solutionToTxt(options_.generationParameters_.outdir_);

        (*pLogStream_) << "# New best is schedule n°" << bestSchedule_
                       << " (score: " << bestScore_ << ")" << std::endl;
        (*pLogStream_) << "# The new best solution was written in "
                       << options_.generationParameters_.outdir_ << std::endl;

      } else {
        (*pLogStream_) << "# Best schedule did not change and is no. "
                       << bestSchedule_ << " (score: " << bestScore_ << ")"
                       << std::endl;
      }
    }
    options_.generationParameters_.printEverySolution_ = printOption;
  }
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
      options_.totalTimeLimitSeconds_ - timerTotal_.dSinceInit();
  Tools::Timer timerSolve("StochasticSolver with increasing demand");
  double timeLastSolve = 0.0;
  int maxNbAddedWeeks = pScenario_->nWeeks() - (pScenario_->thisWeek() + 1);
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
    timerSolve.start();
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Solving week no. " << pScenario_->thisWeek()
                   << " with"
                   << (!options_.generationCostPerturbation_ ? "out" : "")
                   << " PERTURBATIONS." << std::endl;
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
      timeLeft = options_.totalTimeLimitSeconds_ - timerTotal_.dSinceInit();
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
        solutionToTxt(options_.generationParameters_.outdir_);
        previousSolution = solution_;
      }

      // Delete and popback the last generation demand to be consistent
      // with the implementation of solveOneWeekNoGenerationEvaluation
      pGenerationDemands_.pop_back();
    } else if (status_ == INFEASIBLE || status_ == UNSOLVED) {
      Tools::throwError("# solveIterativelyWithIncreasingDemand: "
                        "no solution was found for this instance!");
    } else {
      (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                     << "] The demand is feasible, "
                        "write the schedule based on one week"
                     << std::endl;
      previousSolution = solution_;
      solutionToTxt(options_.generationParameters_.outdir_);
    }

    timerSolve.stop();
    timeLastSolve = timerSolve.dSinceStart();
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
  if (!generateNewSchedule())
    return false;

  if (pEvaluationDemands_.empty())
    generateAllEvaluationDemands();

  return evaluateLastGeneratedSchedule();
}


//----------------------------------------------------------------------------
//
// GENERATION OF DEMANDS FOR THE CURRENT WEEK (=FOR SCHEDULE GENERATION)
// Note that these demands share a common first week which is the week
// we currently try to solve.
//
//----------------------------------------------------------------------------

// Generate a new demand for generation
PDemand StochasticSolver::generateSingleGenerationDemand() {
  int nDaysInDemand = options_.nExtraDaysGenerationDemands_;
  bool isFeasible = false;
  PDemand pCompleteDemand;

  (*pLogStream_) << "# Generating new generation demand..." << std::endl;

  // use the real demand instead of generating a random one
  if (options_.withRealDemand_) {
    pCompleteDemand = std::make_shared<Demand>(*(pScenario_->pDemand()));
    PDemand pFutureDemand = std::make_shared<Demand>(*(demandHistory_[0]));
    pFutureDemand->keepFirstNDays(nDaysInDemand);
    pCompleteDemand = pCompleteDemand->append(pFutureDemand);
  } else {
    PDemand pSingleDemand;
    while (!isFeasible) {
      DemandGenerator dg(1, nDaysInDemand, demandHistory_, pScenario_);
      // no feasibility check here
      pSingleDemand = dg.generateSinglePerturbatedDemand(false);
      pCompleteDemand = pScenario_->pDemand()->append(pSingleDemand);
      isFeasible = dg.checkDemandFeasibility(pCompleteDemand);
      if (!isFeasible) {
        (*pLogStream_) << "# Demand has been regenerated because "
                          "it was infeasible." << std::endl;
      }
    }
  }

  pGenerationDemands_.push_back(pCompleteDemand);
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Generation demand no. "
                 << (pGenerationDemands_.size() - 1)
                 << " created (over "
                 << pGenerationDemands_.back()->nDays_
                 << " days)." << std::endl;

  return pCompleteDemand;
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
  int n = options_.nEvaluationDemands_;
  schedulesFromObjectiveByEvaluationDemand_.resize(n);
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] "
                   << n << " evaluation demands have been created (over "
                   << options_.nDaysEvaluation_ << " days)." << std::endl;
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
  PScenario pScenario = std::make_shared<Scenario>(*pScenario_);
  pScenario->linkWithDemand(std::move(pDemand));
  Solver *pSolver = newSolver(pScenario, options_.generationAlgorithm_,
                              options_.generationParameters_.spType_,
                              options_.mySolverType_);
  pSolver->computeBoundsAccordingToInitialState();
  return pSolver;
}

// Generate a new schedule
bool StochasticSolver::generateNewSchedule() {
  // A. Generate a demand that will be concatenated to the real demand
  // for the generation
  PDemand newDemand = generateSingleGenerationDemand();

  // B. Solve this schedule (in a way that should be defined)
  // Create a new solver if first schedule || if RE-solve is forbidden
  if (nGeneratedSolutions_ == 0 || !options_.withResolveForGeneration_) {
    WeightStrategy previous_strat = pGenerationSolver_ ?
        pGenerationSolver_->getDynamicWeights().getStrategy() : NO_STRAT;
    pGenerationSolver_ = std::unique_ptr<Solver>(
            setGenerationSolverWithInputAlgorithm(newDemand));
    pGenerationSolver_->setDynamicWeightsStrategy(previous_strat);
  }

  int endSchedule = 7 * pScenario_->thisWeek() + newDemand->nDays_;
  if (options_.generationCostPerturbation_ &&
      endSchedule < 7 * pScenario_->nWeeks())
    pGenerationSolver_->boundsAndWeights(
        options_.generationParameters_.weightStrategy_);

  // If first || no RE-solve, solve normally.
  // Otherwise, re-solve with a new demand
  if (pGenerationSolver_->status() == UNSOLVED)
    pGenerationSolver_->solve(options_.generationParameters_);
  else
    pGenerationSolver_->resolve(newDemand, options_.generationParameters_);

  if (pGenerationSolver_->status(true) != FEASIBLE
      && pGenerationSolver_->status() != OPTIMAL) {
    pGenerationDemands_.pop_back();
    return false;
  }

  // C. Store the solution at the end of the first week
  theBaseCosts_.push_back(pGenerationSolver_->computeSolutionCost(7));
  nGeneratedSolutions_++;

  // D. Display
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Candidate schedule no. " << (nGeneratedSolutions_ - 1)
                 << " generated: (length: "
                 << pGenerationSolver_->nDays() << " days)"
                 << std::endl;

  return true;
}



//----------------------------------------------------------------------------
//
// EVALUATION OF SCHEDULES
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule EVALUATION
Solver *StochasticSolver::setEvaluationWithInputAlgorithm(
    PDemand pDemand, const vector<State> &stateEndOfSchedule) {
  PScenario pScen = std::make_shared<Scenario>(*pScenario_);

  // update the scenario to treat next week
  PPreferences pEmptyPref = std::make_shared<Preferences>(
      pScenario_->nNurses(),
      options_.nDaysEvaluation_,
      pScenario_->nShifts());
  pScen->updateNewWeek(std::move(pDemand), pEmptyPref, stateEndOfSchedule);

  Solver *pSolver = newSolver(pScen, options_.evaluationAlgorithm_,
                              options_.evaluationParameters_.spType_,
                              options_.mySolverType_);

  if (options_.evaluationCostPerturbation_)
    pSolver->computeBoundsAccordingToDemandSize();
  else
    pSolver->computeBoundsAccordingToInitialState();

  return pSolver;
}

// Evaluate 1 schedule on all evaluation instances
bool StochasticSolver::evaluateLastGeneratedSchedule() {
  int sched = nGeneratedSolutions_ - 1;
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Evaluation of the schedule no. " << sched << std::endl;

  // use as initial state the state at the end of the first week of the
  // generation solution
  vector<State> initialStates = pGenerationSolver_->statesOfDay(7);
  for (int i = 0; i < pScenario_->nNurses(); i++)
    initialStates[i].dayId_ = 0;

  // create evaluation solver
  auto pSolver = std::unique_ptr<Solver>(
      setEvaluationWithInputAlgorithm(pEvaluationDemands_[0], initialStates));

  for (int j = 0; j < options_.nEvaluationDemands_; j++) {
    double timeLeft =
        options_.totalTimeLimitSeconds_ - timerTotal_.dSinceInit();
    if (timeLeft < 1.0) {
      std::cout << "# Time has run out when evaluating schedule no."
                << sched << std::endl;
      std::cout << options_.totalTimeLimitSeconds_ << "; "
                << timerTotal_.dSinceInit() << std::endl;
      return false;
    }

    (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                   << "] Starting evaluation of schedule no. " << sched
                   << " over evaluation demand no. " << j << std::endl;

    // Only perform the evaluation if the schedule is feasible and
    // there is time for more than one schedule
    double currentCost = costPreviousWeeks_ + theBaseCosts_.back();

    if (pGenerationSolver_->status(true) == INFEASIBLE) {
      currentCost += XLARGE_SCORE;
    } else {
      // Perform the actual evaluation on demand j
      // by running the chosen algorithm
      if (pSolver->status() == UNSOLVED) {
        currentCost += pSolver->solve(options_.evaluationParameters_);
      } else {
        currentCost += pSolver->resolve(pEvaluationDemands_[j],
                                        options_.evaluationParameters_);
      }
      std::cout << "Sol cost: " << currentCost << ", base cost: "
                << theBaseCosts_.back() << ", partial cost: "
                << costPreviousWeeks_ << std::endl;
    }

    // Display
    (*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Schedule no. "
                   << sched << " evaluated over evaluation demand no. " << j
                   << " (solution cost: " << currentCost << ")." << std::endl;

    // Insert the solution cost and solution
    int roundC = static_cast<int>(round(currentCost));
    schedulesFromObjectiveByEvaluationDemand_[j][roundC].insert(sched);
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
  vector<double> theNewScores(nGeneratedSolutions_, 0);
  if (options_.nEvaluationDemands_ == 0)
    theNewScores = theBaseCosts_;

  switch (strategy) {
    case RK_SCORE:
      for (int j = 0; j < options_.nEvaluationDemands_; j++) {
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Solution costs for demand no. " << j << std::endl;
        int localRank = 1;
        const auto &localCosts = schedulesFromObjectiveByEvaluationDemand_[j];
        for (const auto &p : localCosts) {
          double shared_score =
              localRank + (p.second.size() - 1.0) / p.second.size();
          localRank += p.second.size();
          for (int sched : p.second) {
            theNewScores[sched] += shared_score;
            (*pLogStream_) << "#     | sched " << sched << " -> " << p.first
                           << " (score += " << shared_score << ")" << std::endl;
          }
        }
      }
      break;
    case RK_MEAN:
      for (int j = 0; j < options_.nEvaluationDemands_; j++) {
        (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                       << "] Solution costs for demand no. " << j << std::endl;
        const auto &localCosts = schedulesFromObjectiveByEvaluationDemand_[j];
        for (const auto &p : localCosts) {
          double shared_score =
              static_cast<int>(p.first / options_.nEvaluationDemands_);
          for (int sched : p.second) {
            theNewScores[sched] += shared_score;
            (*pLogStream_) << "#     | sched " << sched << " -> " << p.first
                           << " (score += " << shared_score << ")" << std::endl;
          }
        }
      }
      break;
    case RK_NONE:Tools::throwError("Ranking strategy set to NONE.");
      break;
    default:Tools::throwError("Ranking strategy not defined.");
  }

  theScores_ = theNewScores;
  (*pLogStream_) << "# [week=" << pScenario_->thisWeek()
                 << "] Update of the scores and ranking done!" << std::endl;
}
