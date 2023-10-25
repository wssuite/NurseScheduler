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

#ifndef SRC_SOLVERS_STOCHASTICSOLVER_H_
#define SRC_SOLVERS_STOCHASTICSOLVER_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "solvers/Solver.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.),
//  can compute a solution.
//
//-----------------------------------------------------------------------------

class StochasticSolver : public Solver {
 public:
  StochasticSolver(PScenario pScenario,
                   StochasticSolverOptions options,
                   vector2D<PDemand> demandHistory,
                   double costPreviousWeeks = 0);

  ~StochasticSolver();

  //----------------------------------------------------------------------------
  //
  // SOLVE FUNCTIONS
  // The one is general for the whole process
  //
  //----------------------------------------------------------------------------

  // Main function
  double solve(const std::vector<Roster> &initialSolution = {}) override;

  // get the number of generated schedules
  //
  int getNGeneratedSolutions() { return nGeneratedSolutions_; }

 protected:
  // Options that characterize the execution of the stochastic solver
  StochasticSolverOptions options_;

  // Log file that can be useful when calling the solver through simulator
  std::unique_ptr<Tools::LogOutput> pLogStream_;

  //----------------------------------------------------------------------------
  //
  // SUBSOLVE FUNCTIONS
  //
  //----------------------------------------------------------------------------
  // Solves the problem by generation + evaluation of scenarios
  void solveOneWeekGenerationEvaluation();
  // Does everything for one schedule (for one week): Includes generation,
  // evaluation of the score, and update of the rankings and data.
  // Returns false if time has run out
  bool addAndSolveNewSchedule();
  // Iterative solution process in which the week is first solved by itself,
  // before adding one perturbebd week demand and solving the new extended
  // demand demand until no time is left
  void solveIterativelyWithIncreasingDemand();
  // Solves the problem by generating a schedule + using cost penalties
  void solveOneWeekNoGenerationEvaluation();
  // Special case of the last week
  void solveOneWeekWithoutPenalties();


  //----------------------------------------------------------------------------
  //
  // GENERATION OF DEMANDS FOR THE CURRENT WEEK (=FOR SCHEDULE GENERATION)
  // Note that these demands share a common first week which is
  // the week we currently try to solve.
  //
  //----------------------------------------------------------------------------

  // History
  vector2D<PDemand> pDemandsHistory_;
  // Vector of random demands that are used to GENERATE the schedules
  vector2D<PDemand> pGenerationDemands_;
  // Generate a new demand for generation
  vector<PDemand> generateSingleGenerationDemands();



  //----------------------------------------------------------------------------
  //
  // GENERATION OF SCENARIOS FOR THE FUTURE (=FOR SCHEDULE EVALUATION)
  //
  //----------------------------------------------------------------------------

  // Vector of random demands that are used to EVAULATE the generated schedules
  vector2D<PDemand> pEvaluationDemands_;
  // Generate the schedules that are used for evaluation
  void generateAllEvaluationDemands();



  //----------------------------------------------------------------------------
  //
  // GENERATION OF SCHEDULES
  // A solution is a potential candidate to be the chosen schedule
  // for the week we are solving.
  // A result is, given a solution and a potential future, the value obtained
  // for that couple (solution,demand) [i.e. LP bound for instance]
  // A score is, given a solution, the average score it obtains,
  // compared to the other solutions
  // (the precise meaning of "score" should be better defined)
  //
  //----------------------------------------------------------------------------

  // Generation Schedule
  std::unique_ptr<Solver> pGenerationSolver_;
  int nGeneratedSolutions_;

  // Return a solver with the algorithm specified for schedule GENERATION
  Solver *setGenerationSolverWithInputAlgorithm(vector<PDemand> pDemands);
  // Generate a new schedule
  bool generateNewSchedule();

  //----------------------------------------------------------------------------
  //
  // EVALUATION OF SCHEDULES
  //
  //----------------------------------------------------------------------------

  // Empty preferences -> only 1 to avoid multiplying them
  Preferences *pEmptyPreferencesForEvaluation_;
  // Evaluation
  std::vector<std::map<int, std::set<int> > >
      schedulesFromObjectiveByEvaluationDemand_;
  // Scores
  std::vector<double> theScores_;

  int bestSchedule_;
  double bestScore_;
  double costPreviousWeeks_;
  std::vector<double> theBaseCosts_;

  // Return a solver with the algorithm specified for schedule EVALUATION
  Solver *setEvaluationWithInputAlgorithm(
          vector<PDemand> pDemands, const vector<State> &stateEndOfSchedule);

  // Evaluate the last schedule and store the corresponding detailed results
  // (returns false if time has run out)
  bool evaluateLastGeneratedSchedule();

  // Recompute all scores after one schedule evaluation
  void updateRankingsAndScores(RankingStrategy strategy);

 protected:
  // Update the weights
  // For now, the update depends only on the initial states and on the contract
  // of the nurses, on the number of days on the demand, on the number of weeks
  // already treated and on the number of weeks left
  //
  void computeWeightsTotalShifts();
};

#endif  // SRC_SOLVERS_STOCHASTICSOLVER_H_
