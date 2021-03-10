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

#ifndef SRC_SOLVERS_DETERMINISTICSOLVER_H_
#define SRC_SOLVERS_DETERMINISTICSOLVER_H_

#include <string>
#include <vector>

#include "solvers/Solver.h"
#include "tools/InputPaths.h"
#include "tools/GlobalStats.h"

enum NursesSelectionOperator {
  NURSES_RANDOM,
  NURSES_POSITION,
  NURSES_CONTRACT
};
enum DaysSelectionOperator { TWO_WEEKS, FOUR_WEEKS, ALL_WEEKS };
enum RepairOperator {
  REPAIR_TWO_DIVES,
  REPAIR_REPEATED_DIVES,
  REPAIR_OPTIMALITY
};

class DeterministicSolverOptions {
 public:
  DeterministicSolverOptions() {}
  ~DeterministicSolverOptions() {}

  // True -> decompose the process to treat nurses with non connected positions
  // separately
  bool divideIntoConnectedPositions_ = true;

  // True -> solves the problem with a increasing horizon
  // False -> solves the whole horizon directly
  bool withRollingHorizon_ = false;

  // Parameters of the rolling horizon
  int rollingSamplePeriod_ = 7;
  int rollingControlHorizon_ = 14;
  int rollingPredictionHorizon_ = 56;

  // True -> find an initial solution with primal-dual procedure
  // False -> do it otherwise
  bool withPrimalDual_ = false;

  // True -> improve an initial solution with an LNS
  // False -> return directly the best feasible solution
  bool withLNS_ = false;

  // Parameters of the LNS
  int lnsMaxItWithoutImprovement_ = 10;
  bool lnsNursesRandomDestroy_ = true;
  bool lnsNursesPositionDestroy_ = true;
  bool lnsNursesContractDestroy_ = true;
  bool lnsDestroyOverTwoWeeks_ = true;
  bool lnsDestroyOverFourWeeks_ = true;
  bool lnsDestroyOverAllWeeks_ = true;
  int lnsNbNursesDestroyOverTwoWeeks_ = 30;
  int lnsNbNursesDestroyOverFourWeeks_ = 10;
  int lnsNbNursesDestroyOverAllWeeks_ = 5;

  // parameters of column generation
  bool isStabilization_ = false;
  bool isStabUpdateBoxRadius_ = false;
  bool isStabUpdatePenalty_ = false;
  int stopAfterXDegenerateIt_ = 1;

  // Algorithm that we use to solve the problem
  Algorithm solutionAlgorithm_ = GENCOL;

  // Solver used for the LPs
  SolverType MySolverType_ = S_CLP;

  // Level of optimality that is requested to the solvers
  // 0 -> return the first feasible solution
  // 1 -> only get a small number of feasible solutions
  // 2 -> try harder
  // 3 -> go to optimality
  int optimalityLevel_ = 1;

  // Time limit
  int totalTimeLimitSeconds_ = LARGE_TIME;

  // Output options
  std::string logfile_ = "";
  int verbose_ = 1;

  // Initial random seed for each solution of a new deterministic solver
  int randomSeed_ = 0;

  // number of available threads
  int nThreads_ = 1;
};


//-----------------------------------------------------------------------------
//
//  C l a s s   D e t e r m i n i s t i c S o l v e r
//
//  Solves the problem with deterministic demand
//  From a given problem (number of weeks, nurses, etc.), can compute a
//  solution.
//
//-----------------------------------------------------------------------------

class DeterministicSolver : public Solver {
 public:
  explicit DeterministicSolver(PScenario pScenario);
  DeterministicSolver(PScenario pScenario, const InputPaths &inputPaths);

  ~DeterministicSolver();

//----------------------------------------------------------------------------
//
// PARAMETERS FUNCTIONS
// Set the parameters of every solver
//
//----------------------------------------------------------------------------

 public:
  // Getter/setter
  //
  const DeterministicSolverOptions &getOptions() const {
    return options_;
  }
  const SolverParam &getRollingParameters() const {
    return rollingParameters_;
  }
  const SolverParam &getLnsParameters() const {
    return lnsParameters_;
  }
  const SolverParam &getCompleteParameters() const {
    return completeParameters_;
  }

  void copyParameters(DeterministicSolver *pSolver) {
    options_ = pSolver->getOptions();
    rollingParameters_ = pSolver->getRollingParameters();
    lnsParameters_ = pSolver->getLnsParameters();
    completeParameters_ = pSolver->getCompleteParameters();
    generatePResourcesFunc_ = pSolver->generatePResourcesFunc_;
  }

  // Initialize deterministic options with default values
  //
  void initializeOptions(const InputPaths &inputPaths);

  // Read deterministic options from a file
  //
  void readOptionsFromFile(const InputPaths &inputPaths);

  // Set total cpu time available to the solution
  //
  void setTotalTimeLimit(double t) {
    options_.totalTimeLimitSeconds_ = t;
    rollingParameters_.maxSolvingTimeSeconds_ = t;
    completeParameters_.maxSolvingTimeSeconds_ = t;
    lnsParameters_.maxSolvingTimeSeconds_ = t;
  }

 protected:
  // Options that characterize the execution of the stochastic solver
  DeterministicSolverOptions options_;

  //----------------------------------------------------------------------------
  //
  // STATISTICS OF THE OVERALL SOLUTION PROCESS
  //
  //----------------------------------------------------------------------------

 public:
  // getter
  //
  GlobalStats getGlobalStat() { return stats_; }

  // update functions for the most relevant statistics
  //
  void updateInitialStats(Solver *pSolver);
  void updateImproveStats(Solver *pSolver);

 protected:
  GlobalStats stats_;

  //----------------------------------------------------------------------------
  //
  // SOLVE FUNCTIONS
  // The one is general for the whole process
  //
  //----------------------------------------------------------------------------

 public:
  // Main function
  double solve(const std::vector<Roster> &solution = {});

  // Solve the problem using a decomposition of the set nurses by connected
  // components of the rcspp of positions
  double solveByConnectedPositions();

 protected:
  // Ready the solver for the solution process
  void init();

  // After the end of a solution process: retrieve status, solution, etc.
  double treatResults(Solver *pSolver);

  //----------------------------------------------------------------------------
  //
  // SOLVE THE PROBLEM DIRECTLY
  //
  //----------------------------------------------------------------------------

 public:
  // Solve a deterministic input demand with the input algorithm
  //
  double solveCompleteHorizon(const std::vector<Roster> &solution = {});

 private:
  // Solver that will be called to solve each sampling period in the rolling
  // horizon
  //
  Solver *pCompleteSolver_;

  // Parameters of the complete solution
  //
  SolverParam completeParameters_;

  //----------------------------------------------------------------------------
  //
  // SOLUTION WITH ROLLING HORIZON
  //
  //----------------------------------------------------------------------------

 public:
  // Solve the problem with a increasing horizon algorithm
  // The sample period is the number of days for which column variables must
  // be integer.
  //
  double solveWithRollingHorizon(const std::vector<Roster> &solution = {});

 private:
  // Solver that will be called to solve each sampling period in the rolling
  // horizon.
  Solver *pRollingSolver_;

  // Parameters of the rolling horizon solver.
  SolverParam rollingParameters_;

  // Set the optimality level of the rolling horizon solver
  // This function needs to be called before each new solution, and the behavior
  // depends on the first day of the horizon
  //
//  void rollingSetOptimalityLevel(int firstDay);

  //----------------------------------------------------------------------------
  //
  // SOLUTION WITH LARGE NEIGHBORHOOD SEARCH
  //
  //----------------------------------------------------------------------------

 public:
  // Iteratively fix the schedule of complete weeks or of a set of nurses and
  // solve the rest to improve a initial feasible solution
  //
  double solveWithLNS(const std::vector<Roster> &solution = {});

  // Print to a string the statistics of the lns
  //
  std::string lnsStatsToString();

  // Return a solver with the algorithm specified in the options_
  //
  Solver *setSolverWithInputAlgorithm(PDemand pDemand,
                                      Algorithm algorithm,
                                      const SolverParam &param);

 protected:
  // Prepare data structures for LNS
  //
  void initializeLNS();

  // Application of the destroy operator
  //
  void adaptiveDestroy(NursesSelectionOperator nurseOp,
                       DaysSelectionOperator dayOp);

  // Initialize the organized vectors of live nurses
  //
  void organizeTheLiveNursesByPosition();
  void organizeTheLiveNursesByContract();

  // Solver that will be called to solve each subproblem
  //
  Solver *pLNSSolver_;

  // Parameters of the LNS
  //
  SolverParam lnsParameters_;

  // Vector of destroy/repair operators used in the LNS
  //
  std::vector<NursesSelectionOperator> nursesSelectionOperators_;
  std::vector<DaysSelectionOperator> daysSelectionOperators_;
  std::vector<RepairOperator> repairOperators_;

  // Organized live nurses according to their contracts or positions
  //
  int nbPositions_ = 0;
  int nbContracts_ = 0;
  std::vector<std::vector<PLiveNurse> > theLiveNursesByPosition_;
  std::vector<std::vector<PLiveNurse> > theLiveNursesByContract_;

  // Weights of the position or contract when drawing the destroy operator
  //
  std::vector<double> positionWeights_;
  std::vector<double> contractWeights_;
};

#endif  // SRC_SOLVERS_DETERMINISTICSOLVER_H_
