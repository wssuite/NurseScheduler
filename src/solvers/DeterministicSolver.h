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

#include <memory>
#include <string>
#include <vector>

#include "solvers/Solver.h"
#include "tools/InputPaths.h"
#include "tools/GlobalStats.h"


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
  explicit DeterministicSolver(
      const PScenario& pScenario, const InputPaths &inputPaths = InputPaths());
  DeterministicSolver(
      const PScenario& pScenario, const DeterministicSolver &solver);

  ~DeterministicSolver();

  // stop any solver which has been paused
  void stop(bool wait = true);

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
  // return true if ignore rolling horizon, LNS etc.
  bool useCompleteSolverAnyway() const {
    return (pScenario_->nDays() <= 28 && pScenario_->nNurses() <= 8)
        || (pScenario_->nDays() <= 56 && pScenario_->nNurses() <= 5);
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
  double solve(const std::vector<Roster> &solution = {}) override;

  // Solve the problem using a decomposition of the set nurses by connected
  // components of the rcspp of positions
  double solveByConnectedPositions() override;
  double solveOneComponent(const vector<Roster> &solution);

 protected:
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
  std::unique_ptr<Solver> pCompleteSolver_;

  // Parameters of the complete solution
  SolverParam completeParameters_;

  // pool to be able to pause the solving process of the complete solver
  Tools::PThreadsPool pThreadsPool_ = Tools::ThreadsPool::newThreadsPool(1);


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
  std::unique_ptr<Solver> pRollingSolver_;

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
  Solver *newSolverWithInputAlgorithm(Algorithm algorithm,
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
