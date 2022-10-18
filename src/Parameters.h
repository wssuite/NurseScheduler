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

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <climits>
#include <map>
#include <string>
#include <vector>

#include "data/Scenario.h"
#include "tools/Tools.h"

using std::string;

//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r P a r a m
//
//  Structure that gather parameters for a solver. Can be given as an input of
//  the solve function of any solver
//
//-----------------------------------------------------------------------------
enum OptimalityLevel { UNTIL_FEASIBLE, TWO_DIVES, REPEATED_DIVES, OPTIMALITY };

static const std::map<std::string, OptimalityLevel> optimalityLevelByName = {
    {"UNTIL_FEASIBLE", UNTIL_FEASIBLE},
    {"TWO_DIVES", TWO_DIVES},
    {"REPEATED_DIVES", REPEATED_DIVES},
    {"OPTIMALITY", OPTIMALITY}
};

static OptimalityLevel strToOptimalityLevel(const string& strOpt) {
  std::string str = Tools::toUpperCase(strOpt);
  try {
    return optimalityLevelByName.at(strOpt);
  } catch (const std::out_of_range &e) {
    std::cerr << e.what() << std::endl;
    Tools::throwError("Unknown optimality level");
    return UNTIL_FEASIBLE;
  }
}

// Algorithms for the overall solution
//
enum Algorithm { GENCOL, STOCHASTIC_GENCOL, NONE };
static const std::map<std::string, Algorithm> algorithmsByName = {
    {"GENCOL", GENCOL},
    {"STOCHASTIC_GENCOL", STOCHASTIC_GENCOL},
    {"NONE", NONE}
};

enum SolverType { Gurobi, Cplex, CBC, SCIP, CLP, FirstAvailable };
static std::map<std::string, SolverType> solverTypesByName = {
    {"CLP", CLP},
    {"GUROBI", Gurobi},
    {"CPLEX", Cplex},
    {"CBC", CBC},
    {"SCIP", SCIP},
    {"FIRST_AVAIL", FirstAvailable}
};
static std::map<SolverType, std::string> namesBySolverType =
    Tools::buildNamesByType(solverTypesByName);

static bool isGurobiAvailable() {
#ifdef USE_GUROBI
  return true;
#endif
  return false;
}

static bool isCplexAvailable() {
#ifdef USE_CPLEX
  return true;
#endif
  return false;
}

static bool isCBCAvailable() {
#ifdef USE_CBC
  return true;
#endif
  return false;
}

static bool isSolverTypeAvailable(SolverType type) {
  switch (type) {
    case CLP: return true;
    case Gurobi: return isGurobiAvailable();
    case Cplex: return isCplexAvailable();
    case CBC: return isCBCAvailable();
    default: return false;
  }
}

static SolverType getFirstSolverTypeAvailable() {
  for (int t=0; t < CLP; ++t)
    if (isSolverTypeAvailable((SolverType) t))
      return (SolverType) t;
  return CLP;
}

// Solution statuses
//
enum Status { UNSOLVED, FEASIBLE, INFEASIBLE, OPTIMAL, TIME_LIMIT };
static const std::map<Status, std::string> statusToString = {
    {UNSOLVED, "UNSOLVED"},
    {FEASIBLE, "FEASIBLE"},
    {INFEASIBLE, "INFEASIBLE"},
    {OPTIMAL, "OPTIMAL"},
    {TIME_LIMIT, "TIME_LIMIT"}
};

// Subproblem (and also master) type
enum SPType { LONG_ROTATION = 0, ALL_ROTATION = 1, ROSTER = 2 };
static const std::map<std::string, SPType> spTypesByName = {
    {"LONG", LONG_ROTATION},
    {"ALL", ALL_ROTATION},
    {"ROSTER", ROSTER}
};

// RCSPP solver type
enum RCSPPType { BOOST_LABEL_SETTING = 0, LABEL_SETTING = 1 };

static const std::map<std::string, RCSPPType> rcsppTypesByName = {
    {"BOOST", BOOST_LABEL_SETTING},
    {"DEFAULT", LABEL_SETTING}
};

// Strategy for the dynamic scheduler
enum WeightStrategy { MAX, MEAN, RANDOMMEANMAX, NO_STRAT };

static const std::map<std::string, WeightStrategy> weightStrategiesByName =
    {{"MAX", MAX}, {"MEAN", MEAN},
     {"RANDOMMEANMAX", RANDOMMEANMAX}, {"NO_STRAT", NO_STRAT}};

// enum to define how to search the rc graph when solving
enum SPSearchStrategy {
  SP_BREADTH_FIRST, SP_DEPTH_FIRST, SP_BEST_FIRST, SP_DOMINANT_FIRST
};

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

// Parameters of the subproblem (used in the solve function of the subproblem)
class SolverParam;
const double EPSILON = 1e-3;
struct SubProblemParam {
  SubProblemParam() = default;
  explicit SubProblemParam(const SolverParam& param);

  ~SubProblemParam() = default;

  bool setParameter(const string &field, std::fstream *file,
                    const std::vector<SubProblemParam*> &params = {});

  // *** PARAMETERS ***
  int verbose_ = 0;
  double epsilon_ = EPSILON;

  int nbMaxColumnsToAdd_ = 10;
  int spNbRotationsPerNurse_ = 20;
  int spNbNursesToPrice_ = 15;
  double spMaxReducedCostBound_ = 0.0;
  // when the subproblems is looking for all the non dominated solutions,
  // this ratio is used to set the limit on the number of non dominated
  // solutions to find.
  // It's a ratio expressed as a function on the max number of columns to add
  // to the master (nbMaxColumnsToAdd_)
  double spColumnsRatioForNumberPaths_ = 2;
  // when the number of columns generated by a strategy is lower than the ratio,
  // the pricer will use the next strategy
  double spMinColumnsRatioForIncrease_ = .5;

  // Parameters nb_max_paths, shortRotationsStrategy_, maxRotationLength_ and
  // violateConsecutiveConstraints_ are set according to the value of the
  // strategy argument in initSubProblemParam to adopt a more or less
  // agressive attitude
  SPSearchStrategy search_strategy_ = SP_BREADTH_FIRST;

  // true -> authorize the violation of the consecutive constraints
  // false -> forbid any arc that authorized to violate the consecutive
  // constraints
  bool violateConsecutiveConstraints_ = true;

  int strategyLevel_ = 0;

  // 0 -> no short rotations
  // 1 -> day-0 short rotations only
  // 2 -> day-0 and last-day short rotations only
  // 3 -> all short rotations
  int shortRotationsStrategy_ = 3;

  // TODO(JO): we never questioned this option, I think that we could remove the
  //  option and the alternative code that corresponds to the false value
  // true  -> one sink node per day
  // false -> one single sink node for the network
  bool oneSinkNodePerLastDay_ = true;

  // Parameters below are set in the parameter file once and for all

  // reset the parameters of the subproblem at each iteration of column
  // generation. If false, at each new node.
  bool rcsppResetParamAtEachIteration_ = false;

  // true if the RCSPP is solved to optimality. If false, can be satisfied
  // with any negative cost solution
  bool rcsppToOptimality_ = true;

  // if positive, the RC SPP stops as soon as it has generated enough
  // non-dominated path
  int rcsppMaxNegativeLabels_ = -1;

  // if positive, the RC SPP changes its strategy if it hasn't produced
  // this number of labels
  int rcsppMinNegativeLabels_ = -1;

  // true if labels are sorted by increasing costs before domination is checked
  bool rcsppSortLabels_ = true;

  // true if the shortest path from each node to each sink is computed to
  // delete the labels that cannot be extended into a negative cost path
  bool rcsppMinCostToSinks_ = false;

  // true if the domination rule is improved by considering the worst costs that
  // can result from the violation of soft constraints
  bool rcsppImprovedDomination_ = true;

  // true if we enumerate subpaths in the RCSPP graph in order to reduce the
  // number of resources
  bool rcsppEnumSubpaths_ = false;

  // tue  if we enumerate subpaths in a copy of the RCSPP graph, but only use
  // this enumeration to compute a better bound with minimum costs to sinks
  // minimum cost to sinks need to be computed for this option to have an impact
  bool rcsppEnumSubpathsForMinCostToSinks_ = false;

  // true if we apply a decremental state-space relaxation method, where we
  // first ignore the lower bounds of every constraint on any work shift and
  // then consider them only if needed, i.e., if no sufficiently negative
  // reduced cost column was generated in the subproblem
  bool rcsppDssr_ = false;

  // if different from 0, a heuristic search is done in the RCSPP where only
  // the specified number of labels is expanded from each node
  // we take only the labels with smallest (most negative) cost
  int rcsppNbToExpand_ = 0;

  // true if we solve the RCSPP with a bidirectional label-setting algorithm
  bool rcsppBidirectional_ = false;

  // Getters for the class fields
  int shortRotationsStrategy() const { return shortRotationsStrategy_; }
  bool oneSinkNodePerLastDay() const { return oneSinkNodePerLastDay_; }

  // Setters
  void shortRotationsStrategy(int value) { shortRotationsStrategy_ = value; }
  void oneSinkNodePerLastDay(bool value) { oneSinkNodePerLastDay_ = value; }
};


//  Class Print Solution
//  Allow to display a solution
struct PrintSolution {
  PrintSolution() = default;
  virtual ~PrintSolution() = default;
  virtual void saveSolution() = 0;
  virtual std::string currentSolToString() const = 0;
};

// All the parameters for a solver
class SolverParam {
 public:
  SolverParam() = default;
  explicit SolverParam(int verbose,
                       OptimalityLevel level = TWO_DIVES,
                       std::string outfile = "outfiles/SolverResult.txt",
                       std::string logfile = "outfiles/SolverLog.txt");

  bool setParameter(const string &field, std::fstream *file,
                    const std::vector<SolverParam*> &params = {});

  void read(const string &strFile);

  // maximal solving time in s
  int maxSolvingTimeSeconds_ = LARGE_TIME;

  // tolerance
  double epsilon_ = EPSILON;  // precision for the solver

  // print parameters
  int verbose_ = 0;
  bool printEverySolution_ = true;
  std::string outdir_ = "outfiles/";
  std::string logfile_;
  std::vector<int> weekIndices_;
  PrintSolution *saveFunction_ = nullptr;

  bool printRelaxationSol_ = false;
  bool printIntermediarySol_ = false;
  bool printFractionOfInteger_ = false;
  bool printBcpSummary_ = false;
  bool printNodes_ = false;
  bool printFinalSol_ = false;
  bool printBranchStats_ = false;
  bool printRelaxationLp_ = false;

  /* GENERAL SOLUTION PARAMETERS */
  // the solver contains only nurses that belong to the same connected component
  // of positions
  bool isConnectedComponentOfPositions_ = true;

  /* PARAMETERS OF THE BRANCH AND BOUND */
  // relative and absolute gap (with the current LB)
  // if sol below optimalAbsoluteGap_, we stop immediately
  // (the difference between two solution costs is at least 5) -> optimal
  // if sol below absoluteGap_, we stop immediately;
  // if sol below minRelativeGap_, we stop immediately;
  // if sol over relativeGap_, we stop after nbDiveIfRelGap_*dive without new
  // incumbent.
  // the two last conditions are respected just for certain strategy.
  // If we look for optimality, we also use the last feature
  // maxRelativeLPGapToKeepChild_:
  // if the gap between the tree best lb and the node lb is higher than
  // this gap -> backtrack
  double optimalAbsoluteGap_ = 5;
  double absoluteGap_ = 5;
  double minRelativeGap_ = .05;
  double relativeGap_ = .1;
  double maxRelativeLPGapToKeepChild_ = .05;

  int nbDiveIfMinGap_ = 1;
  int nbDiveIfRelGap_ = 2;
  // maximum number of levels between the current node and
  // the most top level active node. If bigger than the max, stop the dive
  // When equals to -1, use dive length instead.
  int maxLevelDifference_ = -1;
  // maximum number of times that the BcpModeler decides to dive on a child
  // node without improvement of the node LB
  int maxDivingWithoutLBImprovements_ = 10e6;
  // number of consecutive dive if branching on columns
  // if 0, the solver won't branch on columns
  int nbDiveIfBranchOnColumns_ = 1;
  // branch on columns with a disjoint argument
  int branchColumnDisjoint_ = true;
  // Set the limit on the total rounded value when branching on columns
  // Default: let the master problem choose it
  double branchColumnUntilValue_ = -1;
  // compute a base cost for each nurse based on the cost of her columns as
  // their number. This base score will be used when choosing the potential
  // branching candidates
  double branchBaseScore_ = 1;

  bool solveToOptimality_ = false;

  // stop the algorithm after finding X solutions
  // if 0, the algorithm computes the relaxation if the algorithm is a column
  // generation procedure
  int stopAfterXSolution_ = 5;

  // number of candidates to create for strong branching
  int nCandidates_ = 1;
  // weekendAdvantage_ gives an advantage to the weekend when choosing
  // the branching variables
  double weekendAdvantage_ = .1;

  // generate variables just after solving LP.
  // The main disadvantage is that the node could have been pruned before
  // if a better solution has been found
  // However, this could simplify the process (less steps between the storage
  // of the LP solution and it's use)
  bool generateColumnsASAP_ = false;

  // perform the heuristic searching an integer solution after visiting X nodes
  // in the tree.
  // if -1, do not perform the heuristic
  int performHeuristicAfterXNode_ = -1;
  double heuristicMinIntegerPercent_ = 50;
  bool performDiveHeuristic_ = false,
        performMIPHeuristic_ = false,
        performLNSHeuristic_ = false;
  // maximum number of iteration that the MIP solver in the heuristic
  // can perform to find a better solution than the current one
  int MIPHeuristicMaxIteration_ = 10e6;
  // stop as soon as the MIP finds a solution better than the limit
  bool MIPHeuristicObjLimit_ = true;
  // if MIPHeuristicObjLimit_, the limit cannot be greater than (1 + gap) * LB
  // If greater, (1 + gap) * LB is used
  double MIPHeuristicGapLimit_ = .5;
  int MIPHeuristicNThreads_ = 1;
  // solver to use for the MIP
  SolverType MIPHeuristicSolver_ = FirstAvailable;
  // if true, transform roster model to a rotation model
  bool MIPHeuristicUseRotations_ = false;
  // number of rosters to keep separate cutting them into rotations
  int MIPHeuristicNColumnsPerNurse_ = 10;
  int MIPHeuristicVerbose_ = 0;

  // parameters of the stabilization : initial costs and bounds of the
  // stabilization variables
  bool isStabilization_ = false;
  // update duals' box radius / primals' cost
  bool isStabUpdateBoxRadius_ = true;
  double stabBoxRadiusIni_ = 1E-4;
  double stabBoxRadiusMax_ = 100;
  double stabBoxBoundMax_ = LARGE_SCORE;
  double stabBoxRadiusFactor_ = 2;
  // update duals' penalty / primals' bounds
  bool isStabUpdatePenalty_ = true;
  double stabPenaltyIni_ = 1;
  double stabPenaltyMax_ = LARGE_SCORE;
  double stabPenaltyFactor_ = 2;

  // other technique against degeneracy: simply stop the solution process when
  // it starts begin degenerate
  int stopAfterXDegenerateIt_ = 10e6;

  // fathom a node is upper bound is smaller than the lagrangian bound
  bool isLagrangianFathom_ = true;

  // Parameters to remove column from the master (only one needs to be violated)
  // max number of consecutive iterations that a column can stay outside of
  // the basis
  int maxInactiveIterations_ = 50;
  // min activity rate (frequency of presence in the basis since creation)
  // that a column must have to remain in the master
  double minActivityRate_ = .1;

  /* PARAMETERS OF THE PRICER */

  // only add disjoint columns if active
  // by default false, but true for roster generation
  bool isColumnDisjoint_ = false;
  // minimum reduced cost that the best solution must have to be considered
  // when disjoint column generation is active
  double minReducedCostDisjoint_ = -1e2;
  int nSubProblemsToSolve_ = 20;

  // primal-dual strategy
  WeightStrategy weightStrategy_ = NO_STRAT;

  SubProblemParam spParam_;
  // type of the subproblem used
  SPType sp_type_ = ROSTER;
  // type of the RCSPP solver used in the subproblem
  RCSPPType rcspp_type_ = LABEL_SETTING;  // BOOST_LABEL_SETTING

  // Initialize all the parameters according to a small number of
  // verbose options
  void verbose(int verbose);

  // Set the parameters relative to the optimality level
  void optimalityLevel(OptimalityLevel level);
};

class DeterministicSolverOptions {
 public:
  DeterministicSolverOptions() = default;
  ~DeterministicSolverOptions() = default;

  bool setParameter(const string &field, std::fstream *file);

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
  SolverType mySolverType_ = CLP;

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

  // if true, pause the thread where the resolution process happens when timeout
  bool pauseSolveOnTimeout_ = false;
};

enum RankingStrategy { RK_MEAN, RK_SCORE, RK_NONE };

static const std::map<std::string, RankingStrategy> rankingsByName =
    {{"MEAN", RK_MEAN}, {"SCORE", RK_SCORE}, {"NONE", RK_NONE}};

class StochasticSolverOptions {
 public:
  StochasticSolverOptions() {
    SolverParam gp;
    generationParameters_ = gp;
    generationParameters_.weightStrategy_ = RANDOMMEANMAX;

    SolverParam ep;
    evaluationParameters_ = ep;
    evaluationParameters_.stopAfterXSolution_ = 0;
//    evaluationParameters_.weightStrategy_ = BOUNDRATIO;
  }

  ~StochasticSolverOptions() = default;

  // Set the options of the stochastic solver
// The solution time depends on the number of nurses
  void setStochasticSolverOptions(
      PScenario pScenario,
      std::string solPath,
      std::string logPathIni,
      double timeout = 10000,
      bool useRotation = true);

  void setStochasticSolverOptions(
      std::string instanceName,
      std::string solPath,
      std::string logPathIni,
      std::string stochasticOptionsFile,
      std::string generationOptionsFile,
      std::string evaluationOptionsFile);

  bool setParameter(const string &field, std::fstream *file);

  void read(const string &strFile);

  // True -> generate several schedules and chose the "best" one
  // (according to ranking strategy)
  // False -> generate only one schedule
  bool withEvaluation_ = true;

  // True -> generate schedules from random demands of increasing size
  // (1 week more each time). Keep the last one.
  // WARNING: Does not work if withEvaluation_=true
  bool withIterativeDemandIncrease_ = false;

  // True -> Perturb the costs when generating the schedules
  // The type of perturbation is set in generationParameters_ (weightStrategy_)
  bool generationCostPerturbation_ = true;

  // True -> When generating a second, third, etc. schedule,
  // warm-start with previously generated columns
  // WARNING: should remain false
  // (if true, no diversity in the generated schedules)
  bool withResolveForGeneration_ = false;

  // True -> use the real demand for the nExtraDaysGenerationDemands_
  // put the real demand in demandHistory_[0]
  bool withRealDemand_ = false;

  Algorithm generationAlgorithm_ = GENCOL;

  // cf. generation
  // withResolve is useful here, particularly when evaluating
  // with LP lowest bound
  bool evaluationCostPerturbation_ = true;
  bool withResolveForEvaluation_ = true;
  Algorithm evaluationAlgorithm_ = GENCOL;

  // Choice of ranking strategy:
  // RK_SCORE: same ranking as for the competition
  // RK_MEAN: keep the schedule with minimum expected cost over
  // the generated evaluation demands
  RankingStrategy rankingStrategy_ = RK_SCORE;
  bool demandingEvaluation_ = true;
  int totalTimeLimitSeconds_ = LARGE_TIME;

  // Number of evaluation demands generated
  // WARNING: if =0 and withEvaluation_=true, ranks the schedules according to
  // their baseCost (i.e. the "real" cost of the 1-week schedule
  // [without min/max costs])
  int nEvaluationDemands_ = 2;
  int nExtraDaysGenerationDemands_ = 7;
  int nDaysEvaluation_ = 14;
  int nGenerationDemandsMax_ = 100;

  std::string logfile_ = "";

  SolverParam generationParameters_;
  SolverParam evaluationParameters_;

#ifdef DBG
  int verbose_ = 1;
#else
  int verbose_ = 0;
#endif
};

#endif  // SRC_PARAMETERS_H_
