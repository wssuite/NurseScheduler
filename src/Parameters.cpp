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

#include <climits>
#include <string>
#include <utility>

#include "Parameters.h"

using std::string;

//-----------------------------------------------------------------------------
//  C l a s s   S o l v e r P a r a m
//  Structure that gather parameters for a solver. Can be given as an input of
//  the solve function of any solver
//-----------------------------------------------------------------------------
// Specific constructor that sets the verbosity level
SolverParam::SolverParam(int v, OptimalityLevel level,
                         std::string outfile, std::string logfile) :
    verbose_(v), outdir_(std::move(outfile)), logfile_(std::move(logfile)) {
  verbose(v);
  optimalityLevel(level);
}

#define ALLPARAMS(f) { *file >> f; for (auto sp : params) sp->f = f; }

bool SolverParam::setParameter(const string &field, std::fstream *file,
                               const std::vector<SolverParam*> &params) {
  try {
    if (Tools::strEndsWith(field, "maxSolvingTimeSeconds")) {
      ALLPARAMS(maxSolvingTimeSeconds_)
    } else if (Tools::strEndsWith(field, "printEverySolution")) {
      ALLPARAMS(printEverySolution_)
    } else if (Tools::strEndsWith(field, "absoluteGap")) {
      ALLPARAMS(absoluteGap_)
    } else if (Tools::strEndsWith(field, "minRelativeGap")) {
      ALLPARAMS(minRelativeGap_)
    } else if (Tools::strEndsWith(field, "relativeGap")) {
      ALLPARAMS(relativeGap_)
    } else if (Tools::strEndsWith(field, "nbDiveIfMinGap")) {
      ALLPARAMS(nbDiveIfMinGap_)
    } else if (Tools::strEndsWith(field, "nbDiveIfRelGap")) {
      ALLPARAMS(nbDiveIfRelGap_)
    } else if (Tools::strEndsWith(field, "maxLevelDifference")) {
      ALLPARAMS(maxLevelDifference_)
    } else if (Tools::strEndsWith(field, "maxDivingWithoutLBImprovements")) {
      ALLPARAMS(maxDivingWithoutLBImprovements_)
    } else if (Tools::strEndsWith(field, "solveToOptimality")) {
      ALLPARAMS(solveToOptimality_)
    } else if (Tools::strEndsWith(field, "solveRelaxationToOptimality")) {
      ALLPARAMS(solveRelaxationToOptimality_)
    } else if (Tools::strEndsWith(field, "maxSolvingTimeRatioForSubproblems")) {
      ALLPARAMS(maxSolvingTimeRatioForSubproblems_)
    } else if (Tools::strEndsWith(field, "stopAfterXSolution")) {
      ALLPARAMS(stopAfterXSolution_)
    } else if (Tools::strEndsWith(field, "printRelaxationSol")) {
      ALLPARAMS(printRelaxationSol_)
    } else if (Tools::strEndsWith(field, "isStabilization")) {
      ALLPARAMS(isStabilization_)
    } else if (Tools::strEndsWith(field, "isStabUpdateBoxRadius")) {
      ALLPARAMS(isStabUpdateBoxRadius_)
    } else if (Tools::strEndsWith(field, "isStabUpdatePenalty")) {
      ALLPARAMS(isStabUpdatePenalty_)
    } else if (Tools::strEndsWith(field, "maxInactiveIterations")) {
      ALLPARAMS(maxInactiveIterations_)
    } else if (Tools::strEndsWith(field, "minActivityRate")) {
      ALLPARAMS(minActivityRate_)
    } else if (Tools::strEndsWith(field, "nbDiveIfBranchOnColumns")) {
      ALLPARAMS(nbDiveIfBranchOnColumns_)
    } else if (Tools::strEndsWith(field, "branchColumnDisjoint")) {
      ALLPARAMS(branchColumnDisjoint_)
    } else if (Tools::strEndsWith(field, "branchColumnUntilValue")) {
      ALLPARAMS(branchColumnUntilValue_)
    } else if (Tools::strEndsWith(field, "branchBaseScore")) {
      ALLPARAMS(branchBaseScore_)
    } else if (Tools::strEndsWith(field, "stopAfterXDegenerateIt")) {
      ALLPARAMS(stopAfterXDegenerateIt_)
    } else if (Tools::strEndsWith(field, "nCandidates")) {
      ALLPARAMS(nCandidates_)
    } else if (Tools::strEndsWith(field, "weekendAdvantage")) {
      ALLPARAMS(weekendAdvantage_)
    } else if (Tools::strEndsWith(field, "heuristicMinIntegerPercent")) {
      ALLPARAMS(heuristicMinIntegerPercent_)
    } else if (Tools::strEndsWith(field, "performHeuristicAfterXNode")) {
      ALLPARAMS(performHeuristicAfterXNode_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicMaxIteration")) {
      ALLPARAMS(MIPHeuristicMaxIteration_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicObjLimit")) {
      ALLPARAMS(MIPHeuristicObjLimit_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicGapLimit")) {
      ALLPARAMS(MIPHeuristicGapLimit_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicNThreads")) {
      ALLPARAMS(MIPHeuristicNThreads_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicSolver")) {
      std::string stype;
      *file >> stype;
      stype = Tools::toUpperCase(stype);
      MIPHeuristicSolver_ = solverTypesByName.at(stype);
    } else if (Tools::strEndsWith(field, "MIPHeuristicUseRotations")) {
      ALLPARAMS(MIPHeuristicUseRotations_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicVerbose")) {
      ALLPARAMS(MIPHeuristicVerbose_)
    } else if (Tools::strEndsWith(field, "MIPHeuristicNRostersPerNurse")) {
      ALLPARAMS(MIPHeuristicNColumnsPerNurse_)
    } else if (Tools::strEndsWith(field, "performDiveHeuristic")) {
      ALLPARAMS(performDiveHeuristic_)
    } else if (Tools::strEndsWith(field, "performMIPHeuristic")) {
      ALLPARAMS(performMIPHeuristic_)
    } else if (Tools::strEndsWith(field, "performLNSHeuristic")) {
      ALLPARAMS(performLNSHeuristic_)
    } else if (Tools::strEndsWith(field, "weightStrategy")) {
      std::string strat;
      *file >> strat;
      strat = Tools::toUpperCase(strat);
      weightStrategy_ = weightStrategiesByName.at(strat);
    } else if (Tools::strEndsWith(field, "spType")) {
      std::string spType;
      *file >> spType;
      spType = Tools::toUpperCase(spType);
      spType_ = spTypesByName.at(spType);
    } else if (Tools::strEndsWith(field, "rcsppType")) {
      std::string rcsppType;
      *file >> rcsppType;
      rcsppType = Tools::toUpperCase(rcsppType);
      rcsppType_ = rcsppTypesByName.at(rcsppType);
    } else {
      return false;
    }

    return true;
  } catch (...) {
    std::cerr << "Exception thrown while processing field " << field
              << std::endl;
    rethrow_exception(std::current_exception());
  }
}

void SolverParam::read(const string &strFile) {
    string content = Tools::loadOptions(strFile,
              [this](const std::string& field, std::fstream *file) {
    return setParameter(field, file);
  });

  Tools::LogOutput log(logfile_, true, true);
  log << "===================================================" << std::endl;
  log.addCurrentTime() << "Solver options : " << strFile << std::endl;
  log << content;
  log << "===================================================" << std::endl;
}

// Initialize all the parameters according to a small number of options that
// represent the strategies we want to test
void SolverParam::verbose(int v) {
  // Default display values are all set to false,
  // so new values are given only if set to true
  //
  verbose_ = v;
  spParam_.verbose_ = v;
  switch (v) {
    case 0: printRelaxationSol_ = false;
      printIntermediarySol_ = false;
      printFractionOfInteger_ = false;
      printBcpSummary_ = false;
      printNodes_ = false;
      printFinalSol_ = false;
      printBranchStats_ = false;
      printRelaxationLp_ = false;
      break;
    case 1: printRelaxationSol_ = false;
      printIntermediarySol_ = false;
      printFractionOfInteger_ = false;
      printBcpSummary_ = true;
      printNodes_ = false;
      printFinalSol_ = false;
      printBranchStats_ = false;
      printRelaxationLp_ = false;
      break;
    case 2: printRelaxationSol_ = false;
      printIntermediarySol_ = true;
      printFractionOfInteger_ = true;
      printBcpSummary_ = true;
      printNodes_ = false;
      printFinalSol_ = false;
      printBranchStats_ = true;
      printRelaxationLp_ = false;
      break;
    case 3: printRelaxationSol_ = true;
      printIntermediarySol_ = true;
      printBcpSummary_ = true;
      printNodes_ = true;
      printFinalSol_ = true;
      printBranchStats_ = true;
      printRelaxationLp_ = true;
      break;
    default: break;
  }
}

SubProblemParam::SubProblemParam(const SolverParam& param):
    SubProblemParam(param.spParam_) {
  verbose_ = param.verbose_;
  epsilon_ = param.epsilon_;
}

// Set the parameters relative to the optimality level
void SolverParam::optimalityLevel(OptimalityLevel level) {
  switch (level) {
    case UNTIL_FEASIBLE: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = .1;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 2;
      solveToOptimality_ = false;
      stopAfterXSolution_ = 1;
      break;
    case TWO_DIVES: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = .1;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 2;
      solveToOptimality_ = false;
      stopAfterXSolution_ = LARGE_SCORE;
      break;
    case REPEATED_DIVES: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = .02;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 4;
      solveToOptimality_ = false;
      stopAfterXSolution_ = LARGE_SCORE;
      break;
    case OPTIMALITY: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = 1e-4;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 2;
      solveToOptimality_ = true;
      stopAfterXSolution_ = LARGE_SCORE;
      solveRelaxationToOptimality_ = true;
      break;
  }
}

bool SubProblemParam::setParameter(
    const string &field,
    std::fstream *file,
    const std::vector<SubProblemParam*> &params) {
  try {
    if (Tools::strEndsWith(field, "spNbRotationsPerNurse")) {
      ALLPARAMS(spNbRotationsPerNurse_)
    } else if (Tools::strEndsWith(field, "spNbNursesToPrice")) {
      ALLPARAMS(spNbNursesToPrice_)
    } else if (Tools::strEndsWith(field, "spMaxReducedCostBound")) {
      ALLPARAMS(spMaxReducedCostBound_)
    } else if (Tools::strEndsWith(field, "spMaxSolvingTimeSeconds")) {
      ALLPARAMS(spMaxSolvingTimeSeconds_)
    } else if (Tools::strEndsWith(
        field, "spMaxSolvingTimeRatioForRelaxation")) {
      ALLPARAMS(spMaxSolvingTimeRatioForRelaxation_)
    } else if (Tools::strEndsWith(field, "spComputeLB")) {
      ALLPARAMS(spComputeLB_)
    } else if (Tools::strEndsWith(field, "rcsppResetParamAtEachIteration")) {
      ALLPARAMS(rcsppResetParamAtEachIteration_)
    } else if (Tools::strEndsWith(
        field, "rcsspWaitBeforeStartingNextExecution")) {
      ALLPARAMS(rcsspWaitBeforeStartingNextExecution_)
    } else if (Tools::strEndsWith(field, "rcsppToOptimality")) {
      ALLPARAMS(rcsppToOptimality_)
    } else if (Tools::strEndsWith(field, "rcsppSortLabels")) {
      ALLPARAMS(rcsppSortLabels_)
    } else if (Tools::strEndsWith(field, "rcsppMinCostToSinks")) {
      ALLPARAMS(rcsppMinCostToSinks_)
    } else if (Tools::strEndsWith(field, "rcsppImprovedDomination")) {
      ALLPARAMS(rcsppImprovedDomination_)
    } else if (Tools::strEndsWith(field, "rcsppEnumSubpaths")) {
      ALLPARAMS(rcsppEnumSubpaths_)
    } else if (Tools::strEndsWith(
        field, "rcsppEnumSubpathsForMinCostToSinks")) {
      ALLPARAMS(rcsppEnumSubpathsForMinCostToSinks_)
    } else if (Tools::strEndsWith(field, "rcsppDssr")) {
      ALLPARAMS(rcsppDssr_)
    } else if (Tools::strEndsWith(field, "rcsppIncrementalDssr")) {
      ALLPARAMS(rcsppIncrementalDssr_)
    } else if (Tools::strEndsWith(field, "rcsppNbToExpand")) {
      ALLPARAMS(rcsppNbToExpand_)
    } else if (Tools::strEndsWith(field, "rcsppBidirectional")) {
      ALLPARAMS(rcsppBidirectional_)
    } else if (Tools::strEndsWith(field, "rcsppRandomStartDay")) {
      ALLPARAMS(rcsppRandomStartDay_)
    } else if (Tools::strEndsWith(field, "rcsppUseGlobalMaxLevel")) {
      ALLPARAMS(rcsppUseGlobalMaxLevel_)
    } else if (Tools::strEndsWith(field, "spDefaultStrategy")) {
      ALLPARAMS(strategyLevel_)
    } else if (Tools::strEndsWith(field, "spSortNursesBasedOnDuals")) {
      ALLPARAMS(spSortNursesBasedOnDuals_)
    } else {
      return false;
    }
    return true;
  } catch (...) {
    std::cerr << "Exception thrown while processing field " << field
              << std::endl;
    rethrow_exception(std::current_exception());
  }
}


bool DeterministicSolverOptions::setParameter(
    const string &field, std::fstream *file) {
  try {
    if (Tools::strEndsWith(field, "divideIntoConnectedPositions")) {
      *file >> divideIntoConnectedPositions_;
    } else if (Tools::strEndsWith(field, "withPrimalDual")) {
      *file >> withPrimalDual_;
    } else if (Tools::strEndsWith(field, "withRollingHorizon")) {
      *file >> withRollingHorizon_;
    } else if (Tools::strEndsWith(field, "rollingSamplePeriod")) {
      *file >> rollingSamplePeriod_;
    } else if (Tools::strEndsWith(field, "rollingControlHorizon")) {
      *file >> rollingControlHorizon_;
    } else if (Tools::strEndsWith(field, "rollingPredictionHorizon")) {
      *file >> rollingPredictionHorizon_;
    } else if (Tools::strEndsWith(field, "withLNS")) {
      *file >> withLNS_;
    } else if (Tools::strEndsWith(field, "lnsMaxItWithoutImprovement")) {
      *file >> lnsMaxItWithoutImprovement_;
    } else if (Tools::strEndsWith(field, "lnsNursesRandomDestroy")) {
      *file >> lnsNursesRandomDestroy_;
    } else if (Tools::strEndsWith(field, "lnsNursesPositionDestroy")) {
      *file >> lnsNursesPositionDestroy_;
    } else if (Tools::strEndsWith(field, "lnsNursesContractDestroy")) {
      *file >> lnsNursesContractDestroy_;
    } else if (Tools::strEndsWith(field, "lnsDestroyOverTwoWeeks")) {
      *file >> lnsDestroyOverTwoWeeks_;
    } else if (Tools::strEndsWith(field, "lnsDestroyOverFourWeeks")) {
      *file >> lnsDestroyOverFourWeeks_;
    } else if (Tools::strEndsWith(field, "lnsDestroyOverAllWeeks")) {
      *file >> lnsDestroyOverAllWeeks_;
    } else if (Tools::strEndsWith(field, "lnsNbNursesDestroyOverTwoWeeks")) {
      *file >> lnsNbNursesDestroyOverTwoWeeks_;
    } else if (Tools::strEndsWith(field, "lnsNbNursesDestroyOverFourWeeks")) {
      *file >> lnsNbNursesDestroyOverFourWeeks_;
    } else if (Tools::strEndsWith(field, "lnsNbNursesDestroyOverAllWeeks")) {
      *file >> lnsNbNursesDestroyOverAllWeeks_;
    } else if (Tools::strEndsWith(field, "verbose")) {
      *file >> verbose_;
    } else if (Tools::strEndsWith(field, "solutionAlgorithm")) {
      std::string algoName;
      *file >> algoName;
      algoName = Tools::toUpperCase(algoName);
      solutionAlgorithm_ = algorithmsByName.at(algoName);
    } else if (Tools::strEndsWith(field, "solverType")) {
      std::string solverName;
      *file >> solverName;
      solverName = Tools::toUpperCase(solverName);
      mySolverType_ = solverTypesByName.at(solverName);
#ifdef NS_DEBUG
      std::cout << "LP solver :" << solverName << std::endl;
#endif
    } else {
      return false;
    }
    return true;
  } catch (...) {
    std::cerr << "Exception thrown while processing field " << field
              << std::endl;
    rethrow_exception(std::current_exception());
  }
}

void StochasticSolverOptions::setStochasticSolverOptions(
    PScenario pScenario, string solPath, string logPathIni,
    double timeout, bool useRotation) {
  string logStochastic =
      logPathIni.empty() ? "" : logPathIni + "LogStochastic.txt";
  string logSolver = logPathIni.empty() ? "" : logPathIni + "LogSolver.txt";

  withIterativeDemandIncrease_ = false;
  withEvaluation_ = true;
  generationCostPerturbation_ = true;
  evaluationCostPerturbation_ = true;
  withResolveForGeneration_ = false;
  generationAlgorithm_ = GENCOL;
  withResolveForEvaluation_ = true;
  evaluationAlgorithm_ = GENCOL;
  totalTimeLimitSeconds_ = timeout;
  nExtraDaysGenerationDemands_ = 7;
  nEvaluationDemands_ = 2;
  nDaysEvaluation_ = 14;
  logfile_ = logStochastic;
  rankingStrategy_ = RK_MEAN;
  demandingEvaluation_ = true;
  verbose_ = 0;
#ifdef NS_DEBUG
  verbose_ = 1;
#endif

  generationParameters_ = SolverParam();
  generationParameters_.verbose(verbose_);
  generationParameters_.maxSolvingTimeSeconds_ = totalTimeLimitSeconds_ - 1;
  generationParameters_.printEverySolution_ = false;
  generationParameters_.outdir_ = solPath;
  generationParameters_.logfile_ = logSolver;
  generationParameters_.absoluteGap_ = 1;
  generationParameters_.minRelativeGap_ = 05;
  generationParameters_.relativeGap_ = .1;
  // ensure a two-dives tree exploration
  generationParameters_.nbDiveIfBranchOnColumns_ = 2;
  generationParameters_.nbDiveIfMinGap_ = 1;
  generationParameters_.nbDiveIfRelGap_ = 2;
  generationParameters_.solveToOptimality_ = false;
  generationParameters_.branchColumnUntilValue_ = .99;
  generationParameters_.weightStrategy_ = RANDOMMEANMAX;
  generationParameters_.spType_ = useRotation ? ALL_ROTATION : ROSTER;

  evaluationParameters_ = SolverParam();
  evaluationParameters_.verbose(verbose_);
  evaluationParameters_.maxSolvingTimeSeconds_ = totalTimeLimitSeconds_ - 1;
  evaluationParameters_.printEverySolution_ = false;
  evaluationParameters_.logfile_ = logSolver;
  evaluationParameters_.absoluteGap_ = 5;
  evaluationParameters_.minRelativeGap_ = .05;
  evaluationParameters_.relativeGap_ = .1;
  evaluationParameters_.nbDiveIfMinGap_ = 1;
  evaluationParameters_.nbDiveIfRelGap_ = 2;
  evaluationParameters_.solveToOptimality_ = false;
  evaluationParameters_.stopAfterXSolution_ = 0;
  evaluationParameters_.spType_ = useRotation ? ALL_ROTATION : ROSTER;
}

void StochasticSolverOptions::setStochasticSolverOptions(
    string instanceName,
    string solPath,
    string logPathIni,
    string stochasticOptionsFile,
    string generationOptionsFile,
    string evaluationOptionsFile) {
  string logStochastic =
      logPathIni.empty() ? "" : logPathIni + "LogStochastic.txt";
  string logSolver = logPathIni.empty() ? "" : logPathIni + "LogSolver.txt";

  read(stochasticOptionsFile);
  logfile_ = logStochastic;

  generationParameters_ = SolverParam();
  generationParameters_.read(generationOptionsFile);
  generationParameters_.outdir_ = solPath;
  generationParameters_.verbose_ = verbose_;

  evaluationParameters_ = SolverParam();
  evaluationParameters_.read(evaluationOptionsFile);
  evaluationParameters_.logfile_ = logSolver;
  evaluationParameters_.verbose_ = verbose_;
}

bool StochasticSolverOptions::setParameter(
    const string &field, std::fstream *file) {
  try {
    if (Tools::strEndsWith(field, "withEvaluation")) {
      *file >> withEvaluation_;
    } else if (Tools::strEndsWith(field, "withIterativeDemandIncrease")) {
      *file >> withIterativeDemandIncrease_;
    } else if (Tools::strEndsWith(field, "generationCostPerturbation")) {
      *file >> generationCostPerturbation_;
    } else if (Tools::strEndsWith(field, "evaluationCostPerturbation")) {
      *file >> evaluationCostPerturbation_;
    } else if (Tools::strEndsWith(field, "generationAlgorithm")) {
      std::string algoName;
      *file >> algoName;
      algoName = Tools::toUpperCase(algoName);
      generationAlgorithm_ = algorithmsByName.at(algoName);
    } else if (Tools::strEndsWith(field, "solverType")) {
      std::string solverName;
      *file >> solverName;
      solverName = Tools::toUpperCase(solverName);
      mySolverType_ = solverTypesByName.at(solverName);
#ifdef NS_DEBUG
      std::cout << "LP solver :" << solverName << std::endl;
#endif
    } else if (Tools::strEndsWith(field, "evaluationAlgorithm")) {
      std::string algoName;
      *file >> algoName;
      algoName = Tools::toUpperCase(algoName);
      evaluationAlgorithm_ = algorithmsByName.at(algoName);
    } else if (Tools::strEndsWith(field, "rankingStrategy")) {
      std::string stratName;
      *file >> stratName;
      stratName = Tools::toUpperCase(stratName);
      rankingStrategy_ = rankingsByName.at(stratName);
    } else if (Tools::strEndsWith(field, "nExtraDaysGenerationDemands")) {
      *file >> nExtraDaysGenerationDemands_;
    } else if (Tools::strEndsWith(field, "nEvaluationDemands")) {
      *file >> nEvaluationDemands_;
    } else if (Tools::strEndsWith(field, "nDaysEvaluation")) {
      *file >> nDaysEvaluation_;
    } else {
      return false;
    }
    return true;
  } catch (...) {
    std::cerr << "Exception thrown while processing field " << field
              << std::endl;
    rethrow_exception(std::current_exception());
  }
}

void StochasticSolverOptions::read(const string &strFile) {
  string content = Tools::loadOptions(
      strFile, [this](const std::string &field, std::fstream *file) {
        return setParameter(field, file);
      });

  Tools::LogOutput log(logfile_, true, true);
  log << "===================================================" << std::endl;
  log.addCurrentTime() << "Stochastic options : " << strFile << std::endl;
  log << content;
  log << "===================================================" << std::endl;
}
