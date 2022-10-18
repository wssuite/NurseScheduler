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
  } else {
    return false;
  }
  return true;
}

void SolverParam::read(const string &strFile) {
    Tools::loadOptions(strFile,
              [this](const std::string& field, std::fstream *file) {
    return setParameter(field, file);
  });
}

// Initialize all the parameters according to a small number of options that
// represent the strategies we want to test
void SolverParam::verbose(int v) {
  // Default display values are all set to false,
  // so new values are given only if set to true
  //
  verbose_ = v;
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
      stopAfterXSolution_ = 999999;
      break;
    case REPEATED_DIVES: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = .02;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 4;
      solveToOptimality_ = false;
      stopAfterXSolution_ = 9999999;
      break;
    case OPTIMALITY: absoluteGap_ = optimalAbsoluteGap_;
      minRelativeGap_ = 1e-4;
      relativeGap_ = 1e-4;
      nbDiveIfMinGap_ = 1;
      nbDiveIfRelGap_ = 2;
      solveToOptimality_ = true;
      stopAfterXSolution_ = 9999999;
      break;
  }
}

bool SubProblemParam::setParameter(
    const string &field,
    std::fstream *file,
    const std::vector<SubProblemParam*> &params) {
  if (Tools::strEndsWith(field, "spNbRotationsPerNurse")) {
    ALLPARAMS(spNbRotationsPerNurse_)
  } else if (Tools::strEndsWith(field, "spNbNursesToPrice")) {
    ALLPARAMS(spNbNursesToPrice_)
  } else if (Tools::strEndsWith(field, "spMaxReducedCostBound")) {
    ALLPARAMS(spMaxReducedCostBound_)
  } else if (Tools::strEndsWith(field, "rcsppResetParamAtEachIteration")) {
    ALLPARAMS(rcsppResetParamAtEachIteration_)
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
  } else if (Tools::strEndsWith(field, "rcsppEnumSubpathsForMinCostToSinks")) {
    ALLPARAMS(rcsppEnumSubpathsForMinCostToSinks_)
  } else if (Tools::strEndsWith(field, "rcsppDssr")) {
    ALLPARAMS(rcsppDssr_)
  } else if (Tools::strEndsWith(field, "rcsppNbToExpand")) {
    ALLPARAMS(rcsppNbToExpand_)
  } else if (Tools::strEndsWith(field, "rcsppBidirectional")) {
    ALLPARAMS(rcsppBidirectional_)
  } else if (Tools::strEndsWith(field, "spDefaultStrategy")) {
    ALLPARAMS(strategyLevel_)
  } else {
    return false;
  }
  return true;
}


bool DeterministicSolverOptions::setParameter(
    const string &field, std::fstream *file) {
  if (Tools::strEndsWith(field, "divideIntoConnectedPositions")) {
    *file >> divideIntoConnectedPositions_;
  }  else if (Tools::strEndsWith(field, "withPrimalDual")) {
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
  }     else if (Tools::strEndsWith(field, "solutionAlgorithm")) {
    std::string algoName;
    *file >> algoName;
    algoName = Tools::toUpperCase(algoName);
    solutionAlgorithm_ = algorithmsByName.at(algoName);
  } else if (Tools::strEndsWith(field, "solverType")) {
    std::string solverName;
    *file >> solverName;
    solverName = Tools::toUpperCase(solverName);
    mySolverType_ = solverTypesByName.at(solverName);
#ifdef DBG
    std::cout << "LP solver :" << solverName << std::endl;
#endif
  } else {
    return false;
  }
  return true;
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
  rankingStrategy_ = RK_SCORE;
  nExtraDaysGenerationDemands_ = 7;
  nEvaluationDemands_ = 2;
  nDaysEvaluation_ = 14;
  nGenerationDemandsMax_ = 100;
  logfile_ = logStochastic;
  rankingStrategy_ = RK_SCORE;
  demandingEvaluation_ = false;
  verbose_ = 1;
#ifdef DBG
  verbose_ = 1;
#endif

  SolverParam generationParameters;
  generationParameters.verbose(verbose_);
  generationParameters.maxSolvingTimeSeconds_ = totalTimeLimitSeconds_ - 1;
  generationParameters.printEverySolution_ = false;
  generationParameters.outdir_ = solPath;
  generationParameters.logfile_ = logSolver;
  generationParameters.absoluteGap_ = 5;
  generationParameters.minRelativeGap_ = .05;
  generationParameters.relativeGap_ = .1;
  generationParameters.nbDiveIfMinGap_ = 1;
  generationParameters.nbDiveIfRelGap_ = 2;
  generationParameters.solveToOptimality_ = false;
  generationParameters.weightStrategy_ = RANDOMMEANMAX;
  generationParameters.sp_type_ = useRotation ? ALL_ROTATION : ROSTER;

  generationParameters_ = generationParameters;

  SolverParam evaluationParameters;
  evaluationParameters.verbose(verbose_);
  evaluationParameters.maxSolvingTimeSeconds_ = totalTimeLimitSeconds_ - 1;
  evaluationParameters.printEverySolution_ = false;
  // evaluationParameters.outdir_ = "";
  evaluationParameters.logfile_ = logSolver;
  evaluationParameters.absoluteGap_ = 5;
  evaluationParameters.minRelativeGap_ = .05;
  evaluationParameters.relativeGap_ = .1;
  evaluationParameters.nbDiveIfMinGap_ = 1;
  evaluationParameters.nbDiveIfRelGap_ = 2;
  evaluationParameters.solveToOptimality_ = false;
//  evaluationParameters.weightStrategy_ = BOUNDRATIO;
  evaluationParameters.stopAfterXSolution_ = 0;
  evaluationParameters.sp_type_ = useRotation ? ALL_ROTATION : ROSTER;

  evaluationParameters_ = evaluationParameters;
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

  SolverParam generationParameters;
  generationParameters.read(generationOptionsFile);
  generationParameters.outdir_ = solPath;
  generationParameters_ = generationParameters;
  generationParameters_.verbose_ = verbose_;

  SolverParam evaluationParameters;
  evaluationParameters.read(evaluationOptionsFile);
  evaluationParameters.logfile_ = logSolver;
  evaluationParameters_ = evaluationParameters;
  evaluationParameters_.verbose_ = verbose_;
}

bool StochasticSolverOptions::setParameter(
    const string &field, std::fstream *file) {
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
  } else if (Tools::strEndsWith(field, "nGenerationDemandsMax")) {
    *file >> nGenerationDemandsMax_;
  } else {
    return false;
  }
  return true;
}

void StochasticSolverOptions::read(const string &strFile) {
  Tools::loadOptions(
      strFile, [this](const std::string &field, std::fstream *file) {
        return setParameter(field, file);
      });
}
