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

#include "solvers/mp/sp/RosterSP.h"
#include "solvers/mp/sp/RotationSP.h"
#include "solvers/mp/sp/rcspp/boost/RosterSP.h"
#include "solvers/mp/sp/rcspp/boost/RotationSP.h"
#include "solvers/InitializeSolver.h"
#include "solvers/mp/sp/SubProblem.h"
#include "solvers/mp/RosterMP.h"
#include "solvers/mp/RotationMP.h"
#include "data/Shift.h"

std::pair<float, float> comparePricing(MasterProblem *pMaster,
                                       SubProblem *mSP,
                                       SubProblem *pSP2,
                                       SubProblemParam spParam,
                                       bool *errorFound) {
  PDualCosts pDualCosts = std::make_shared<DualCosts>(pMaster);
  pDualCosts->randomUpdateDuals(true);
  int verbose = spParam.verbose_;

  if (verbose)
    std::cout << "\n      SOLVE WITH MY RCSPP SOLVER : \n" <<std::endl;

  // solve it
  mSP->solve(pDualCosts);

  if (verbose) {
    std::cout << "## Optimal value = " << mSP->bestReducedCost() <<
              std::endl;
    std::cout << "## Cpu time = " << mSP->cpuInLastSolve() << std::endl;
    std::cout << "\n-----------------------------------\n" << std::endl;

    std::cout << "\n      SOLVE WITH 2nd RCSPP SOLVER : \n" << std::endl;
  }

  // solve it
  pSP2->solve(pDualCosts);

  if (verbose) {
    std::cout << "## Optimal value = " << pSP2->bestReducedCost() << std::endl;
    std::cout << "## Cpu time = " << pSP2->cpuInLastSolve() << std::endl;
  }

  if (std::fabs(pSP2->bestReducedCost() - mSP->bestReducedCost()) > EPSILON) {
    std::cerr << "\nNew pricer value = " << mSP->bestReducedCost() <<
              "; 2nd value = " << pSP2->bestReducedCost() << std::endl;
    /*Tools::throwError("The new pricer does not find the same optimal value "
                      "as Boost");*/
    vector<RCSolution> bSols = pSP2->getSolutions(),
        mSols = mSP->getSolutions();
    RCSolution::sort(&bSols);
    RCSolution::sort(&mSols);
#ifdef DBG
    std::cout << bSols.front().toString();
    std::cout << mSols.front().toString();
    if (pMaster->parameters().sp_type_ == ROSTER) {
      RosterColumn bPat(bSols.front(), 0);
      bPat.checkReducedCost(*pDualCosts, true);
      RosterColumn mPat(mSols.front(), 0);
      mPat.checkReducedCost(*pDualCosts, true);
    } else {
      RotationColumn bPat(bSols.front(), 0);
      bPat.checkReducedCost(*pDualCosts, true);
      RotationColumn mPat(mSols.front(), 0);
      mPat.checkReducedCost(*pDualCosts, true);
    }
#endif
    *errorFound = true;
  }

  return {mSP->cpuInLastSolve(), pSP2->cpuInLastSolve()};
}

float comparePricing(
    MasterProblem *pMaster, bool *errorFound, bool compareToBoost, int nTest) {
  // solve the pricing problems
  PScenario pScenario = pMaster->pScenario();
  float totalCPU = 0;
  std::cout << "CPU reduction factors = ";
  for (const PLiveNurse& pNurse : pMaster->pLiveNurses()) {
    // Solve with boost RCSPP solver first
    if (pMaster->parameters().verbose_)
      std::cout << "Solve subproblem for nurse " << pNurse->name_ << std::endl;
    // create param
    SubProblemParam spParam(pMaster->parameters());
    spParam.strategyLevel_ =
        boostRCSPP::SubProblem::maxSubproblemStrategyLevel_;
    // retrieve a boost subproblem
    SubProblem *mSP, *pSP2;
    SubProblemParam spParam2 = spParam;
    if (!compareToBoost)
      spParam2.rcsppImprovedDomination_ = false;
    if (pMaster->parameters().sp_type_ == ROSTER) {
      mSP = new RosterSP(pScenario, 0, pScenario->nDays(), pNurse,
                         pMaster->getSPResources(pNurse), spParam);
      if (compareToBoost)
        pSP2 = new boostRCSPP::RosterSP(
            pScenario, pScenario->nDays(), pNurse, spParam2);
      else
        pSP2 = new RosterSP(pScenario, 0, pScenario->nDays(), pNurse,
                            pMaster->getSPResources(pNurse), spParam2);
    } else {
      mSP = new RotationSP(pScenario, 0, pScenario->nDays(), pNurse,
                           pMaster->getSPResources(pNurse), spParam);
      if (compareToBoost)
        pSP2 = new boostRCSPP::RotationSP(
            pScenario, pScenario->nDays(), pNurse, spParam2);
      else
        pSP2 = new RotationSP(pScenario, 0, pScenario->nDays(), pNurse,
                              pMaster->getSPResources(pNurse), spParam2);
    }
    mSP->build();
    pSP2->build();

    // build random dual costs
    float cpu = 0;
    for (int i=0; i < nTest; i++) {
      auto p = comparePricing(pMaster, mSP, pSP2, spParam, errorFound);
      cpu += p.second/p.first;
    }
    cpu /= static_cast<float>(nTest);
    totalCPU += cpu;

    // delete
    delete mSP;
    delete pSP2;

    std::cout << cpu;
    if ((pNurse->num_+1)%10 == 0)
      std::cout << std::endl << "                        ";
    else
      std::cout << ", ";
    std::flush(std::cout);
  }

  totalCPU /= static_cast<float>(pMaster->nNurses());

  if (pMaster->parameters().verbose_)
    std::cout << "\nCPU reduction factor = " << totalCPU << std::endl;

  return totalCPU;
}

float test_pricer(const std::string &instance,
                  bool *errorFound,
                  bool compareToBoost = true,
                  int verbose = 0,
                  int nTest = 1,
                  const SubProblemParam& spParam = SubProblemParam()) {
  // set input path
  std::vector<std::string> p = Tools::tokenize<std::string>(instance, '_');
  InputPaths inputPaths(
      "datasets/INRC2/",
      p[0],  // instance name
      p.size() > 2 ? std::stoi(p[1]) :0,  // history
      Tools::tokenize<int>(p[p.size()-1], '-'),  // week indices
      "", "", "paramfiles/default.txt", LARGE_TIME, verbose, 1, "ROSTER", 1);
  // set the scenario
  if (verbose)
    std::cout << "# INITIALIZE THE SCENARIO" << std::endl;
  PScenario pScenario;
  if (inputPaths.nbWeeks() > 1)
    pScenario = initializeMultipleWeeksINRC2(inputPaths);
  else
    pScenario = initializeScenarioINRC2(inputPaths);
  std::cout << std::endl;
  // initialize random of tools
  Tools::initializeRandomGenerator(inputPaths.randSeed());

  // initialize the solver, build a master problem,
  // then solve the all the pricing problems
  if (inputPaths.verbose())
    std::cout << "# SOLVE THE INSTANCE" << std::endl;
  auto pSolver = new DeterministicSolver(pScenario, inputPaths);
  SolverParam param = pSolver->getCompleteParameters();
  param.spParam_ = spParam;
  param.verbose_ = inputPaths.verbose();
  param.sp_type_ = ROSTER;
  auto pMaster =
      dynamic_cast<MasterProblem *>(
          pSolver->setSolverWithInputAlgorithm(
              pSolver->getOptions().solutionAlgorithm_,
              param));
  pMaster->initialize(param, {});
  // make the tests
  float cpu = comparePricing(pMaster, errorFound, compareToBoost, nTest);
  delete pMaster;
  std::cout << "\nComparison finished for " << instance << "." << std::endl;
  std::cout << "CPU reduction factor = " << cpu  << std::endl;
  return cpu;
}

// TODO(AL): some comment would be welcome about the format of input
int main(int argc, char **argv) {
  std::cout << "# TEST THE NEW PRICER" << std::endl;

  std::vector<string> insts;
  if (argc <= 1)
    insts = {"n005w4_2-0-2-1", "n030w4_1-2-3-4", "n012w8_3-5-0-2-0-4-5-2"};
  else
    insts = Tools::tokenize<string>(string(argv[1]), ',');
  int ntests = argc <= 2 ? 30 : std::stoi(argv[2]);
  bool errorFound = false;

  auto test = [&](const string& name, const SubProblemParam& spParam) {
    std::cout << "###############################" << std::endl;
    std::vector<float> cpus;
    cpus.reserve(insts.size());
    for (const string &inst : insts)
      cpus.push_back(test_pricer(inst, &errorFound, true, 0, ntests, spParam));
    for (size_t i = 0; i < insts.size(); ++i)
      std::cout << name << " CPU reduction factor for "
                << insts[i] << " = " << cpus[i] << std::endl;
  };

  SubProblemParam spParam;
  test("[Default]", spParam);

//  spParam.rcsppImprovedDomination_ = true;
//  spParam.rcsppEnumSubpaths_ = true;
//  test("[Sub path]", spParam);
//
//  spParam.rcsppEnumSubpaths_ = false;
//  spParam.rcsppMinCostToSinks_ = true;
//  test("[Min cost sink]", spParam);
//
//  spParam.rcsppEnumSubpathsForMinCostToSinks_ = true;
//  test("[Min sink sub path cost]", spParam);

  return errorFound;
}
