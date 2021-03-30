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
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/sp/SubProblem.h"
#include "data/Shift.h"

#define DBG_AG

std::pair<float, float> comparePricing(PDualCosts pDualCosts,
                                       PLiveNurse pNurse,
                                       SubProblem *bSP,
                                       SubProblem *mSP,
                                       const SubProblemParam &spParam,
                                       bool *errorFound) {
  if (spParam.verbose_)
    std::cout << "\n      SOLVE WITH BOOST RCSPP SOLVER : \n" <<std::endl;

  // solve it
  bSP->solve(pDualCosts);

  if (spParam.verbose_) {
    std::cout << "## Optimal value = " << bSP->bestReducedCost() << std::endl;
    std::cout << "## Cpu time = " << bSP->cpuInLastSolve() << std::endl;

    std::cout << "\n      SOLVE WITH MY RCSPP SOLVER : \n" <<std::endl;
  }

  // solve it
  mSP->solve(pDualCosts);

  if (spParam.verbose_) {
    std::cout << "## Optimal value = " << mSP->bestReducedCost() <<
              std::endl;
    std::cout << "## Cpu time = " << mSP->cpuInLastSolve() << std::endl;
    std::cout << "\n-----------------------------------\n" << std::endl;
  }

  if (spParam.rcsppToOptimality_ &&
      std::fabs(bSP->bestReducedCost() - mSP->bestReducedCost()) > 1.0e-4 ) {
    std::cerr << "\nBoost value = " << bSP->bestReducedCost() <<
              "; New pricer value = " << mSP->bestReducedCost() << std::endl;
    /*Tools::throwError("The new pricer does not find the same optimal value "
                      "as Boost");*/
    *errorFound = true;
  }

  return {bSP->cpuInLastSolve(), mSP->cpuInLastSolve()};
}

float comparePricingToBoost(
    MasterProblem *pMaster, bool *errorFound, int nTest) {
  // solve the pricing problems
  PScenario pScenario = pMaster->pScenario();
  float totalCPU = 0;
  std::cout << "CPU reduction factors = ";
  for (const PLiveNurse pNurse : pMaster->liveNurses()) {
    // Solve with boost RCSPP solver first
    if (pMaster->parameters().verbose_)
      std::cout << "Solve subproblem for nurse " << pNurse->name_ << std::endl;
    // create param
    SubProblemParam spParam(pNurse, pMaster->parameters());
    spParam.strategyLevel_ =
        boostRCSPP::SubProblem::maxSubproblemStrategyLevel_;
    // retrieve a boost subproblem
    SubProblem *bSP, *mSP;
    if (pMaster->parameters().sp_type_ == ROSTER) {
      bSP = new boostRCSPP::RosterSP(
          pScenario, pScenario->nDays(), pNurse, spParam);
      mSP = new RosterSP(pScenario, pScenario->nDays(), pNurse,
                         pMaster->createPResources(pNurse), spParam);
    } else {
      bSP = new boostRCSPP::RotationSP(
          pScenario, pScenario->nDays(), pNurse, spParam);
      mSP = new RotationSP(pScenario, pScenario->nDays(), pNurse,
                           pMaster->createPResources(pNurse), spParam);
    }
    bSP->build();
    mSP->build();

    // build random dual costs
    float cpu = 0;
    for (int i=0; i < nTest; i++) {
      PDualCosts pDualCosts = pMaster->buildRandomDualCosts(true);
      auto p =
          comparePricing(pDualCosts, pNurse, bSP, mSP, spParam, errorFound);
      cpu += p.first/p.second;
    }
    cpu /= nTest;
    totalCPU += cpu;

    // delete
    delete bSP;
    delete mSP;

    std::cout << cpu;
    if ((pNurse->num_+1)%10 == 0)
      std::cout << std::endl << "                        ";
    else
      std::cout << ", ";
    std::flush(std::cout);
  }

  totalCPU /= pMaster->nNurses();

  if (pMaster->parameters().verbose_)
    std::cout << "\nCPU reduction factor = " << totalCPU << std::endl;

  return totalCPU;
}

float test_pricer(const std::string &instance,
                  bool *errorFound,
                  int verbose = 0,
                  int nTest = 1,
                  SubProblemParam spParam = SubProblemParam()) {
  // set input path
  std::vector<std::string> p = Tools::tokenize<std::string>(instance, '_');
  InputPaths inputPaths(
      "datasets/",
      p[0],  // instance name
      p.size() > 2 ? std::stoi(p[1]) :0,  // history
      Tools::tokenize<int>(p[p.size()-1], '-'),  // week indices
      "", "", "paramfiles/default.txt", LARGE_TIME, verbose, 1, "ROSTER", 1);
  // set the scenario
  if (verbose)
    std::cout << "# INITIALIZE THE SCENARIO" << std::endl;
  PScenario pScenario;
  if (inputPaths.nbWeeks() > 1)
    pScenario = initializeMultipleWeeks(inputPaths);
  else
    pScenario = initializeScenario(inputPaths);
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
  pMaster->initialize(param);
  // make the tests
  float cpu = comparePricingToBoost(pMaster, errorFound, nTest);
  delete pMaster;
  std::cout << "\nComparison finished for " << instance << "." << std::endl;
  std::cout << "CPU reduction factor = " << cpu  << std::endl;
  return cpu;
}

int main(int argc, char **argv) {
  std::cout << "# TEST THE NEW PRICER" << std::endl;

  std::vector<string> insts;
  if (argc <= 1)
    insts = {"n005w4_2-0-2-1", "n030w4_1-2-3-4", "n012w8_3-5-0-2-0-4-5-2"};
  else
    insts = Tools::tokenize<string>(string(argv[1]), ',');
  int n = argc <= 2 ? 30 : std::stoi(argv[2]);
  bool errorFound = false;

  auto test = [&](string name, SubProblemParam spParam) {
    std::cout << "###############################" << std::endl;
    std::vector<float> cpus;
    for (const string &inst : insts)
      cpus.push_back(test_pricer(inst, &errorFound, 0, n, spParam));
    for (size_t i = 0; i < insts.size(); ++i)
      std::cout << name << " CPU reduction factor for "
                << insts[i] << " = " << cpus[i] << std::endl;
  };

  SubProblemParam spParam;
  test("[Default]", spParam);

  spParam.rcsppEnumSubpaths_ = true;
  test("[Sub path]", spParam);

  spParam.rcsppEnumSubpaths_ = false;
  spParam.rcsppMinCostToSinks_ = true;
  test("[Min cost sink]", spParam);

  spParam.rcsppEnumSubpathsForMinCostToSinks_ = true;
  test("[Min sink sub path cost]", spParam);

  return errorFound;
}
