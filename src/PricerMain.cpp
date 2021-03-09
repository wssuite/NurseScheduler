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

#include <solvers/mp/sp/RosterSP.h>
#include <solvers/mp/sp/rcspp/boost/RosterSP.h>
#include "solvers/InitializeSolver.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/sp/SubProblem.h"
#include "data/Shift.h"

#define DBG_AG

std::pair<float, float> comparePricing(PDualCosts pDualCosts,
                                       PLiveNurse pNurse,
                                       SubProblem *bSP,
                                       SubProblem *mSP,
                                       const SubproblemParam &spParam) {
  if (spParam.verbose_)
    std::cout << "\n      SOLVE WITH BOOST RCSPP SOLVER : \n" <<std::endl;

  // solve it
  bSP->solve(pNurse, pDualCosts, spParam);

  if (spParam.verbose_) {
    std::cout << "## Optimal value = " << bSP->bestReducedCost() << std::endl;
    std::cout << "## Cpu time = " << bSP->cpuInLastSolve() << std::endl;

    std::cout << "\n      SOLVE WITH MY RCSPP SOLVER : \n" <<std::endl;
  }

  // solve it
  mSP->solve(pNurse, pDualCosts, spParam);

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
  }

  return {bSP->cpuInLastSolve(), mSP->cpuInLastSolve()};
}

float comparePricingToBoost(const MasterProblem *pMaster) {
  // solve the pricing problems
  PScenario pScenario = pMaster->pScenario();
  double cpuBoost = 0.0;
  double cpuNew = 0.0;
  for (const PLiveNurse pNurse : pMaster->liveNurses()) {
    // Solve with boost RCSPP solver first
    if (pMaster->parameters().verbose_)
      std::cout << "Solve subproblem for nurse " << pNurse->name_ << std::endl;
    // create param
    SubproblemParam spParam(SubproblemParam::maxSubproblemStrategyLevel_,
                            pNurse,
                            pMaster->parameters());
    // retrieve a boost subproblem
    SubProblem *bSP =
        new boostRCSPP::RosterSP(pScenario,
                                 pScenario->nbDays(),
                                 pNurse->pContract_,
                                 pScenario->pInitialState());
    bSP->build();
    // Solve with my solver under development
    SubProblem *mSP =
        new RosterSP(pScenario, pScenario->nbDays(), pNurse,
                     pMaster->createResources(pNurse), spParam);
    mSP->build();

    // build random dual costs
    PDualCosts pDualCosts = pMaster->buildRandomDualCosts(true);
    auto p = comparePricing(pDualCosts, pNurse, bSP, mSP, spParam);
    cpuBoost += p.first;
    cpuNew += p.second;

    // resolve with new random dual costs
    pDualCosts = pMaster->buildRandomDualCosts(true);
    p = comparePricing(pDualCosts, pNurse, bSP, mSP, spParam);
    cpuBoost += p.first;
    cpuNew += p.second;

    // delete
    delete bSP;
    delete mSP;
  }

  if (pMaster->parameters().verbose_)
    std::cout << "\nCPU reduction factor = " << cpuBoost/cpuNew << std::endl;

  return cpuBoost/cpuNew;
}

float test_pricer(const std::string &instance, int verbose = 0, int nTest = 1) {
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
  param.verbose_ = inputPaths.verbose();  // 4: print also the solutions
  auto pMaster =
      dynamic_cast<MasterProblem *>(
          pSolver->setSolverWithInputAlgorithm(
              pScenario->pWeekDemand(),
              pSolver->getOptions().solutionAlgorithm_,
              param));
  pMaster->initialize(param);
  // make the tests
  float cpu = 0;
  std::cout << "CPU reduction factors = ";
  for (int i=0; i < nTest; i++) {
    float c = comparePricingToBoost(pMaster);
    std::cout << c;
    if ((i+1)%10 == 0) std::cout << std::endl << "                        ";
    else
      std::cout << ", ";
    std::flush(std::cout);
    cpu += c;
  }
  cpu /= nTest;
  delete pMaster;
  std::cout << "\nComparison finished for " << instance << "." << std::endl;
  std::cout << "CPU reduction factor = " << cpu  << std::endl;
  std::cout << "\n###############################"  << std::endl;
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
  std::vector<float> cpus;
  for (const string& inst : insts)
    cpus.push_back(test_pricer(inst, 0, n));

  std::cout << "\n###############################" << std::endl;
  std::cout << "###############################\n" << std::endl;

  for (size_t i = 0; i < insts.size(); ++i)
    std::cout << "\nCPU reduction factor for " << insts[i] << " = "
              << cpus[i] << std::endl;

  return 0;
}
