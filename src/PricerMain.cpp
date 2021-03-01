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

#include <solvers/mp/sp/MyRosterSP.h>
#include <solvers/mp/sp/BoostRosterSP.h>
#include "solvers/InitializeSolver.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/sp/SubProblem.h"
#include "data/Shift.h"

#define DBG_AG

std::pair<float, float> comparePricing(DualCosts *dualCosts,
                                       PLiveNurse pNurse,
                                       SubProblem *bSP,
                                       SubProblem *mSP,
                                       const SubproblemParam &spParam) {
  if (spParam.verbose_)
    std::cout << "\n      SOLVE WITH BOOST RCSPP SOLVER : \n" <<std::endl;

  // solve it
  bSP->solve(pNurse, dualCosts, spParam);

  if (spParam.verbose_) {
    std::cout << "## Optimal value = " << bSP->bestReducedCost() << std::endl;
    std::cout << "## Cpu time = " << bSP->cpuInLastSolve() << std::endl;

    std::cout << "\n      SOLVE WITH MY RCSPP SOLVER : \n" <<std::endl;
  }

  // solve it
  mSP->solve(pNurse, dualCosts, spParam);

  if (spParam.verbose_) {
    std::cout << "## Optimal value = " << mSP->bestReducedCost() <<
              std::endl;
    std::cout << "## Cpu time = " << mSP->cpuInLastSolve() << std::endl;
    std::cout << "\n-----------------------------------\n" << std::endl;
  }

  if (std::abs(bSP->bestReducedCost() - mSP->bestReducedCost()) > 1.0e-4 )
    Tools::throwError("The new pricer does not find the same optimal value "
                      "as Boost");

  return {bSP->cpuInLastSolve(), mSP->cpuInLastSolve()};
}

float comparePricingToBoost(const InputPaths &inputPaths) {
  // set the scenario
  if (inputPaths.verbose())
    std::cout << "# INITIALIZE THE SCENARIO" << std::endl;
  PScenario pScenario;
  if (inputPaths.nbWeeks() > 1) {
    pScenario = initializeMultipleWeeks(inputPaths);
  } else {
    pScenario = initializeScenario(inputPaths);
  }
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
  // solve the pricing problems
  double cpuBoost = 0.0;
  double cpuNew = 0.0;
  for (const PLiveNurse pNurse : pMaster->getLiveNurses()) {
    // Solve with boost RCSPP solver first
    if (inputPaths.verbose())
      std::cout << "Solve subproblem for nurse " << pNurse->name_ << std::endl;
    // create param
    SubproblemParam spParam(SubproblemParam::maxSubproblemStrategyLevel_,
                            pNurse,
                            pMaster);
    // retrieve a boost subproblem
    SubProblem *bSP = new BoostRosterSP(pScenario,
                                        pScenario->nbDays(),
                                        pNurse->pContract_,
                                        pScenario->pInitialState());
    bSP->build();
    // Solve with my solver under development
    SubProblem *mSP =
        new MyRosterSP(pScenario, pScenario->nbDays(), pNurse, spParam);
    mSP->build();

    // build random dual costs
    DualCosts dualCosts = pMaster->buildRandomDualCosts(true);
    auto p = comparePricing(&dualCosts, pNurse, bSP, mSP, spParam);
    cpuBoost += p.first;
    cpuNew += p.second;

    // resolve with new random dual costs
    dualCosts = pMaster->buildRandomDualCosts(true);
    p = comparePricing(&dualCosts, pNurse, bSP, mSP, spParam);
    cpuBoost += p.first;
    cpuNew += p.second;

    // delete
    delete bSP;
    delete mSP;
  }

  if (inputPaths.verbose())
    std::cout << "\nRelative CPU reduction = " << 100.0* (cpuBoost - cpuNew)
        /cpuBoost << "%" << std::endl;

  return 100.0 * (cpuBoost - cpuNew) / cpuBoost;
}

float test_pricer(const std::string &instance, int verbose) {
  std::vector<std::string> p = Tools::tokenize<std::string>(instance, '_');
  InputPaths *pInputPaths = new InputPaths(
      "datasets/",
      p[0],  // instance name
      p.size() > 2 ? std::stoi(p[1]) :0,  // history
      Tools::tokenize<int>(p[p.size()-1], '-'),  // week indices
      "", "", "paramfiles/default.txt", LARGE_TIME, verbose, 1, "ROSTER", 1);
  float cpu = comparePricingToBoost(*pInputPaths);
  std::cout << "\nComparison finished for " << instance << "." << std::endl;
  std::cout << "\n###############################"  << std::endl;
  delete pInputPaths;
  return cpu;
}

int main(int argc, char **argv) {
  std::cout << "# TEST THE NEW PRICER" << std::endl;

  std::vector<string> insts;
  if (argc <= 1)
    insts = {"n005w4_2-0-2-1", "n030w4_1-2-3-4", "n012w8_3-5-0-2-0-4-5-2"};
  else
    insts = std::vector<string>(argv+1, argv+argc);

  std::vector<float> cpus;
  for (const string& inst : insts)
    cpus.push_back(test_pricer(inst, 0));

  std::cout << "\n###############################" << std::endl;
  std::cout << "###############################\n" << std::endl;

  for (size_t i = 0; i < insts.size(); ++i)
    std::cout << "\nRelative CPU reduction for " << insts[i] << " = "
              << cpus[i] << "%" << std::endl;

  return 0;
}
