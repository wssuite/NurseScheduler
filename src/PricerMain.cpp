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

#include "solvers/InitializeSolver.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/sp/SubProblem.h"

int solvePricing(const InputPaths &inputPaths) {
  // set the scenario
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
  std::cout << "# SOLVE THE INSTANCE" << std::endl;
  DeterministicSolver *pSolver = new DeterministicSolver(pScenario, inputPaths);
  SolverParam param = pSolver->getCompleteParameters();
  param.verbose_ = 3;  // 4: print also the solutions
  MasterProblem *pMaster =
      dynamic_cast<MasterProblem *>(
          pSolver->setSolverWithInputAlgorithm(
              pScenario->pWeekDemand(),
              pSolver->getOptions().solutionAlgorithm_,
              param));
  pMaster->initialize(param);
  // solve the pricing problems
  // fetch pricer
  RCPricer *pPricer = dynamic_cast<RCPricer *>(pMaster->getPricer());
  for (PLiveNurse pNurse : pMaster->getLiveNurses()) {
    std::cout << "Solve subproblem for nurse " << pNurse->name_ << std::endl;
    // build random dual costs
    DualCosts dualCosts = pMaster->buildRandomDualCosts();
    // retreive a subproblem
    SubProblem *pSP = pPricer->retriveSubproblem(pNurse);
    // create param
    SubproblemParam spParam(SubproblemParam::maxSubproblemStrategyLevel_,
                            pNurse,
                            pMaster);
    // solve it
    pSP->solve(pNurse, &dualCosts, spParam);
    // release subproblem
    pPricer->releaseSubproblem(pNurse, pSP);
  }

  return 0;
}

int main(int argc, char **argv) {
  std::cout << "# SOLVE THE PRICING PROBLEM" << std::endl;
  std::cout << "Number of arguments= " << argc << std::endl;

  // Detect errors in the number of arguments
  //
  if (argc % 2 != 1)
    Tools::throwError("main: There should be an even number of arguments!");

  // Read the arguments and store them in pInputPaths
  // If in non compact format, each week is input,
  // so there are at least 19 arguments.
  // In compact format, the number of arguments is smaller than that
  InputPaths *pInputPaths;
  if (argc >= 21) {
    pInputPaths = readNonCompactArguments(argc, argv);
  } else {
    pInputPaths = readCompactArguments(argc, argv);
  }

  // Initialize the random seed
  srand(pInputPaths->randSeed());

  solvePricing(*pInputPaths);

  delete pInputPaths;
}
