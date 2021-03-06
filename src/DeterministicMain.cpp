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

#include <exception>

#include "solvers/InitializeSolver.h"
#include "tools/Tools.h"
#include "tools/ReadWrite.h"
#include "solvers/Solver.h"
#include "solvers/DeterministicSolver.h"
#include "DeterministicMain_test.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

/******************************************************************************
* Solve the complete planning horizon with the deterministic solver
******************************************************************************/

int solveDeterministic(const InputPaths &inputPaths, double timeout) {
  // set the scenario
  //
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

  // initialize the solver and call the generic solution where the
  // specific solution processes are called
  //
  std::cout << "# SOLVE THE INSTANCE" << std::endl;
  DeterministicSolver *pSolver = new DeterministicSolver(pScenario, inputPaths);
  double objValue = pSolver->solve();
  std::cout << std::endl;

  // Display the solution and write the files for the validator
  //
  std::cout << "# FINAL SOLUTION" << std::endl;
  std::cout << "# Solution status = "
            << statusToString.at(pSolver->status()) << std::endl;
  std::cout << "# Objective value = ";
  bool noSolution =
      pSolver->status() == INFEASIBLE || pSolver->status() == UNSOLVED;
  if (noSolution)
    std::cout << "  -  ";
  else
    std::cout << objValue;
  std::cout << std::endl;
  pSolver->displaySolutionMultipleWeeks(inputPaths);

  // Write the final statistics
  //
  string statPath =
      inputPaths.solutionPath().empty() ? "" : inputPaths.solutionPath()
          + "/stat.txt";
  Tools::LogOutput statStream(statPath);
  statStream << pSolver->getGlobalStat().toString() << std::endl;

  if (pSolver->getOptions().withLNS_) {
    string lnsStatPath =
        inputPaths.solutionPath().empty() ? "" : inputPaths.solutionPath()
            + "/lns_stat.txt";
    Tools::LogOutput lnsStatStream(lnsStatPath);
    lnsStatStream << pSolver->getGlobalStat().lnsStatsToString() << std::endl;
    lnsStatStream.close();
  }

  //  release memory
  delete pSolver;
  statStream.close();

  // 1 if no solution found, 0 otherwise
  return noSolution;
}

/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char **argv) {
  std::cout << "# SOLVE THE PROBLEM WITH DETERMINISTIC DEMAND" << std::endl;
  std::cout << "Number of arguments= " << argc << std::endl;

  // Detect errors in the number of arguments
  //
  if (argc % 2 != 1) {
    Tools::throwError("main: There should be an even number of arguments!");
  }

  // Retrieve the file names in arguments
  //
  int narg = 1;
  InputPaths *pInputPaths = 0;
  string solutionFile = "";
  double timeout = 100.0;

  // On se limite à trois arguments pour les tests
  //
  if (argc == 3 && !strcmp(argv[1], "--test")) {
    std::cout << "arg = " << argv[narg] << " " << argv[narg + 1] << std::endl;

    // Procédures de test
    if (!strcmp(argv[2], "divide")) {
      testDivideIntoConnectedComponents();
    }

    return 0;
  }

  // Read the arguments and store them in pInputPaths
  pInputPaths = readArguments(argc, argv);


  // Initialize the random seed
  //
  srand(pInputPaths->randSeed());

  // Solve the problem
  //
  int r = solveDeterministic(*pInputPaths, timeout);

  // Release memory
  //
  delete pInputPaths;

  return r;
}
