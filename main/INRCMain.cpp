/*
 * Copyright (C) 2021 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "solvers/mp/RotationMP.h"
#include "solvers/mp/RosterMP.h"
#include "solvers/mp/sp/RosterSP.h"
#include "solvers/mp/sp/RotationSP.h"
#include "solvers/mp/sp/rcspp/boost/RosterSP.h"
#include "solvers/mp/sp/rcspp/boost/RotationSP.h"
#include "solvers/InitializeSolver.h"
#include "solvers/mp/sp/SubProblem.h"
#include "data/Shift.h"
#include "ReadWrite.h"
#include "ParseArguments.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

/******************************************************************************
* Solve the complete planning horizon with the deterministic solver
******************************************************************************/

int solveINRC(InputPaths *pInputPaths) {
  // set the scenario
  // TODO(JO): there must be an error due to alternative skills that would be
  //  counted twice
  std::cout << "# INITIALIZE THE SCENARIO" << std::endl;
  PScenario pScenario = ReadWrite::readINRCInstance(pInputPaths->scenario());
  std::cout << std::endl;

  // initialize the solver and call the generic solution where the
  // specific solution processes are called
  std::cout << "# SOLVE THE INSTANCE" << std::endl;
  auto *pSolver = new DeterministicSolver(pScenario, *pInputPaths);
  double objValue = pSolver->solve();
  std::cout << std::endl;

  // Display the solution and write the files for the validator
  bool noSolution =
      pSolver->status() == INFEASIBLE || pSolver->status() == UNSOLVED;
  if (!noSolution && !pSolver->solution().empty()) {
    std::cout << "# FINAL SOLUTION" << std::endl;
    std::cout << pSolver->writeResourceCostsPerNurse() << std::endl;
  }
  std::cout << "# Solution status = "
            << statusToString.at(pSolver->status()) << std::endl;
  std::cout << "# Objective value = ";
  if (noSolution)
    std::cout << "  -  ";
  else
    std::cout << objValue;
  std::cout << std::endl;

  // Write the final statistics
  string statPath = pInputPaths->solutionPath()+"stat.txt";
  Tools::LogOutput statStream(statPath, true);
  statStream << pScenario->name() << "    "
             << pSolver->getGlobalStat().toString() << std::endl;

  // Write the solution in an xml file with INRC format
  if (!noSolution && !pSolver->solution().empty())
    pSolver->solutionToXmlINRC(pInputPaths->solutionPath());

  //  release memory
  delete pSolver;

  // 1 if no solution found, 0 otherwise
  return noSolution;
}


/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char **argv) {
  InputPaths *pInputPaths;
  if (argc == 3) {
    // Use basic arguments: datadir then instance file
    pInputPaths = new InputPaths(argv[1], argv[2]);
  } else {
    // Read the arguments and store them in pInputPaths
    pInputPaths = readArguments(argc, argv);
  }

  std::cout << "# SOLVE THE INRC STATIC INSTANCE: "
            << pInputPaths->scenario() << std::endl;

  // Initialize the random seed
  Tools::initializeRandomGenerator(0);

  // Solve the problem
  int r = solveINRC(pInputPaths);

  // Release memory
  delete pInputPaths;

  return r;
}


