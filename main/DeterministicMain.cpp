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


#include "ParseArguments.h"
#include "InitializeInstance.h"
#include "solvers/Solver.h"
#include "solvers/DeterministicSolver.h"
#include "tools/Tools.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

/******************************************************************************
* Solve the complete planning horizon with the deterministic solver
******************************************************************************/

int solveDeterministic(const InputPaths &inputPaths) {
  // initialize random of tools
  Tools::initializeRandomGenerator(inputPaths.randSeed());

  // set the scenario
  std::cout << "# INITIALIZE THE SCENARIO: " << inputPaths.scenario()
            << std::endl;
  PScenario pScenario = buildInstance(inputPaths);
  std::cout << std::endl;

  // initialize the solver and call the generic solution where the
  // specific solution processes are called
  //
  std::cout << "# SOLVE THE INSTANCE" << std::endl;

  auto *pSolver = new DeterministicSolver(pScenario, inputPaths);

  double objValue = pSolver->solve();
  std::cout << std::endl;

  // Display the solution and write the files for the validator
  //
  bool noSolution =
      pSolver->status() == INFEASIBLE || pSolver->status() == UNSOLVED;
  std::cout << "# FINAL SOLUTION" << std::endl;
  std::cout << pSolver->coverageToString() << std::endl;
  std::cout << pSolver->solutionToLogString() << std::endl;
  std::cout << pSolver->writeResourceCosts() << std::endl;
  if (!noSolution && !pSolver->solution().empty())
    std::cout << pSolver->writeResourceCostsPerNurse() << std::endl;
  std::cout << "# Solution status = "
            << statusToString.at(pSolver->status()) << std::endl;
  std::cout << "# Objective value = ";
  if (noSolution)
    std::cout << "  -  ";
  else
    std::cout << objValue;
  std::cout << std::endl;
  std::cout << "# Running time = " << std::setprecision(1) << std::fixed
            << pSolver->getGlobalStat().timeTotal() << " sec." << std::endl;
  if (!noSolution) pSolver->displaySolutionMultipleWeeks();

  // Write the final statistics
  string statPath = inputPaths.solutionPath().empty() ?
                    "" : inputPaths.solutionPath() + "stat.txt";
  Tools::LogOutput statStream(statPath);
  statStream.printnl(pSolver->getGlobalStat().toString());

  if (pSolver->getOptions().withLNS_) {
    string lnsStatPath =
        inputPaths.solutionPath().empty() ? "" : inputPaths.solutionPath()
            + "lns_stat.txt";
    Tools::LogOutput lnsStatStream(lnsStatPath);
    statStream.printnl(pSolver->getGlobalStat().lnsStatsToString());
  }
  // Write the outputfile
  if (inputPaths.ui()) {
    pSolver->solutionToUI(inputPaths.solutionPath());
  } else if (inputPaths.inrc() && !noSolution && !pSolver->solution().empty()) {
    // Write the solution in an xml file with INRC format
    pSolver->solutionToXmlINRC(inputPaths.solutionPath());
  }

  // display memory used
  double memGB = Tools::getResidentMemoryGB();
  std::cout << std::endl << "The program has consumed "
            << std::setprecision(3) << memGB << " GB of memory." << std::endl;

  //  release memory
  delete pSolver;

  // 1 if no solution found, 0 otherwise
  return noSolution;
}

/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char **argv) {
  std::cout << "# SOLVE THE PROBLEM WITH DETERMINISTIC DEMAND" << std::endl;
  std::cout << "Number of arguments= " << argc << std::endl;

  // Read the arguments and store them in pInputPaths
  InputPaths *pInputPaths = readArguments(argc, argv);

  if (pInputPaths->inrc2()) PrintSolution::writeMultiWeeks = true;
  if (pInputPaths->inrc()) PrintSolution::writeXML = true;

  // Solve the problem
  int r = solveDeterministic(*pInputPaths);

  // Release memory
  delete pInputPaths;

  return r;
}
