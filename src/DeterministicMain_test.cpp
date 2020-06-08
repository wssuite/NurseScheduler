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

# include "DeterministicMain_test.h"

// some include files to go through the files of an input directory
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "solvers/InitializeSolver.h"
#include "tools/ReadWrite.h"
#include "solvers/DeterministicSolver.h"
#include "tools/MyTools.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

// Test the result of the method that divides the scenario
// according to the connected components of positions
//
bool testDivideIntoConnectedComponents() {
  // Test directory and instance, the same is used in every test
  string dataDir = "datasets/";  // testdatasets datasets
  string instanceName = "n030w4";  // n100w4 n030w4 n005w4
  string outDir = "outfiles/" + instanceName + "/";
  int historyIndex = 1;
  vector<int> weekIndices = {6, 2, 9, 1};

  // build the paths of the input files
  InputPaths inputPaths(dataDir, instanceName, historyIndex, weekIndices);

  // set the scenario
  PScenario pScenarioInitial;
  pScenarioInitial = initializeMultipleWeeks(inputPaths);

  std::cout << "INITIAL SCENARIO" << std::endl;
  std::cout << pScenarioInitial->toString() << std::endl;

  DeterministicSolverOptions optionsInitial;
  DeterministicSolver *pSolverInitial =
      new DeterministicSolver(pScenarioInitial, inputPaths);

  pSolverInitial->solveByConnectedPositions();

  // Display the complete solution
  pSolverInitial->displaySolutionMultipleWeeks(inputPaths);

  delete pSolverInitial;

  return true;
}
