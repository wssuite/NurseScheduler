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

#ifndef SRC_SOLVERS_INITIALIZESOLVER_H_
#define SRC_SOLVERS_INITIALIZESOLVER_H_

#include <string>
#include <vector>

#include "solvers/StochasticSolver.h"
#include "tools/InputPaths.h"
#include "solvers/DeterministicSolver.h"
#include "solvers/Solver.h"


// Initialize the resources of all the nurses
void initializeResourcesINRC2(const PScenario& pScenario);

// Initialize the week scenario by reading the input files
PScenario initializeScenarioINRC2(const InputPaths &inputPaths,
                                  const std::string& logPath = "");

// Initialize the scenario for multiple weeks
// When calling this function, the intent is to solve all the weeks at once
PScenario initializeMultipleWeeksINRC2(const std::string& dataDir,
                                       const std::string& instanceName,
                                       int historyIndex,
                                       std::vector<int> weekIndices,
                                       const std::string& logPath = "");
PScenario initializeMultipleWeeksINRC2(const InputPaths &inputPaths,
                                       const std::string& logPath = "");

// Separate the scenario into multiple scenarios that only focus on the nurses
// whose positions are in the same connected component of positions
std::vector<PScenario> divideScenarioIntoConnectedPositions(
    const PScenario& pScenario);

// When a solution of multiple consecutive weeks is available, load it in a
// solver for all the weeks and  display the results
void displaySolutionMultipleWeeks(const std::string& dataDir,
                                  const std::string& instanceName,
                                  int historyIndex,
                                  const std::vector<int>& weekIndices,
                                  const std::vector<Roster> &solution,
                                  Status status,
                                  const std::string& outPath = "");
void displaySolutionMultipleWeeks(const InputPaths& inputPaths,
                                  const std::vector<Roster> &solution,
                                  Status status,
                                  const std::string& outDir = "");

// Compute and record stats on all the demand files of all the instances in the
// input directory
void computeStatsOnTheDemandsOfAllInstances(const std::string& inputDir);

#endif  // SRC_SOLVERS_INITIALIZESOLVER_H_
