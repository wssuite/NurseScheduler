//
//  InitializeSolver.h
//  RosterDesNurses
//
//  Created by Jeremy Omer on January 21, 2016
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//

#include "solvers/StochasticSolver.h"
#include "tools/InputPaths.h"
#include "solvers/DeterministicSolver.h"
#include "solvers/Solver.h"


// Read the arguments in noncompact/compact format
InputPaths* readNonCompactArguments(int argc, char** argv);
InputPaths* readCompactArguments(int argc, char** argv);

//Initialize the week scenario by reading the input files
PScenario initializeScenario(std::string scenFile, std::string demandFile, std::string historyFile, std::string logFile="");
PScenario initializeScenario(const InputPaths & inputPaths, std::string logPath="");

// Initialize the scenario for multiple weeks
// When calling this function, the intent is to solve all the weeks at once
PScenario initializeMultipleWeeks(std::string dataDir, std::string instanceName,
  int historyIndex, std::vector<int> weekIndices, std::string logPath="");
PScenario initializeMultipleWeeks(const InputPaths & inputPaths, std::string logPath="");

// Separate the scenario into multiple scenarios that only focus on the nurses
// whose positions are in the same connex component of positions
std::vector<PScenario> divideScenarioIntoConnexPositions(PScenario pScenario);

// Solve the complete planning horizon with the deterministic solver
//void solveDeterministic(InputPaths inputPaths, string solPath, string logPathIni, double timeout);

// Create a solver of the class specified by the input algorithm type
Solver* setSolverWithInputAlgorithm(PScenario pScen, Algorithm algorithm);

// When a solution of multiple consecutive weeks is available, load it in a
// solver for all the weeks and  display the results
void displaySolutionMultipleWeeks(std::string dataDir, std::string instanceName,
	int historyIndex, std::vector<int> weekIndices, std::vector<Roster> &solution, Status status, std::string outPath="");
void displaySolutionMultipleWeeks(InputPaths inputPaths, std::vector<Roster> &solution, Status status, std::string outDir="");

// Compute and record stats on all the demand files of all the instances in the
// input directory
void computeStatsOnTheDemandsOfAllInstances(std::string inputDir);
