//
//  InitializeSolver.h
//  RosterDesNurses
//
//  Created by Jeremy Omer on January 21, 2016
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//

#include "StochasticSolver.h"
#include "InputPaths.h"
#include "DeterministicSolver.h"
#include "Solver.h"


// Read the arguments in noncompact/compact format
InputPaths* readNonCompactArguments(int argc, char** argv);
InputPaths* readCompactArguments(int argc, char** argv);

//Initialize the week scenario by reading the input files
Scenario* initializeScenario(string scenFile, string demandFile, string historyFile, string logFile="");
Scenario* initializeScenario(const InputPaths & inputPaths, string logPath="");

// Initialize the scenario for multiple weeks
// When calling this function, the intent is to solve all the weeks at once
Scenario* initializeMultipleWeeks(string dataDir, string instanceName,
  int historyIndex, vector<int> weekIndices, string logPath="");
Scenario* initializeMultipleWeeks(const InputPaths & inputPaths, string logPath="");

// Separate the scenario into multiple scenarios that only focus on the nurses
// whose positions are in the same connex component of positions
vector<Scenario*> divideScenarioIntoConnexPositions(Scenario* pScenario);

// Solve the complete planning horizon with the deterministic solver
//void solveDeterministic(InputPaths inputPaths, string solPath, string logPathIni, double timeout);

// Create a solver of the class specified by the input algorithm type
Solver* setSolverWithInputAlgorithm(Scenario* pScen, Algorithm algorithm);

// When a solution of multiple consecutive weeks is available, load it in a
// solver for all the weeks and  display the results
void displaySolutionMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, vector<Roster> &solution, Status status, string outPath="");
void displaySolutionMultipleWeeks(InputPaths inputPaths, vector<Roster> &solution, Status status, string outDir="");

// Compute and record stats on all the demand files of all the instances in the
// input directory
void computeStatsOnTheDemandsOfAllInstances(string inputDir);
