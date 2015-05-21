//
//  main_test.h
//  RosterDesNurses
//
//  Created by J��r��my Omer on 04/03/2015.
//  Copyright (c) 2015 J��r��my Omer. All rights reserved.
//

#include "StochasticSolver.h"
#include "Solver.h"


// The instances of this class contain the paths of the input files of the
// problem
class InputPaths{
public:
  InputPaths(string dataDir, string instanceName,int historyIndex, vector<int> weekIndices);

protected:
  string scenario_;
  string history_;
  vector<string> weeks_;

public:
  string scenario() {return scenario_;}
  string history() {return history_;}
  vector<string> weeks() {return weeks_;}
  string week(int w) {return weeks_[w];}
};

// Function for testing parts of the code (Antoine)
void testFunction_Antoine();

// Function for testing parts of the code (Jeremy)
void testFunction_Jeremy();

// Function for testing parts of the code (Samuel)
void testFunction_Samuel();


//Initialize the week scenario by reading the input files
Scenario* initializeScenario(string scenFile, string demandFile, string historyFile, string logFile="");

// Initialize the scenario for multiple weeks
// When calling this function, the intent is to solve all the weeks at once
Scenario* initializeMultipleWeeks(string dataDir, string instanceName,
  int historyIndex, vector<int> weekIndices, string logPath="");

// Solve one week inside the stochastic process
void solveOneWeek(string scenPath, string demandPath, string historyPath, string solPath, string logPath);

// Set the options of the stochastic solver
// This is not automated, so the options need to be changed inside the code 
// during the tests
// The solution time depends on the number of nurses and on the computed
void setStochasticSolverOptions(StochasticSolverOptions& options, Scenario* pScenario, string solPath);

// Solve a deterministic input demand with the input algorithm
// In this method, we assume that all the demands are knwon in advance
// (the method can also treat only one week)
void testMultipleWeeksDeterministic(string dataDir, string instanceName,
   int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outPath);
void testMultipleWeeksDeterministic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outPath, SolverParam param);

// Test a solution on multiple weeks
// In this method, the weeks are solved sequentially without knowledge of future
// demand
void testMultipleWeeksStochastic(string dataDir, string instanceName,
		int historyIndex, vector<int> weekIndices, StochasticSolverOptions stochasticSolverOptions, string outDir = "");

// Create a solver of the class specified by the input algorithm type
Solver* setSolverWithInputAlgorithm(Scenario* pScen, Algorithm algorithm);

// When a solution of multiple consecutive weeks is available, load it in a
// solver for all the weeks and  display the results
void displaySolutionMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, vector<Roster> &solution, Status status, string outPath="");

// Compute and record stats on all the demand files of all the instances in the
// input directory
void computeStatsOnTheDemandsOfAllInstances(string inputDir);

// Test the random demand generator
void testRandomDemandGenerator(int nbDemands,string logFile, Scenario* pScen);

// Test the cbc modeler
void testCbc(Scenario* pScen);
