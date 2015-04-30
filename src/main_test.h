//
//  main_test.h
//  RosterDesNurses
//
//  Created by Jérémy Omer on 04/03/2015.
//  Copyright (c) 2015 Jérémy Omer. All rights reserved.
//

#include "Solver.h"

// Main test function directly called in main.cpp
void main_test();

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

// Test the cbc modeler
void testCbc(Scenario* pScen);

// Test the random demand generator
void testRandomDemandGenerator(int nbDemands,string logFile, Scenario* pScen);

// Print the main characteristics of all the demands of an input directory
// This is done to find some invariant properties among demands
void compareDemands(std::string inputDir);

// Solve a deterministic input demand with the input algorithm
// In this method, we assume that all the demands are knwon in advance
// (the method can also treat only one week)
void testMultipleWeeksDeterministic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outPath);

// Test a solution on multiple weeks
// In this method, the weeks are solved sequentially without knowledge of future
// demand
void testMultipleWeeksStochastic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algo, string logPath="");

// When a solution of multiple consecutive weeks is available, load it in a
// solver for all the weeks and  display the results
void displaySolutionMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, vector<Roster> &solution, string outPath="");
