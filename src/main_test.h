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

// Test the cbc modeler
void testCbc(Scenario* pScen, Demand* pDemand, Preferences* pPref,
  std::vector<State>* pStateIni, std::vector<Roster>& solIni);

// Test the random demand generator
void testRandomDemandGenerator(int nbDemands,string logFile, Scenario* pScen);

// Print the main characteristics of all the demands of an input directory
// This is done to find some invariant properties among demands
void compareDemands(std::string inputDir);
