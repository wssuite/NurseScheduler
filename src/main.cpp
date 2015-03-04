//
//  main.cpp
//  RosterDesNurses
//
//  Created by Jérémy Omer on 16/11/2014.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//

#include "main_test.h"
#include "MyTools.h"
#include "ReadWrite.h"
#include "Solver.h"
#include "Greedy.h"

int main(int argc, char** argv)
{
	std::cout << "Number of arguments= " << argc << std::endl;

	// Detect errors in the number of arguments
	//
	if (argc%2 != 1) {
		Tools::throwError("main: There should be an even number of arguments!");
	}
	else if (argc > 1 && (argc < 9 || argc > 15)) {
		Tools::throwError("main: There is either too many or not enough arguments!");
	}

	// During tests, only run some test methods
	//
	if (argc == 1) {
		std::cout << "Running the test methods..."<< std::endl;
		// Tests functions to check the functions one by one
		main_test();
	}
	// Nominal behavior of the executable, as required by INRCII
	//
	else {
		// Retrieve the file names in arguments
		//
		int narg = 1;
		string scenarioFile, initialHistoryFile, weekDataFile, solutionFile;
		string customInputFile, customOutputFile, randSeed;
		while (narg < argc) {
			std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;
			if (!strcmp(argv[narg],"--sce")) {
				scenarioFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--his")) {
				initialHistoryFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--week")) {
				weekDataFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--out")) {
				solutionFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--cusIn")) {
				customInputFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--cusOut")) {
				customOutputFile = argv[narg+1];
				narg += 2;
			}
			else if (!strcmp(argv[narg],"--rand")) {
				randSeed = argv[narg+1];
				narg += 2;
			}
			else {
				Tools::throwError("main: the argument does not match the expected list!");
			}
		}
		// Throw an error if a necessary input file is missing
		if ( scenarioFile.empty() || initialHistoryFile.empty()
				|| weekDataFile.empty() || solutionFile.empty() ) {
			Tools::throwError("A necessary file name is missing!");
		}

		// Read the input files
		//
		Scenario* pScen = ReadWrite::readScenario(scenarioFile);
		Demand* pWeekDemand = ReadWrite::readWeek(weekDataFile, pScen);
		ReadWrite::readHistory(initialHistoryFile,pScen);
		if (!customInputFile.empty()) {
			ReadWrite::readCustom(customInputFile, pScen);
		}

		// Instantiate the solver class as a test
		//
		Greedy* pSolverTest =
		new Greedy(pScen, pWeekDemand,	pScen->pWeekPreferences(), pScen->pInitialState());
		pSolverTest->constructiveGreedy();

		// Write the solution in the required output format
		//
		Tools::LogOutput outStream(solutionFile);
		outStream << pSolverTest->solutionToString();
		if (!customOutputFile.empty()) {
			// Todo: the method that writes custom outputs. Which outputs ?
		}
		// Todo: the method that writes the history file corresponding to the
		// solution
		string outputHistoryFile("history-week");
		outputHistoryFile += std::to_string(pScen->thisWeek()) + ".txt";
		std::cout << "Output history file: " << outputHistoryFile << std::endl;
	}

}
