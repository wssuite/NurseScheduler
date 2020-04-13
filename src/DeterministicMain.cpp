//
//  main.cpp
//  RosterDesNurses
//
//  Created by Jeremy Omer on 16/11/2014.
//  Copyright (c) 2014 Jeremy Omer. All rights reserved.
//

#include "solvers/InitializeSolver.h"
#include "tools/MyTools.h"
#include "tools/ReadWrite.h"
#include "solvers/Solver.h"
//#include "Greedy.h"
#include <exception>
#include "solvers/mp/modeler/Modeler.h"
#include "solvers/MasterProblem.h"
#include "solvers/DeterministicSolver.h"
#include "DeterministicMain_test.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

/******************************************************************************
* Solve the complete planning horizon with the deterministic solver
******************************************************************************/

void solveDeterministic(InputPaths inputPaths, string solPath, string logPathIni, double timeout) {

	// set the scenario
	//
	std::cout << "# INITIALIZE THE SCENARIO" << std::endl;
	Scenario* pScenario;
	if (inputPaths.nbWeeks() > 1) {
		pScenario = initializeMultipleWeeks(inputPaths);
	}
	else {
		pScenario = initializeScenario(inputPaths);
	}
	std::cout << std::endl;

	//initialize random of tools
	Tools::initializeRandomGenerator(inputPaths.randSeed());

	// DBG
	std::cout << "Next random : " << Tools::randomInt(0, RAND_MAX) << std::endl;
	std::cout << "Next random : " << Tools::randomInt(0, RAND_MAX) << std::endl;


	// initiialize the solver and call the generic solution where the
	// specific solution processes are called
	//
	std::cout << "# SOLVE THE INSTANCE" << std::endl;
	DeterministicSolver* pSolver = new DeterministicSolver(pScenario,inputPaths);
	double objValue = pSolver->solve();
	std::cout << std::endl;

	// Display the solution and write the files for the validator
	//
	std::cout << "# FINAL SOLUTION" << std::endl;
	std::string solutionStatus = statusToString[pSolver->getStatus()];
	std::cout << "# Solution status = " << solutionStatus <<  std::endl;
	std::cout << "# Objective value = " << objValue <<  std::endl;
	pSolver->displaySolutionMultipleWeeks(inputPaths);

	// Write the final statistics
	//
	string statPath = inputPaths.solutionPath().empty() ? "" : inputPaths.solutionPath()+"/stat.txt";
	Tools::LogOutput statStream(statPath);
	statStream << pSolver->getGlobalStat().toString() << std::endl;

	if (pSolver->getOptions().withLNS_) {
		string lnsStatPath = inputPaths.solutionPath().empty() ? "" : inputPaths.solutionPath()+"/lns_stat.txt";
		Tools::LogOutput lnsStatStream(lnsStatPath);
		lnsStatStream << pSolver->getGlobalStat().lnsStatsToString() << std::endl;
		lnsStatStream.close();
	}


	//  release memory
	if (pSolver) delete pSolver;
	if (pScenario) delete pScenario;
	statStream.close();
}


/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char** argv)
{
	std::cout << "# SOLVE THE PROBLEM WITH DETERMINISTIC DEMAND" << std::endl;
	std::cout << "Number of arguments= " << argc << std::endl;

	// Detect errors in the number of arguments
	//
	if (argc%2 != 1) {
		Tools::throwError("main: There should be an even number of arguments!");
	}

	// Retrieve the file names in arguments
	//
	int narg = 1;
	InputPaths* pInputPaths=0;
	string solutionFile="";
	double timeout = 100.0;

	// On se limite à trois arguments pour les tests
	//
	if (argc == 3 && !strcmp(argv[1],"--test") ) {
		std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;

		// Procédures de test
		if (!strcmp(argv[2], "divide")) {
			testDivideIntoConnexComponents();
		}

		return 0;
	}

	// Read the arguments and store them in pInputPaths
	// If in non compact format, each week is input, so there are at least 19 arguments
	// In compact format, the number of arguments is smaller than that
	//
	if (argc >= 21) {
		pInputPaths = readNonCompactArguments(argc,argv);
	}
	else {
		pInputPaths = readCompactArguments(argc,argv);
	}

	// Initialize the random seed
	//
	srand(pInputPaths->randSeed());

	// Solve the problem
	//
	solveDeterministic(*pInputPaths, pInputPaths->solutionPath(), pInputPaths->logPath(), timeout);

	// Release memory
	//
	if (pInputPaths) delete pInputPaths;

}
