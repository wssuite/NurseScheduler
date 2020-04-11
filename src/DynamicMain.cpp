//
//  main.cpp
//  RosterDesNurses
//
//  Created by J������r������my Omer on 16/11/2014.
//  Copyright (c) 2014 J������r������my Omer. All rights reserved.
//

#include <exception>

#include "tools/MyTools.h"
#include "tools/ReadWrite.h"
#include "tools/DemandGenerator.h"
#include "solvers/Greedy.h"
#include "solvers/MasterProblem.h"
#include "solvers/StochasticSolver.h"
#include "solvers/mp/SubProblem.h"
//#include "CbcModeler.h"
#include "tools/MyTools.h"
#include "solvers/InitializeSolver.h"

/******************************************************************************
* Solve one week inside the stochastic process
******************************************************************************/
void solveOneWeek(string scenPath, string demandPath, string historyPath, string customInputFile,
	string solPath, double timeout) {

  unsigned found = solPath.find_last_of(".");
	string logPathIni = solPath.substr(0,found),
    logPath = logPathIni+"Log.txt";
	Tools::LogOutput logStream(logPath);


	// set the scenario
	logStream << "# Initialize the scenario" << std::endl;
	Scenario* pScen = initializeScenario(scenPath,demandPath,historyPath);

	// set the options of the stochastic solver
	// (the corresponding method needs to be change manually for tests)
	logStream << "# Set the options" << std::endl;
	StochasticSolverOptions options;
	setStochasticSolverOptions(options, pScen, solPath, logPathIni, timeout);

	// check if a parameter file in the the in the log path ini
  string pathIni = "";
  found = solPath.find_last_of("/");
  if(found != string::npos)
    pathIni = solPath.substr(0,found+1);
	string stochasticOptions = pathIni+"solverOptions.txt",
					generationOptions = pathIni+"generationOptions.txt",
					evaluationOptions = pathIni+"evaluationOptions.txt";
	try {
		logStream << "Stochastic options:" << endl <<
      ReadWrite::readStochasticSolverOptions(stochasticOptions, options) << endl;
	} catch(const std::string& ex) {}
	try {
    logStream << "Generation options:" << endl <<
      ReadWrite::readSolverOptions(generationOptions, options.generationParameters_) << endl;
	} catch(const std::string& ex) {}
	try {
    logStream << "Evaluation options:" << endl <<
      ReadWrite::readSolverOptions(evaluationOptions, options.evaluationParameters_) << endl;
	} catch(const std::string& ex) {}

	// get history demands by reading the custom file
	//
	vector<Demand*> demandHistory;
	demandHistory.push_back(new Demand (*(pScen->pWeekDemand())) );
	if (!customInputFile.empty()) {
		// int coWeek = ReadWrite::readCustom(customInputFile, pScen, demandHistory);
	 ReadWrite::readCustom(customInputFile, pScen, demandHistory);
	}

	Solver* pSolver = new StochasticSolver(pScen, options, demandHistory);

	logStream << "# Solve the week" << std::endl;
	pSolver->solve();
	int solutionStatus = pSolver->getStatus();
	logStream << "# Solution status = " << solutionStatus <<  std::endl;

	Tools::LogOutput solStream(solPath);
	solStream << pSolver->solutionToString() << std::endl;

	//  release memory
	if (pScen) delete pScen;
	if (pSolver) delete pSolver;
	while (!demandHistory.empty()) {
		delete demandHistory.back();
		demandHistory.pop_back();
	}
	logStream.close();
}

/******************************************************************************
* Test a solution on multiple weeks
* In this method, the weeks are solved sequentially without knowledge of future
* demand
******************************************************************************/

pair<double, int> testMultipleWeeksStochastic(string dataDir, string instanceName, int historyIndex,
		vector<int> weekIndices, StochasticSolverOptions stochasticSolverOptions, string outdir, int seed) {

	// build the paths of the input files
	InputPaths inputPaths(dataDir, instanceName,historyIndex,weekIndices);

	// initialize the scenario object of the first week
	Scenario* pScen = initializeScenario(inputPaths.scenario(),inputPaths.week(0),inputPaths.history(),"");

	// solve the problem for each week and store the solution in the vector below
	vector<Roster> solution;
	int nbWeeks = weekIndices.size();
	Status solutionStatus;

	vector<Demand*> demandHistory;
	double currentCost = 0;
	int nbSched = 0;

	vector<int> seeds;

	for (int week = 0; week < nbWeeks; week++) {
	   if(week==0 and seed != -1)
	      seeds.push_back(seed);
	   else
	      seeds.push_back(std::rand());
	   std::srand(seeds[week]);

		demandHistory.push_back(new Demand (*(pScen->pWeekDemand())) );

		StochasticSolver* pSolver = new StochasticSolver(pScen, stochasticSolverOptions, demandHistory, currentCost);

		currentCost += pSolver->solve();
		nbSched += pSolver->getNbSchedules();
		printf( "Current cost = %.2f \n", currentCost);
		solutionStatus = pSolver->getStatus();
		if (solutionStatus == INFEASIBLE) {
			delete pSolver;
			break;
		}

		// update the overall solution with the solution of the week that was just
		// treated
		// warning: we must make sure that the main solution in pSolver applies only
		// to the demand of one week
		vector<Roster> weekSolution = pSolver->getSolutionAtDay(6);
		if (solution.empty()) solution = weekSolution;
		else {
			for (int n = 0; n < pScen->nbNurses_; n++) {
				solution[n].push_back(weekSolution[n]);
			}
		}

		std::cout << pSolver->solutionToLogString() << std::endl;

		// prepare the scenario for next week if we did not reach the last week yet
		if (week < nbWeeks-1) {

			// Initialize demand and preferences
			Demand* pDemand(0);
			Preferences* pPref(0);

			// Read the demand and preferences and link them with the scenario
			ReadWrite::readWeek(inputPaths.week(week+1), pScen, &pDemand, &pPref);

			// read the initial state of the new week from the last state of the
			// last week
			// modify the dayId_ to show that this is the first day of the new week
			vector<State> initialStates = pSolver->getStatesOfDay(6);
			for (int i = 0; i < pScen->nbNurses_; i++) {
				initialStates[i].dayId_ = 0;
			}

			// update the scenario to treat next week
			pScen->updateNewWeek(pDemand, *pPref, initialStates);
		}

		delete pSolver;
	}

	printf( "Total cost = %.2f \n", currentCost);

   string seedOutfile = outdir+"seeds.txt";
   Tools::LogOutput seedStream(seedOutfile, true);
   char str[50];
   sprintf(str, "Cost %.2f; NbGene %d; Seeds", currentCost, nbSched);
   seedStream << str;
   for(int s: seeds)
   seedStream << " " << s;
   seedStream << std::endl;

	// Display the solution
	// displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices, solution, solutionStatus, outDir);

	delete pScen;

	return pair<double, int>(currentCost, nbSched);
}


/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char** argv)
{
  std::cout << "Number of arguments= " << argc << std::endl;

  // Detect errors in the number of arguments
  //
  if (argc%2 != 1) {
    Tools::throwError("main: There should be an even number of arguments!");
  }
  else if (argc > 1 && (argc < 9 || argc > 17)) {
    Tools::throwError("main: There is either too many or not enough arguments!");
  }

  // Simulate default behavior for a test instance
  //
  if (argc == 1) {
     std::cout << "Running the default method..."<< std::endl;
     string dataDir = "datasets/";
     string instanceName = "n030w4";
     int historyIndex = 0;
     vector<int> weekIndices = {3, 5, 0, 2};
     StochasticSolverOptions stochasticSolverOptions;
     stochasticSolverOptions.totalTimeLimitSeconds_ = 40;
     string outdir = "outfiles/" + instanceName + "/";
     int seed = 1;
     testMultipleWeeksStochastic(dataDir, instanceName, historyIndex,
     		weekIndices, stochasticSolverOptions, outdir, seed);
  }
  // Nominal behavior of the executable, as required by INRCII
  //
  else {
    // Retrieve the file names in arguments
    //
    int narg = 1;
    string scenarioFile="", initialHistoryFile="", weekDataFile="", solutionFile="";
    string customInputFile="", customOutputFile="";
    int randSeed=0;
    double timeout =0.0;

    while (narg < argc) {
     std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;
     // Attention usine ������ gaz: problem with the java simulator that add quote
     // marks in the arguments, which fucks up the open file methods
     // the code below is here to remove these quote marks
     //
     string str(argv[narg+1]);
     std::size_t found = str.find("\"");
     while (found!=std::string::npos) {
        str.erase(found,1);
        found = str.find("\"");
     }

     if (!strcmp(argv[narg],"--sce")) {
        scenarioFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--his")) {
        initialHistoryFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--week")) {
        weekDataFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--sol")) {
        solutionFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--cusIn")) {
        customInputFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--cusOut")) {
        customOutputFile = str;
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--timeout")) {
        timeout = std::stod(str);
        narg += 2;
     }
     else if (!strcmp(argv[narg],"--rand")) {
        randSeed = std::stoi(str);
        narg += 2;
     }
     else {
        Tools::throwError("main: the argument does not match the expected list!");
     }
   }

    // Throw an error if a necessary input file is missing
    if ( scenarioFile.empty() || initialHistoryFile.empty()
          || weekDataFile.empty() || solutionFile.empty() ) {
       throw Tools::myException("A necessary file name is missing!",__LINE__);
    }

    srand(randSeed);

    // Solve the week
    solveOneWeek(scenarioFile, weekDataFile, initialHistoryFile, customInputFile, solutionFile, timeout);

    // Write the solution in the required output format
    //
    if (!customOutputFile.empty()) {
       ReadWrite::writeCustom(customOutputFile,weekDataFile,customInputFile);
    }
    cout << "Custom output file : " << customOutputFile << endl;
    // Todo: the method that writes the history file corresponding to the
    // solution
    // string outputHistoryFile("history-week");
    // outputHistoryFile += std::to_string(pScen->thisWeek()) + ".txt";
    // std::cout << "Output history file: " << outputHistoryFile << std::endl;
  }
}
