//
//  main.cpp
//  RosterDesNurses
//
//  Created by Jérémy Omer on 16/11/2014.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//

#include "main_test.h"
#include "MyTools.h"
<<<<<<< HEAD
#include "ReadWrite.h"
#include "Solver.h"
#include "Greedy.h"
#include <exception>
=======
#include "Vrp.h"

// Function for testing parts of the code (Antoine)
void testFunction_Antoine(){

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

   // Create a log file
   string logFile = "logs/test.log";
   Tools::LogOutput logStream(logFile);

   //Create an output file
   string outFile = "outfiles/test.out";
   Tools::LogOutput outStream(outFile);

   string inst = "n100w4";//n030w4

   string scenarPath = "datasets/" + inst + "/Sc-" + inst + ".txt";
   string firstWeekPath = "datasets/" + inst + "/WD-" + inst + "-1.txt";
   string firstHistoryPath = "datasets/" + inst + "/H0-" + inst + "-0.txt";

   // Read the input data from files
   Scenario* pScen = ReadWrite::readScenario(scenarPath);
   Demand* pWeekDemand = ReadWrite::readWeek(firstWeekPath, pScen);
   ReadWrite::readHistory(firstHistoryPath,pScen);

   // Check that the scenario was read properly
   //
   // logStream << *pScen << std::endl;
   logStream << pScen->toString() << std::endl;
   logStream << pWeekDemand->toString(true) << std::endl;

   // Write the aggregate information on the demand
   //


   // Write aggregate information on the cover capacity of the staff
   // (TBD)

   //Compute initial solution
   //
   Greedy* pGreedy =
      new Greedy(pScen, pWeekDemand,   pScen->pWeekPreferences(), pScen->pInitialState());
   pGreedy->constructiveGreedy();
   outStream << pGreedy->solutionToString();

   // Instantiate the solver class as a test
   //
   MasterProblem* pSolverTest =
      new MasterProblem(pScen, pWeekDemand,   pScen->pWeekPreferences(), pScen->pInitialState(), pGreedy->getSolution());
   pSolverTest->solve();

   // Write the solution in an output file
   outStream << pSolverTest->solutionToString();

   // Display the total time spent in the algorithm
   timertotal->stop();
   logStream.print("Total time spent in the algorithm : ");
   logStream.print(timertotal->dSinceInit());
   logStream.print("\n");


   //test vrp example of scip
   //   string dataFile = "datasets/vrp/eil22.vrp";
   //   Vrp* vrp = new Vrp(dataFile);

   // free the allocated pointers
   //
   //   delete vrp;
   delete timertotal;
   delete pWeekDemand;
   delete pScen;
   delete pGreedy;
   delete pSolverTest;
}

// Function for testing parts of the code (Jeremy)
void testFunction_Jeremy(){

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

   // Create a log file
   string logFile = "../logfiles/test.log";
   Tools::LogOutput logStream(logFile);


   // Read the input data from files
   Scenario* pScen = ReadWrite::readScenario("../datasets/n030w4/Sc-n030w4.txt");
   Demand* pWeekDemand = ReadWrite::readWeek("../datasets/n030w4/WD-n030w4-1.txt", pScen);
   ReadWrite::readHistory("../datasets/n030w4/H0-n030w4-0.txt",pScen);

   // Check that the scenario was read properly
   //
   // logStream << *pScen << std::endl;
   logStream << pScen->toString() << std::endl;
   logStream << pWeekDemand->toString(true) << std::endl;

   // Write the aggregate information on the demand
   //


   // Write aggregate information on the cover capacity of the staff
   // (TBD)

   // Instantiate the solver class as a test
   //
   Greedy* pSolverTest =
      new Greedy(pScen, pWeekDemand,	pScen->pWeekPreferences(), pScen->pInitialState());
   pSolverTest->constructiveGreedy();

   // Write the solution in an output file
   string outFile = "../outfiles/test.out";
   Tools::LogOutput outStream(outFile);
   outStream << pSolverTest->solutionToString();

   // Display the total time spent in the algorithm
   timertotal->stop();
   logStream.print("Total time spent in the algorithm : ");
   logStream.print(timertotal->dSinceInit());
   logStream.print("\n");

   // free the allocated pointers
   //
   delete timertotal;
   delete pWeekDemand;
   delete pScen;
   delete pSolverTest;


}

// Function for testing parts of the code (Samuel)
void testFunction_Samuel(){

   Tools::Timer* timertest = new Tools::Timer();
   timertest->init();
   timertest->start();

   string inst = "n100w4";

   string scenar = "/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/" + inst + "/Sc-" + inst + ".txt";
   string firstWeek = "/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/" + inst + "/WD-" + inst + "-1.txt";
   string firstHistory = "/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/" + inst + "/H0-" + inst + "-0.txt";

   Scenario * s = ReadWrite::readScenario(scenar);
   Demand* pWeekDemand = ReadWrite::readWeek(firstWeek,s);
   ReadWrite::readHistory(firstHistory,s);

   cout << s->toString() << endl;

   timertest->stop();

   string logFile = "../logfiles/samuel_test.log";
   Tools::LogOutput logStream(logFile);

   // RqJO : attention, j'ai enlevé ta surcharge de << parce qu'elle me faisait
   // des segfaults
   logStream << s->toString() << std::endl;
   logStream.print("Total time spent in the algorithm : ");
   logStream.print(timertest->dSinceInit());
   logStream.print("\n");

   for(map<string,Contract*>::const_iterator it = s->contracts_.begin(); it != s->contracts_.end(); ++it){
      Contract * c = it->second;
      cout << "# " << endl;
      cout << "# " << endl;
      cout << "# +----------------------------------------------------------------------------------------" << endl;
      cout << "# CONTRACT : " << c->toString() << endl;
      cout << "# +----------------------------------------------------------------------------------------" << endl;
      SubProblem sp (s, c);
      cout << "# +----------------------------------------------------------------------------------------" << endl;
      cout << "# " << endl;
      cout << "# " << endl;
   }

   //SubProblem sp;
   //sp.testGraph_spprc();


}
>>>>>>> refs/remotes/origin/AL-SCIP.0

int main(int argc, char** argv)
{
	std::cout << "Number of arguments= " << argc << std::endl;

<<<<<<< HEAD
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
			// Attention usine à gaz: problem with the java simulator that add quote
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
			else if (!strcmp(argv[narg],"--rand")) {
				randSeed = str;
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

		// Read the input files
		//
		Scenario* pScen = ReadWrite::readScenario(scenarioFile);
		std::cout << "Scenario file is read\n";
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
=======
   // Tests functions to check the functions one by one
   testFunction_Antoine();
   // testFunction_Jeremy();
   //	testFunction_Samuel();
>>>>>>> refs/remotes/origin/AL-SCIP.0

}
