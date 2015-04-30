//
//  main_test.cpp
//  RosterDesNurses
//
//  Created by Jeremy Omer on 04/03/2015.
//  Copyright (c) 2015 Jeremy Omer. All rights reserved.
//

#include "main_test.h"
#include "ReadWrite.h"
#include "DemandGenerator.h"
#include "Greedy.h"
#include "MasterProblem.h"
#include "StochasticSolver.h"
#include "SubProblem.h"
#include "CbcModeler.h"
#include "MyTools.h"


void main_test()
{
	testFunction_Antoine();
	//testFunction_Jeremy();
	//testFunction_Samuel();
}

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

   string data = "testdatasets/";// testdatasets datasets
   const char* inst = "n005w4";// n100w4 n030w4 n005w4

   string scenarPath = data + inst + "/Sc-" + inst + ".txt";
   //n005w4: {1, 2, 3, 3}
   //n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
   vector<int> numberWeek = {1, 2, 3, 3};
   vector<string> weekPaths(numberWeek.size());
   for(int i=0; i<numberWeek.size(); ++i){
      string path = data + inst + "/WD-" + inst + "-"+std::to_string(numberWeek[i])+".txt";
      weekPaths[i] = path;
   }
   string firstHistoryPath = data + inst + "/H0-" + inst + "-0.txt";

   // Read the input data from files
   Scenario* pScen = ReadWrite::readScenario(scenarPath);
   Preferences preferences;
   Demand* pWeekDemand = ReadWrite::readWeeks(weekPaths, pScen);
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
   outStream << pGreedy->solutionToLogString();

   // Instantiate the solver class as a test
   //
   MasterProblem* pBCP =
      new MasterProblem(pScen, pWeekDemand,   pScen->pWeekPreferences(), pScen->pInitialState(), S_BCP, pGreedy->getSolution());
   pBCP->solve();

   // Write the solution in the required output format
   vector<string> solutions = pBCP->solutionToString(pScen->nbWeeks());
   for(int w=0; w<pScen->nbWeeks(); ++w){
      int thisWeek = w+pScen->thisWeek();
      char solutionFile[30];
      snprintf ( solutionFile, 30, "outfiles/Sol-%s-%d-%d.txt", inst, numberWeek[w], thisWeek );
      Tools::LogOutput solutionStream(solutionFile);
      solutionStream << solutions[w];
   }

   // Write the solution in an output file
   outStream << pBCP->solutionToLogString();

   // Display the total time spent in the algorithm
   timertotal->stop();
   logStream.print("Total time spent in the algorithm : ");
   logStream.print(timertotal->dSinceInit());
   logStream.print("\n");


   // free the allocated pointers
   //
   //   delete vrp;
   delete timertotal;
   delete pWeekDemand;
   delete pScen;
   delete pGreedy;
   delete pBCP;
}

// Function for testing parts of the code (Jeremy)
void testFunction_Jeremy(){

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

	/************************************************************************
	* Go through the demands of the directory to find invariants in the demand
	*************************************************************************/

	// ReadWrite::compareDemands("testdatasets/n005w4","outfiles/comparedemands_n005w4.log");

	/************************************************************************
	* Initialize the week scenario by reading the input files
	*************************************************************************/

	Scenario* pScen(0);
	pScen = initializeScenario("datasets/n030w4/Sc-n030w4.txt",
		"datasets/n030w4/WD-n030w4-1.txt", "datasets/n030w4/H0-n030w4-0.txt","outfiles/inputdata.log");

	/************************************************************************
	* Test the random demand generator
	*************************************************************************/
	testRandomDemandGenerator(1,"outfiles/randomdemands.out",pScen);

	/****************************************
	* Run the greedy to get an initial solution
	*****************************************/

   Greedy* pGreedy =
      new Greedy(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
   pGreedy->constructiveGreedy();

	// Write the solution in the required output format
	string greedyFile = "outfiles/greedy.out";
	Tools::LogOutput greedyStream(greedyFile);
	greedyStream << pGreedy->solutionToString();

	// Write the solution and advanced information in a more convenient format
	string greedyLog = "outfiles/greedylog.out";
	Tools::LogOutput greedyLogStream(greedyLog);
	greedyLogStream << pGreedy->solutionToLogString();

	/****************************************
	* Test the CBC modeler
	*****************************************/
	std::vector<Roster> solIni = pGreedy->getSolution();
	testCbc(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState(),solIni);


   // Display the total time spent in the tests
	//
   timertotal->stop();
	 Tools::LogOutput logStream("outfiles/execution.log");
   logStream.print("Total time spent in the tests : ");
   logStream.print(timertotal->dSinceInit());
   logStream.print("\n");

   // free the allocated pointers
   //
   delete timertotal;
   delete pScen;
   delete pGreedy;
}

// Function for testing parts of the code (Samuel)
void testFunction_Samuel(){

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

	   string data = "datasets/";// testdatasets datasets
	   const char* inst = "n030w4";// n100w4 n030w4 n005w4 n012w8

	   string scenarPath = data + inst + "/Sc-" + inst + ".txt";
	   //n005w4: {1, 2, 3, 3}
	   //n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
	   vector<int> numberWeek = {1, 2, 3, 3};
	   vector<string> weekPaths(numberWeek.size());
	   for(int i=0; i<numberWeek.size(); ++i){
	      string path = data + inst + "/WD-" + inst + "-"+std::to_string(numberWeek[i])+".txt";
	      weekPaths[i] = path;
	   }
	   string firstHistoryPath = data + inst + "/H0-" + inst + "-0.txt";

	   // Read the input data from files
	   Scenario* pScen = ReadWrite::readScenario(scenarPath);
	   cout << pScen->toString() << endl;
	   Preferences preferences;
	   Demand* pWeekDemand = ReadWrite::readWeeks(weekPaths, pScen);
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
	   outStream << pGreedy->solutionToLogString();

	   // Instantiate the solver class as a test
	   //
	   MasterProblem* pBCP =
	      new MasterProblem(pScen, pWeekDemand,   pScen->pWeekPreferences(), pScen->pInitialState(), S_BCP, pGreedy->getSolution());
	   pBCP->solve();

	   // Write the solution in the required output format
	   vector<string> solutions = pBCP->solutionToString(pScen->nbWeeks());
	   for(int w=0; w<pScen->nbWeeks(); ++w){
	      int thisWeek = w+pScen->thisWeek();
	      char solutionFile[30];
	      snprintf ( solutionFile, 30, "outfiles/Sol-%s-%d-%d.txt", inst, numberWeek[w], thisWeek );
	      Tools::LogOutput solutionStream(solutionFile);
	      solutionStream << solutions[w];
	   }

	   // Write the solution in an output file
	   outStream << pBCP->solutionToLogString();

	   // Display the total time spent in the algorithm
	   timertotal->stop();
	   logStream.print("Total time spent in the algorithm : ");
	   logStream.print(timertotal->dSinceInit());
	   logStream.print("\n");


	   // free the allocated pointers
	   //
	   //   delete vrp;
	   delete timertotal;
	   delete pWeekDemand;
	   delete pScen;
	   delete pGreedy;
	   delete pBCP;
}

/************************************************************************
* Initialize the week scenario by reading the input files
*************************************************************************/

Scenario* initializeScenario(string scenFile, string demandFile, string historyFile, string logFile) {

	// Create a log file
	Tools::LogOutput logStream(logFile);

	// Initialize demand and preferences
	Demand* pDemand(0);
	Preferences* pPref(0);

	// Read the scenario
	Scenario* pScen = ReadWrite::readScenario(scenFile);

	// Read the demand and preferences and link them with the scenario
	ReadWrite::readWeek(demandFile, pScen,&pDemand,&pPref);
	pScen->linkWithDemand(pDemand);
	pScen->linkWithPreferences(*pPref);

	// Read the history
	ReadWrite::readHistory(historyFile, pScen);

	// Check that the scenario was read properly if logfile specified in input
	logStream << pScen->toString() << std::endl;
	logStream << pScen->pWeekDemand()->toString(true) << std::endl;


	return pScen;
}

/****************************************
* Test the CBC modeler
*****************************************/
void testCbc(Scenario* pScen, Demand* pDemand, Preferences* pPref,
	std::vector<State>* pStateIni, std::vector<Roster>& solIni) {
		// First method: directly instantiate a master problem equipped with a Cbc
		// modeler as input
		//
		MasterProblem* pMPCbc;
		pMPCbc = new MasterProblem(pScen, pDemand, pPref, pStateIni, S_CBC, solIni);
		pMPCbc->solve();

		// Write the solution in the required output format
		string outFile = "outfiles/cbctest1.out";
		Tools::LogOutput outStream(outFile);
		outStream << pMPCbc->solutionToString();

		// Second method, load the Cbc modeler from the model of the MP
		//
		CoinModeler* coinModel = (CoinModeler*) pMPCbc->getModel();
		CbcModeler* cbcModel =
			new CbcModeler(coinModel->getCoreVars(),coinModel->getColumns(),coinModel->getCons());
		cbcModel->solve();

		// a new method is needed to get the solution in the proper format from this
		// external Cbc model
}

/************************************************************************
* Test the random demand generator
*************************************************************************/

void testRandomDemandGenerator(int nbDemands,string logFile, Scenario* pScen) {

	Tools::LogOutput logStream(logFile);

	vector<Demand*> demandHistory;
	demandHistory.push_back(pScen->pWeekDemand());
	DemandGenerator generator(nbDemands,demandHistory,pScen);
	vector<Demand*> randomDemands = generator.generatePerturbedDemands();

	while (!randomDemands.empty()) {
		logStream << randomDemands.back()->toString(true) << std::endl;
		if (randomDemands.back()) delete randomDemands.back();
		randomDemands.pop_back();
	}
}
