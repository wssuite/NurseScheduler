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
<<<<<<< HEAD
	testFunction_Antoine();
	//testFunction_Jeremy();
=======
	//testFunction_Antoine();
	testFunction_Jeremy();
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses
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
//   while(!pGreedy->constructiveGreedy())
//      cout << "No solution found by the greedy." << endl;
//   outStream << pGreedy->solutionToLogString();

   // Instantiate the solver class as a test
   //
   MasterProblem* pBCP =
      new MasterProblem(pScen, pWeekDemand,   pScen->pWeekPreferences(), pScen->pInitialState(), S_BCP);//, pGreedy->getSolution());
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

	// Test directory and instance, the same is used in every test
	string dataDir = "testdatasets/";// testdatasets datasets
	string instanceName = "n005w4";// n100w4 n030w4 n005w4
	string outDir = "outfiles/" + instanceName + "/";

   /************************************************************************
   * Go through the demands of the directory to find invariants in the demand
   *************************************************************************/

   // ReadWrite::compareDemands("testdatasets/n005w4","outfiles/comparedemands_n005w4.log");

	/***************************************************************************
	* Test the solution of only one week
	* Here, we don't want to compare deterministic to stochastic, only greedy to
	* column generation
	****************************************************************************/
	vector<int> weekIndex = {1};
	testMultipleWeeksDeterministic(dataDir, instanceName, 0, weekIndex, GREEDY,
		(string)(outDir+"greedyW"+std::to_string(weekIndex[0])+".log"));


	/******************************************************************************
	* Test a solution on multiple weeks
	* In this method, the weeks are solved sequentially without knowledge of future
	* demand
	******************************************************************************/
	vector<int> weekIndices = {1, 2, 3, 3};
	//n005w4: {1, 2, 3, 3}
	//n012w8: {3, 5, 0, 2, 0, 4, 5, 2}

	testMultipleWeeksStochastic(dataDir, instanceName, 0, weekIndices, GREEDY, (string)(outDir+"greedyStochastic.log"));
	testMultipleWeeksDeterministic(dataDir, instanceName, 0, weekIndices, GREEDY,  (string)(outDir+"greedyDeterministic.log"));

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
   * Test the CBC modeler
   *****************************************/

   testCbc(pScen);


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

<<<<<<< HEAD
	   string data = "datasets/";// testdatasets datasets
	   const char* inst = "n030w4";// n100w4 n030w4 n005w4 n012w8
=======
      string data = "datasets/";// testdatasets datasets
      const char* inst = "n030w4";// n100w4 n030w4 n005w4
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses

<<<<<<< HEAD
      string scenarPath = data + inst + "/Sc-" + inst + ".txt";
      //n005w4: {1, 2, 3, 3}
      //n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
      vector<int> numberWeek = {1,2,5,0};
      vector<string> weekPaths(numberWeek.size());
      for(int i=0; i<numberWeek.size(); ++i){
         string path = data + inst + "/WD-" + inst + "-"+std::to_string(numberWeek[i])+".txt";
         weekPaths[i] = path;
      }
      string firstHistoryPath = data + inst + "/H0-" + inst + "-0.txt";
=======
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
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses

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

Scenario* initializeScenario(string scenPath, string demandPath, string historyPath, string logPath) {

	// Initialize demand and preferences
	Demand* pDemand(0);
	Preferences* pPref(0);

	// Read the scenario
	Scenario* pScen = ReadWrite::readScenario(scenPath);

	// Read the demand and preferences and link them with the scenario
	ReadWrite::readWeek(demandPath, pScen,&pDemand,&pPref);
	pScen->linkWithDemand(pDemand);
	pScen->linkWithPreferences(*pPref);

	// Read the history
	ReadWrite::readHistory(historyPath, pScen);

	// Check that the scenario was read properly if logfile specified in input
	if (!logPath.empty()) {
		Tools::LogOutput logStream(logPath);
		logStream << pScen->toString() << std::endl;
		logStream << pScen->pWeekDemand()->toString(true) << std::endl;
	}


	return pScen;
}

/*****************************************************************************
* Initialize the scenario for multiple weeks
* When calling this function, the intent is to solve all the weeks at once
******************************************************************************/

Scenario* initializeMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, string logPath) {

	int nbWeeks = weekIndices.size();
	string instanceDir = dataDir + instanceName + "/";

	// initialize the scenario and history file names
	string scenPath = instanceDir + "Sc-" + instanceName + ".txt";
	string historyPath = instanceDir + "H0" + "-" + instanceName + "-" + std::to_string(historyIndex) + ".txt";

	// initialize the file names for each week demand
	vector<string> weekPaths;
	for(int week: weekIndices){
		string path = instanceDir + "WD-" + instanceName + "-" + std::to_string(week) + ".txt";
		weekPaths.push_back(path);
	}

	// Read the scenario
	Scenario* pScen = ReadWrite::readScenario(scenPath);

	// Read the demand and preferences and link them with the scenario
	ReadWrite::readWeeks(weekPaths, pScen);

	// Read the history
	ReadWrite::readHistory(historyPath, pScen);

	// Check that the scenario was read properly if logfile specified in input
	if (!logPath.empty()) {
		Tools::LogOutput logStream(logPath);
		logStream << pScen->toString() << std::endl;
		logStream << pScen->pWeekDemand()->toString(true) << std::endl;
	}

	return pScen;
}

/****************************************
* Test the CBC modeler
*****************************************/
void testCbc(Scenario* pScen) {

	Demand* pDemand = pScen->pWeekDemand();
	Preferences* pPref = pScen->pWeekPreferences();
	vector<State>* pStateIni = pScen->pInitialState();

	// first, get an initial set of columns to input to Cbc by calling the greedy
	Greedy* pGreedy = new Greedy(pScen, pDemand, pPref, pStateIni);
	pGreedy->solve();

  // First method: directly instantiate a master problem equipped with a Cbc
  // modeler as input
  //
  MasterProblem* pMPCbc;
  pMPCbc = new MasterProblem(pScen, pDemand, pPref, pStateIni, S_CBC, pGreedy->getSolution());
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

/******************************************************************************
* Solve a deterministic input demand with the input algorithm
* In this method, we assume that all the demands are knwon in advance
* (the method can also treat only one week)
******************************************************************************/

void testMultipleWeeksDeterministic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outPath) {

	Scenario* pScen = initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);

	Solver* pSolver;
	switch(algorithm){
	case GREEDY:
		pSolver = new Greedy(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
		break;
	case GENCOL:
		pSolver = new MasterProblem(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState(), S_BCP);
		break;
	default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}

	pSolver->solve();

	// display the solution if outfile specified in input
	if (!outPath.empty()) {
		Tools::LogOutput outStream(outPath);
		outStream << pSolver->solutionToLogString();
	}

	delete pSolver;
	delete pScen;
}

/******************************************************************************
* Test a solution on multiple weeks
* In this method, the weeks are solved sequentially without knowledge of future
* demand
******************************************************************************/

void testMultipleWeeksStochastic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outPath) {

	int nbWeeks = weekIndices.size();
	string instanceDir = dataDir + instanceName + "/";

	// initialize the scenario and history file names
	string scenPath = instanceDir + "Sc-" + instanceName + ".txt";
	string historyPath = instanceDir + "H0" + "-" + instanceName + "-" + std::to_string(historyIndex) + ".txt";

	// initialize the file names for each week demand
	vector<string> weekPaths;
	for(int week: weekIndices){
		string path = instanceDir + "WD-" + instanceName + "-" + std::to_string(week) + ".txt";
		weekPaths.push_back(path);
	}

	// initialize the scenario object of the first week
	Scenario* pScen = initializeScenario(scenPath,weekPaths[0],historyPath,"");
	vector<Roster> solution;

	for (int week = 0; week < nbWeeks; week++) {
		Solver* pSolver;
		switch(algorithm){
		case GREEDY:
			pSolver = new Greedy(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
			break;
		case GENCOL:
			pSolver = new MasterProblem(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState(), S_BCP);
			break;
		default:
			Tools::throwError("The algorithm is not handled yet");
			break;
		}

		pSolver->solve();

		// update the overall solution with the solution of the week that was just
		// treated
		// warning: we must make sure that the main solution in pSolver applies only
		// to the demand of one week
		vector<Roster> weekSolution = pSolver->getSolution();
		if (solution.empty()) solution = weekSolution;
		else {
			for (int n = 0; n < pScen->nbNurses_; n++) {
				solution[n].push_back(weekSolution[n]);
			}
		}

		// prepare the scenario for next week if we did not reach the last week yet
		if (week < nbWeeks-1) {

			// Initialize demand and preferences
			Demand* pDemand(0);
			Preferences* pPref(0);

			// Read the demand and preferences and link them with the scenario
			ReadWrite::readWeek(weekPaths[week+1], pScen, &pDemand, &pPref);

			// read the initial state of the new week from the last state of the
			// last week
			vector<State> initialStates = pSolver->getFinalStates();

			// update the scenario to treat next week
			pScen->updateNewWeek(pDemand, *pPref, initialStates);
		}

		delete pSolver;
	}

	// Display the solution if outfile specified in input
	if (!outPath.empty()) {
		displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices, solution, outPath);
	}

	delete pScen;
}

/******************************************************************************
* When a solution of multiple consecutive weeks is available, load it in a
* solver for all the weeks and  display the results
******************************************************************************/
void displaySolutionMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, vector<Roster> &solution, string outPath) {

	// load the solution in a new solver
	Scenario* pScen = initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);
	Solver* pSolver = new Solver(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
	pSolver->loadSolution(solution);

	// display the solution if outfile specified in input
	if (!outPath.empty()) {
		Tools::LogOutput outStream(outPath);
		outStream << pSolver->solutionToLogString();
	}

	delete pSolver;
	delete pScen;
}
