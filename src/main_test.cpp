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

// some include files to go through the files of an input directory
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

// Function for testing parts of the code (Antoine)
void testFunction_Antoine(){

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

   string data = "datasets/";// testdatasets datasets
   const char* inst = "n120w8";// n100w4 n030w4 n005w4

   string scenarPath = data + inst + "/Sc-" + inst + ".txt";
   //n005w4: {1, 2, 3, 3}
   //n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
   //n021w4:
   //n120w8: {3, 2}
   vector<int> numberWeek = {3, 5, 0, 2, 0, 4, 5, 2};


//   testMultipleWeeksDeterministic(data, inst, 0, numberWeek, GENCOL, "outfiles/");
   testMultipleWeeksStochastic(data, inst, 0, numberWeek, STOCHASTIC_GENCOL, "outfiles/");

   // Display the total time spent in the algorithm
   timertotal->stop();
   cout << "Total time spent in the algorithm : " << timertotal->dSinceInit() << endl;


   // free the allocated pointers
   //
   //   delete vrp;
   delete timertotal;
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
	// computeStatsOnTheDemandsOfAllInstances("testdatasets/");
	// computeStatsOnTheDemandsOfAllInstances("datasets/");

	/***************************************************************************
	* Test the solution of only one week
	* Here, we don't want to compare deterministic to stochastic, only greedy to
	* column generation
	****************************************************************************/
	vector<int> weekIndex = {1};
	testMultipleWeeksDeterministic(dataDir, instanceName, 0, weekIndex, STOCHASTIC_GREEDY,
		(string)(outDir+"GreedyStochastic/"));



	/******************************************************************************
	* Test a solution on multiple weeks
	* In this method, the weeks are solved sequentially without knowledge of future
	* demand
	******************************************************************************/
	vector<int> weekIndices = {1, 2, 3, 3};
	//n005w4: {1, 2, 3, 3}
	//n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
  testMultipleWeeksDeterministic(dataDir, instanceName, 0, weekIndices, GENCOL,
    (string)(outDir+"BCP/"));

	testMultipleWeeksStochastic(dataDir, instanceName, 0, weekIndices, GREEDY, (string)(outDir+"GreedyStochastic/"));
	testMultipleWeeksDeterministic(dataDir, instanceName, 0, weekIndices, GREEDY,  (string)(outDir+"Greedy/"));

   /************************************************************************
   * Test the random demand generator
   *************************************************************************/
  //  Scenario* pScen(0);
  //  pScen = initializeScenario("datasets/n030w4/Sc-n030w4.txt",
  //     "datasets/n030w4/WD-n030w4-1.txt", "datasets/n030w4/H0-n030w4-0.txt","outfiles/inputdata.log");
  //  testRandomDemandGenerator(1,"outfiles/randomdemands.out",pScen);
  // delete pScen;

   /****************************************
   * Test the CBC modeler
   *****************************************/
  //  testCbc(pScen);


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
}

// Function for testing parts of the code (Samuel)
void testFunction_Samuel(){

	// Time the complete execution of the algorithm
	Tools::Timer* timertotal = new Tools::Timer();
	timertotal->init();
	timertotal->start();


	string data = "datasets/";// testdatasets datasets userdatasets
	const char* inst = "n030w4";// n100w4 n030w4 n005w4 n005w1

	string scenarPath = data + inst + "/Sc-" + inst + ".txt";
	//n005w4: {1, 2, 3, 3}
	//n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
	//n021w4:
	//n120w8: {3, 2}
	vector<int> numberWeek = {4,6};


	testMultipleWeeksDeterministic(data, inst, 0, numberWeek, GENCOL, "outfiles/");

	// Display the total time spent in the algorithm
	timertotal->stop();
	cout << "Total time spent in the algorithm : " << timertotal->dSinceInit() << endl;


	// free the allocated pointers
	//
	//   delete vrp;
	delete timertotal;
}


/******************************************************************************
* The instances of InputPaths contain the paths of the input files of the
* problem
*******************************************************************************/

InputPaths::InputPaths(string dataDir, string instanceName,int historyIndex, vector<int> weekIndices) {

	int nbWeeks = weekIndices.size();
	string instanceDir = dataDir + instanceName + "/";

	// initialize the scenario and history file names
	scenario_ = instanceDir + "Sc-" + instanceName + ".txt";
	history_ = instanceDir + "H0" + "-" + instanceName + "-" + std::to_string(historyIndex) + ".txt";

	// initialize the file names for each week demand
	for(int week: weekIndices){
		string path = instanceDir + "WD-" + instanceName + "-" + std::to_string(week) + ".txt";
		weeks_.push_back(path);
	}
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

	// build the paths of the input files
	InputPaths inputPaths(dataDir, instanceName,historyIndex,weekIndices);

	// Read the scenario
	Scenario* pScen = ReadWrite::readScenario(inputPaths.scenario());

	// Read the demand and preferences and link them with the scenario
	ReadWrite::readWeeks(inputPaths.weeks(), pScen);

	// Read the history
	ReadWrite::readHistory(inputPaths.history(), pScen);

	// Check that the scenario was read properly if logfile specified in input
	if (!logPath.empty()) {
		Tools::LogOutput logStream(logPath);
		logStream << pScen->toString() << std::endl;
		logStream << pScen->pWeekDemand()->toString(true) << std::endl;
	}

	return pScen;
}

/******************************************************************************
* Solve a deterministic input demand with the input algorithm
* In this method, we assume that all the demands are knwon in advance
* (the method can also treat only one week)
******************************************************************************/

void testMultipleWeeksDeterministic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outDir) {

	Scenario* pScen = initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);

	Solver* pSolver = setSolverWithInputAlgorithm(pScen, algorithm);
	pSolver->evaluate();

	// Display the solution
	vector<Roster> solution = pSolver->getSolution();
  Status status = pSolver->getStatus();
	displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices, solution,status, outDir);

	delete pSolver;
	delete pScen;
}

/******************************************************************************
* Test a solution on multiple weeks
* In this method, the weeks are solved sequentially without knowledge of future
* demand
******************************************************************************/

void testMultipleWeeksStochastic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outDir) {

	// build the paths of the input files
	InputPaths inputPaths(dataDir, instanceName,historyIndex,weekIndices);

	// initialize the scenario object of the first week
	Scenario* pScen = initializeScenario(inputPaths.scenario(),inputPaths.week(0),inputPaths.history(),"");

	// solve the problem for each week and store the solution in the vector below
	vector<Roster> solution;
	int nbWeeks = weekIndices.size();
	Status solutionStatus;
	for (int week = 0; week < nbWeeks; week++) {

		Solver* pSolver = setSolverWithInputAlgorithm(pScen, algorithm);

		pSolver->solve();
		solutionStatus = pSolver->getStatus();
		if (solutionStatus == INFEASIBLE) {
			delete pSolver;
			break;
		}

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
			ReadWrite::readWeek(inputPaths.week(week+1), pScen, &pDemand, &pPref);

			// read the initial state of the new week from the last state of the
			// last week
			// modify the dayId_ to show that this is the first day of the new week
			vector<State> initialStates = pSolver->getFinalStates();
			for (int i = 0; i < pScen->nbNurses_; i++) {
				initialStates[i].dayId_ = 0;
			}

			// update the scenario to treat next week
			pScen->updateNewWeek(pDemand, *pPref, initialStates);
		}

		delete pSolver;
	}

	// Display the solution
	displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices, solution, solutionStatus, outDir);

	delete pScen;
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
  pMPCbc = new MasterProblem(pScen, pDemand, pPref, pStateIni, S_CBC);
  pMPCbc->solve(pGreedy->getSolution());

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
   DemandGenerator generator(nbDemands,12,demandHistory,pScen);
   vector<Demand*> randomDemands = generator.generatePerturbedDemands();

   while (!randomDemands.empty()) {
      logStream << randomDemands.back()->toString(true) << std::endl;
      if (randomDemands.back()) delete randomDemands.back();
      randomDemands.pop_back();
   }
}


/*****************************************************************************
* Create a solver of the class specified by the input algorithm type
******************************************************************************/

Solver* setSolverWithInputAlgorithm(Scenario* pScen, Algorithm algorithm) {

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
	return pSolver;
}

/*******************************************************************************
* Create a stochastic solver of the class specified by the input algorithm type
* (in a way, this function is not necessary, but it fits well with the rest)
********************************************************************************/

Solver* setStochasticSolverWithInputAlgorithm(Scenario* pScen, Algorithm generationAlgorithm, Algorithm evaluationAlgorithm,
		int nExtraDaysGenerationDemands, int nEvaluationDemands, int nDaysEvaluation, int nGenerationDemands, vector<Demand*> demandHistory) {
	Solver* pSolver;
	pSolver = new StochasticSolver(pScen, generationAlgorithm, evaluationAlgorithm, nExtraDaysGenerationDemands, nEvaluationDemands,
			nDaysEvaluation, nGenerationDemands, demandHistory);
	return pSolver;
}

/******************************************************************************
* When a solution of multiple consecutive weeks is available, load it in a
* solver for all the weeks and  display the results
******************************************************************************/
void displaySolutionMultipleWeeks(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, vector<Roster> &solution, Status status, string outDir) {

	if (outDir.empty()) return;

  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = weekIndices.size();
  string catWeeks;
  for (int w=0; w< nbWeeks; w++) catWeeks += std::to_string(weekIndices[w]);
  string logPath = outDir+"Log-"+catWeeks+".txt";
  Tools::LogOutput outStream(logPath);

  // treat the case where the solver was unable to find a feasible solution
  if (status == INFEASIBLE) {
    outStream << "The solver was not able to find a solution\n";
    return;
  }

	// load the solution in a new solver
	Scenario* pScen = initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);
	Solver* pSolver = new Solver(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
	pSolver->loadSolution(solution);

	// write the log file for all the weeks
	outStream << pSolver->solutionToLogString();

	// write separately the solutions of each week in the required output format
	vector<string> solutions = pSolver->solutionToString(nbWeeks);
	for(int w=0; w < nbWeeks; ++w){
		string solutionFile = outDir+"Sol-"+catWeeks+"-"+std::to_string(weekIndices[w])+"-"+std::to_string(w)+".txt";
		Tools::LogOutput solutionStream(solutionFile);
		solutionStream << solutions[w];
	}

	delete pSolver;
	delete pScen;
}

/******************************************************************************
* Compute and record stats on all the demand files of all the instances in the
* input directory
******************************************************************************/

void computeStatsOnTheDemandsOfAllInstances(string inputDir) {
	struct dirent *dirp;
	struct stat filestat;

	// Open the input directory
	DIR* dp = opendir( inputDir.c_str() );
	if (dp == NULL) {
		Tools::throwError("Error while opening ");
	}
	else{
		std::cout << "Reading from directory " << inputDir << std::endl;
	}
	while ((dirp = readdir( dp )))
	{
		std::string filename(dirp->d_name);

		// The instance names start with "WD"
		if (filename[0] != 'n') continue;
		ReadWrite::compareDemands((string) (inputDir+filename),(string) ("outfiles/statDemands/"+filename+".txt"));
	}

}
