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
//#include "CbcModeler.h"
#include "MyTools.h"

// some include files to go through the files of an input directory
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

//initialize the counter of object
unsigned int MyObject::s_count = 0;
unsigned int Rotation::s_count = 0;

// Function for testing parts of the code (Antoine)
void testFunction_Antoine(){

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

   string data = "datasets/";// testdatasets datasets
   const char* inst = "n030w4";// n100w4 n030w4 n005w4

   string scenarPath = data + inst + "/Sc-" + inst + ".txt";
   //n005w4: {1, 2, 3, 3}
   //n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
   //n021w4:
   //n120w8: {3, 2}
   vector<int> numberWeek = {3, 5, 0, 2}; // , 0, 4, 5, 2};


   SolverParam optParam;
   optParam.stopAfterXSolution_ = 2;
//   optParam.nbDiveIfMinGap_ = 0;
//   optParam.nbDiveIfRelGap_ = 0;
//   testMultipleWeeksDeterministic(data, inst, 0, numberWeek, GENCOL, "outfiles/MyTests/", optParam);
//   StochasticSolverOptions stochasticSolverOptions;
//   testMultipleWeeksStochastic(data, inst, 0, numberWeek, stochasticSolverOptions, "outfiles/MyTests/");

   Scenario* pScen = initializeMultipleWeeks(data, inst, 0, {0});

   Solver* pSolver = setSolverWithInputAlgorithm(pScen, GENCOL);
   pSolver->solve(optParam);

   for(int i=1; i<10; ++i){
      printf("**********************************************\n"
         "Demand %d"
         "\n**********************************************\n", i);

      Tools::Timer timer;
      timer.init();

      //build the corresponding scenario
      vector<int> demand = {i};
      InputPaths inputPaths(data, inst,0,demand);
      Scenario* pScen2 = initializeMultipleWeeks(data, inst, 0, demand);
      Solver* pSolver0 = setSolverWithInputAlgorithm(pScen2, GREEDY);
      pSolver0->solve();

      //read the new demand
      Demand* newDemand(0);
      Preferences* newPref(0);
      ReadWrite::readWeek(inputPaths.week(0), pScen, &newDemand, &newPref);
      delete newPref;
      //resolve
      timer.start();
      pSolver->resolve(newDemand, optParam,pSolver0->getSolution());
      timer.stop();
      cout << "Total time spent in the algorithm : " << timer.dSinceInit() << endl;

      //solve
      timer.init();
      Solver* pSolver2 = setSolverWithInputAlgorithm(pScen2, GENCOL);
      timer.start();
      pSolver2->solve(optParam, pSolver0->getSolution());
      timer.stop();
      cout << "Total time spent in the algorithm : " << timer.dSinceInit() << endl;

      delete pSolver0;
      delete pSolver2;
      delete pScen2;
   }

   delete pSolver;
   delete pScen;
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

//	testMultipleWeeksStochastic(dataDir, instanceName, 0, weekIndices, GREEDY, (string)(outDir+"GreedyStochastic/"));
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

	std::cout << "# Test Samuel" << endl;

	// Time the complete execution of the algorithm
	Tools::Timer* timertotal = new Tools::Timer();
	timertotal->init();
	timertotal->start();

	string data = "datasets/";// testdatasets datasets userdataset
	string inst = "n030w4";// n100w4 n030w4 n012w8 n005w4 n005w1

	int maxTimeAllowed = allowedTime(inst,SAM);

	string scenarPath = data + inst + "/Sc-" + inst + ".txt";
	//n005w4: {1, 2, 3, 3}
	//n012w8: {3, 5, 0, 2, 0, 4, 5, 2}
	//n021w4:
	//n120w8: {3, 2}
	vector<int> numberWeek = {6,2,9,1};
	int historyId = 1;

	string catWeek;
	for(int w: numberWeek) catWeek += std::to_string(w);

	StochasticSolverOptions stochasticSolverOptions;
	stochasticSolverOptions.withIterativeDemandIncrease_ = false;
	stochasticSolverOptions.withEvaluation_ = true;
	stochasticSolverOptions.generationCostPerturbation_ = true;
	stochasticSolverOptions.evaluationCostPerturbation_ = true;
	stochasticSolverOptions.withResolveForGeneration_ = false;
	stochasticSolverOptions.generationAlgorithm_ = GENCOL;
	stochasticSolverOptions.withResolveForEvaluation_ = true;
	stochasticSolverOptions.evaluationAlgorithm_ = GENCOL;
	stochasticSolverOptions.totalTimeLimitSeconds_ = maxTimeAllowed;
	stochasticSolverOptions.nExtraDaysGenerationDemands_ = 7;
	stochasticSolverOptions.nEvaluationDemands_ = 4;
	stochasticSolverOptions.nDaysEvaluation_ = 14;
	stochasticSolverOptions.nGenerationDemandsMax_ = 100;

	SolverParam generationParameters;
	generationParameters.maxSolvingTimeSeconds_ = 3000;
	generationParameters.printEverySolution_ = false;
	string of = "outfiles/" + inst + "-" + catWeek + "-sol-week";
	generationParameters.outfile_ = of;
//	generationParameters.logfile_ = generationParameters.outfile_;
	generationParameters.absoluteGap_ = 5;
	generationParameters.minRelativeGap_ = .05;
	generationParameters.relativeGap_ = .1;
	generationParameters.nbDiveIfMinGap_ = 1;
	generationParameters.nbDiveIfRelGap_ = 2;
	generationParameters.solveToOptimality_ = false;
	generationParameters.weightStrategy_ = RANDOMMEANMAX;

	stochasticSolverOptions.generationParameters_ = generationParameters;

	SolverParam evaluationParameters;
	evaluationParameters.maxSolvingTimeSeconds_ = 3000;
	evaluationParameters.printEverySolution_ = false;
//	evaluationParameters.outfile_ = "";
//	evaluationParameters.logfile_ = evaluationParameters.outfile_;
	evaluationParameters.absoluteGap_ = 5;
	evaluationParameters.minRelativeGap_ = .05;
	evaluationParameters.relativeGap_ = .1;
	evaluationParameters.nbDiveIfMinGap_ = 1;
	evaluationParameters.nbDiveIfRelGap_ = 2;
	evaluationParameters.solveToOptimality_ = false;
	evaluationParameters.weightStrategy_ = BOUNDRATIO;
	evaluationParameters.stopAfterXSolution_ = 0;

	stochasticSolverOptions.evaluationParameters_ = evaluationParameters;

	testMultipleWeeksStochastic(data, inst, historyId, numberWeek, stochasticSolverOptions, "outfiles/");


	// Display the total time spent in the algorithm
	timertotal->stop();
	cout << "Total time spent in the algorithm : " << timertotal->dSinceInit() << endl;


	// free the allocated pointers
	//
	//   delete vrp;
	delete timertotal;
}

/******************************************************************************
* Solve one week inside the stochastic process
******************************************************************************/
void solveOneWeek(string scenPath, string demandPath, string historyPath, string customInputFile, 
	string solPath, string logPathIni, double timeout) {

	string logPath = logPathIni+"LogStochastic.txt";
	Tools::LogOutput logStream(logPath);


	// set the scenario
	logStream << "# Initialize the scenario" << std::endl; 
	Scenario* pScen = initializeScenario(scenPath,demandPath,historyPath);

	// set the options of the stochastic solver 
	// (the corresponding method needs to be change manually for tests)
	logStream << "# Set the options" << std::endl; 
	StochasticSolverOptions options;
	setStochasticSolverOptions(options, pScen, solPath, logPathIni, timeout); 

	// get history demands by reading the custom file
	// 
	vector<Demand*> demandHistory;
	demandHistory.push_back(new Demand (*(pScen->pWeekDemand())) );
	if (!customInputFile.empty()) {
		int coWeek = ReadWrite::readCustom(customInputFile, pScen, demandHistory);
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
* Set the options of the stochastic solver
* This is not automated, so the options need to be changed inside the code 
* during the tests
* The solution time depends on the number of nurses and on the computed
******************************************************************************/

void setStochasticSolverOptions(StochasticSolverOptions& options, Scenario* pScenario, string solPath, string logPathIni, double timeout) {
   #ifdef __MACH__
   double cpuMaxFor30Nurses = 60.0;
   double cpuMaxPer10Nurses = 45.0;
   #else
   double cpuMaxFor30Nurses = 45.0;
   double cpuMaxPer10Nurses = 35.0;
   #endif

   string logStochastic = logPathIni.empty() ? "":logPathIni+"LogStochastic.txt";
   string logSolver = logPathIni.empty() ? "":logPathIni+"LogSolver.txt";

   options.withEvaluation_ = false;
   options.generationCostPerturbation_ = true;
   options.evaluationCostPerturbation_ = true;
   options.generationAlgorithm_ = GENCOL;
   options.evaluationAlgorithm_ = GENCOL;
   options.totalTimeLimitSeconds_ = timeout;
   options.nExtraDaysGenerationDemands_ = 7;
   options.nEvaluationDemands_ = 4;
   options.nDaysEvaluation_ = 14;
   options.nGenerationDemandsMax_ = 100;
   options.logfile_ = logStochastic;

   SolverParam generationParameters;
   generationParameters.maxSolvingTimeSeconds_ = options.totalTimeLimitSeconds_-1.0;
   generationParameters.printEverySolution_ = false;
   generationParameters.outfile_ = solPath;
   generationParameters.logfile_ = logSolver;
   generationParameters.absoluteGap_ = 5;
   generationParameters.minRelativeGap_ = .05;
   generationParameters.relativeGap_ = .1;
   generationParameters.nbDiveIfMinGap_ = 1;
   generationParameters.nbDiveIfRelGap_ = 2;
   generationParameters.solveToOptimality_ = false;

   options.generationParameters_ = generationParameters;

   SolverParam evaluationParameters;
   evaluationParameters.maxSolvingTimeSeconds_ = (options.totalTimeLimitSeconds_-1.0)/(2.0*options.nEvaluationDemands_);
   evaluationParameters.printEverySolution_ = false;
   evaluationParameters.outfile_ = "outdir/";
   evaluationParameters.logfile_ = logSolver;
   evaluationParameters.absoluteGap_ = 5;
   evaluationParameters.minRelativeGap_ = .05;
   evaluationParameters.relativeGap_ = .1;
   evaluationParameters.nbDiveIfMinGap_ = 1;
   evaluationParameters.nbDiveIfRelGap_ = 2;
   evaluationParameters.solveToOptimality_ = false;
   evaluationParameters.stopAfterXSolution_ = 0;

   options.evaluationParameters_ = evaluationParameters;
}

void setStochasticSolverOptions(StochasticSolverOptions& stochasticSolverOptions, Computer computer, string instanceName,
   string solPath, string logPathIni) {

	string logStochastic = logPathIni.empty() ? "":logPathIni+"LogStochastic.txt";
	string logSolver = logPathIni.empty() ? "":logPathIni+"LogSolver.txt";


   int maxTimeAllowed = allowedTime(instanceName,computer);


   stochasticSolverOptions.withIterativeDemandIncrease_ = false;
   stochasticSolverOptions.withEvaluation_ = true;
   stochasticSolverOptions.generationCostPerturbation_ = true;
   stochasticSolverOptions.evaluationCostPerturbation_ = true;
   stochasticSolverOptions.withResolveForGeneration_ = false;
   stochasticSolverOptions.generationAlgorithm_ = GENCOL;
   stochasticSolverOptions.withResolveForEvaluation_ = true;
   stochasticSolverOptions.evaluationAlgorithm_ = GENCOL;
   stochasticSolverOptions.rankingStrategy_ = RK_SCORE;
   stochasticSolverOptions.totalTimeLimitSeconds_ = maxTimeAllowed;
   stochasticSolverOptions.nExtraDaysGenerationDemands_ = 7;
   stochasticSolverOptions.nEvaluationDemands_ = 4;
   stochasticSolverOptions.nDaysEvaluation_ = 14;
   stochasticSolverOptions.nGenerationDemandsMax_ = 100;
   stochasticSolverOptions.logfile_ = logStochastic;

   SolverParam generationParameters;
   generationParameters.maxSolvingTimeSeconds_ = 3000;
   generationParameters.printEverySolution_ = false;
   generationParameters.outfile_ = solPath;
// generationParameters.logfile_ = generationParameters.outfile_;
   generationParameters.absoluteGap_ = 5;
   generationParameters.minRelativeGap_ = .05;
   generationParameters.relativeGap_ = .1;
   generationParameters.nbDiveIfMinGap_ = 1;
   generationParameters.nbDiveIfRelGap_ = 2;
   generationParameters.solveToOptimality_ = false;
   generationParameters.weightStrategy_ = RANDOMMEANMAX;

   stochasticSolverOptions.generationParameters_ = generationParameters;

   SolverParam evaluationParameters;
   evaluationParameters.maxSolvingTimeSeconds_ = 3000;
   evaluationParameters.printEverySolution_ = false;
// evaluationParameters.outfile_ = "";
   evaluationParameters.logfile_ = logSolver;
   evaluationParameters.absoluteGap_ = 5;
   evaluationParameters.minRelativeGap_ = .05;
   evaluationParameters.relativeGap_ = .1;
   evaluationParameters.nbDiveIfMinGap_ = 1;
   evaluationParameters.nbDiveIfRelGap_ = 2;
   evaluationParameters.solveToOptimality_ = false;
   evaluationParameters.weightStrategy_ = BOUNDRATIO;
   evaluationParameters.stopAfterXSolution_ = 0;

   stochasticSolverOptions.evaluationParameters_ = evaluationParameters;
}

void setStochasticSolverOptions(StochasticSolverOptions& stochasticSolverOptions, Computer computer, string instanceName,
   string solPath, string logPathIni, string stochasticOptionsFile, string generationOptionsFile, string evaluationOptionsFile) {

   string logStochastic = logPathIni.empty() ? "":logPathIni+"LogStochastic.txt";
   string logSolver = logPathIni.empty() ? "":logPathIni+"LogSolver.txt";

   ReadWrite::readStochasticSolverOptions(stochasticOptionsFile, stochasticSolverOptions);
   int maxTimeAllowed = allowedTime(instanceName,computer);
   stochasticSolverOptions.totalTimeLimitSeconds_ = maxTimeAllowed;
   stochasticSolverOptions.logfile_ = logStochastic;

   SolverParam generationParameters;
   ReadWrite::readSolverOptions(generationOptionsFile, generationParameters);
   generationParameters.outfile_ = solPath;
   stochasticSolverOptions.generationParameters_ = generationParameters;
   stochasticSolverOptions.generationParameters_.verbose_ = stochasticSolverOptions.verbose_;

   SolverParam evaluationParameters;
   ReadWrite::readSolverOptions(evaluationOptionsFile, evaluationParameters);
   evaluationParameters.logfile_ = logSolver;
   stochasticSolverOptions.evaluationParameters_ = evaluationParameters;
   stochasticSolverOptions.evaluationParameters_.verbose_ = stochasticSolverOptions.verbose_;
}


/******************************************************************************
* The instances of InputPaths contain the paths of the input files of the
* problem
*******************************************************************************/

int allowedTime(string instance, Computer computer){

	switch(computer){
	case SAM:
		if(instance.at(2) == '3')
			return 45;
		if(instance.at(2) == '4')
			return 79;
		if(instance.at(2) == '5')
			return 112;
		if(instance.at(2) == '6')
			return 146;
		if(instance.at(2) == '8')
			return 212;
		if(instance.at(2) == '0')
			return 279;
		if(instance.at(2) == '2')
			return 346;
		break;
	case BUCAREST:
	case SUNGRID:
		if(instance.at(2) == '3')
			return 40;
		if(instance.at(2) == '4')
			return 70;
		if(instance.at(2) == '5')
			return 100;
		if(instance.at(2) == '6')
			return 130;
		if(instance.at(2) == '8')
			return 190;
		if(instance.at(2) == '0')
			return 250;
		if(instance.at(2) == '2')
			return 310;
		break;
	case VALGRIND:
		return 20 * allowedTime(instance, SAM);
		break;
	default:
		Tools::throwError("Computer not known.");
	}

	std::cout << "# Input problem for the max time function..." << endl;
	return -1;
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
double testMultipleWeeksDeterministic(string dataDir, string instanceName,
   int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outDir) {
   SolverParam param;
   return testMultipleWeeksDeterministic(dataDir, instanceName, historyIndex, weekIndices, algorithm, outDir, param);
}

double testMultipleWeeksDeterministic(string dataDir, string instanceName,
	int historyIndex, vector<int> weekIndices, Algorithm algorithm, string outDir, SolverParam current_param) {

	Scenario* pScen = initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);

	Solver* pSolver = setSolverWithInputAlgorithm(pScen, algorithm);
	pSolver->solve(current_param);

	// Display the solution
	vector<Roster> solution = pSolver->getSolution();
  Status status = pSolver->getStatus();

	displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices, solution,status, outDir);

	delete pSolver;
	delete pScen;

	return pSolver->solutionCost();
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
		string solutionFile = outDir+"Sol-"+instanceName+"-"+catWeeks+"-"+std::to_string(weekIndices[w])+"-"+std::to_string(w)+".txt";
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
}

