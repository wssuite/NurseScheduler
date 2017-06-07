//
//  InitializeSolver.cpp
//  RosterDesNurses
//
//  Created by Jeremy Omer on January 21, 2016
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//

#include "InitializeSolver.h"
//#include "main_test.h"
#include "ReadWrite.h"
#include "Greedy.h"
#include "MasterProblem.h"
#include "SubProblem.h"
#include "MyTools.h"

// some include files to go through the files of an input directory
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

//initialize the counter of objects
unsigned int MyObject::s_count = 0;
unsigned int Rotation::s_count = 0;



/******************************************************************************
* Read the arguments in non compact format
******************************************************************************/

InputPaths* readNonCompactArguments(int argc, char** argv) {
	InputPaths* pInputPaths;

	// Read the arguments and store them in inputPaths
	//
	int narg = 1;
	while (narg < argc) {
		std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;
		// Attention usine a gaz: problem with the java simulator that add quote
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
			pInputPaths->scenario(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--his")) {
			pInputPaths->history(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--week")) {
			pInputPaths->addWeek(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--sol")) {
			pInputPaths->solutionPath(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--param")) {
			pInputPaths->paramFile(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--timeout")) {
			pInputPaths->timeOut(std::stod(str));
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--rand")) {
			pInputPaths->randSeed(std::stoi(str));
			narg += 2;
		}
		else {
			Tools::throwError("main: the argument does not match the expected list!");
		}
	}
	// Throw an error if a necessary input file is missing
	if ( pInputPaths->scenario().empty() || pInputPaths->history().empty()
			|| pInputPaths->weeks().empty() ) {
		throw Tools::myException("readNonCompactArguments: A necessary file name is missing!",__LINE__);
	}

	// Non compact format is only for release versions, so no log file is required
	pInputPaths->logPath("");

	return pInputPaths;
}

/******************************************************************************
* Read the arguments in compact format
******************************************************************************/

InputPaths* readCompactArguments(int argc, char** argv) {

	// Default arguments are set to enable simple call to the function without argument
	//
	std::string dataDir = "datasets/",instanceName = "n005w4",solutionPath="",logPath="",paramFile="";
	int historyIndex = 1, randSeed=0;
	std::vector<int> weekIndices = {6, 2, 9, 1};
	double timeOut = LARGE_TIME;

	// Read the arguments and store them in inputPaths
	//
	int narg = 1;
	while (narg < argc) {
		std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;
		// Attention usine a gaz: problem with the java simulator that add quote
		// marks in the arguments, which fucks up the open file methods
		// the code below is here to remove these quote marks
		//
		string str(argv[narg+1]);
		std::size_t found = str.find("\"");
		while (found!=std::string::npos) {
			str.erase(found,1);
			found = str.find("\"");
		}

		if (!strcmp(argv[narg],"--dir")) {
			dataDir = str;
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--instance")) {
			instanceName = str;
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--his")) {
			historyIndex = std::stoi(str);
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--weeks")) {
			weekIndices = Tools::parseList(str,'-');
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--sol")) {
			solutionPath = str;
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--log")) {
			logPath = str;
			narg+= 2;
		}
		else if (!strcmp(argv[narg],"--param")) {
			paramFile = str;
			narg += 2;
		}
		else if (!strcmp(argv[narg],"--timeout")) {
			timeOut = std::stod(str);
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

	// Initialize the input paths
	//
	InputPaths* pInputPaths =
	new InputPaths(dataDir, instanceName, historyIndex,weekIndices,solutionPath,logPath,paramFile,timeOut,randSeed);

	return pInputPaths;
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

	delete pPref;

	return pScen;
}

Scenario* initializeScenario(const InputPaths & inputPaths, string logPath) {

	// Initialize demand and preferences
	Demand* pDemand(0);
	Preferences* pPref(0);

	// Read the scenario
	Scenario* pScen = ReadWrite::readScenario(inputPaths.scenario());

	// Read the demand and preferences and link them with the scenario
	ReadWrite::readWeek(inputPaths.week(0), pScen,&pDemand,&pPref);
	pScen->linkWithDemand(pDemand);
	pScen->linkWithPreferences(*pPref);

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

Scenario* initializeMultipleWeeks(const InputPaths & inputPaths, string logPath) {

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

/*****************************************************************************
* Separate the scenario into multiple scenarios that only focus on the nurses
* whose positions are in the same connex component of positions
******************************************************************************/

vector<Scenario*> divideScenarioIntoConnexPositions(Scenario* pScenario) {

	vector<Scenario*> scenariosPerComponent;

	// First, identify the connex components of the graph of positions
	pScenario->computeConnexPositions();


	for (int c = 0; c < pScenario->nbOfConnexComponentsOfPositions(); c++) {
		vector<Position*> positionsInTheComponent = pScenario->componentOfConnexPositions(c);
		vector<Nurse> nursesInTheComponent = pScenario->nursesInConnexComponentOfPositions(c);

		// retrieve a vector containing the indices of the skills contained in the component (decreasing order)
		// use a set first, because it manages duplicate skills automatically
		std::set<int> skillsInTheComponent;
		for (Position* pPosition: positionsInTheComponent) {
			for (int skill: pPosition->skills()) {
				skillsInTheComponent.insert(skill);
			}
		}
		vector<int> skillsVector;
		for (set<int>::iterator it = skillsInTheComponent.begin(); it != skillsInTheComponent.end(); it++) {
			skillsVector.push_back(*it);
		}
		std::stable_sort(skillsVector.begin(),skillsVector.end(), Tools::compareDecreasing);

		// build a vector containinf the skills that need to be removed from the scenario
		vector<int> skillsToRemove;
		for (int skill=0; skill < pScenario->nbSkills_; skill++) {
			skillsToRemove.push_back(skill);
		}
		for (int skill:skillsVector) {
			skillsToRemove.erase(skillsToRemove.begin()+skill);
		}
		std::stable_sort(skillsToRemove.begin(),skillsToRemove.end(), Tools::compareDecreasing);

		// shorten the vectors intToSkill and skillToInt to match the list of skills
		vector<string> intToSkill(pScenario->intToSkill_);
		map<string,int> skillToInt(pScenario->skillToInt_);

		for (int skill:skillsToRemove) {
			skillToInt.erase(pScenario->intToSkill_[skill]);
			intToSkill.erase(intToSkill.begin()+skill);
		}

		// create the demand that relates only to input skills
		//
		Demand* pDemand = pScenario->pWeekDemand();

		// erase the skills to remove from the minimum and optimal demands
		vector3D minDemand = pDemand->minDemand_;
		vector3D optDemand = pDemand->optDemand_;
		for (int day = 0; day < pDemand->nbDays_; day++) {
    		for (int shift = 0; shift < pDemand->nbShifts_; shift++) {
				for (int skill:skillsToRemove) {
					minDemand[day][shift][skill] = 0;
					optDemand[day][shift][skill] = 0;
				}
			}
		}
		Demand* pDemandInTheComponent = new Demand(pDemand->nbDays_, pDemand->firstDay_, pDemand->nbShifts_,
		 pDemand->nbSkills_, pDemand->name_, minDemand, optDemand);

		// create the preferences that relate only to the nurses of the connex component
		//
		Preferences* pPreferencesInTheComponent = new Preferences(nursesInTheComponent,pDemand->nbDays_,pScenario->nbShifts_);
		Preferences* pPreferences = pScenario->pWeekPreferences();

		// only keep the demand of the nurses in the component
		for (int i = 0; i < nursesInTheComponent.size(); i++) {
			Nurse nurse = nursesInTheComponent[i];
			map<int,std::set<int> >::iterator itDay;
			for (itDay = pPreferences->nurseWishesOff(nurse.id_)->begin(); itDay != pPreferences->nurseWishesOff(nurse.id_)->end(); itDay++) {
				set<int>::iterator itShift;
				for (itShift = (*itDay).second.begin(); itShift != (*itDay).second.end(); itShift++) {
					pPreferencesInTheComponent->addShiftOff(i,(*itDay).first, *itShift);
				}
			}
		}

		// Create the new scenario
		//
		Scenario* pScenarioInTheConnexComponent = new Scenario(pScenario,nursesInTheComponent,pDemandInTheComponent,pPreferencesInTheComponent);

		// create the initial states that relate only to the nurses of the connex component
		//
		vector<State> intialStatesInTheComponent;
		vector<State>* pInitialState = pScenario->pInitialState();

		for (Nurse nurse: nursesInTheComponent) {
			intialStatesInTheComponent.push_back(pInitialState->at(nurse.id_));
		}
		pScenarioInTheConnexComponent->setInitialState(intialStatesInTheComponent);

		// Push back in the vector of scenarios
		//
		scenariosPerComponent.push_back(pScenarioInTheConnexComponent);

		//delete pDemandInTheComponent;
		delete pPreferencesInTheComponent;
	}

	return scenariosPerComponent;
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
		// DBG: ICI, ON CHOISIT SI ON UTILISE CLP OU GUROBI, A TERME IL CHOISIR EN FONCTION D'UNE OPTION...
//	   		pSolver = new MasterProblem(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState(), S_Gurobi);
		pSolver = new MasterProblem(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState(), S_CLP);
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

void displaySolutionMultipleWeeks(InputPaths inputPaths, vector<Roster> &solution, Status status, string outDir) {

	if (outDir.empty()) return;

	// initialize the log stream
	// first, concatenate the week numbers
	int nbWeeks = inputPaths.nbWeeks();
	string catWeeks;
	for (int w=0; w< nbWeeks; w++) catWeeks += inputPaths.week(w);
	string logPath = outDir+"Log-"+catWeeks+".txt";
	Tools::LogOutput outStream(logPath);

	// treat the case where the solver was unable to find a feasible solution
	if (status == INFEASIBLE) {
		outStream << "The solver was not able to find a solution\n";
		return;
	}

	// load the solution in a new solver
	Scenario* pScen = initializeMultipleWeeks(inputPaths);
	Solver* pSolver = new Solver(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
	pSolver->loadSolution(solution);

	// write the log file for all the weeks
	outStream << pSolver->solutionToLogString();

	// write separately the solutions of each week in the required output format
	vector<string> solutions = pSolver->solutionToString(nbWeeks);
	for(int w=0; w < nbWeeks; ++w){
		string solutionFile = outDir+"Sol-"+inputPaths.instance()+"-"+catWeeks+"-"+inputPaths.week(w)+"-"+std::to_string(w)+".txt";
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
