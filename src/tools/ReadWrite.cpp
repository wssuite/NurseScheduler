#include "tools/ReadWrite.h"
#include "tools/MyTools.h"
#include "data/Scenario.h"
#include "solvers/Solver.h"
#include "solvers/StochasticSolver.h"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

// some include files to go through the files of an input directory
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/assign/list_of.hpp>

using std::string;
using std::vector;
using std::map;
using std::pair;

std::map<std::string, Algorithm> stringToAlgorithm =
	boost::assign::map_list_of("GREEDY", GREEDY)("GENCOL", GENCOL)("STOCHASTIC_GREEDY",STOCHASTIC_GREEDY)("STOCHASTIC_GENCOL",STOCHASTIC_GENCOL)("NONE",NONE);
std::map<std::string, WeightStrategy> stringToWeightStrategy =
	boost::assign::map_list_of("MAX", MAX)("MEAN", MEAN)("RANDOMMEANMAX",RANDOMMEANMAX)("BOUNDRATIO",BOUNDRATIO)("NO_STRAT",NO_STRAT);
std::map<std::string, RankingStrategy> stringToRankingStrategy =
	boost::assign::map_list_of("SCORE", RK_SCORE)("MEAN", RK_MEAN);


//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
Scenario* ReadWrite::readScenario(string fileName) {
	// open the file
	std::fstream file;
	std::cout << "Reading " << fileName << std::endl;
	file.open(fileName.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read the file " << fileName << std::endl;
		std::cout << "The input file was not opened properly!" << std::endl;

		throw Tools::myException("The input file was not opened properly!",__LINE__);
	}

	string title;
	string strTmp;
	int intTmp;
	// declare the attributes that will initialize the Scenario instance
	//
	string name;
	int nbWeeks=-1, nbSkills=-1, nbShifts=-1, nbShiftsType=-1, nbContracts=-1, nbNurses=-1;
	vector<string> intToSkill, intToShift, intToShiftType, intToContract;
	map<string,int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
	vector<int> minConsShiftType, maxConsShiftType, nbForbiddenSuccessors, hoursInShift, shiftIDToShiftTypeID;
        vector2D<int> shiftTypeIDToShiftID;
	vector2D<int> forbiddenSuccessors;
	map<string,Contract*> contracts;
	vector<Nurse> theNurses;

	bool  foundShift = false;

	// fill the attributes of the scenario structure
	//
	while(file.good()){
		readUntilChar(&file, '=', &title);

		// Read the name of the scenario
		//
		if(strEndsWith(title, "SCENARIO ")){
			file >> name;
		}

		// Read the number of weeks in scenario
		//
		else if (strEndsWith(title, "WEEKS ")) {
			file >> nbWeeks;
		}

		// Read the number of weeks in scenario
		//
		else if (strEndsWith(title, "SKILLS ")) {
			file >> nbSkills;
			for(int i=0; i<nbSkills; i++){
				file >> strTmp;
				intToSkill.push_back(strTmp);
				skillToInt.insert(pair<string,int>(strTmp,i));
			}
		}

		// Read the different shift types and forbidden successions
		//
		else if (strEndsWith(title, "SHIFT_TYPES ")) {

			// Number of shifts : Given number + REST_SHIFT
			file >> intTmp;
			nbShiftsType = intTmp+1;

			// IMPORTANT : INSERT REST SHIFT !!!!!!
			// It is given 0 and 99 as bounds so that they never perturbate the cost
			//
			intToShiftType.push_back(REST_SHIFT);
			shiftTypeToInt.insert(pair<string,int>(REST_SHIFT,0));
			minConsShiftType.push_back(0);
			maxConsShiftType.push_back(99);

			// Other shift types
			//
			for(int i=1; i<nbShiftsType; i++){
				// Name
				file >> strTmp;
				intToShiftType.push_back(strTmp);
				shiftTypeToInt.insert(pair<string,int>(strTmp,i));
				readUntilChar(&file,'(',&strTmp);
				// Min consecutive
				file >> intTmp;
				minConsShiftType.push_back(intTmp);
				readUntilChar(&file, ',', &strTmp);
				// Max consecutive
				file >> intTmp;
				maxConsShiftType.push_back(intTmp);
				readUntilChar(&file,'\n',&strTmp);
			}


			// Forbidden successions
			//
			for(int i=0; i<nbShiftsType; i++){
				vector<int> v;
				forbiddenSuccessors.push_back(v);
				nbForbiddenSuccessors.push_back(0);
			}
			while(!strEndsWith(title,"FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
				file >> title;
			// Reading all lines
			for(int i=1; i<nbShiftsType; i++){
				// Which current shift type ?
				string currentShiftType;
				file >> currentShiftType;
				int currentShiftTypeId = shiftTypeToInt.at(currentShiftType);
				// How many forbidden after it ?
				file >> intTmp;
				nbForbiddenSuccessors[currentShiftTypeId] = intTmp;
				// Which ones are forbidden ?
				for(int j=0; j<nbForbiddenSuccessors[currentShiftTypeId]; j++){
					file >> strTmp;
					forbiddenSuccessors[currentShiftTypeId].push_back(shiftTypeToInt.at(strTmp));
				}
				readUntilChar(&file,'\n',&strTmp);

			}
		}

		// Read the different shifts
		//
		else if (strEndsWith(title, "SHIFTS ")) {
		  foundShift = true;

			// IMPORTANT : INSERT REST SHIFT !!!!!!
			//
			intToShift.push_back(REST_SHIFT);
			shiftToInt.insert(pair<string,int>(REST_SHIFT,0));
			hoursInShift.push_back(0);
			int currentShiftTypeId = shiftTypeToInt.at(REST_SHIFT);
			shiftIDToShiftTypeID.push_back(currentShiftTypeId);

			shiftTypeIDToShiftID.resize(nbShiftsType);
			shiftTypeIDToShiftID[0].push_back(0);
	
			// Number of shifts 
			file >> intTmp;
			nbShifts = intTmp+1;

			// Shifts
			//
			for(int i=1; i<nbShifts; i++){
				// Name
				file >> strTmp;
				intToShift.push_back(strTmp);
				shiftToInt.insert(pair<string,int>(strTmp,i));

				int hours;
				file >> hours;
				hoursInShift.push_back(hours);

				string  currentShiftType;
				file >> currentShiftType;

				int currentShiftTypeId = shiftTypeToInt.at(currentShiftType);
				shiftIDToShiftTypeID.push_back(currentShiftTypeId);
				shiftTypeIDToShiftID[currentShiftTypeId].push_back(i);
			}
		}

		// Read the different contracts type
		//
		else if (strEndsWith(title, "CONTRACTS ")) {
			file >> intTmp;
			nbContracts = intTmp;
			// Read each contract type
			for(int i=0; i<nbContracts; i++){
				string contractName;
				int minDays, maxDays, minConsWork, maxConsWork, minConsRest, maxConsRest, maxWeekends, isTotalWeekend;
				file >> contractName;
				readUntilChar(&file,'(',&strTmp);
				file >> minDays;
				readUntilChar(&file,',',&strTmp);
				file >> maxDays;
				readUntilChar(&file,'(',&strTmp);
				file >> minConsWork;
				readUntilChar(&file,',',&strTmp);
				file >> maxConsWork;
				readUntilChar(&file,'(',&strTmp);
				file >> minConsRest;
				readUntilChar(&file,',',&strTmp);
				file >> maxConsRest;
				readUntilChar(&file,' ',&strTmp);
				file >> maxWeekends;
				readUntilChar(&file,' ',&strTmp);
				file >> isTotalWeekend;
				readUntilChar(&file,'\n',&strTmp);

				Contract * pContract = new Contract (i, contractName, minDays, maxDays, minConsWork, maxConsWork, minConsRest, maxConsRest, maxWeekends, isTotalWeekend);
				contracts.insert(pair<string,Contract*>(contractName,pContract));
				intToContract.push_back(contractName);
			}
		}

		// Read all nurses
		//
		else if (strEndsWith(title, "NURSES ")) {
			file >> nbNurses;
			for(int i=0; i<nbNurses; i++){
				string nurseName, contractName;
				int nbSkills;
				vector<int> skills;
				// Read everything on the line
				file >> nurseName;
				file >> contractName;
				file >> nbSkills;
				for(int j=0; j<nbSkills; j++){
					file >> strTmp;
					skills.push_back(skillToInt.at(strTmp));
				}
				// sort the skill indices before initializing the nurse
				std::sort (skills.begin(), skills.end());

				Nurse nurse (i, nurseName, nbSkills, skills, contracts.at(contractName));
				theNurses.push_back(nurse);
				nurseNameToInt.insert(pair<string,int>(nurseName,i));
			}

		}
	}

	//  to be backward compatible with old style of input file (without the SHIFTS section)
	//  set default shifts
	
	if (!foundShift) {
	  nbShifts=nbShiftsType;
    intToShift = intToShiftType;
    shiftToInt = shiftTypeToInt;
    shiftTypeIDToShiftID.resize(nbShiftsType);
    for(int i=0; i<nbShiftsType; i++){
	    hoursInShift.push_back(1); // 1 as default (could be days, hours, ...)
	    shiftIDToShiftTypeID.push_back(i);
      shiftTypeIDToShiftID[i].push_back(i);
	  }
	}

	// Check that all fields were initialized before initializing the scenario
	//
	if ( nbWeeks==-1 || nbSkills==-1 || nbShifts==-1 || nbContracts==-1 || nbNurses==-1 ) {
		Tools::throwError("In readScenario: missing fields in the initialization");
	}

	return new Scenario(name, nbWeeks, nbSkills, intToSkill, skillToInt, nbShifts,
			    intToShift, shiftToInt, hoursInShift, shiftIDToShiftTypeID,
			    nbShiftsType, intToShiftType, shiftTypeToInt, shiftTypeIDToShiftID, 
			    minConsShiftType, maxConsShiftType, nbForbiddenSuccessors,forbiddenSuccessors,
			    nbContracts, intToContract, contracts, nbNurses, theNurses, nurseNameToInt) ;
}

Demand* ReadWrite::readWeeks(std::vector<std::string> strWeekFiles, Scenario* pScenario)
{
	//initialize pDemand
	Demand* pDemand(nullptr);
	Preferences* pPref(nullptr);

	for(string strWeekFile: strWeekFiles)
		if(pDemand == 0){
			ReadWrite::readWeek(strWeekFile, pScenario,&pDemand,&pPref);
      std::cout << pPref->toString(pScenario) << std::endl;
		}
		else{
			//load the next week
			Demand* nextDemand(nullptr);
			Preferences* nextPref(nullptr);
			ReadWrite::readWeek(strWeekFile, pScenario, &nextDemand, &nextPref);
      std::cout << nextPref->toString(pScenario) << std::endl;
			//update the current weeks
			pDemand->push_back(nextDemand);
			pPref->push_back(nextPref);
			pScenario->addAWeek();
      std::cout << pPref->toString(pScenario) << std::endl;
			//delete the demand and the preferences which we have created
			delete nextDemand;
			delete nextPref;
		}

	//link the scenario to the current demand and preferences
	pScenario->linkWithDemand(pDemand);
	pScenario->linkWithPreferences(pPref);

	return pDemand;
}

// Read the Week file and store the content in a Scenario instance
//
void ReadWrite::readWeek(std::string strWeekFile, Scenario* pScenario,
  Demand** pDemand, Preferences** pPref){
	// open the file
	std::fstream file;
	std::cout << "Reading " << strWeekFile << std::endl;
	file.open(strWeekFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read the file " << strWeekFile << std::endl;
		std::cout << "The input file was not opened properly!" << std::endl;
		throw Tools::myException("The input file was not opened properly!",__LINE__);
	}

	string title;
	string strTmp;
	int intTmp;

	// declare the attributes to be updated in the Scenario*
	//
	string weekName;
	vector3D<int> minWeekDemand;
	vector3D<int> optWeekDemand;
	if (*pPref) delete *pPref;
	if (*pDemand) delete *pDemand;


	// fill the attributes when reading the week file
	//
	while(file.good()){
		readUntilOneOfTwoChar(&file, '\n', '=', &title);

		// Read the name of the week
		//
		if(strEndsWith(title, "WEEK_DATA")){
			file >> weekName;
		}

		// Read the requirements
		//
		else if (strEndsWith(title, "REQUIREMENTS")) {
			string shiftName, skillName;
			int shiftId, skillId;
			// init the vectors
			Tools::initVector3D(minWeekDemand, 7, pScenario->nbShifts_, pScenario->nbSkills_, 0);
			Tools::initVector3D(optWeekDemand, 7, pScenario->nbShifts_, pScenario->nbSkills_, 0);

			// Do not take the rest shift into account here (by initialization, requirements already at 0
			for(int i=1; i<pScenario->nbShifts_; i++){
				for(int j=0; j<pScenario->nbSkills_; j++){
					// Read shift and skill
					file >> shiftName;
					file >> skillName;
					shiftId = pScenario->shiftToInt_.at(shiftName);
					skillId = pScenario->skillToInt_.at(skillName);
					// For every day in the week, read min and opt values
					for (int day = 0; day<7; day++){
						readUntilChar(&file,'(',&strTmp);
						file >> intTmp;
						minWeekDemand[day][shiftId][skillId] = intTmp;
						readUntilChar(&file,',',&strTmp);
						file >> intTmp;
						optWeekDemand[day][shiftId][skillId] = intTmp;
					}
					readUntilChar(&file,')',&strTmp);
				}
			}
		}

		// Read the shift off requests
		//
		else if(strEndsWith(title,"SHIFT_OFF_REQUESTS ")){
      if(!*pPref)
        *pPref = new Preferences(pScenario->nbNurses_, 7, pScenario->nbShifts_);
			// Temporary vars
			string nurseName, shift, day;
			int nbShiftOffRequests, nurseId, shiftId, dayId, level;
			file >> nbShiftOffRequests;
			for (int i=0; i<nbShiftOffRequests; i++){
				file >> nurseName;
				file >> shift;
				file >> day;
				file >> level;
				nurseId = pScenario->nurseNameToInt_.at(nurseName);
				dayId = Tools::dayToInt(day);

				if(shift == "Any")
				  (*pPref)->addDayOff(nurseId, dayId, level);
				else {
					// shiftId = pScenario->shiftTypeToInt_.at(shift);
					shiftId = pScenario->shiftToInt_.at(shift);
					(*pPref)->addShiftOff(nurseId, dayId, shiftId, level);
				}
			}
		}

		// Read the shift on requests
		//
		else if(strEndsWith(title,"SHIFT_ON_REQUESTS ")){
			if(!*pPref)
			  *pPref = new Preferences(pScenario->nbNurses_, 7, pScenario->nbShifts_);
			// Temporary vars
			string nurseName, shift, day;
			int nbShiftOnRequests, nurseId, shiftId, dayId, level;
			file >> nbShiftOnRequests;
			for (int i=0; i<nbShiftOnRequests; i++){
				file >> nurseName;
				file >> shift;
				file >> day;
				file >> level;
				nurseId = pScenario->nurseNameToInt_.at(nurseName);
				dayId = Tools::dayToInt(day);

				if(shift == "Any")
				  (*pPref)->addDayOn(nurseId, dayId, level);
				else {
					// shiftId = pScenario->shiftTypeToInt_.at(shift);
					shiftId = pScenario->shiftToInt_.at(shift);
					(*pPref)->addShiftOn(nurseId, dayId, shiftId, level);
				}
			}
		}
	}

	// Define a new instance of demand
	*pDemand = new Demand(7, 0, pScenario->nbShifts_,pScenario->nbSkills_, weekName,
		minWeekDemand, optWeekDemand);
  std::cout << "Demand created" << std::endl;

	// Now, add all these objects to the Scenario
  //  pScenario->linkWithDemand(pDemand);
  //  pScenario->linkWithPreferences(weekPreferences);

	// return the demand
  //  return pDemand;

}

// Read the history file
//
void ReadWrite::readHistory(std::string strHistoryFile, Scenario* pScenario){
	// open the file
	std::fstream file;
	std::cout << "Reading " << strHistoryFile << std::endl;
	file.open(strHistoryFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strHistoryFile << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;
	string strTmp;

	// declare the attributes to be updated in the Scenario*
	//
	int thisWeek;
	string weekName;
	vector<State> initialState;


	// fill the attributes of the week structure
	//
	while(file.good()){
		readUntilChar(&file,'\n', &title);

		// Read the index and name of the week
		//
		if(!strcmp(title.c_str(), "HISTORY")){
			file >> thisWeek;
			file >> weekName;
			// Raise exception if it does not match the week previously read !
			if (strcmp(weekName.c_str(),(pScenario->weekName()).c_str())) {
				std::cout << "The given history file requires week " << weekName << std::endl;
				std::cout << " but a different one (" << pScenario->weekName() << ") has been given!" << std::endl;
				Tools::throwError("History file and week data file do not match!");
			}
		}

		// Read each nurse's initial state
		//
		else if (strEndsWith(title, "NURSE_HISTORY")) {
			for(int n=0; n<pScenario->nbNurses_; n++){
				string nurseName, shiftTypeName;
				// int nurseId;
				int shiftTypeId, totalTimeWorked, totalWeekendsWorked, consDaysWorked, consShiftWorked, consRest, consShifts;
				file >> nurseName;
				// nurseId = pScenario->nurseNameToInt_.at(nurseName);
				file >> totalTimeWorked;
				file >> totalWeekendsWorked;
				file >> shiftTypeName;
				shiftTypeId = pScenario->shiftTypeToInt_.at(shiftTypeName);
				file >> consShiftWorked;
				file >> consDaysWorked;
				file >> consRest;

				if(consRest == 0 && consDaysWorked == 0)
                    Tools::throwError("History of nurse " + nurseName + " is invalid as one must either work or rest.");

				int  shiftID = pScenario->shiftTypeIDToShiftID_[shiftTypeId].front();
				consShifts = (shiftTypeId == 0) ? consRest : consShiftWorked;
				State nurseState (0, totalTimeWorked, totalWeekendsWorked,
						  consDaysWorked, consShifts, consRest, shiftTypeId, shiftID);
				initialState.push_back(nurseState);
			}
		}
	}
	pScenario->setThisWeek(thisWeek);
	pScenario->setInitialState(initialState);
}


// Read the input custom file
// Store the result in a vector of historical demands and return the number of treated weeks
//
int ReadWrite::readCustom(string strCustomInputFile, Scenario* pScenario, vector<Demand*>& demandHistory) {
	// open the file
	std::fstream file;
	std::cout << "Reading " << strCustomInputFile << std::endl;
	file.open(strCustomInputFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strCustomInputFile << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;
	int nbWeeks;

	// get the custom information
	//
	while(file.good()){
		readUntilOneOfTwoChar(&file, '\n', '=', &title);

		// Read the file names of the past demand
		//
		if(!strcmp(title.c_str(), "PAST_DEMAND_FILES")){
			file >> nbWeeks;
			if (!nbWeeks) continue;

			string strDemandFile;
			for (int i = 0; i < nbWeeks; i++) {
				file >> strDemandFile;
				Demand* pDemand = nullptr;
				Preferences* pPref = nullptr;
				readWeek(strDemandFile,pScenario,&pDemand,&pPref);
				demandHistory.push_back(pDemand);
			}
		}
	}
	return nbWeeks;
}

void ReadWrite::writeCustom(string strCustomOutputFile, string strWeekFile, string strCustomInputFile) {

	Tools::LogOutput outStream(strCustomOutputFile);

	// if there is no custom input file, this is the first week
	if (strCustomInputFile.empty()) {
		outStream << "PAST_DEMAND_FILES= " << 1 << std::endl;
		outStream << strWeekFile<< std::endl;
		return;
	}

	// open the custom input file
	// we want the content of the input custom file in the custom output file
	std::fstream file;
	std::cout << "Reading " << strCustomInputFile << std::endl;
	file.open(strCustomInputFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strCustomInputFile << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;
	int nbWeeks;

	// fill the attributes of the week structure
	//
	while(file.good()){
		readUntilOneOfTwoChar(&file, '\n', '=', &title);

		// Read the file names of the past demand
		//
		if(!strcmp(title.c_str(), "PAST_DEMAND_FILES")){
			file >> nbWeeks;
			outStream << "PAST_DEMAND_FILES= " << nbWeeks+1 << std::endl;
			if (!nbWeeks) continue;

			string strDemandFile;
			for (int i = 0; i < nbWeeks; i++) {
				file >> strDemandFile;
				outStream << strDemandFile << std::endl;
			}
			outStream << strWeekFile << std::endl;
		}
	}
}


/************************************************************************
* Read the options of the stochastic and ot the other solvers
*************************************************************************/
std::string ReadWrite::readStochasticSolverOptions(string strOptionFile, StochasticSolverOptions& options) {

	// open the file
	std::fstream file;
	std::cout << "Reading " << strOptionFile << std::endl;
	file.open(strOptionFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strOptionFile << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;

	// fill the attributes of the options structure
	//
	while(file.good()){
		readUntilOneOfTwoChar(&file, '\n', '=', &title);

		if(!strcmp(title.c_str(), "withEvaluation")){
			file >> options.withEvaluation_;
		}
		if(!strcmp(title.c_str(), "withIterativeDemandIncrease")){
			file >> options.withIterativeDemandIncrease_;
		}
		if(!strcmp(title.c_str(), "generationCostPerturbation")){
			file >> options.generationCostPerturbation_;
		}
		if(!strcmp(title.c_str(), "evaluationCostPerturbation")){
			file >> options.evaluationCostPerturbation_;
		}
		if(!strcmp(title.c_str(), "generationAlgorithm")){
			string strtmp;
			file >> strtmp;
			options.generationAlgorithm_ = stringToAlgorithm[strtmp];
		}
		if(!strcmp(title.c_str(), "evaluationAlgorithm")){
			string strtmp;
			file >> strtmp;
			options.evaluationAlgorithm_ = stringToAlgorithm[strtmp];
		}
		if(!strcmp(title.c_str(), "rankingStrategy")){
			string strtmp;
			file >> strtmp;
			options.rankingStrategy_ = stringToRankingStrategy[strtmp];
		}
		if(!strcmp(title.c_str(), "nExtraDaysGenerationDemands")){
			file >> options.nExtraDaysGenerationDemands_;
		}
		if(!strcmp(title.c_str(), "nEvaluationDemands")){
			file >> options.nEvaluationDemands_;
		}
		if(!strcmp(title.c_str(), "nDaysEvaluation")){
			file >> options.nDaysEvaluation_;
		}
		if(!strcmp(title.c_str(), "nGenerationDemandsMax")){
			file >> options.nGenerationDemandsMax_;
		}
	}

  std::ifstream fin(strOptionFile.c_str());
  std::ostringstream sout;
  while(fin.good()) copy(std::istreambuf_iterator<char>(fin),
       std::istreambuf_iterator<char>(),
       std::ostreambuf_iterator<char>(sout));
   return sout.str();
}

std::string ReadWrite::readSolverOptions(string strOptionFile, SolverParam& options) {
	// open the file
	std::fstream file;
	std::cout << "Reading " << strOptionFile << std::endl;
	file.open(strOptionFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strOptionFile << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;

	// fill the attributes of the options structure
	//
	while(file.good()){
		readUntilOneOfTwoChar(&file, '\n', '=', &title);

		if(!strcmp(title.c_str(), "maxSolvingTimeSeconds")){
			file >> options.maxSolvingTimeSeconds_;
		}
		if(!strcmp(title.c_str(), "printEverySolution")){
			file >> options.printEverySolution_;
		}
		if(!strcmp(title.c_str(), "absoluteGap")){
			file >> options.absoluteGap_;
		}
		if(!strcmp(title.c_str(), "minRelativeGap")){
			file >> options.minRelativeGap_;
		}
		if(!strcmp(title.c_str(), "relativeGap")){
			file >> options.relativeGap_;
		}
		if(!strcmp(title.c_str(), "nbDiveIfMinGap")){
			file >> options.nbDiveIfMinGap_;
		}
		if(!strcmp(title.c_str(), "nbDiveIfRelGap")){
			file >> options.nbDiveIfRelGap_;
		}
		if(!strcmp(title.c_str(), "solveToOptimality")){
			file >> options.solveToOptimality_;
		}
		if(!strcmp(title.c_str(), "weightStrategy")){
			string strtmp;
			file >> strtmp;
			options.weightStrategy_ = stringToWeightStrategy[strtmp];
		}
		if(!strcmp(title.c_str(), "stopAfterXSolution")){
			file >> options.stopAfterXSolution_;
		}
	}

  std::ifstream fin(strOptionFile.c_str());
  std::ostringstream sout;
  while(fin.good()) copy(std::istreambuf_iterator<char>(fin),
       std::istreambuf_iterator<char>(),
       std::ostreambuf_iterator<char>(sout));
   return sout.str();
}

// Read the solution from multiple week solution files
//
vector<Roster> ReadWrite::readSolutionMultipleWeeks(vector<std::string> strWeekSolFiles, Scenario* pScenario){

	int nbWeeks = strWeekSolFiles.size();
	string  strNurse, strDay, strShift, strSkill;
	int nurse, day, shift, skill;


	// initialize the solution
	vector<Roster> solution;
	for (int n=0; n < pScenario->nbNurses(); n++) {
		vector<int> shifts(7*nbWeeks,0);
		vector<int> skills(7*nbWeeks,-1);

		Roster roster(7*nbWeeks, 0, shifts, skills);
		solution.push_back(roster);
	}

	// read the solution of each week
	for (int w = 0; w < nbWeeks; w++) {
		std::string strSolFile = strWeekSolFiles[w];
		int firstDay = 7*w;

		// open the week solution file
		std::fstream file;
		std::cout << "Reading " << strSolFile << std::endl;
		file.open(strSolFile.c_str(), std::fstream::in);
		if (!file.is_open()) {
			std::cout << "While trying to read " << strSolFile << std::endl;
			Tools::throwError("The input file was not opened properly!");
		}
		string title;

		// parse the file until reaching the number of assignments
		std::size_t found = title.find("ASSIGNMENTS");
		while (found!=std::string::npos) {
			readUntilOneOfTwoChar(&file, '\n', '=', &title);
			found = title.find("ASSIGNMENTS");
			if (!file.good()) {
				Tools::throwError("The solution file is not formatted as it should!");
			}
		}
		readUntilOneOfTwoChar(&file, '\n', '\n', &title);

		// parse the assignments
		while (file.good()) {
			file >> strNurse >> strDay >> strShift >> strSkill;
			nurse = pScenario->nurseNameToInt_.at(strNurse);
			day = firstDay+Tools::dayToInt(strDay);
			shift = pScenario->shiftToInt_.at(strShift);
			skill = pScenario->skillToInt_.at(strSkill);

			solution[nurse].assignTask(day,shift,skill);
		}
	}

	return solution;
}

/************************************************************************
* Print the main characteristics of all the demands of an input directory
* This is done to find some invariant properties among demands
*************************************************************************/

void ReadWrite::compareDemands(string inputDir, string logFile) {

	struct dirent *dirp;
	Tools::LogOutput logStream(logFile,8);

	vector2D<int> minPerShift, optPerShift, minPerSkill,optPerSkill;
	vector2D<int> minHighestPerSkill, optHighestPerSkill;
  vector<int> minTotal, optTotal;

	// Open the input directory
	DIR* dp = opendir( inputDir.c_str() );
	if (dp == NULL) {
		Tools::throwError("Error while opening ");
	}
	else{
		std::cout << "Reading from directory " << inputDir << std::endl;
	}

	// Read the scenario that appears in the directory
	unsigned found = inputDir.find_last_of("/");
  string instanceName = inputDir.substr(found+1);
	string scenFile = inputDir+"/Sc-"+instanceName+".txt";
  string historyFile = inputDir+"/H0-"+instanceName+"-0.txt";

	Scenario* pScen = ReadWrite::readScenario(scenFile);

	// Go through all the demand files of the directory
  int coDemand = 0;
	while ((dirp = readdir( dp )))
	{
		std::string filename(dirp->d_name);

		// The file names of week demands start with "WD"
		std::size_t found= filename.find("WD");
		if (found > 0) continue;

		string filepath = inputDir + "/" + filename;

		Demand* pDemand(nullptr);
		Preferences* pPref(nullptr);
		ReadWrite::readWeek(filepath, pScen, &pDemand, &pPref);

		logStream << "#####################################\n";
		logStream << "# DEMAND FILE: " << filepath << std::endl;
		logStream << "#####################################\n\n";
		logStream << pDemand->toString(true) << std::endl << std::endl;

		// record the advanced data on the demand
	 minTotal.push_back(pDemand->minTotal_);
	 optTotal.push_back(pDemand->optTotal_);
		minPerShift.push_back(pDemand->minPerShift_);
		optPerShift.push_back(pDemand->optPerShift_);
		minPerSkill.push_back(pDemand->minPerSkill_);
		optPerSkill.push_back(pDemand->optPerSkill_);
		minHighestPerSkill.push_back(pDemand->minHighestPerSkill_);
		optHighestPerSkill.push_back(pDemand->optHighestPerSkill_);

	 // link the scenario with the first demand and preferences to be able to
	 // retrieve information about the nurses
	 if (!coDemand) {
		pScen->linkWithDemand(pDemand);
		pScen->linkWithPreferences(pPref);
	 }
	 else {
		delete pDemand;
		delete pPref;
		coDemand++;
	 }
	}

  // Also preprocess the nurses to get statistics on the capacity of the nurses
  // to cover the demand
  //
  ReadWrite::readHistory(historyFile, pScen);
  Solver* pSolver = new Solver(pScen, pScen->pWeekDemand(), pScen->pWeekPreferences(), pScen->pInitialState());
  pSolver->preprocessTheNurses();

	// Write a summary  of the advanced data computed for the demands
	logStream << "#####################################" << std::endl;
	logStream << "# SUMMARY OF THE STATISTICS" << std::endl;
	logStream << "#####################################\n\n";

  logStream  << "# Total capacity of the nurses (without unavoidable penalty/with average work): ";
  logStream << Tools::itoa(pSolver->maxTotalStaffNoPenalty_)+"/"+Tools::itoa(pSolver->maxTotalStaffAvgWork_);
  logStream.endl();
  logStream.endl();


  logStream  << "# Total capacity of the nurses per skill (without unavoidable penalty/with average work):\n";
  for (int i = 0; i < pScen->nbSkills_; i++)	{
	 logStream << "# SK"+Tools::itoa(i)+":";
	 logStream << Tools::itoa(pSolver->maxStaffPerSkillNoPenalty_[i]) +"/"+Tools::itoa(pSolver->maxStaffPerSkillAvgWork_[i]);
	 logStream.endl();
  }
  logStream.endl();

  logStream  << "# Total demand (min/opt): \n" ;
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
	 logStream << "WD"+Tools::itoa(d);
  }
  logStream.endl();
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
	 logStream << Tools::itoa(minTotal[d])+"/"+Tools::itoa(optTotal[d]);
  }
  logStream.endl();
  logStream.endl();


	logStream  << "# Demand per shift\n";
	logStream << "#";
	for (unsigned int d = 0; d < minPerShift.size(); d++) {
	 logStream << "WD"+Tools::itoa(d);
	}
	logStream.endl();
	for (int i = 1; i < pScen->nbShifts_; i++)	{
		logStream << "# SH"+Tools::itoa(i)+":";
		for (unsigned int d = 0; d < minPerShift.size(); d++) {
			logStream << Tools::itoa(minPerShift[d][i])+"/"+Tools::itoa(optPerShift[d][i]);
		}
		logStream.endl();
	}
  logStream.endl();


	logStream  << "# Demand per skill\n";
  logStream << "#";
	for (unsigned int d = 0; d < minPerSkill.size(); d++) {
	 logStream << "WD"+Tools::itoa(d);
	}
	logStream.endl();
	for (int i = 0; i < pScen->nbSkills_; i++)	{
	 logStream.setWidth(12);
	 logStream << "# " + pScen->intToSkill_[i] +":";
	 logStream.setWidth(10);
	 for (unsigned int d = 0; d < minPerShift.size(); d++) {
			logStream << Tools::itoa(minPerSkill[d][i])+"/"+Tools::itoa(optPerSkill[d][i]);
		}
		logStream.endl();
	}
  logStream.endl();


	logStream <<  "# Highest demand per skill for one shift" << std::endl;
  logStream << "#";
	for (unsigned int d = 0; d < minPerSkill.size(); d++) {
	 logStream << "WD"+Tools::itoa(d);
	}
	logStream.endl();
	for (int i = 0; i < pScen->nbSkills_; i++)	{
	 logStream << "# SK"+Tools::itoa(i)+":";
		for (unsigned int d = 0; d < minPerShift.size(); d++) {
			logStream << Tools::itoa(minHighestPerSkill[d][i])+"/"+Tools::itoa(optHighestPerSkill[d][i]);
		}
		logStream.endl();
	}
  logStream.endl();

logStream.setWidth(12);
logStream <<  "# Agregate indicators" << std::endl;
logStream <<  "# Total understaffing without unavoidable penalty" << std::endl;
logStream << "#";
for (unsigned int d = 0; d < minPerShift.size(); d++) {
  logStream << "WD"+Tools::itoa(d);
}
logStream << "Average" << "Std dev" << std::endl;
logStream << "#";
double averageMin=0, averageOpt=0, stdDevMin=0, stdDevOpt=0;
for (unsigned int d = 0; d < minPerShift.size(); d++) {
  logStream << Tools::itoa(minTotal[d]-pSolver->maxTotalStaffNoPenalty_)+"/"+Tools::itoa(optTotal[d]-pSolver->maxTotalStaffNoPenalty_);
  averageMin += (minTotal[d]-pSolver->maxTotalStaffNoPenalty_)/(double)minPerShift.size();
  averageOpt += (optTotal[d]-pSolver->maxTotalStaffNoPenalty_)/(double)minPerShift.size();
  stdDevMin += pow(minTotal[d]-pSolver->maxTotalStaffNoPenalty_,2);
  stdDevOpt += pow(optTotal[d]-pSolver->maxTotalStaffNoPenalty_,2);
}
stdDevMin = sqrt(stdDevMin/(double)minPerShift.size()-pow(averageMin,2));
stdDevOpt = sqrt(stdDevOpt/(double)minPerShift.size()-pow(averageOpt,2));
logStream << Tools::itoa(averageMin)+"/"+Tools::itoa(averageOpt);
logStream << Tools::itoa(stdDevMin)+"/"+Tools::itoa(stdDevOpt);
logStream.endl();
logStream.endl();

logStream <<  "# Total understaffing with average work" << std::endl;
logStream << "#";
for (unsigned int d = 0; d < minPerShift.size(); d++) {
  logStream << "WD"+Tools::itoa(d);
}
logStream << "Average" << "Std dev" << std::endl;
logStream << "#";
averageMin=0, averageOpt=0, stdDevMin=0, stdDevOpt=0;
for (unsigned int d = 0; d < minPerShift.size(); d++) {
  logStream << Tools::itoa(minTotal[d]-pSolver->maxTotalStaffAvgWork_)+"/"+Tools::itoa(optTotal[d]-pSolver->maxTotalStaffAvgWork_);
  averageMin += (minTotal[d]-pSolver->maxTotalStaffAvgWork_)/(double)minPerShift.size();
  averageOpt += (optTotal[d]-pSolver->maxTotalStaffAvgWork_)/(double)minPerShift.size();
  stdDevMin += pow(minTotal[d]-pSolver->maxTotalStaffAvgWork_,2);
  stdDevOpt += pow(optTotal[d]-pSolver->maxTotalStaffAvgWork_,2);
}
stdDevMin = sqrt(stdDevMin/(double)minPerShift.size()-pow(averageMin,2));
stdDevOpt = sqrt(stdDevOpt/(double)minPerShift.size()-pow(averageOpt,2));
logStream << Tools::itoa(averageMin)+"/"+Tools::itoa(averageOpt);
logStream << Tools::itoa(stdDevMin)+"/"+Tools::itoa(stdDevOpt);
logStream.endl();
logStream.endl();

logStream.setPrecision(1);
logStream <<  "# Understaffing per skill without unavoidable penalty" << std::endl;
logStream << "#";
for (unsigned int d = 0; d < minPerSkill.size(); d++) {
  logStream << "WD"+Tools::itoa(d);
}
logStream << "Average" << "Std dev" << std::endl;
for (int i = 0; i < pScen->nbSkills_; i++)	{
  logStream << "# " + pScen->intToSkill_[i] +":";
  averageMin=0, averageOpt=0, stdDevMin=0, stdDevOpt=0;
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
	 logStream << Tools::itoa(minPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i])+"/"
		+Tools::itoa(optPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i]);
	 averageMin += (minPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i])/(double)minPerShift.size();
	 averageOpt += (optPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i])/(double)minPerShift.size();
	 stdDevMin += pow(minPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i],2);
	 stdDevOpt += pow(optPerSkill[d][i]-pSolver->maxStaffPerSkillNoPenalty_[i],2);
  }
  stdDevMin = sqrt(stdDevMin/(double)minPerShift.size()-pow(averageMin,2));
  stdDevOpt = sqrt(stdDevOpt/(double)minPerShift.size()-pow(averageOpt,2));
  logStream << Tools::itoa(averageMin)+"/"+Tools::itoa(averageOpt);
  logStream << Tools::itoa(stdDevMin)+"/"+Tools::itoa(stdDevOpt);
  logStream.endl();
}
logStream.endl();

logStream <<  "# Understaffing per skill with average work" << std::endl;
logStream << "#";
for (unsigned int d = 0; d < minPerSkill.size(); d++) {
  logStream << "WD"+Tools::itoa(d);
}
logStream << "Average" << "Std dev" << std::endl;
for (int i = 0; i < pScen->nbSkills_; i++)	{
  logStream << "# " + pScen->intToSkill_[i] +":";
  averageMin=0, averageOpt=0, stdDevMin=0, stdDevOpt=0;
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
	 logStream << Tools::itoa(minPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i])+"/"
		+Tools::itoa(optPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i]);
	 averageMin += (minPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i])/(double)minPerShift.size();
	 averageOpt += (optPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i])/(double)minPerShift.size();
	 stdDevMin += pow(minPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i],2);
	 stdDevOpt += pow(optPerSkill[d][i]-pSolver->maxStaffPerSkillAvgWork_[i],2);
  }
  stdDevMin = sqrt(stdDevMin/(double)minPerShift.size()-pow(averageMin,2));
  stdDevOpt = sqrt(stdDevOpt/(double)minPerShift.size()-pow(averageOpt,2));
  logStream << Tools::itoa(averageMin)+"/"+Tools::itoa(averageOpt);
  logStream << Tools::itoa(stdDevMin)+"/"+Tools::itoa(stdDevOpt);
  logStream.endl();
}
logStream.endl();





	delete pScen;

}

//--------------------------------------------------------------------------
// Useful parsing functions
//

// Read a file stream until the separating character is met
// Store the characters read until the separating character in pStrRead
//
bool ReadWrite::readUntilChar(std::fstream *pFile, char separater, std::string *pStrRead) {
	char cTmp = 'A';

	// empty the title string if it is not
	//
	if (!pStrRead->empty())
		pStrRead->erase();

	// go through the file until the delimiter is met
	//
	if (pFile->good()) {
		cTmp = pFile->get();
	}
	while (cTmp != separater && pFile->good() )  {
		pStrRead->push_back(cTmp);
		cTmp = pFile->get();
	}

	if (!pFile->good())
		return false;

	return true;
}

// Read a file stream until the separating character is met
// Store the characters read until the separating character in pStrRead
//
bool ReadWrite::readUntilOneOfTwoChar(std::fstream *pFile, char separater1, char separater2, std::string *pStrRead) {
	char cTmp = 'A';

	// empty the title string if it is not
	//
	if (!pStrRead->empty())
		pStrRead->erase();

	// go through the file until the delimiter is met
	//
	if (pFile->good()) {
		cTmp = pFile->get();
	}
	while (cTmp != separater1 && cTmp != separater2 && pFile->good() )  {
		pStrRead->push_back(cTmp);
		cTmp = pFile->get();
	}

	if (!pFile->good())
		return false;

	return true;
}

// Reads a file stream until the end of the line is met
// Stores the line in pStrRead
//



// Checks if the string (sentence) ends with the given substring (word)
//
bool ReadWrite::strEndsWith(string sentence, string word){
	int lWord = word.length();
	int lSentence = sentence.length();
	if(lWord > lSentence)
		return false;
	else{
		string endOfSentence = sentence.substr(lSentence-lWord, lWord);
		return (!strcmp(word.c_str(),endOfSentence.c_str()));
	}
}

//--------------------------------------------------------------------------
