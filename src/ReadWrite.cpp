#include "ReadWrite.h"
#include "MyTools.h"
#include "Scenario.h"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>


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
		std::cout << "While trying to read fucking " << fileName << std::endl;
		std::cout << "The input file was not opened properly!" << std::endl;

		throw Tools::myException("The input file was not opened properly!",__LINE__);
	}

	string title;
	char charTmp[256];
	string strTmp;
	int intTmp;
	// declare the attributes that will intialize the Scenario instance
	//
	string name;
	int nbWeeks=-1, nbSkills=-1, nbShifts=-1, nbContracts=-1, nbNurses=-1;
	vector<string> intToSkill, intToShift, intToContract;
	map<string,int> skillToInt, shiftToInt, nurseNameToInt;
	vector<int> minConsShifts, maxConsShifts, nbForbiddenSuccessors;
	vector2D forbiddenSuccessors;
	map<string,Contract*> contracts;
	vector<Nurse> theNurses;


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
			nbShifts = intTmp+1;

			// IMPORTANT : INSERT REST SHIFT !!!!!!
			// It is given 0 and 99 as bounds so that they never perturbate the cost
			//
			intToShift.push_back(REST_SHIFT);
			shiftToInt.insert(pair<string,int>(REST_SHIFT,0));
			minConsShifts.push_back(0);
			maxConsShifts.push_back(99);

			// Other shift types
			//
			for(int i=1; i<nbShifts; i++){
				// Name
				file >> strTmp;
				intToShift.push_back(strTmp);
				shiftToInt.insert(pair<string,int>(strTmp,i));
				readUntilChar(&file,'(',&strTmp);
				// Min consecutive
				file >> intTmp;
				minConsShifts.push_back(intTmp);
				readUntilChar(&file, ',', &strTmp);
				// Max consecutive
				file >> intTmp;
				maxConsShifts.push_back(intTmp);
				readUntilChar(&file,'\n',&strTmp);
			}


			// Forbidden successions
			//
			for(int i=0; i<nbShifts; i++){
				vector<int> v;
				forbiddenSuccessors.push_back(v);
				nbForbiddenSuccessors.push_back(0);
			}
			while(!strEndsWith(title,"FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
				file >> title;
			// Reading all lines
			for(int i=1; i<nbShifts; i++){
				// Which current shift ?
				string currentShift;
				file >> currentShift;
				int currentShiftId = shiftToInt.at(currentShift);
				// How many forbidden after it ?
				file >> intTmp;
				nbForbiddenSuccessors[currentShiftId] = intTmp;
				// Which ones are forbidden ?
				for(int j=0; j<nbForbiddenSuccessors[currentShiftId]; j++){
					file >> strTmp;
					forbiddenSuccessors[currentShiftId].push_back(shiftToInt.at(strTmp));
				}
				readUntilChar(&file,'\n',&strTmp);

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

				Contract * pContract = new Contract (contractName, minDays, maxDays, minConsWork, maxConsWork, minConsRest, maxConsRest, maxWeekends, isTotalWeekend);
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

	// Check that all fields were initialized before initializing the scenario
	//
	if ( nbWeeks==-1 || nbSkills==-1 || nbShifts==-1 || nbContracts==-1 || nbNurses==-1 ) {
		Tools::throwError("In readScenario: missing fields in the initialization");
	}

	return new Scenario(name, nbWeeks, nbSkills, intToSkill, skillToInt, nbShifts,
		intToShift, shiftToInt, minConsShifts, maxConsShifts,
		nbForbiddenSuccessors,forbiddenSuccessors,
		nbContracts, intToContract, contracts, nbNurses, theNurses, nurseNameToInt) ;
}

// Read the Week file and store the content in a Scenario instance
//
Demand* ReadWrite::readWeek(std::string strWeekFile, Scenario* pScenario){
	// open the file
	std::fstream file;
	std::cout << "Reading " << strWeekFile << std::endl;
	file.open(strWeekFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read fucking " << strWeekFile << std::endl;
		std::cout << "The input file was not opened properly!" << std::endl;
		throw Tools::myException("The input file was not opened properly!",__LINE__);
	}

	string title;
	string strTmp;
	int intTmp;

	// declare the attributes to be updated in the Scenario*
	//
	string weekName;
	vector3D minWeekDemand;
	vector3D optWeekDemand;
	int nbShiftOffRequests;
	Preferences weekPreferences;


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
			Tools::initVector3D(&minWeekDemand, 7, pScenario->nbShifts_, pScenario->nbSkills_);
			Tools::initVector3D(&optWeekDemand, 7, pScenario->nbShifts_, pScenario->nbSkills_);

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
			Preferences pref (pScenario->nbNurses_, 7, pScenario->nbShifts_);
			// Temporary vars
			string nurseName, shift, day;
			int nurseId, shiftId, dayId;
			file >> nbShiftOffRequests;
			for (int i=0; i<nbShiftOffRequests; i++){
				file >> nurseName;
				file >> shift;
				file >> day;
				nurseId = pScenario->nurseNameToInt_.at(nurseName);
				dayId = Tools::dayToInt(day);

				if(shift == "Any")
					pref.addDayOff(nurseId, dayId);
				else {
					shiftId = pScenario->shiftToInt_.at(shift);
					pref.addShiftOff(nurseId, dayId, shiftId);
				}
				weekPreferences = pref;
			}
		}
	}

	// Define a new instance of demand
	Demand* pDemand = new Demand(7, 0, pScenario->nbShifts_,pScenario->nbSkills_,
	minWeekDemand, optWeekDemand);

	// Now, add all these objects to the Scenario
	pScenario->setWeekName(weekName);
	pScenario->setWeekDemand(pDemand);
	pScenario->setTNbShiftOffRequests(nbShiftOffRequests);
	pScenario->setWeekPreferences(weekPreferences);

	// return the demand
	return pDemand;

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
	int intTmp;

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
				string nurseName, shiftName;
				int nurseId, shiftId, totalDaysWorked, totalWeekendsWorked, consDaysWorked, consShiftWorked, consRest, consShifts;
				file >> nurseName;
				nurseId = pScenario->nurseNameToInt_.at(nurseName);
				file >> totalDaysWorked;
				file >> totalWeekendsWorked;
				file >> shiftName;
				shiftId = pScenario->shiftToInt_.at(shiftName);
				file >> consShiftWorked;
				file >> consDaysWorked;
				file >> consRest;

				consShifts = (shiftId == 0) ? consRest : consShiftWorked;
				State nurseState (0, totalDaysWorked, totalWeekendsWorked,
					consDaysWorked, consShifts, consRest, shiftId);
				initialState.push_back(nurseState);
			}
		}
	}
	pScenario->setThisWeek(thisWeek);
	pScenario->setInitialState(initialState);
}


// Read the input custom file and store the content in a Scenario instance
//
void ReadWrite::readCustom(std::string strCustomInputFile, Scenario* pScenario) {

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
