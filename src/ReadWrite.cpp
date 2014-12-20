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
		std::cout << "While trying to read " << fileName << std::endl;
		Tools::throwError("The input file was not opened properly!");
	}

	string title;
	char charTmp[256];
	string strTmp;
	int intTmp;
	// declare the attributes that will intialize the Scenario instance
	//
	string name;
	int nbWeeks, nbSkills, nbShifts, nbContracts, nbNurses;
	vector<string> intToSkill, intToShift;
	map<string,int> skillToInt, shiftToInt;
	vector<int> minConsShifts, maxConsShifts, nbForbiddenSuccessors;
	vector2D forbiddenSuccessors;
	map<string,Contract> contracts;
	vector<Nurse> theNurses;

	// Booleans that say if the field has already been completely treated
	//
	bool isScenarioRead = false;
	bool isWeeksRead = false;
	bool isSkillsRead = false;
	bool isShifTypesRead = false;
	bool isForbiddenShiftSuccessionRead = false;
	bool isContractsRead = false;
	bool isNursesRead = false;


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

			// Shift types
			//
			file >> nbShifts;
			for(int i=0; i<nbShifts; i++){
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
			for(int i=0; i<nbShifts; i++){ vector<int> v; forbiddenSuccessors.push_back(v);}
			while(!strEndsWith(title,"FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
				readUntilChar(&file,'\n',&title);
			// Reading all lines
			for(int i=0; i<nbShifts; i++){
				// Which current shift ?
				string currentShift;
				file >> currentShift;
				// How many forbidden after it ?
				file >> intTmp;
				nbForbiddenSuccessors.push_back(intTmp);
				// Which ones are forbidden ?
				for(int j=0; j<nbForbiddenSuccessors[shiftToInt.at(currentShift)]; j++){
					file >> strTmp;
					forbiddenSuccessors[i].push_back(shiftToInt.at(strTmp));
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

				Contract c (contractName, minDays, maxDays, minConsWork, maxConsWork, minConsRest, maxConsRest, maxWeekends, isTotalWeekend);
				contracts.insert(pair<string,Contract>(contractName,c));

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
				Contract * c;
				c = &(contracts.at(contractName));
				Nurse nurse (i, nurseName, nbSkills, skills, c);
				theNurses.push_back(Nurse(i, nurseName, nbSkills, skills, c));
			}

		}

		//getchar();



	}

	cout << endl << endl << endl;
	cout << "# Scenario_name = [" << name << "]" << endl;
	cout << "# Number_of_weeks = [" << nbWeeks << "]" << endl;
	cout << "# Number_of_skills = [" << nbSkills << "]" << endl;
	cout << "# List of skills = "; for(int i=0; i<nbSkills; i++) {cout << "[" << i << ":" << intToSkill[i] << "] ";} cout << endl;
	cout << "# Number_of_shifts = [" << nbShifts << "]" << endl;
	cout << "# List of shifts = ["; for(int i=0; i<nbShifts; i++) {cout << "[" << intToShift[i] << "|" << minConsShifts[i] << "<" << maxConsShifts[i] << "]";} cout << endl;
	cout << "# Forbidden successions : " << endl;
	for(int i=0; i<nbShifts; i++){
		cout << "#   | " << intToShift[i] << " [" << nbForbiddenSuccessors[i] << "]  ->  ";
		for(int j=0; j<nbForbiddenSuccessors[i]; j++){
			cout << "[" << intToShift[j] << "][" << j << "]  ";
		}
		cout << endl;
	}
	cout << "# Contracts : [" << nbContracts << "]" << endl;
	for(map<string,Contract>::iterator itC = contracts.begin(); itC != contracts.end(); ++itC){
		cout << "#   | " << (itC->second) << endl;
	}
	cout << "# Nurses : [" << nbNurses << "]" << endl;
	for(int i=0; i<nbNurses; i++){
		cout << theNurses[i] << endl;
	}

	cout << endl << endl << endl;

	// VERIFICATION :
	Scenario * scenario  = new Scenario(
			name, nbWeeks, nbSkills, intToSkill, skillToInt, nbShifts, intToShift, shiftToInt, minConsShifts, maxConsShifts,
			nbForbiddenSuccessors, forbiddenSuccessors, nbContracts, contracts, nbNurses, theNurses) ;
	std::cout << endl << *scenario << endl;
	return scenario;
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
