#include "ReadWrite.h"
#include "MyTools.h"
#include <iostream>
#include <fstream>
#include <streambuf>
#include <fstream>
#include <math.h>
#include <time.h>


//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
static Scenario* ReadWrite::readScenario(string fileName) {
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
	int nbWeeks, nbSkills, nbShifts;
	vector<string> intToSkill, intToShift;
	map<string,int> skillToInt, shiftToInt;
	vector<int> minConsShifts, maxConsShifts, nbForbiddenSuccessors;
	vector2D pForbiddenSuccessors;

	// fill the attributes of the scenario structure
	//
	Tools::readUntilChar(&file, '=', &title);
	while ( file.good() ) {
		if (!strcmp(title.c_str(), "SCENARIO")) {
			file >> name;
		}
		else if (!strcmp(title.c_str(), "WEEKS")) {
			file >> nbWeeks;
		}
		else if (!strcmp(title.c_str(), "SKILLS"))   {
			file >> nbSkills;
			for (int i = 0; i < nbSkills; i++) {
				file >> strTmp;
				intToSkill.push_back(strTmp);
				skillToInt[strTmp] = i;
			}
		}
		else if (!strcmp(title.c_str(), "SHIFT_TYPES"))   {
			file >> nbShifts;
			for (int i = 0; i < nbShifts; i++) {
				file >> strTmp;
				intToShift.push_back(strTmp);
				shiftToInt[strTmp] = i;
				Tools::readUntilChar(&file, '(', &title);
				file >> intTmp;
				minConsShifts.push_back(intTmp);
				Tools::readUntilChar(&file, ',', &title);
				file >> intTmp;
				maxConsShifts.push_back(intTmp);
				file.getline(charTmp, 256);
			}
			file.getline(charTmp, 256);
			file >> title;
			if (!strcmp(title.c_str(), "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS")) {
				std::cout << "While trying to read " << fileName << std::endl;
				std::cout << "It is presenty reading " << title << std::endl;
				Tools::throwError("The reader has no idea what it is reading!");
			}
			for (int i =0; i < nbShifts; i++) {
				file >> strTmp;
				file >> intTmp;
				nbForbiddenSuccessors.push_back(intTmp);
				// terminer l'affectation du nombre de successeurs
			}

			Scenario* pScenario =
				new Scenario(name, nbWeeks, nbSkills,  intToSkill, skillToInt,
				nbShifts, intToShift, shiftToInt,
				minConsShifts,  maxConsShifts,
				nbForbiddenSuccessors, pForbiddenSuccessors)

		}
	}
	return pScenario;
}


//--------------------------------------------------------------------------
// Useful parsing functions
//

// Read a file stream until the separating character is met
// Store the characters read until the separating character in pStrRead
//
static bool ReadWrite::readUntilChar(std::fstream *pFile, char separater, std::string *pStrRead) {
	char cTmp = 'A';

	// empty the title string if it is not
	//
	if (!pStrRead->empty())
		pStrRead->erase();

	// go through the file until the delimiter is met
	//
	if (file->good()) {
		cTmp = file->get();
	}
	while (cTmp != separater && file->good() )  {
		pStrRead->push_back(cTmp);
		cTmp = file->get();
	}

	if (!file->good())
		return false;

	return true;
}

//--------------------------------------------------------------------------
