
#ifndef _ReadWrite_h
#define _ReadWrite_h

#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <limits.h>

#include "Scenario.h"

using std::string;
using std::cout;
using std::endl;

//BADVALDBL used as a flag value (not initialised parameter) when there's no ambiguity.
#define BADVALDBL -666.0

#ifdef WIN32
	#ifndef NAN
		const unsigned long nan[2]={0xffffffff, 0x7fffffff};
		#define NAN (*(const double *) nan)
	#endif
#endif


//--------------------------------------------------------------------------
//
//  C l a s s   R e a d W r i t e
//
//  Contains the (static) functions to read the input and write the output
//
//--------------------------------------------------------------------------

class ReadWrite{

// All functions in this class shall be public
public:

	//--------------------------------------------------------------------------
	// Methods that read all the input files and store the content in the
	// input scenario instance
	//
	// Read the scenario file and store the content in a Scenario instance
	//
	static Scenario* readScenario(std::string strScenarioFile);

	//Read several week files and strore the content in one demand and one preference
	//
   static Demand* readWeeks(std::vector<std::string> strWeekFiles, Scenario* pScenario);
	// Read the Week file and store the content in a Scenario instance
	//
	static Demand* readWeek(std::string strWeekFile, Scenario* pScenario);

	// Read the history file
	//
	static void readHistory(std::string strHistoryFile, Scenario* pScenario);

	// Read the input custom file and store the content in a Scenario instance
	//
	static void readCustom(std::string strCustomInputFile, Scenario* pScenario);

	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	// Methods that write the ouputs of the solver

	// Write the solution file for the current week
	//
	// void writeSolution(std::string strCustomOutputFile, Solution* pSolution);
	//
	// // Write the output custom file from values in the scenario and the solution
	// // instances
	// //
	// void writeCustom(std::string strCustomOutputFile, Scenario* pScenario, Solution* pSolution);

	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	// Useful parsing functions
	// Read a file stream until the separating character (or one of them) is met
	// Store the characters read until the separating character in pStrRead
	//
	static bool readUntilChar(std::fstream *pFile, char separater, std::string *pStrRead);
	static bool readUntilOneOfTwoChar(std::fstream *pFile, char separater1, char separater2, std::string *pStrRead);

	// Checks if the string (sentence) ends with the given substring (word)
	//
	static bool strEndsWith(string sentence, string word);

	//--------------------------------------------------------------------------
};

#endif
