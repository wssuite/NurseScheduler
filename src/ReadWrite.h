
#ifndef _ReadWrite_h
#define _ReadWrite_h

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <limits.h>

#define TRUE 1
#define FALSE 0


//BADVALDBL used as a flag value (not initialised parameter) when there's no ambiguity.
#define BADVALDBL -666.0

#ifdef WIN32
	#ifndef NAN
		const unsigned long nan[2]={0xffffffff, 0x7fffffff};
		#define NAN (*(const double *) nan)
	#endif
#endif


//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scneario file and store the content in a Scenario instance
//
void readScenario(std::string strWeekFile, Scenario* pScenario);

// Read the Week file and store the content in a Scenario instance
//
void readWeek(std::string strWeekFile, Scenario* pScenario);


// Read the history file
//
void readHistory(std::string strHistoryFile, Scenario* pScenario);

// Read the input custom file and store the content in a Scenario instance
//
void readCustom(std::string strCustomInputFile, Scenario* pScenario);

//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// Methods that write the ouputs of the solver

// Write the solution file for the current week
//
void writeSolution(std::string strCustomOutputFile, Solution* pSolution);

// Write the output custom file from values in the scenario and the solution
// instances
//
void writeCustom(std::string strCustomOutputFile, Scenario* pScenario, Solution* pSolution);

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Useful parsing functions

// Read a file stream until the separating character is met
// Store the characters read until the separating character in pStrRead
//
bool readUntilChar(std::fstream *pFile, char separater, std::string *pStrRead);
//--------------------------------------------------------------------------


#endif
