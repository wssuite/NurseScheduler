
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
#include "StochasticSolver.h"

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
   static Demand* readWeeks(vector<std::string> strWeekFiles, Scenario* pScenario);
	// Read the Week file and store the content in a Scenario instance
	//
	static void readWeek(string strWeekFile, Scenario* pScenario, Demand** pDemand, Preferences** pPref);

	// Read the history file
	//
	static void readHistory(string strHistoryFile, Scenario* pScenario);

	// Read the input custom file
	// Store the result in a vector of historical demands and return the number of treated weeks
	//
	static int readCustom(string strCustomInputFile, Scenario* pScenario, vector<Demand*>& demandHistory);

	// Read the options of the stochastic and ot the other solvers
	//
	static string readStochasticSolverOptions(string strOptionFile, StochasticSolverOptions& options);
	static string readSolverOptions(string strOptionFile, SolverParam& options);


	// Print the main characteristics of all the demands of an input directory
	// This is done to find some invariant properties among demands
	static void compareDemands(string inputDir, string logFile);

	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	// Methods that write the ouputs of the solver

	// Write the solution file for the current week
	//
	// void writeSolution(std::string strCustomOutputFile, Solution* pSolution);
	//

	// Write the output custom file from values in the scenario and the solution
	// instances
	//
	static void writeCustom(string stdCustomOutputFile, string strWeekFile, string strCustomInputFile="");

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

	// writes a string in a stream with a constant number of character
	//
	template<typename T>
	static void writeConstantWidth(std::ostream &out, int width, const T output);


	//--------------------------------------------------------------------------
};

// writes a string in a stream with a constant number of character
//
template<typename T>
void ReadWrite::writeConstantWidth(std::ostream &out, int width, const T output) {
  char buffer[100];
  int cx;

  cx = snprintf(buffer, 100, "%s", output);
  if (cx > width){
    std::cout << "Warning: the string is larger than the width!" << std::endl;
    width = cx;
  }

  std::fill_n(buffer+cx, width-cx, ' ');

  out.write(buffer, width);
}

#endif
