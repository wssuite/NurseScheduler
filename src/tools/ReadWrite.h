
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

#include "data/Scenario.h"
#include "solvers/StochasticSolver.h"


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
	static PScenario readScenario(std::string strScenarioFile);

	//Read several week files and strore the content in one demand and one preference
	//
   static PDemand readWeeks(std::vector<std::string> strWeekFiles, PScenario pScenario);
	// Read the Week file and store the content in a Scenario instance
	//
	static void readWeek(std::string strWeekFile, PScenario pScenario, PDemand* pDemand, PPreferences* pPref);

	// Read the history file
	//
	static void readHistory(std::string strHistoryFile, PScenario pScenario);

	// Read the input custom file
	// Store the result in a vector of historical demands and return the number of treated weeks
	//
	static int readCustom(std::string strCustomInputFile, PScenario pScenario, std::vector<PDemand>& demandHistory);

	// Read the options of the stochastic and ot the other solvers
	//
	static std::string readStochasticSolverOptions(std::string strOptionFile, StochasticSolverOptions& options);
	static std::string readSolverOptions(std::string strOptionFile, SolverParam& options);

	// Read the solution from multiple week solution files
	//
	static std::vector<Roster> readSolutionMultipleWeeks(std::vector<std::string> strWeekSolFiles, PScenario pScenario);

	// Print the main characteristics of all the demands of an input directory
	// This is done to find some invariant properties among demands
	static void compareDemands(std::string inputDir, std::string logFile);

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
	static void writeCustom(std::string stdCustomOutputFile, std::string strWeekFile, std::string strCustomInputFile="");

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
	static bool strEndsWith(std::string sentence, std::string word);

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
