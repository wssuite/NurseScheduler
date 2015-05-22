//
//  main.cpp
//  RosterDesNurses
//
//  Created by J������r������my Omer on 16/11/2014.
//  Copyright (c) 2014 J������r������my Omer. All rights reserved.
//

#include "main_test.h"
#include "MyTools.h"
#include "ReadWrite.h"
#include "Solver.h"
#include "Greedy.h"
#include <exception>
#include "Modeler.h"
#include "MasterProblem.h"
#include "StochasticSolver.h"


void main_test()
{
   //testFunction_Antoine();
   //testFunction_Jeremy();
	testFunction_Samuel();
}

int main(int argc, char** argv)
{
   std::cout << "Number of arguments= " << argc << std::endl;

   // Detect errors in the number of arguments
   //
   if (argc%2 != 1) {
      Tools::throwError("main: There should be an even number of arguments!");
   }
   else if (argc > 1 && (argc < 9 || argc > 15)) {
      Tools::throwError("main: There is either too many or not enough arguments!");
   }

   // During tests, only run some test methods
   //
   if (argc == 1) {
      std::cout << "Running the test methods..."<< std::endl;
      // Tests functions to check the functions one by one
      main_test();
   }
   // Nominal behavior of the executable, as required by INRCII
   //
   else {
      // Retrieve the file names in arguments
      //
      int narg = 1;
      string scenarioFile="", initialHistoryFile="", weekDataFile="", solutionFile="";
      string customInputFile="", customOutputFile="", randSeed="";

      while (narg < argc) {
         std::cout << "arg = " << argv[narg] << " " << argv[narg+1] << std::endl;
         // Attention usine ������ gaz: problem with the java simulator that add quote
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
            scenarioFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--his")) {
            initialHistoryFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--week")) {
            weekDataFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--sol")) {
            solutionFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--cusIn")) {
            customInputFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--cusOut")) {
            customOutputFile = str;
            narg += 2;
         }
         else if (!strcmp(argv[narg],"--rand")) {
            randSeed = str;
            narg += 2;
         }
         else {
            Tools::throwError("main: the argument does not match the expected list!");
         }
      }
      // Throw an error if a necessary input file is missing
      if ( scenarioFile.empty() || initialHistoryFile.empty()
            || weekDataFile.empty() || solutionFile.empty() ) {
         throw Tools::myException("A necessary file name is missing!",__LINE__);
      }

      // Read the input files
      //
      // if (!customInputFile.empty()) {
      //    vector<Demand*> demandHistory;
      //    demandHistory.push_back(pScen->pWeekDemand());
      //    int coWeek = ReadWrite::readCustom(customInputFile, pScen, demandHistory);
      // }

		unsigned found = solutionFile.find_last_of(".");
		string logFile = solutionFile.substr(0,found);

      // Solve the week
      solveOneWeek(scenarioFile, weekDataFile, initialHistoryFile, solutionFile, logFile);


      // Write the solution in the required output format
      //
      // if (!customOutputFile.empty()) {
      //    ReadWrite::writeCustom(customOutputFile,weekDataFile,customInputFile);
      // }
      // Todo: the method that writes the history file corresponding to the
      // solution
      // string outputHistoryFile("history-week");
      // outputHistoryFile += std::to_string(pScen->thisWeek()) + ".txt";
      // std::cout << "Output history file: " << outputHistoryFile << std::endl;
   }
}
