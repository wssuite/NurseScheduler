//
//  main.cpp
//  CftSolver
//
//  Created by Jérémy Omer on 16/01/2014.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//

//#include "Scenario.h"
//#include "Nurse.h"
#include "ReadWrite.h"
#include "Solver.h"
#include "MyTools.h"

// Function for testing parts of the code (Antoine)
void testFunction_Antoine(){

}

// Function for testing parts of the code (Jeremy)
void testFunction_Jeremy(){

	// Time the complete execution of the algorithm
	Tools::Timer* timertotal = new Tools::Timer();
	timertotal->init();
	timertotal->start();

	// Create a log file
	string logFile = "../logfiles/test.log";
	Tools::LogOutput logStream(logFile);


	// Read the input data from files
	Scenario* pScen = ReadWrite::readScenario("../datasets/n030w4/Sc-n030w4.txt");
	Demand* pWeekDemand = ReadWrite::readWeek("../datasets/n030w4/WD-n030w4-1.txt", pScen);
	ReadWrite::readHistory("../datasets/n030w4/H0-n030w4-0.txt",pScen);

	// Check that the scenario was read properly
	//
	// logStream << *pScen << std::endl;
	logStream << pScen->toString() << std::endl;
	logStream << pWeekDemand->toString(true) << std::endl;

	// Write the aggregate information on the demand
	//


	// Write aggregate information on the cover capacity of the staff
	// (TBD)

	// Instantiate the solver class as a test
	//
	Solver* pSolverTest =
	new Solver(pScen, pWeekDemand,	pScen->pWeekPreferences(), pScen->pInitialState());
	pSolverTest->preprocessTheNurses();

	// Display the total time spent in the algorithm
	timertotal->stop();
	logStream.print("Total time spent in the algorithm : ");
	logStream.print(timertotal->dSinceInit());
	logStream.print("\n");

	// free the allocated pointers
	//
	delete timertotal;
	delete pWeekDemand;
	delete pScen;
	delete pSolverTest;


}

// Function for testing parts of the code (Samuel)
void testFunction_Samuel(){

	Tools::Timer* timertest = new Tools::Timer();
	timertest->init();
	timertest->start();

	Scenario * s = ReadWrite::readScenario("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/Sc-n030w4.txt");
	Demand* pWeekDemand = ReadWrite::readWeek("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/WD-n030w4-1.txt",s);
	ReadWrite::readHistory("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/H0-n030w4-0.txt",s);

	timertest->stop();

	string logFile = "../logfiles/samuel_test.log";
	Tools::LogOutput logStream(logFile);

  // RqJO : attention, j'ai enlevé ta surcharge de << parce qu'elle me faisait
	// des segfaults
	logStream << s->toString() << std::endl;
	logStream.print("Total time spent in the algorithm : ");
	logStream.print(timertest->dSinceInit());
	logStream.print("\n");

}

int main(int argc, char** argv)
{

	// Tests functions to check the functions one by one
	// testFunction_Antoine();
	testFunction_Jeremy();
	// testFunction_Samuel();

}
