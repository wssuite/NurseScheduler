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
#include "MyTools.h"

// Function for testing parts of the code (Antoine)
void testFunction_Antoine(){

}

// Function for testing parts of the code (Jeremy)
void testFunction_Jeremy(){
	Tools::Timer* timertest = new Tools::Timer();

	timertest->init();
	timertest->start();

	double x;
	for (int i = 0; i < 100000000;i++) 	{
	 x = 1005*190;
	}

	timertest->stop();

	string logFile = "../logfiles/test.log";
	Tools::LogOutput logStream(logFile);

	logStream << "Total time spent in the algorithm : ";

	logStream.switchToFormatted(25);
	logStream << timertest->dSinceInit();

	logStream.switchToUnformatted();
	logStream << " seconds" << std::endl;

	logStream.print("Total time spent in the algorithm : ");
	logStream.print(timertest->dSinceInit());
	logStream.print("\n");
}

// Function for testing parts of the code (Samuel)
void testFunction_Samuel(){

	Tools::Timer* timertest = new Tools::Timer();
	timertest->init();
	timertest->start();

	Scenario * s = ReadWrite::readScenario("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/Sc-n030w4.txt");
	ReadWrite::readWeek("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/WD-n030w4-1.txt",s);
	ReadWrite::readHistory("/home/samuel/Dropbox/Nurse Rostering Competition/Data/datasets_txt/n030w4/H0-n030w4-0.txt",s);

	timertest->stop();

	string logFile = "../logfiles/samuel_test.log";
	Tools::LogOutput logStream(logFile);

	logStream << *s << std::endl;;
	logStream.print("Total time spent in the algorithm : ");
	logStream.print(timertest->dSinceInit());
	logStream.print("\n");

}

int main(int argc, char** argv)
{

	// Tests functions to check the functions one by one
	testFunction_Antoine();
	testFunction_Jeremy();
	testFunction_Samuel();

}
