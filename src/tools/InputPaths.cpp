//
//  InputPaths.cpp
//  RosterDesNurses
//
//  Created by Jeremy Omer on 16/02/2016.
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//

#include "tools/InputPaths.h"

/******************************************************************************
* The instances of InputPaths contain the paths of the input files of the
* problem
*******************************************************************************/
InputPaths::InputPaths() {
	instance_ = "";
	scenario_ = "";
	history_ = "";
}

InputPaths::InputPaths(std::string dataDir, std::string instanceName,int historyIndex, std::vector<int> weekIndices,
		       std::string solutionPath,std::string logPath, std::string paramFile, double timeOut, int randSeed,
                       std::string SPType,int SPStrategy, int nThreads):
	instance_(instanceName), historyIndex_(historyIndex), weekIndices_(weekIndices),
	solutionPath_(solutionPath), logPath_(logPath), paramFile_(paramFile), randSeed_(randSeed), timeOut_(timeOut),
  SPType_(SPType), SPStrategy_(SPStrategy), nThreads_(nThreads) {

	// int nbWeeks = weekIndices.size();

  
	std::string instanceDir = dataDir + instanceName + "/";
	// initialize the scenario and history file names
	scenario_ = instanceDir + "Sc-" + instanceName + ".txt";
	history_ = instanceDir + "H0" + "-" + instanceName + "-" + std::to_string(historyIndex) + ".txt";

	// initialize the file names for each week demand
	for(int week: weekIndices){
		std::string path = instanceDir + "WD-" + instanceName + "-" + std::to_string(week) + ".txt";
		weeks_.push_back(path);
	}
}
