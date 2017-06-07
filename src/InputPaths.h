//
//  InputPaths.h
//  RosterDesNurses
//
//  Created by Jeremy Omer on 16/02/2016.
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//


#ifndef __InputPaths__
#define __InputPaths__

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>


// The instances of this class contain the paths of the input paths of the
// problem
class InputPaths{
public:
	InputPaths();
	InputPaths(std::string dataDir, std::string instanceName,int historyIndex, std::vector<int> weekIndices, 
		std::string solutionPath="",std::string logPath="", std::string paramFile="" , double timeOut=3600.0, int randSeed=0);

protected:
	std::string instance_;
	std::string scenario_;
	int historyIndex_;
	std::string history_;
	std::vector<std::string> weeks_;
	std::vector<int> weekIndices_;
	std::string solutionPath_="";
	std::string logPath_="";
	std::string paramFile_="";
	int randSeed_=0;
	double timeOut_=3600;

public:
	// get/set attributes
	std::string instance() const {return instance_;}
	std::string scenario() const {return scenario_;}
	inline void scenario(std::string scenario) {scenario_ = scenario;}

	int historyIndex() const {return historyIndex_;}
	std::string history() const {return history_;}
	inline void history(std::string history) {history_ = history;}

	std::vector<std::string> weeks() const {return weeks_;}
	std::string week(int w) const {return weeks_[w];}
	int weekIndex(int w) const {return weekIndices_[w];}
	int nbWeeks() const {return weeks_.size();}
	inline void addWeek(std::string week) {weeks_.push_back(week);}

	std::string paramFile() {return paramFile_;}
	inline void paramFile(std::string file) {paramFile_=file;}
	std::string solutionPath() {return solutionPath_;}
	inline void solutionPath(std::string path) {solutionPath_=path;}
	std::string logPath() {return logPath_;}
	inline void logPath(std::string path) {logPath_=path;}

	int randSeed() {return randSeed_;}
	inline void randSeed(int seed) {randSeed_ =  seed;}
  double timeOut() {return timeOut_;}
  inline void timeOut(double t) {timeOut_= t;}
};

#endif
