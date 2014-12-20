//
//  Scenario.cpp
//  Project:RosterDesNurses
//
//  Created by Jérémy Omer on 18/12/2013.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//

#include <sstream>
#include <string>

#include "Scenario.h"
#include "MyTools.h"
#include "Nurse.h"

using std::cout;
using std::endl;

Scenario::~Scenario(){

}

// Display methods: toString + override operator<< (easier)
//
string Scenario::toString(){
	//return "\ntoto\n";
	//std::stringstream std::cout;
	std::cout << "# Scenario_name = [" << name_ << "]" << std::endl;
	std::cout << "# Number_of_weeks = [" << nbWeeks_ << "]" << std::endl;
	std::cout << "# Number_of_skills = [" << nbSkills_ << "]" << std::endl;
	std::cout << "# List of skills = "; for(int i=0; i<nbSkills_; i++){
		std::cout << "[" << i << ":" << intToSkill_[i] << "] ";
	}
	std::cout << std::endl;
	std::cout << "# Number_of_shifts = [" << nbShifts_ << "]" << std::endl;

	cout << endl;
	std::cout << intToShift_[0] << endl;
	std::cout << intToShift_[1] << endl;
	std::cout << intToShift_[2] << endl;
	std::cout << intToShift_[3] << endl;

	std::cout << "More to be printed in the end... debug in progress" << endl;

	/*
	cout << endl;
	cout << minConsShifts_[0] << endl;;
	cout << minConsShifts_[1] << endl;;
	cout << minConsShifts_[2] << endl;;
	cout << minConsShifts_[3] << endl;;
	*/

	// @TODO : finish here to debug... probably need to change the data structure either to (nonconst) vector<int> or maybe to const vector<int>* (the former seems better)

	/*
	std::cout << "# List of shifts = [";
	for(int i=0; i<nbShifts_; i++){
		std::cout << "[" << intToShift_[i] << "|" << minConsShifts_[i] << "<" << maxConsShifts_[i] << "]";
	}
	std::cout << std::endl;
	std::cout << "# Forbidden successions : " << std::endl;
	for(int i=0; i<nbShifts_; i++){
		std::cout << "#   | " << intToShift_[i] << " [" << nbForbiddenSuccessors_[i] << "]  ->  ";
		for(int j=0; j<nbForbiddenSuccessors_[i]; j++){
			std::cout << "[" << intToShift_[j] << "][" << j << "]  ";
		}
		std::cout << std::endl;
	}
	std::cout << "# Contracts : [" << nbContracts_ << "]" << std::endl;
	for(map<string,Contract>::const_iterator itC = contracts_.begin(); itC != contracts_.end(); ++itC){
		std::cout << "#   | " << (itC->second) << std::endl;
	}
	std::cout << "# Nurses : [" << nbNurses_ << "]" << std::endl;
	for(int i=0; i<nbNurses_; i++){
		std::cout << theNurses_[i] << std::endl;
	}
	//return std::cout.str();
	*/
	return "";
}

