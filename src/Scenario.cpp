//
//  Scenario.cpp
//  Project:RosterDesNurses
//
//  Created by J��r��my Omer on 18/12/2013.
//  Copyright (c) 2014 J��r��my Omer. All rights reserved.
//

#include <sstream>
#include <string>

#include "Scenario.h"
#include "MyTools.h"
#include "Nurse.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S c e n a r i o
//
//  Class that contains all the attributes describing the scenario
//
//-----------------------------------------------------------------------------

// Constructor and destructor
//
Scenario::Scenario(string name, int nbWeeks,
		int nbSkills, vector<string> intToSkill, map<string,int> skillToInt,
		int nbShifts, vector<string> intToShift, map<string,int> shiftToInt,
		vector<int> minConsShifts, vector<int> maxConsShifts,
		vector<int> nbForbiddenSuccessors, vector2D forbiddenSuccessors,
		int nbContracts, vector<string> intToContract, map<string,Contract*> contracts,
		int nbNurses, vector<Nurse> theNurses, map<string,int> nurseNameToInt) :
		name_(name), nbWeeks_(nbWeeks),
		nbSkills_(nbSkills), intToSkill_(intToSkill), skillToInt_(skillToInt),
		nbShifts_(nbShifts), intToShift_(intToShift), shiftToInt_(shiftToInt),
		minConsShifts_(minConsShifts), maxConsShifts_(maxConsShifts),
		nbForbiddenSuccessors_(nbForbiddenSuccessors), forbiddenSuccessors_(forbiddenSuccessors),
		nbContracts_(nbContracts), intToContract_(intToContract), contracts_(contracts),
		nbNurses_(nbNurses), theNurses_(theNurses), nurseNameToInt_(nurseNameToInt),
		nbPositions_(0), nbShiftOffRequests_(0),
		pWeekDemand_(0){

	// To make sure that it is modified later when reading the history data file
	//
	thisWeek_ = -1;
	nbWeeksLoaded_ = 1;

	// Preprocess the vector of nurses
	// This creates the positions
	//
	this->preprocessTheNurses();
}

Scenario::~Scenario(){
	// delete the contracts
	for(map<string,Contract*>::const_iterator itC = contracts_.begin(); itC != contracts_.end(); ++itC){
		delete (itC->second);
	}
   //delete pPositions_;
	for(Position* position: pPositions_)
	   delete position;

	delete pWeekDemand_;
}

// return true if the shift shNext is a forbidden successor of sh
//
bool Scenario::isForbiddenSuccessor(int shNext, int shLast) {
	if (shLast <= 0) return false;

	for (int i = 0; i < nbForbiddenSuccessors_[shLast]; i++) {
		if (shNext == forbiddenSuccessors_[shLast][i])  {
			return true;
		}
	}
	return false;
}

// update the scenario to treat a new week
//
void Scenario::updateNewWeek(Demand* pDemand, Preferences& preferences, vector<State> &initialStates) {

	// delete the current demand
	delete pWeekDemand_;

	// set the demand, preferences and initial states
	this->linkWithDemand(pDemand);
	this->linkWithPreferences(preferences);
	this->setInitialState(initialStates);

	// update the index of the week
	thisWeek_++;

}

//------------------------------------------------
// Display functions
//------------------------------------------------

// Display methods: toString + override operator<< (easier)
//
string Scenario::toString(){
	std::stringstream rep;
	rep << "############################################################################" << std::endl;
	rep << "##############################    Scenario    ##############################" << std::endl;
	rep << "############################################################################" << std::endl;
	rep << "# " << std::endl;
	rep << "# NAME             \t= " << name_ << std::endl;
	rep << "# NUMBER_OF_WEEKS  \t= " << nbWeeks_ <<std::endl;
	rep << "# " << std::endl;
	rep << "# SKILLS           \t= " << nbSkills_ << std::endl;
	for(int i=0; i<nbSkills_; i++){
		rep << "#                  \t= " << i << ":" << intToSkill_[i] << std::endl;
	}
	rep << "# " << std::endl;
	rep << "# SHIFTS           \t= " << nbShifts_ << std::endl;
	for(int i=0; i<nbShifts_; i++){
		rep << "#                  \t= ";
		rep << i << ":" << intToShift_[i] << " \t(" << minConsShifts_[i] << "," << maxConsShifts_[i] << ")" << std::endl;
	}
	rep << "# " << std::endl;
	rep << "# FORBIDDEN        " << std::endl;
	for(int i=0; i<nbShifts_; i++){
		rep << "#\t\t\t" << intToShift_[i] << "\t-> ";
		for(int j=0; j<nbForbiddenSuccessors_[i]; j++){
			rep << intToShift_[forbiddenSuccessors_[i][j]] << " ";
		}
		rep << std::endl;
	}
	rep << "# CONTRACTS        " << std::endl;
	for(map<string,Contract*>::const_iterator itC = contracts_.begin(); itC != contracts_.end(); ++itC){
		rep << "#\t\t\t" << *(itC->second) << std::endl;
	}
	rep << "# " << std::endl;
	rep << "# NURSES           \t= " << nbNurses_ << std::endl;
	for(int i=0; i<nbNurses_; i++){
		rep << "#\t\t\t" << theNurses_[i].toString() << std::endl;
	}
	rep << "# " << std::endl;
	rep << "# POSITIONS        \t= " << nbPositions_ << std::endl;
	for (int i=0; i<nbPositions_; i++) {
		rep << "#\t\t\t" << pPositions_[i]->toString() << std::endl;
	}
	if (weekName_!=""){
		// write the demand using the member method toString
		// do not write the preprocessed information at this stage
		//
		rep << pWeekDemand_->toString(false) << std::endl;

		// write the preferences
		//
		rep << "# " << std::endl;
		rep << "# WISHED SHIFTS OFF" << std::endl;
		for(int n=0; n<nbNurses_; n++){
			// Display only if the nurse has preferences
			map<int,set<int> > prefNurse = weekPreferences_.wishesOff_[n];
			if(!prefNurse.empty()){
				rep << "#\t\t\t" << n << "\t" << theNurses_[n].name_ << "\t";
				for(map<int,set<int> >::iterator itWishlist = prefNurse.begin(); itWishlist != prefNurse.end(); ++itWishlist){
					rep << Tools::intToDay(itWishlist->first) << ": ";
					set<int> dayList = itWishlist->second;
					bool first = true;
					for(set<int>::iterator itShift = dayList.begin(); itShift != dayList.end(); ++itShift){
						if(first) first = false; else rep << ",";
						rep << intToShift_[*itShift];
					}
					rep << "    ";
				}
				rep << std::endl;
			}
		}
	}
	if(thisWeek_ > -1){
		rep << "# " << std::endl;
		rep << "# INITIAL STATE    \t= WEEK Nb " << thisWeek_ << std::endl;
		for(int n=0; n<nbNurses_; n++){
			rep << "#\t\t\t" << theNurses_[n].name_ << " ";
			State s = initialState_[n];
			rep << s.totalDaysWorked_ << " " << s.totalWeekendsWorked_ << " " << intToShift_[s.shift_] << " ";
			if(s.shift_) rep << s.consShifts_ << " " << s.consDaysWorked_; else	rep << "0 0";
			if(s.shift_) rep << " 0"; else rep << " " << s.consShifts_;
			rep << std::endl;
		}
	}
	rep << "############################################################################" << std::endl;
	return rep.str();
}


//------------------------------------------------
// Preprocess functions
//------------------------------------------------

// preprocess the nurses to get the types
//
void Scenario::preprocessTheNurses() {

	if (nbPositions_) {
		Tools::throwError("The nurse preprocessing is run for the second time!");
	}

	// Go through the nurses, and created their positions when it has not already
	// been done
	//
	vector<Nurse>::const_iterator itNurse = theNurses_.begin();
	for (itNurse = theNurses_.begin(); itNurse != theNurses_.end(); itNurse++)	{
		bool positionExists = nbPositions_? true:false;
		int nbSkills = (*itNurse).nbSkills_;
		vector<int> skills = (*itNurse).skills_;
		vector<Position*>::iterator itPos = pPositions_.begin();

		// go through every existing position to see if the position of this nurse
		// has already been created
		for (itPos = pPositions_.begin(); itPos != pPositions_.end(); itPos++)	{
			positionExists = true;
			if ((*itPos)->nbSkills_ == nbSkills) {
				for (int i = 0; i < nbSkills; i++) {
					if (skills[i] != (*itPos)->skills_[i])	{
						positionExists = false;
						break;
					}
				}
			}
			else positionExists = false;
			if (positionExists) break;
		}

		// create the position if if doesn't exist
		if (!positionExists) {
			pPositions_.push_back(new Position(nbPositions_, nbSkills, skills));
			nbPositions_++;
		}
	}

	// build the list of position dominance
	for (int i=0; i<nbPositions_; i++) {
		for (int j=i+1; j<nbPositions_; j++) {
			if(pPositions_[i]->compare(*pPositions_[j]) == 1) {
				pPositions_[i]->addBelow(pPositions_[j]);
				pPositions_[j]->addAbove(pPositions_[i]);
			}
			if(pPositions_[i]->compare(*pPositions_[j]) == -1) {
				pPositions_[i]->addAbove(pPositions_[j]);
				pPositions_[j]->addBelow(pPositions_[i]);
			}
		}
	}

	// compute the rank of each position
	vector<bool> isRanked;
	bool isAllRanked = false;
	for (int i =0; i < nbPositions_; i++)	{
		isRanked.push_back(false);
	}
	while (!isAllRanked) {
		for (int i=0; i<nbPositions_; i++) {
			if (!isRanked[i]) {
				int rankIni = pPositions_[i]->rank();

				// go through the positions above the considered position to increment
				// its rank
				if (pPositions_[i]->nbAbove()) {
					for (int j = 0; j < pPositions_[i]->nbAbove(); j++) {
						int currentRank = pPositions_[i]->rank();
						int newRank = pPositions_[j]->rank()+1;
						pPositions_[i]->rank(std::max(currentRank, newRank));
					}
				}

				// the position is treated when the rank is not modified by the loop above
				if (pPositions_[i]->rank() == rankIni) isRanked[i] = true;
			}
		}

		// check if all the positions are ranked
		isAllRanked = true;
		for (int i=0; i<nbPositions_; i++) {
			if (!isRanked[i]) {
				isAllRanked = false;
				break;
			}
		}
	}


}
