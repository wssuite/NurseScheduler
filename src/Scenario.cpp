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

Scenario::~Scenario(){

}

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
		rep << "#\t\t\t" << theNurses_[i] << std::endl;
	}
	if (weekName_!=""){
		rep << "# " << std::endl;
		rep << "# DEMAND" << std::endl;

		rep << "#\t\t\t";
		for(int dayId=0; dayId<7; dayId++){
			rep << " " << Tools::intToDay(dayId) << "\t";
		}
		rep << "# " << std::endl;
		for(int sh = 0; sh < nbShifts_; sh ++){
			for (int sk = 0; sk < nbSkills_; sk++){
				string s = "#   " + intToShift_[sh] + " " + intToSkill_[sk] + " ";
				rep << s;
				if(s.length() < 16) rep << "\t";
				for(int day = 0; day < 7; day ++){
					rep << "\t(" << minWeekDemand_[day][sh][sk] << "," << optWeekDemand_[day][sh][sk] << ")";
				}
				rep << std::endl;
			}
			rep << "# " << std::endl;
		}
		rep << "# " << std::endl;
		rep << "# WISHED SHIFTS OFF" << std::endl;
		for(int n=0; n<nbNurses_; n++){
			// cout only if the nurse has preferences
			map<int,set<int> > prefNurse = weekPreferences_.wishesOff_[n];
			if(!prefNurse.empty()){
				rep << "#\t\t\t" << n << "\t" << theNurses_[n].name_ << "\t";
				for(map<int,set<int> >::iterator itWishlist = prefNurse.begin(); itWishlist != prefNurse.end(); ++itWishlist){
					rep << Tools::intToDay(itWishlist->first) << ": ";
					set<int> dayList = itWishlist->second;
					for(set<int>::iterator itShift = dayList.begin(); itShift != dayList.end(); ++itShift){
						rep << *itShift << " ";
						//rep << intToShift_[*itShift] << ",";
					}
					rep << "    ";
				}
				rep << std::endl;
			}
		}
	}
	rep << "############################################################################" << std::endl;
	return rep.str();
}

