//
//  Scenario.cpp
//  Project:RosterDesNurses
//
//  Created by Jérémy Omer on 18/12/2013.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//


#include "Scenario.h"
#include "MyTools.h"
#include "Nurse.h"

#include <iostream>
#include <streambuf>
#include <fstream>
#include <math.h>
#include <time.h>


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//-----------------------------------------------------------------------------

// Destructor
State::~State(){}

// Updates the state if a new day is worked on shift newShift
void State::updateWithNewDay(int newShift){

	// Consecutives : +1 iff it is the same as the previous one
	consShifts_ = (shift_==newShift) ? (consShifts_ + 1) : 1;

	// Consecutive Days Worked : +1 if the new one is worked (!=0), 0 if it is a rest (==0)
	consDaysWorked_ = shift_ ? (consDaysWorked_ + 1) : 0;

	// Current shift worked : updated with the new one
	shift_ = newShift;
}



//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   P r e f e r e n c e
//
//-----------------------------------------------------------------------------

// Destructor
Preferences::~Preferences(){}

// Initialization with an vector or size nNurses with no wished Shift-Off.
Preferences::Preferences(int nbNurses, int nbShifts) :
	nbNurses_(nbNurses), nbShifts_(nbShifts){
	// Wish lists are initialized to empty
	vector<map<int,set<int>>> wishesOff;
	for(int i=0; i<nbNurses_; i++){
		map<int,set<int>> m;
		wishesOff.push_back(m);
	}
	wishesOff_ = wishesOff;
}

// Add a wished day-shift off for a nurse
void Preferences::addShiftOff(int nurse, int day, int shift){
	// If the nurse does not already have a wish for that day -> insert a new emptyset on that day
	if(wishesOff_[nurse].find(day) == wishesOff_[nurse].end()){
		set<int> emptyset;
		wishesOff_[nurse].insert(pair<int,set<int>>(day,emptyset));
	}
	// Insert the wished shift in the set
	wishesOff_[nurse][day].insert(shift);
}

// Adds the whole day to the wish-list
void Preferences::addDayOff(int nurse, int day){
	set<int> wishList;
	for(int s=1; s<nbShifts_; s++){				// Starts from 1 because no nurse wished "not to rest"
		wishList.insert(s);
	}
	wishesOff_[nurse].insert(pair<int,set<int>>(day,wishList));
}

// Returns true if the nurses wants that shift off
bool Preferences::wantsTheShiftOff(int nurse, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,set<int>>::iterator itM = wishesOff_[nurse].find(day);
	if(itM == wishesOff_[nurse].end())
		return false;
	// If the shift is not in the wish-list for that day, return false
	else if(itM->second.find(shift) == itM->second.end())
		return false;
	else
		return true;
}

// True if the nurses wants the whole day off
bool Preferences::wantsTheDayOff(int nurse, int day){
	// If the day is not in the wish-list, return false
	map<int,set<int>>::iterator itM = wishesOff_[nurse].find(day);
	if(itM == wishesOff_[nurse].end())
		return false;
	else if(itM->second.size() < nbShifts_-1)		// Set does not repeat its elements. Wants the day off if and only if she wants all shifts off (-1 because REST does not appear)
		return false;
	else
		return true;
}




