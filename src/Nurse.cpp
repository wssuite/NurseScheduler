#include "MyTools.h"
#include "Nurse.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//-----------------------------------------------------------------------------

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

// Initialization with an vector or size nNurses with no wished Shift-Off.
Preferences::Preferences(int nNurses){
	nNurses_ = nNurses;
	vector<map<int,set<int>>> wishesOff
	for(int i=0; i<nNurses; i++){
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

