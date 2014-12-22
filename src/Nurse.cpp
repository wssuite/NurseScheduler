#include "MyTools.h"
#include "Nurse.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <streambuf>
#include <time.h>


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------

// Print method
//
string Contract::toString(){
	std::stringstream rep;
	rep << name_ << "  -  ";
	rep << "Tot:" << minTotalShifts_ << "<" << maxTotalShifts_ << "  |  ";
	rep << "Work:" << minConsDaysWork_ << "<" << maxConsDaysWork_ << "  |  ";
	rep << "Rest:" << minConsDaysOff_ << "<" << maxConsDaysOff_ << "  |  ";
	rep << "WE:" << maxTotalWeekends_ << "  |  ";
	if(!isCompleteWeekends_) rep << "NOT";
	rep << "complete";
	return rep.str();
};

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   S t a t e
//
//-----------------------------------------------------------------------------

// Destructor
State::~State(){}

// Updates the state if a new day is worked on shift newShift
void State::addNewDay(int newShift){

	// Total shifts worked if it is a worked day
	totalDaysWorked_ += (newShift ? 1 : 0);

	// Total weekends worked :
	// +1 IF : new day is worked AND (new day is saturday OR (new day is sunday AND previous day was not worked) )
	if( newShift and
			( (dayId_/7==4) or ((dayId_/7==5) and !shift_)))
		totalWeekendsWorked_ ++;

	// Consecutives : +1 iff it is the same as the previous one
	consShifts_ = (shift_==newShift) ? (consShifts_ + 1) : 1;

	// Consecutive Days Worked : +1 if the new one is worked (!=0), 0 if it is a rest (==0)
	consDaysWorked_ = shift_ ? (consDaysWorked_ + 1) : 0;

	// Current shift worked : updated with the new one
	shift_ = newShift;

	// Finally, the day index
	dayId_++;
}

// Display method: toString
//
string State::toString(){
	std::stringstream rep;
	rep << totalDaysWorked_ << " " << totalWeekendsWorked_ << " " << shift_ << " ";
	if(shift_) rep << consShifts_ << " " << consDaysWorked_; else rep << "0 0";
	if(shift_) rep << " 0"; else rep << " " << consShifts_;
	rep << std::endl;
	return rep.str();
}



//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   P r e f e r e n c e
//
//-----------------------------------------------------------------------------


// Default constructor and destructor
Preferences::Preferences(){}
Preferences::~Preferences(){}

// Initialization with an vector or size nNurses with no wished Shift-Off.
Preferences::Preferences(int nbNurses, int nbDays, int nbShifts) :
	nbNurses_(nbNurses), nbDays_(nbDays), nbShifts_(nbShifts){
	// Wish lists are initialized to empty
	vector<map<int,set<int> > > wishesOff;
	for(int i=0; i<nbNurses_; i++){
		map<int,set<int> > m;
		wishesOff.push_back(m);
	}
	wishesOff_ = wishesOff;
}

// Add a wished day-shift off for a nurse
void Preferences::addShiftOff(int nurse, int day, int shift){
	// If the nurse does not already have a wish for that day -> insert a new emptyset on that day
	if(wishesOff_[nurse].find(day) == wishesOff_[nurse].end()){
		set<int> emptyset;
		wishesOff_[nurse].insert(pair<int,set<int> >(day,emptyset));
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
	wishesOff_[nurse].insert(pair<int,set<int> >(day,wishList));
}

// Returns true if the nurses wants that shift off
bool Preferences::wantsTheShiftOff(int nurse, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,set<int> >::iterator itM = wishesOff_[nurse].find(day);
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
	map<int,set<int> >::iterator itM = wishesOff_[nurse].find(day);
	if(itM == wishesOff_[nurse].end())
		return false;
	else if(itM->second.size() < nbShifts_-1)		// Set does not repeat its elements. Wants the day off if and only if she wants all shifts off (-1 because REST does not appear)
		return false;
	else
		return true;
}

// Display method: toString()
//
string Preferences::toString(){
	std::stringstream rep;
	rep << "# Preference display not implemented yet..." << std::endl;
	return rep.str();
}


//-----------------------------------------------------------------------------
//
//	C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

// constructor and destructor
//
Demand::Demand(int nbDays, int nbShifts, int nbSkills,
	vector3D minDemand, vector3D optDemand):
		nbDays_(nbDays), nbShifts_(nbShifts), nbSkills_(nbSkills),
		minDemand_(minDemand), optDemand_(optDemand),
		minTotal_(0), optTotal_(0)
		{
			// initialize the preprocessed vectors
			Tools::initVector(&minPerDay_, nbDays_);
			Tools::initVector(&optPerDay_, nbDays_);
			Tools::initVector(&minPerShift_, nbShifts_);
			Tools::initVector(&optPerShift_, nbShifts_);
			Tools::initVector(&minPerSkill_, nbSkills_);
			Tools::initVector(&optPerSkill_, nbSkills_);
			Tools::initVector(&minHighestPerSkill_, nbSkills_);
			Tools::initVector(&optHighestPerSkill_, nbSkills_);

			// run the preprocessing
			this->preprocessDemand();
		}

Demand::~Demand()
{}

// compute all the potentially helpful attributes of a demand
// this includes the total demand per skill, per shift,
//
void Demand::preprocessDemand() {
	for (int day = 0; day < nbDays_; day++)	{
		for (int shift = 0; shift < nbShifts_; shift++) {
			for (int skill = 0; skill < nbSkills_; skill++)	{
				// update the total demand
				minTotal_ += minDemand_[day][shift][skill];
				optTotal_ += optDemand_[day][shift][skill];

				// update the demand per day
				minPerDay_[day] += minDemand_[day][shift][skill];
				optPerDay_[day] += optDemand_[day][shift][skill];

				// update the demand per shift
				minPerShift_[shift] += minDemand_[day][shift][skill];
				optPerShift_[shift] += optDemand_[day][shift][skill];

				// update the demand per skill
				minPerSkill_[skill] += minDemand_[day][shift][skill];
				optPerSkill_[skill] += optDemand_[day][shift][skill];

				// update the demand per day
				minHighestPerSkill_[skill] +=
					std::max(minDemand_[day][shift][skill],minHighestPerSkill_[skill]);
				optHighestPerSkill_[skill] +=
					std::max(optDemand_[day][shift][skill],optHighestPerSkill_[skill]);
			}
		}
	}
}



//-----------------------------------------------------------------------------
//
//  C l a s s   N u r s e
//
//  Class that contains all the attributes describing the characteristics and
//  the planning of each nurse
//
//-----------------------------------------------------------------------------

// Destructor
//
Nurse::~Nurse(){
	// WARNING: Do NOT delete Contract* contract (eventhough it is a pointer.
	//          Contracts are common to all nurses and don't "belong" to them -> should not be deleted.
}

// Print method
//
string Nurse::toString(){
	std::stringstream rep;
	rep << id_ << ": ";
	if(id_<10) rep << " ";
	if(id_<100) rep << " ";
	rep << name_ << "\t" << nbSkills_ << " [ ";
	for(int i=0; i<nbSkills_; i++) rep << skills_[i] << " ";
	rep << "]\t";
	if(nbSkills_==1) rep << "\t";
	rep << pContract_->name_;
	return rep.str();
}

// Assignment operator -> vector request
//
Nurse& Nurse::operator=(const Nurse& n){
	int id = n.id_;
	string name = n.name_;
	int nbSkills = n.nbSkills_;
	vector<int> skills = n.skills_;
	const Contract* contract;
	contract = n.pContract_;
	Nurse * n2 = new Nurse(id, name, nbSkills, skills, contract);
	return *n2;
}
