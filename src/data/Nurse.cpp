#include "data/Nurse.h"

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
	if(!needCompleteWeekends_) rep << "NOT";
	rep << "complete";
	return rep.str();
}


//-----------------------------------------------------------------------------
//
//  C l a s s   P o s i t i o n
//
//  A position (job) is a set of skills that a nurse may possess
//
//-----------------------------------------------------------------------------

// Constructor and Destructor
//
Position::Position(int index, int nbSkills, vector<int> skills):
			id_(index), nbSkills_(nbSkills), skills_(skills), nbBelow_(0), nbAbove_(0), rank_(0)  {
	// Verify that the vecor of skills is sorted
	//
	for (vector<int>::const_iterator it=skills_.begin(); it!=skills_.end()-1; ++it) {
		if (*it >= *(it+1))  {
			Tools::throwError("The skills in a position are not sorted or some skill is repeated!");
		}
	}
	for (int sk = 0; sk < nbSkills_; sk++) {
	  skillRarity_.push_back(1.0);
	}
}

// Set positions above and below
//
void Position::addBelow(Position* pPosition) {
	positionsBelow_.push_back(pPosition);
	nbBelow_++;
}
void Position::addAbove(Position* pPosition) {
	positionsAbove_.push_back(pPosition);
	nbAbove_++;
}

// Print method
//
string Position::toString() const{
	std::stringstream rep;
	rep << id_ << ": ";
	if(id_<10) rep << " ";
	if(id_<100) rep << " ";

	// print the rank and the list of skills
	rep << ": rank = " << rank_;
	rep << " ; skills = [ ";
	for(int i=0; i<nbSkills_; i++) rep << skills_[i] << " ";
	rep << "]\t";

	// print the rank
	// print the list of positions above and below this one
	if (!positionsBelow_.empty()) {
		rep << std::endl;
		rep << "#\t\t\t\t\tdominates positions:      ";
		for (vector<Position*>::const_iterator it=positionsBelow_.begin(); it!=positionsBelow_.end();it++) {
			rep << "\t" << (*it)->id_;
		}
	}
	if (!positionsAbove_.empty()) {
		rep << std::endl;
		rep << "#\t\t\t\t\tis dominated by positions:";
		for (vector<Position*>::const_iterator it=positionsAbove_.begin(); it!=positionsAbove_.end();it++) {
			rep << "\t" << (*it)->id_;
		}
	}
	return rep.str();
}


 // update the rarity of the skills
 // the input is the vector of the rarity of all the skills
 // the vector is sorted without record of the corresponding skill because it
 // is used only to compare two positions with the same rank
 //
 void Position::updateRarities(vector<double> allRarities) {
	for (int sk = 0; sk < this->nbSkills_; sk++) {
	  skillRarity_[sk] = allRarities[skills_[sk]];
	}
	std::sort(skillRarity_.begin(),skillRarity_.end(), std::greater<double>());
 }

// Compare this position with the input position
// The dominance criterion is that a position p1 with skills sk1 dominates p2
// with skills sk2 if and only if (sk1 contains sk2) and sk1 has more skills
// than sk2
// The function returns 1 if this position dominates, -1 if it is dominated
// and 0 if there is no dominance
//
int Position::compare(const Position &p)  {
	vector<int>::const_iterator it1;
	vector<int>::const_iterator it2;

	// no possible dominance if both positions have as many skills
	if (p.nbSkills_ == this->nbSkills_) {
		return 0;
	}
	// only p can dominate if it has more skills
	// the comparison of the two skill lists is based on the fact that they are
	// both sorted
	else if(p.nbSkills_ > this->nbSkills_) {
		it1 = p.skills_.begin();
		for (it2 = this->skills_.begin(); it2 != this->skills_.end(); it2++) {
			while (*it1 < *it2) {
				if (it1 == p.skills_.end()-1) return 0;
				else it1++;
			}
			if (*it1 != *it2) return 0;
		}
		return -1;
	}
	// only this position can dominate if it has more skills
	else if(this->nbSkills_ > p.nbSkills_) {
		it1 = this->skills_.begin();
		for (it2 = p.skills_.begin(); it2 != p.skills_.end(); it2++)   {
			while (*it1 < *it2) {
				if (it1 == this->skills_.end()-1) return 0;
				else it1++;
			}
			if (*it1 != *it2) return 0;
		}
		return 1;
	}

	return 0;
}

// returns true if the position shares at least one skill with the input position
//
bool Position::shareSkill(const Position&p) {
	vector<int>::const_iterator it1;
	vector<int>::const_iterator it2;

	for (it2 = this->skills_.begin(); it2 != this->skills_.end(); it2++) {
		for (it1 = p.skills_.begin(); it1 != p.skills_.end(); it1++)   {
			if (*it1== *it2) return true;
		}
	}

	return false;
}

// reset the list of positions below and above
//
void Position::resetAbove() {
	positionsAbove_.clear();
	nbAbove_ = 0;
}
void Position::resetBelow() {
	positionsBelow_.clear();
	nbBelow_ = 0;
}


//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   P r e f e r e n c e
//
//-----------------------------------------------------------------------------


// Default constructor and destructor
Preferences::Preferences(){}
Preferences::~Preferences(){}

// Initialization with a map of size nNurses with no wished Shift-Off.
Preferences::Preferences(int nbNurses, int nbDays, int nbShifts) :
				nbNurses_(nbNurses), nbDays_(nbDays), nbShifts_(nbShifts){
	// Wish lists are initialized to empty
	map<int, map<int,vector<Wish> > > wishesOff;
	for(int i=0; i<nbNurses_; i++){
		map<int,vector<Wish> > m;
		wishesOff[i] = m;
	}
	wishesOff_ = wishesOff;
	
	map<int, map<int,vector<Wish> > > wishesOn;
	for(int i=0; i<nbNurses_; i++){
		map<int,vector<Wish> > m;
		wishesOn[i] = m;
	}
	wishesOn_ = wishesOn;
}

// Initialization with a map corresponding to the input nurses and no wished Shift-Off.
Preferences::Preferences(vector<Nurse>& nurses, int nbDays, int nbShifts) :
				nbNurses_(nurses.size()), nbDays_(nbDays), nbShifts_(nbShifts){
	// Wish lists are initialized to empty
	map<int, map<int,vector<Wish> > > wishesOff;
	for(Nurse nurse: nurses){
		map<int,vector<Wish> > m;
		wishesOff[nurse.id_] = m;
	}
	wishesOff_ = wishesOff;
}

// Add a wished day-shift off for a nurse
void Preferences::addShiftOff(int nurseId, int day, int shift, int level){
	// If the nurse does not already have a wish for that day -> insert a new emptyset on that day
	if(wishesOff_[nurseId].find(day) == wishesOff_[nurseId].end()){
		vector<Wish> emptyset;
		wishesOff_[nurseId].insert(pair<int,vector<Wish> >(day,emptyset));
	}

	Wish  wish{shift, level};
	// Insert the wished shift in the set
	wishesOff_[nurseId][day].push_back(wish);
	
	// If the nurse does not already have a wish for that day -> insert a new emptyset on that day
	if(wishesOn_[nurseId].find(day) == wishesOn_[nurseId].end()){
		vector<Wish> emptyset;
		wishesOn_[nurseId].insert(pair<int,vector<Wish> >(day,emptyset));
	}

	Wish  wish2{shift, level};
	// Insert the wished shift in the set
	wishesOn_[nurseId][day].push_back(wish2);
}

// Adds the whole day to the wish-list
void Preferences::addDayOff(int nurseId, int day, int level){
	vector<Wish> wishList;
	for(int s=1; s<nbShifts_; s++){           // Starts from 1 because no nurse wished "not to rest"
	  Wish  wish{s, level};
	  wishList.push_back(wish);
	}
	wishesOff_[nurseId].insert(pair<int,vector<Wish> >(day,wishList));
}

// Returns true if the nurse wants that shift off
bool Preferences::wantsTheShiftOff(int nurseId, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOff_[nurseId].find(day);
	if(itM == wishesOff_[nurseId].end())
		return false;
	// If the shift is not in the wish-list for that day, return false
	else {
	  vector<Wish>::iterator itShift;
	  for (itShift = itM->second.begin(); itShift != itM->second.end(); itShift++) {
	    if(itShift->shift == shift)
	      return true;
	  }
	}

	return false;
}

// Returns level if the nurse wants that shift off : -1 otherwise
int Preferences::wantsTheShiftOffLevel(int nurseId, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOff_[nurseId].find(day);
	if(itM == wishesOff_[nurseId].end())
	  return -1;
	// If the shift is not in the wish-list for that day, return false
	else {
	  vector<Wish>::iterator itShift;
	  for (itShift = itM->second.begin(); itShift != itM->second.end(); itShift++) {
	    if(itShift->shift == shift)
	      return itShift->level;
	  }
	}

	return -1;
}

// True if the nurse wants the whole day off
bool Preferences::wantsTheDayOff(int nurseId, int day){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOff_[nurseId].find(day);
	if(itM == wishesOff_[nurseId].end())
		return false;
	else if((int) itM->second.size() < nbShifts_-1)    // Set does not repeat its elements. Wants the day off if and only if she wants all shifts off (-1 because REST does not appear)
		return false;
	else
		return true;
}

// Total number of shifts off that the nurse wants
int Preferences::howManyShiftsOff(int nurseId) {
	int nbShiftsOff = 0;

	// look at every wishes of the nurse
	for( pair<int,vector<Wish> > itDayOff: wishesOff_[nurseId]){
		nbShiftsOff += itDayOff.second.size();
	}
	return nbShiftsOff;
}

// Number of whole days off that the nurse wants
int Preferences::howManyDaysOff(int nurseId, int dayMin, int dayMax){
	int nbDayOff = 0;
	// look at every wishes of the nurse
	for( pair<int,vector<Wish> > itDayOff: wishesOff_[nurseId]){
		nbDayOff += ( (itDayOff.first >= dayMin) && (itDayOff.first <= dayMax) );
	}
	return nbDayOff;
}

// add another week preferences at the end of the current one
//
void Preferences::push_backOff(Preferences* pPreferences){
	// check if same scenario
	if( (nbShifts_ != pPreferences->nbShifts_) || (nbNurses_ != pPreferences->nbNurses_) ){
		string error = "Preferences are not compatible";
		Tools::throwError(error.c_str());
	}

	//update the preferences
	map<int, map<int,std::vector<Wish> > >::iterator it = wishesOff_.begin();
	for(it = wishesOff_.begin(); it!=wishesOff_.end(); it++){
		int nurseId = (*it).first;
		for(pair<int,std::vector<Wish> > pair1: pPreferences->wishesOff_[nurseId]){
			pair<int,std::vector<Wish> > pair2(pair1.first+nbDays_, pair1.second);
			(*it).second.insert(pair2);
		}
	}

	// update the number of days
	nbDays_  += pPreferences->nbDays_;
}

// K the preferences relative to the nbDays first days
Preferences* Preferences::keepOff(int begin, int end) {

   Preferences* pPref = new Preferences();

   for (int i=0; i < nbNurses_; i++) {
      for(pair<int,std::vector<Wish> > pair1: wishesOff_[i]){
         if (pair1.first >= begin && pair1.first < end) {
            pair<int,std::vector<Wish> > pair2(pair1.first-begin, pair1.second);
            pPref->wishesOff_[i].insert(pair2);
         }
      }
   }

   return pPref;
}

// Remove the preferences relative to the nbDays first days
Preferences* Preferences::removeNFirstDayOff(int nbDays) {

	Preferences* pPref = new Preferences();

	for (int i=0; i < nbNurses_; i++) {
		for(pair<int,std::vector<Wish> > pair1: wishesOff_[i]){
			if (pair1.first >= nbDays) {
				pair<int,std::vector<Wish> > pair2(pair1.first-nbDays_, pair1.second);
				pPref->wishesOff_[i].insert(pair2);
			}
		}
	}

	return pPref;
}

//////////////////////////////////////////////////////////////////

// Add a wished day-shift off for a nurse
void Preferences::addShiftOn(int nurseId, int day, int shift, int level){
	// If the nurse does not already have a wish for that day -> insert a new emptyset on that day
	if(wishesOn_[nurseId].find(day) == wishesOn_[nurseId].end()){
		vector<Wish> emptyset;
		wishesOn_[nurseId].insert(pair<int,vector<Wish> >(day,emptyset));
	}

	Wish  wish{shift, level};
	// Insert the wished shift in the set
	wishesOn_[nurseId][day].push_back(wish);
}

// Adds the whole day to the wish-list
void Preferences::addDayOn(int nurseId, int day, int level){
	vector<Wish> wishList;
	for(int s=1; s<nbShifts_; s++){           // Starts from 1 because no nurse wished "not to rest"
	  Wish  wish{s, level};
	  wishList.push_back(wish);
	}
	wishesOn_[nurseId].insert(pair<int,vector<Wish> >(day,wishList));
}

// Returns true if the nurse wants that shift off
bool Preferences::wantsTheShiftOn(int nurseId, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOn_[nurseId].find(day);
	if(itM == wishesOn_[nurseId].end())
		return false;
	// If the shift is not in the wish-list for that day, return false
	else {
	  vector<Wish>::iterator itShift;
	  for (itShift = itM->second.begin(); itShift != itM->second.end(); itShift++) {
	    if(itShift->shift == shift)
	      return true;
	  }
	}

	return false;
}

// Returns level if the nurse wants that shift off : -1 otherwise
int Preferences::wantsTheShiftOnLevel(int nurseId, int day, int shift){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOn_[nurseId].find(day);
	if(itM == wishesOn_[nurseId].end())
	  return -1;
	// If the shift is not in the wish-list for that day, return false
	else {
	  vector<Wish>::iterator itShift;
	  for (itShift = itM->second.begin(); itShift != itM->second.end(); itShift++) {
	    if(itShift->shift == shift)
	      return itShift->level;
	  }
	}

	return -1;
}

// True if the nurse wants the whole day off
bool Preferences::wantsTheDayOn(int nurseId, int day){
	// If the day is not in the wish-list, return false
	map<int,vector<Wish> >::iterator itM = wishesOn_[nurseId].find(day);
	if(itM == wishesOn_[nurseId].end())
		return false;
	else if((int) itM->second.size() < nbShifts_-1)    // Set does not repeat its elements. Wants the day off if and only if she wants all shifts off (-1 because REST does not appear)
		return false;
	else
		return true;
}

// Total number of shifts off that the nurse wants
int Preferences::howManyShiftsOn(int nurseId) {
	int nbShiftsOn = 0;

	// look at every wishes of the nurse
	for( pair<int,vector<Wish> > itDayOn: wishesOn_[nurseId]){
		nbShiftsOn += itDayOn.second.size();
	}
	return nbShiftsOn;
}

// Number of whole days off that the nurse wants
int Preferences::howManyDaysOn(int nurseId, int dayMin, int dayMax){
	int nbDayOn = 0;
	// look at every wishes of the nurse
	for( pair<int,vector<Wish> > itDayOn: wishesOn_[nurseId]){
		nbDayOn += ( (itDayOn.first >= dayMin) && (itDayOn.first <= dayMax) );
	}
	return nbDayOn;
}

// add another week preferences at the end of the current one
//
void Preferences::push_backOn(Preferences* pPreferences){
	// check if same scenario
	if( (nbShifts_ != pPreferences->nbShifts_) || (nbNurses_ != pPreferences->nbNurses_) ){
		string error = "Preferences are not compatible";
		Tools::throwError(error.c_str());
	}

	//update the preferences
	map<int, map<int,std::vector<Wish> > >::iterator it = wishesOn_.begin();
	for(it = wishesOn_.begin(); it!=wishesOn_.end(); it++){
		int nurseId = (*it).first;
		for(pair<int,std::vector<Wish> > pair1: pPreferences->wishesOn_[nurseId]){
			pair<int,std::vector<Wish> > pair2(pair1.first+nbDays_, pair1.second);
			(*it).second.insert(pair2);
		}
	}

	// update the number of days
	nbDays_  += pPreferences->nbDays_;
}

// K the preferences relative to the nbDays first days
Preferences* Preferences::keepOn(int begin, int end) {

   Preferences* pPref = new Preferences();

   for (int i=0; i < nbNurses_; i++) {
      for(pair<int,std::vector<Wish> > pair1: wishesOn_[i]){
         if (pair1.first >= begin && pair1.first < end) {
            pair<int,std::vector<Wish> > pair2(pair1.first-begin, pair1.second);
            pPref->wishesOn_[i].insert(pair2);
         }
      }
   }

   return pPref;
}

// Remove the preferences relative to the nbDays first days
Preferences* Preferences::removeNFirstDayOn(int nbDays) {

	Preferences* pPref = new Preferences();

	for (int i=0; i < nbNurses_; i++) {
		for(pair<int,std::vector<Wish> > pair1: wishesOn_[i]){
			if (pair1.first >= nbDays) {
				pair<int,std::vector<Wish> > pair2(pair1.first-nbDays_, pair1.second);
				pPref->wishesOn_[i].insert(pair2);
			}
		}
	}

	return pPref;
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
//  C l a s s   N u r s e
//
//  Class that contains all the attributes describing the characteristics and
//  the planning of each nurse
//
//-----------------------------------------------------------------------------

// Constructor and destructor
// Note : need both with const Contract and (non-const) Contract because non-const is used in our code,
//        and const is needed so that we can override the operator= and have vector<Nurse>. We need to
//        override it because vector members should have some properties (assignable a.o., which implies
//        non-const)
//
Nurse::Nurse(int id, string name, int nbSkills, vector<int> skills, Contract* contract) :
			id_(id), name_(name), nbSkills_(nbSkills), skills_(skills), pContract_(contract){
	// Verify that the vecor of skills is sorted
	//
	for (vector<int>::const_iterator it=skills_.begin(); it!=skills_.end()-1; ++it) {
		if (*it >= *(it+1))  {
			Tools::throwError("The skills in a nurse are not sorted or some skill is repeated!");
		}
	}
}
Nurse::Nurse(int id, string name, int nbSkills, vector<int> skills, const Contract* contract) :
			id_(id), name_(name), nbSkills_(nbSkills), skills_(skills), pContract_(contract){
	// Verify that the vecor of skills is sorted
	//
	for (vector<int>::const_iterator it=skills_.begin(); it!=skills_.end()-1; ++it) {
		if (*it >= *(it+1))  {
			Tools::throwError("The skills in a nurse are not sorted or some skill is repeated!");
		}
	}
}

Nurse::~Nurse(){
	// WARNING: Do NOT delete Contract* contract (eventhough it is a pointer.
	//          Contracts are common to all nurses and don't "belong" to them -> should not be deleted.
}

// Check that the nurse has a given skill
//
bool Nurse::hasSkill(int skill) const {
	for (int i = 0; i < nbSkills_; i++) {
		if (skills_[i] == skill)   return true;
	}
	return false;
}
// Print method
//
string Nurse::toString() const{
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
