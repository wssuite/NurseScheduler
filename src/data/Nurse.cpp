
#include "data/Nurse.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <streambuf>
#include <time.h>


using std::vector;
using std::map;
using std::string;
using std::pair;

//-----------------------------------------------------------------------------
//
//  S t r u c t u r e   C o n t r a c t
//
//  A contract as defined in the subject
//
//-----------------------------------------------------------------------------

// Cost function for consecutive identical shifts
//
double Contract::consDaysCost(int n) const {
  if(minConsDaysWork_ - n > 0) return (WEIGHT_CONS_DAYS_WORK * ( minConsDaysWork_ - n ) );
  if(n - maxConsDaysWork_ > 0) return (WEIGHT_CONS_DAYS_WORK * ( n - maxConsDaysWork_ ) );
  return 0;
}

// Print method
//
std::string Contract::toString(){
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
Position::Position(int index, int nbSkills, std::vector<int> skills):
			id_(index), nbSkills_(nbSkills), skills_(skills), nbBelow_(0), nbAbove_(0), rank_(0)  {
	// Verify that the vecor of skills is sorted
	//
	for (auto it=skills_.begin(); it!=skills_.end()-1; ++it) {
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
void Position::addBelow(PPosition pPosition) {
	positionsBelow_.push_back(pPosition);
	nbBelow_++;
}
void Position::addAbove(PPosition pPosition) {
	positionsAbove_.push_back(pPosition);
	nbAbove_++;
}

// Print method
//
std::string Position::toString() const{
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
		for (auto it=positionsBelow_.begin(); it!=positionsBelow_.end();it++) {
			rep << "\t" << (*it)->id_;
		}
	}
	if (!positionsAbove_.empty()) {
		rep << std::endl;
		rep << "#\t\t\t\t\tis dominated by positions:";
		for (auto it=positionsAbove_.begin(); it!=positionsAbove_.end();it++) {
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
 void Position::updateRarities(std::vector<double> allRarities) {
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
	// no possible dominance if both positions have as many skills
	if (p.nbSkills_ == this->nbSkills_) {
		return 0;
	}
	// only p can dominate if it has more skills
	// the comparison of the two skill lists is based on the fact that they are
	// both sorted
	else if(p.nbSkills_ > this->nbSkills_) {
		auto it1 = p.skills_.begin();
		for (auto it2 = this->skills_.begin(); it2 != this->skills_.end(); it2++) {
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
		auto it1 = this->skills_.begin();
		for (auto it2 = p.skills_.begin(); it2 != p.skills_.end(); it2++)   {
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
	for (auto it2 = this->skills_.begin(); it2 != this->skills_.end(); it2++) {
		for (auto it1 = p.skills_.begin(); it1 != p.skills_.end(); it1++)   {
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
	for(int i=0; i<nbNurses_; i++){
    wishesOff_[i];
    wishesOn_[i];
	}
}

// Initialization with a map corresponding to the input nurses and no wished Shift-Off.
Preferences::Preferences(const vector<PNurse>& nurses, int nbDays, int nbShifts) :
				nbNurses_(nurses.size()), nbDays_(nbDays), nbShifts_(nbShifts){
	// Wish lists are initialized to empty
	for(PNurse pNurse: nurses){
    wishesOff_[pNurse->id_];
    wishesOn_[pNurse->id_];
	}
}


int Preferences::wishLevel(const std::map<int, std::vector<Wish> > &wishes, int day, int shift) {
  auto itM = wishes.find(day);
  // If the day is not in the wish-list, no possible violation
  if (itM == wishes.end())
    return -1;
  // no preference either in the wish-list for that day
  else {
    for (const Wish &w: itM->second)
      if (w.shift == shift)
        return w.level;
  }
  return -1;
}

// Add a wished day-shift off for a nurse
void Preferences::addShiftOff(int nurseId, int day, int shift, PREF_LEVEL level){
	// Insert the wished shift in the set
	wishesOff_[nurseId][day].push_back({shift, level});
}

// Adds the whole day to the wish-list
void Preferences::addDayOff(int nurseId, int day, PREF_LEVEL level){
	vector<Wish>& wishList = wishesOff_[nurseId][day];
	for(int s=1; s<nbShifts_; s++)           // Starts from 1 because it's a rest wish
	  wishList.push_back({s, level});
}

// Returns true if the nurse wants that shift off
bool Preferences::wantsTheShiftOff(int nurseId, int day, int shift) const{
  return wantsTheShiftOffLevel(nurseId, day, shift) != -1;
}

// Returns level if the nurse wants that shift off : -1 otherwise
int Preferences::wantsTheShiftOffLevel(int nurseId, int day, int shift) const {
	// If the day is not in the wish-list, return -1
	auto itM = wishesOff_.at(nurseId).find(day);
	if(itM == wishesOff_.at(nurseId).end())
	  return -1;
	// If the shift is not in the wish-list for that day, return -1
	const std::vector<Wish> &wishes = itM->second;
	for(const Wish &wish: wishes)
	  if(wish.shift == shift) return wish.level;
  return -1;
}

// True if the nurse wants the whole day off
bool Preferences::wantsTheDayOff(int nurseId, int day) const{
	// If the day is not in the wish-list, return false
	auto itM = wishesOff_.at(nurseId).find(day);
	if(itM == wishesOff_.at(nurseId).end())
		return false;
  // Set does not repeat its elements. Wants the day off if and only if all shifts off (-1 because REST does not appear)
	return itM->second.size() == nbShifts_-1;
}

// Total number of shifts off that the nurse wants
int Preferences::howManyShiftsOff(int nurseId) const {
	int nbShiftsOff = 0;

	// look at every wishes of the nurse
	for( const auto & dayOff: wishesOff_.at(nurseId)){
		nbShiftsOff += dayOff.second.size();
	}
	return nbShiftsOff;
}

// Number of whole days off that the nurse wants
int Preferences::howManyDaysOff(int nurseId, int dayMin, int dayMax) const {
	int nbDayOff = 0;
	// look at every wishes of the nurse
	for( const auto & dayOff: wishesOff_.at(nurseId)){
		nbDayOff += ( (dayOff.first >= dayMin) && (dayOff.first <= dayMax) && (dayOff.second.size() == nbShifts_-1));
	}
	return nbDayOff;
}

//////////////////////////////////////////////////////////////////

// Add a wished day-shift on for a nurse
void Preferences::addShiftOn(int nurseId, int day, int shift, PREF_LEVEL level){
  // Insert the wished shift in the set
  wishesOn_[nurseId][day].push_back({shift, level});
}

// Adds the whole day to the wish-list
void Preferences::addDayOn(int nurseId, int day, PREF_LEVEL level){
  vector<Wish>& wishList = wishesOn_[nurseId][day];
  for(int s=1; s<nbShifts_; s++)           // Starts from 1 because it's a rest wish
    wishList.push_back({s, level});
}

// Returns true if the nurse wants that shift on
bool Preferences::wantsTheShiftOn(int nurseId, int day, int shift) const{
  return wantsTheShiftOnLevel(nurseId, day, shift) != -1;
}

// Returns level if the nurse wants that shift on : -1 otherwise
int Preferences::wantsTheShiftOnLevel(int nurseId, int day, int shift) const{
  // If the day is not in the wish-list, return -1
  auto itM = wishesOn_.at(nurseId).find(day);
  if(itM == wishesOn_.at(nurseId).end())
    return -1;
  // If the shift is not in the wish-list for that day, return -1
  for(const Wish &wish: itM->second)
    if(wish.shift == shift) return wish.level;
  return -1;
}

// True if the nurse wants the whole day on
bool Preferences::wantsTheDayOn(int nurseId, int day) const{
  // If the day is not in the wish-list, return false
  auto itM = wishesOn_.at(nurseId).find(day);
  if(itM == wishesOn_.at(nurseId).end())
    return false;
  // Set does not repeat its elements. Wants the day on if and only if all shifts off (-1 because REST does not appear)
  return itM->second.size() == nbShifts_-1;
}

// Total number of shifts on that the nurse wants
int Preferences::howManyShiftsOn(int nurseId) const {
  int nbShiftsOn = 0;

  // look at every wishes of the nurse
  for( const auto & dayOn: wishesOn_.at(nurseId)){
    nbShiftsOn += dayOn.second.size();
  }
  return nbShiftsOn;
}

// Number of whole days on that the nurse wants
int Preferences::howManyDaysOn(int nurseId, int dayMin, int dayMax) const {
  int nbDayOn = 0;
  // look at every wishes of the nurse
  for( const auto & dayOn: wishesOn_.at(nurseId)){
    nbDayOn += ( (dayOn.first >= dayMin) && (dayOn.first <= dayMax) && (dayOn.second.size() == nbShifts_-1));
  }
  return nbDayOn;
}

// add another week preferences at the end of the current one
//
void Preferences::push_back(PPreferences pPref){
  // check if same scenario
  if( (nbShifts_ != pPref->nbShifts_) || (nbNurses_ != pPref->nbNurses_) ){
    string error = "Preferences are not compatible";
    Tools::throwError(error.c_str());
  }

  //update the whishes off
  for(const auto &pWishes: pPref->wishesOff_)
    for(const auto &pWishes2: pWishes.second) {
      std::vector<Wish>& wishes = wishesOff_[pWishes.first][pWishes2.first+nbDays_];
      wishes.insert(wishes.end(), pWishes2.second.begin(), pWishes2.second.end());
    }

  //update the whishes on
  for(const auto &pWishes: pPref->wishesOn_)
    for(const auto &pWishes2: pWishes.second) {
      std::vector<Wish>& wishes = wishesOn_[pWishes.first][pWishes2.first+nbDays_];
      wishes.insert(wishes.end(), pWishes2.second.begin(), pWishes2.second.end());
    }

  // update the number of days
  nbDays_  += pPref->nbDays_;
}

// K the preferences relative to the nbDays first days
PPreferences Preferences::keep(int begin, int end) {

   PPreferences pPref = std::make_shared<Preferences>();

   for (int i=0; i < nbNurses_; i++) {
     for(pair<int,std::vector<Wish> > pair1: wishesOff_[i]){
       if (pair1.first >= begin && pair1.first < end) {
         pair<int,std::vector<Wish> > pair2(pair1.first-begin, pair1.second);
         pPref->wishesOff_[i].insert(pair2);
       }
     }

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
PPreferences Preferences::removeNFirstDays(int nbDays) {

	PPreferences pPref = std::make_shared<Preferences>();

	for (int i=0; i < nbNurses_; i++) {
    for(pair<int,std::vector<Wish> > pair1: wishesOff_[i]){
      if (pair1.first >= nbDays) {
        pair<int,std::vector<Wish> > pair2(pair1.first-nbDays_, pair1.second);
        pPref->wishesOff_[i].insert(pair2);
      }
    }

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
string Preferences::toString(PScenario pScenario) const{
	std::stringstream rep;
  rep << "# Preferences:" << std::endl;
  rep << "# Wishes off:" << std::endl;
  for(const auto & pWishes: wishesOff_)
    for(const auto & pWishes2: pWishes.second) {
      rep <<  "      | " << pWishes.first << ":" << pWishes2.first << "  ->  ";
      for(const Wish& w : pWishes2.second)
        rep << (pScenario ? pScenario->intToShift_[w.shift] : std::to_string(w.shift))
            << " (" << levelsToString.at(w.level) <<")\t";
      rep << std::endl;
    }
  rep << "# Wished on:" << std::endl;
  for(const auto & pWishes: wishesOn_)
    for(const auto & pWishes2: pWishes.second) {
    rep <<  "      | " << pWishes.first << ":" << pWishes2.first << "  ->  ";
    for(const Wish& w : pWishes2.second)
      rep << (pScenario ? pScenario->intToShift_[w.shift] : std::to_string(w.shift))
          << " (" << levelsToString.at(w.level) <<")\t";
    rep << std::endl;
  }
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
Nurse::Nurse(int id, string name, int nbSkills, vector<int> skills, PConstContract contract) :
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
