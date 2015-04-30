//
//  Roster.h
//  RosterDesNurses
//

#ifndef __Solver__
#define __Solver__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"
#include "Roster.h"

using std::vector;


//-----------------------------------------------------------------------------
// struct task:
// a task is a shift performed on a given day with a given skill
//
//-----------------------------------------------------------------------------
typedef struct {int shift; int skill;} task;

//-----------------------------------------------------------------------------
//
//  C l a s s   R o s t e r
//
//  Schedule of a single nurse
//
//-----------------------------------------------------------------------------

class Roster{

public:

  // Default constructor
  //
  Roster():nbDays_(0), firstDay_(0) {}

  // Constructor form no particular planning
  //
  Roster(int nbDays, int firstDay);

  // Constructor: initialize planning from an input set of shifts for the nurse
  //
  Roster(int nbDays, int firstDay, vector<int> shifts);

  // Constructor: initialize planning from an input set of shifts and skills
  //
  Roster(int nbDays, int firstDay, vector<int> shifts, vector<int> skills);

  // Destructor
  ~Roster();

private:

  // number of days in the roster and index of the first day
  //
  int nbDays_, firstDay_;

  // pointer to the scenario, the nurse under consideration and her wishes in
  // terms of days off
  // (the key of the map is the day and the value is the set of wishes)
  //
  const Scenario* pScenario_;
  Nurse* pNurse_;
  std::map< int,std::set<int> >* pWishesOff_;

  // vector containing for each day the shift assigned to the nurse
  // the vector contains exactly one element per day
  // the shift 0 corresponds to a rest
  //
  vector<int> shifts_;

  // vector containing for each day the shift assigned to the nurse
  // the vector contains exactly one element per day
  // if the nurse is resting, the skill has no importance
  //
  vector<int> skills_;

public:
  // Basic getters
  //
  int firstDay() {return firstDay_;}
  int nbDays() {return nbDays_;}
  int shift(int day) const {return shifts_[day];}
  int skill(int day) const {return skills_[day];}

  // initialize the roster
  //
  void init(int nbDays, int firstDay, int shiftDefault=0);

  // get a vector of consecutive states that will result from applying the
  // the roster from a given initial state
  //
  vector<State> getStates(const State& pStateIni);

  // assign a task at on a given day
  //
  void assignTask(int day, int shift, int skill=0);

  // add a roster at the end of the roster
  //
  void push_back(Roster& roster);

  // copy the input roster
  //
  void copy(Roster& roster);

};


#endif /* defined(__ATCSolver__CftSolver__)*/
