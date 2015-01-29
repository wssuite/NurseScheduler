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
  Roster() {}
  
  // Constructor form no particular planning
  //
  Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map< int,std::set<int> >* pWishesOff, const State& initialState);

  // Constructor: initialize planning from an input set of shifts for the nurse
  //
  Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map< int,std::set<int> >* pWishesOff, const State& initialState,
  vector<int> shifts);

  // Constructor: initialize planning from an input set of shifts and skills
  //
  Roster(int nbDays, int firstDay, Scenario* pScenario, Nurse* pNurse,
  std::map< int,std::set<int> >* pWishesOff, const State& initialState,
  vector<int> shifts, vector<int> skills);

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

  // vector containing for each day the state of the nurse
  // the size is the number of days of the roster plus one, since the initial
  // and the final states are of importance
  //
  vector<State> states_;

  // vectors of booleans that for each day, is equal to 1 if the nurse is about
  // to go from a working day to a day off or from a day off to a working day
  //
  vector<bool> switchOff_;

  // vectors of booleans that for each day, is equal to 1 if the nurse is about
  // to take another type of shift
  //
  vector<bool> switchShift_;

  // costs for the violation of soft constraints
  //
  vector<int> costConsDays_; // the same vector also accounts for consecutive days off
  vector<int> costConsShifts_;
  vector<int> costPreferences_;
  vector<int> costCompleteWeekEnd_;

  // vector of booleans equal to true if the corresponding hard contraint is
  // violated on each day
  //
  vector<bool> violationSuccShifts_; // forbidden successive shifts
  vector<bool> violationSkill_; // missing required skill

public:
  // assign a task at on a given day and update the states of the nurse
  //
  void assignTask(task t, int day);

  // check the satisfaction of the hard constraints and record the violations
  //
  void checkHardConstraints();

  // check the soft constraints and record the costs of the violations and the
  // remaining margin for the satisfied ones.
  //
  void checkSoftConstraints();

};


#endif /* defined(__ATCSolver__CftSolver__)*/
