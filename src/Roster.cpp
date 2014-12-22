#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"
#include "Roster.h"


// Constructor
//
Roster::Roster(int nbDays):nbDays_(nbDays){

  // initialize the roster with a rest week
  for (int i = 0; i < nbDays; i++) {
    tasks_.push_back(std::pair<int,int> (0,0));
  }
};

// Destructor
Roster::~Roster(){}
