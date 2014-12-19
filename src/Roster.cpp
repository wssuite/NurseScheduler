#include "MyTools.h"
#include "Scenario.h"
#include "Nurse.h"
#include "Roster.h"


// Constructor
Roster::Roster(Scenario* pScenario, vector<Nurse>* pTheNurses):
        pScenario_(pScenario), pTheNurses_(pTheNurses){

	// Known attributes initializede from args
	//
	nbShifts_ = pScenario->nbShifts_;

	// Others initialized to 0
	//
	nbDays_ = 0;
	totalCostUnderStaffing_ = 0;

}

// Destructor
Roster::~Roster(){}
