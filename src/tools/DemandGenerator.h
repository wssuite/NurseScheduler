/*

* DemandGenerator.h
*
*  Created on: April 27, 2015
*      Author: jeremy omer
*/

#ifndef __DemandGenerator__
#define __DemandGenerator__

#include <map>
#include <string>
#include <vector>

#include "data/Scenario.h"
#include "data/Demand.h"


//-----------------------------------------------------------------------------
//
//	C l a s s  D e m a n d G e n e r a t o r
//
// Generate and store random demand scenarios for stochastic solution of the
// problem
//
//-----------------------------------------------------------------------------

class DemandGenerator{
public:
	// default constructor and destructor
	DemandGenerator(int nbDemands, int nbDays, std::vector<PDemand> demands, PScenario pScenario):
		nbDemandsToGenerate_(nbDemands), nbDaysInGeneratedDemands_(nbDays),demandHistory_(demands), pScenario_(pScenario),
	   rdm_(Tools::getANewRandomGenerator()) {
	}
	~DemandGenerator();

public:

	// check the feasibility of a demand scenario
	bool checkDemandFeasibility(PDemand pDemand);

	// generate nbScenarios_ through perturbations of the demand history
	std::vector<PDemand> generatePerturbedDemands();

	// generate 1 demand through perturbations of the demand history
	PDemand generateSinglePerturbatedDemand(bool checkFeasibility = true);

protected:
	// number of demand scenarios that should be generated
	int nbDemandsToGenerate_;

	// number of days that must be considered in each generated demand
	int nbDaysInGeneratedDemands_;

	// demand history from which the random scenarios should be generated
	std::vector<PDemand> demandHistory_;

	// nurse rostering scenario under study
	// this attribute is necessary to check the feasibility of the generated demands
	PScenario pScenario_;

	  //random generator
	  std::minstd_rand rdm_;
};

 #endif
