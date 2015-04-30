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

#include "Scenario.h"
#include "Demand.h"

/* namespace usage */
using std::map;
using std::pair;
using std::string;
using std::vector;

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
	DemandGenerator(int nbDemands, vector<Demand*> demands, Scenario* pScenario):
		nbDemands_(nbDemands), demandHistory_(demands), pScenario_(pScenario) {
	}
	~DemandGenerator();

public:
	// basic getters
	int nbDemands() {return nbDemands_;}

	// check the feasibility of a demand scenario
	bool checkDemandFeasibility(Demand* pDemand);

	// generate nbScenarios_ through perturbations of the demand history
	vector<Demand*> generatePerturbedDemands();

protected:
	// number of demand scenarios that should be generated
	int nbDemands_;

	// demand history from which the random scenarios should be generated
	vector<Demand*> demandHistory_;

	// nurse rostering scenario under study
	// this attribute is necessary to check the feasibility of the generated demands
	Scenario* pScenario_;
};

 #endif
