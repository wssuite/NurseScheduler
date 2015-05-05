/*
* SotchasticSolver.cpp
*
*  Created on: April 29, 2015
*      Author: jeremy omer
*/

#include "StochasticSolver.h"
#include "DemandGenerator.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

StochasticSolver::StochasticSolver(Scenario* pScenario, int nbRandomDemands, int nbDaysRandom,vector<Demand*> demandHistory):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
nbRandomDemands_(nbRandomDemands),nbDaysRandom_(nbDaysRandom),

maxTotalWeekendsAvg_(pScenario->nbNurses()), weightTotalWeekendsAvg_(pScenario->nbNurses()),
minTotalShifts_(pScenario->nbNurses()), maxTotalShifts_(pScenario->nbNurses()),
minTotalShiftsAvg_(pScenario->nbNurses()), maxTotalShiftsAvg_(pScenario->nbNurses()), weightTotalShiftsAvg_(pScenario->nbNurses())
{

	// Generate the random demand
	// Keep in mind that the generated demand are for one week
	DemandGenerator generator(nbRandomDemands,demandHistory,pScenario);
	vector<Demand*> pRandomDemands_ = generator.generatePerturbedDemands();

	// Compute the weights and bound for the total numbers of working days and
	// week-ends
	this->computeWeightsTotalShifts();


}
StochasticSolver::~StochasticSolver() {
	// delete the random demands
	while (!pRandomDemands_.empty()) {
		if (pRandomDemands_.back()) delete pRandomDemands_.back();
		pRandomDemands_.pop_back();
	}

	// delete the solvers
	while (!pRandomSolvers_.empty()) {
		if (pRandomSolvers_.back()) delete pRandomSolvers_.back();
		pRandomSolvers_.pop_back();
	}
}

//-----------------------------------------------------------------------------
// Generate extended demand that includes the demand for the current week and
// randomly generated demand for the 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Solve the problem
//-----------------------------------------------------------------------------

void StochasticSolver::solve() {

	// The weights must the computed for the total demand including the random demand
	// We thus initialize a solver with the augmented demand to preprocess the nursess
	// with the proper demand
	Demand demandExtended = *(pScenario_->pWeekDemand());
	demandExtended->push_back(*pRandomDemands_.back());
	Solver* pSolver = new Solver(pScenario_,&demandExtended,pScenario->pWeekPreferences(), pScenario->pInitialState());

	for (Solver* pSolver:pRandomSolvers_) {
  	pSolver->solve();
  }

}
