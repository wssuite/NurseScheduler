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
// Compute the weights o the violation of the min/max number of working days
// For now, the update depends only on the initial states and on the contract
// of the nurses, on the number of days on the demand, on the number of weeks
// already treated and on the number of weeks left
// The required data on the nurses is mostly computed in preprocessTheNurses
//-----------------------------------------------------------------------------

void StochasticSolver::computeWeightsTotalShifts() {

	// The weights must the computed for the total demand including the random demand
	// We thus initialize a solver with the augmented demand to preprocess the nursess
	// with the proper demand
	// Demand* pDemand = new Demand(pScenario_->pWeekDemand());
	// pDemand->push_back(*pRandomDemands_.back());
	// Solver* pSolver = new Solver(pScenario,pDemand,pScenario->pWeekPreferences(), pScenario->pInitialState());

	// The nurses must be preprocessed to retrieve the information relative to the
	// past activity of the nurses and to their capacity to work more in the future
	if (!isPreprocessedNurses_) this->preprocessTheNurses();

	// The important value to infer the importance of respecting the strict constraints 
	// on the total number of working days/week-ends is the remaining number of days/week-ends
	// after the demand currently treated
	int remainingDays = 7*pScenario_->nbWeeks()-7*(pScenario_->thisWeek()+1)-nbDaysRandom_;
	double factorRemainingDays = (double) remainingDays/(double)(7*pScenario_->nbWeeks());
	int remainingWeekends = pScenario_->nbWeeks()-(pScenario_->thisWeek()+1)-nbDaysRandom_/7;
	double factorRemainingWeekends = (double)remainingWeekends/(double)pScenario_->nbWeeks();

	// Compute the non-penalized intervals and the associated penalties 
	for (int n = 0; n < pScenario_->nbNurses(); n++) {
		LiveNurse* pNurse =  theLiveNurses_[n]; 

		// first compute the values relative to the average number of working days
		// the interval is larger for the first weeks and the associated penalty is smaller
		minTotalShiftsAvg_[n] = (1.0-factorRemainingDays)*pNurse->minAvgWorkDaysNoPenaltyTotalDays_;
		maxTotalShiftsAvg_[n] = (1.0+factorRemainingDays)*pNurse->maxAvgWorkDaysNoPenaltyTotalDays_;
		weightTotalShiftsAvg_[n] = (1.0-factorRemainingDays)*(double)WEIGHT_TOTAL_SHIFTS;

		// compute the interval that must be respected to have a chance of not paying 
		// penalties in the future
		minTotalShifts_[n] = pNurse->minWorkDaysNoPenaltyTotalDays_;
		maxTotalShifts_[n] = pNurse->maxWorkDaysNoPenaltyTotalDays_;

		// Number of worked week-ends below which there is no penalty for the 
	  // total number of working week-ends
	  // This interval is computed from the max number of working week-ends averaged
	  // over the number of remaining weeks
	  maxTotalWeekendsAvg_[n] = factorRemainingWeekends*(double)pNurse->maxTotalWeekends();

	}
}

//-----------------------------------------------------------------------------
// Solve the problem
//-----------------------------------------------------------------------------

void StochasticSolver::solve() {

	for (Solver* pSolver:pRandomSolvers_) {
  	pSolver->solve();
  }

}
