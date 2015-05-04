/*
* SotchasticSolver.cpp
*
*  Created on: April 29, 2015
*      Author: jeremy omer
*/

#include "StochasticSolver.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

StochasticSolver(Scenario* pScenario, int nbRandomDemands_, int nbDaysRandom_,vector<Demand*> demandHistory):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
nbRandomDemands_(nbRandomDemands),nbDaysRandom_(nbDaysRandom) {

	// Generate the random demand
	// Keep in mind that the generated demand are for one week
	DemandGenerator generator(nbDemands,demandHistory,pScen);
	vector<Demand*> pRandomDemands_ = generator.generatePerturbedDemands();

	// Initialize the attributes relative to the penalties added in the model
	// to take the uncertainties into account
	for (int n = 0; n < pScenario_->nbNurses(); n++) {
		maxTotalWeekEndsAvg_.push_back(0.0);
		weightTotalWeekEndsAvg_.push_back(0.0);
		minTotalShifts_.push_back(0.0);
		maxTotalShifts_.push_back(0.0);
		minTotalShiftsAvg_.push_back(0.0);
		maxTotalShiftsAvg_.push_back(0.0);
		weightTotalShiftsAvg_.push_back(0.0);
	}

	// Compute the weights and bound for the total numbers of working days and
	// week-ends
	this->computeWeightsTotalShifts();


}
~StochasticSolver() {
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

	// The nurses must be preprocessed to retrieve the information relative to the
	// past activity of the nurses and to their capacity to work more in the future
	if (!isPrepreprocessedNurses_) this->preprocessTheNurses();

	// The important value to infer the importance of respecting the strict constraints 
	// on the total number of working days/week-ends is the remaining number of days/week-ends
	// after the demand currently treated
	int remainingDays = 7*pScenario_->nbWeeks()-7*(pScenario_->thisWeek()+1)-nbDaysRandom_;
	double factorRemainingDays = (double) remainingDays/(double)(7*pScenario_->nbWeeks());
	int remainingWeekEnds = pScenario_->nbWeeks()-(pScenario_->thisWeek()+1)-nbDaysRandom_/7;
	double factorRemainingWeekends = (double)remainingWeekEnds/(double)pScenario_->nbWeeks();

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
	  maxTotalWeekEndsAvg_[n] = factorRemainingWeekends*(double)pNurse->maxTotalWeekends();

	}

  // // short names
  // int nbNurses = pScenario_->nbNurses_;
  // int nbContracts = pScenario_->nbContracts_;
  //
  // // number of days that will be covered after the current demand
  // int nbDaysFuture = 7*(pScenario_->nbWeeks()-pScenario_->thisWeek())-pDemand_->nbDays_;
  //
  //
  // // get the average target number of working days for this week
  // double avgMinDays, avgMaxDays, avgMaxWeekends, minDays, maxDays, maxWeekends;
  // double factorWeek = (double)pScenario_->thisWeek()/(double) pScenario_->nbWeeks_;
  // avgMinDays = factorWeek * (double) nurse.minTotalShifts();
  // avgMaxDays = factorWeek * (double) nurse.maxTotalShifts();
  // avgMaxWeekends = factorWeek*(double) nurse.maxTotalWeekends();
  // minDays = (int) (avgMinDays-(1.0-factorWeek)*((double)nurse.minTotalShifts()/(double) pScenario_->nbWeeks_-(double)(day%7)));
  // maxDays = (int) (avgMaxDays+(1.0-factorWeek)*((double)nurse.maxTotalShifts()/(double) pScenario_->nbWeeks_-(double)(day%7))) +1;
  // maxWeekends = (int) avgMaxWeekends;
  //
  // if (state.totalDaysWorked_ >= maxDays) {
  //   cost += factorWeek*WEIGHT_TOTAL_SHIFTS;
  // }
  // if ( (day%7==5 || (day%7==6 && lastShift==0)) && state.totalWeekendsWorked_ >= maxWeekends) {
  //   cost += factorWeek*WEIGHT_TOTAL_WEEKENDS;
  // }
  // if the nurse is below the number of total assignments, getting a shift is
  // rewarded with a negative cost
  // RqJO: not sure yet

  // Add a small penalty when the considered skill is not the rarest in the
  // skill list of the nurse
  // double maxRarity = skillRarity_[skill];
  // double arbitraryWeight = 5;
  // for (int sk = 0; sk < nurse.nbSkills_; sk++) {
  //   maxRarity = std::max(skillRarity_[sk],maxRarity);
  // }
  // cost += (maxRarity/skillRarity_[skill]-1)*arbitraryWeight;
}

//-----------------------------------------------------------------------------
// Solve the problem
//-----------------------------------------------------------------------------

void StochasticSolver::solve() {

  pSolver_->solve();

}
