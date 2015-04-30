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

StochasticSolver::StochasticSolver() {}
StochasticSolver::~StochasticSolver() {}

//-----------------------------------------------------------------------------
// Compute the weights o the violation of the min/max number of working days
// For now, the update depends only on the initial states and on the contract
// of the nurses, on the number of days on the demand, on the number of weeks
// already treated and on the number of weeks left
// The required data on the nurses is mostly computed in preprocessTheNurses
//-----------------------------------------------------------------------------

void StochasticSolver::computeWeightsTotalShifts() {

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
