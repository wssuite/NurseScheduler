/*

* StochasticSolver.h
*
*  Created on: April 29, 2015
*      Author: jeremy omer
*/

#ifndef __StochasticSolver__
#define __StochasticSolver__

#include "Solver.h"
#include "MasterProblem.h"


//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

class StochasticSolver:public Solver {

public:
  StochasticSolver(Scenario* pScenario, int nbRandomDemands_, int nbDaysRandom_);
  ~StochasticSolver();

protected:

  // Number of worked week-ends below which there is no penalty for the 
  // total number of working week-ends
  // This interval is computed from the max number of working week-ends averaged
  // over the number of remaining weeks
  vector<double> maxTotalWeekEndsAvg_;

  // Penalties for the number of working weekends on the current period
  // (for each nurse)
  vector<double> weightTotalWeekEndsAvg_;

  // Interval inside of which there is no penalty for the total number of
  // working days (for each nurse)
  // This interval is computed from the max/min number of working days averaged
  // over the number of remaining weeks
  vector<double> minTotalShiftsAvg_;
  vector<double> maxTotalShiftsAvg_;

  // Penalties for values outside of [minTotalShiftsAvg_,maxTotalShiftsAvg_]
  vector<double> weightTotalShiftsAvg_;

  // Interval outside of which the violation of the total number of working
  // days is penalized with the complete weight
  vector<double> minTotalShifts_;
  vector<double> maxTotalShifts_;

  // Random scenarios that extrapolate the future weeks
  int nbRandomDemands_;
  vector<Demand*> pRandomDemands_;

  // Number of days in the random demands that should be taken into account
  // in the solver
  int nbDaysRandom_;

  // Vector of solvers that will be used to solve the random instances
  vector<Solver*> pRandomsSolvers_;

  // Algorithm that is used to solve the stochastic demands
  Algorithm algo_;

protected:
  // Update the weights
  // For now, the update depends only on the initial states and on the contract
  // of the nurses, on the number of days on the demand, on the number of weeks
  // already treated and on the number of weeks left
  //
  void computeWeightsTotalShifts();

  // Solve the problem
  //
  virtual void solve();
};

#endif
