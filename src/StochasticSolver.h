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
  StochasticSolver();
  ~StochasticSolver();

protected:
  // Penalties for the total number of working days on the current period
  // (for each nurse)
  vector<double> weightTotalShiftsAvg_;

  // Penalties for the number of working weekends on the current period
  // (for each nurse)
  vector<double> weightTotalWeekEnds_;

  // Interval inside of which there is no penalty for the total number of
  // working days (for each nurse)
  vector<int> minTotalShiftsAvg_;
  vector<int> maxTotalShiftsAvg_;

  Solver* pSolver_;

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
