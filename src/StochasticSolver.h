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
  StochasticSolver(Scenario* pScenario, Algorithm algo);
  StochasticSolver(Scenario* pScenario, int nbRandomDemands, int nbDaysRandom, Algorithm algo, vector<Demand*> demandHistory);
  ~StochasticSolver();

protected:

  // Random scenarios that extrapolate the future weeks
  int nbRandomDemands_;
  vector<Demand*> pRandomDemands_;

  // Number of days in the random demands that should be taken into account
  // in the solver
  int nbDaysRandom_;

  // Vector of solvers that will be used to solve the random instances
  vector<Solver*> pRandomSolvers_;

  // Algorithm that is used to solve the stochastic demands
  Algorithm algorithm_;

protected:
  // Update the weights
  // For now, the update depends only on the initial states and on the contract
  // of the nurses, on the number of days on the demand, on the number of weeks
  // already treated and on the number of weeks left
  //
  void computeWeightsTotalShifts();

  // Return a solver with the algorithm specified in the attributes and the
  // in input
  //
  Solver* setSubSolverWithInputAlgorithm(Demand* pDemand);

  // Solve the problem with the algorithm in input and modified penalties for
  // the min/max of total working days or week-ends
  // There is no random scenario involved in this basic method
  //
  void solveOneWeekWithPenalties();

  // Solve the problem
  //
  virtual double solve(vector<Roster> solution = {});
};

#endif
