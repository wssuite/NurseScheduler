/*
* SotchasticSolver.cpp
*
*  Created on: April 29, 2015
*      Author: jeremy omer
*/

#include "StochasticSolver.h"
#include "DemandGenerator.h"
#include "Greedy.h"
#include "MasterProblem.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------
StochasticSolver::StochasticSolver(Scenario* pScenario, Algorithm algo):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
nbRandomDemands_(0),nbDaysRandom_(0), algorithm_(algo) {
}

StochasticSolver::StochasticSolver(Scenario* pScenario, int nbRandomDemands, int nbDaysRandom,Algorithm algo,vector<Demand*> demandHistory):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
nbRandomDemands_(nbRandomDemands),nbDaysRandom_(nbDaysRandom), algorithm_(algo)
{
	// Generate the random demand
	// Keep in mind that the generated demand are for one week
	DemandGenerator generator(nbRandomDemands,nbDaysRandom,demandHistory,pScenario);
	vector<Demand*> pRandomDemands_ = generator.generatePerturbedDemands();
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

Solver* StochasticSolver::setSubSolverWithInputAlgorithm(Demand* pDemand) {

	Solver* pSolver;
	switch(algorithm_){
	case GREEDY:
		pSolver = new Greedy(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState());
		break;
	case GENCOL:
		pSolver = new MasterProblem(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState(), S_BCP);
		break;
	default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}
	return pSolver;
}


//-----------------------------------------------------------------------------
// Solve the problem
//-----------------------------------------------------------------------------

void StochasticSolver::solve() {

	// create solvers with the extended demands obtained by appending the random
	// demands to the current weekly demand
	// for (Demand* pRandomDemand: pRandomDemands_) {
	// 	Demand* pExtendedDemand = new Demand(*(pScenario_->pWeekDemand()));
	// 	pExtendedDemand->push_back(pRandomDemand);
	//
	// 	Solver* pSolver = NULL;
	// 	switch(algorithm_){
	// 	case GREEDY:
	// 		pSolver = new Greedy(pScenario_, pExtendedDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState());
	// 		break;
	// 	case GENCOL:
	// 		pSolver = new MasterProblem(pScenario_, pExtendedDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState(), S_BCP);
	// 		break;
	// 	default:
	// 		Tools::throwError("The algorithm is not handled yet");
	// 		break;
	// 	}
	//
	// 	pSolver->computeWeightsTotalShiftsForStochastic();
	// 	pRandomSolvers_.push_back(pSolver);
	// }
	//
	// for (Solver* pSolver:pRandomSolvers_) {
  // 	pSolver->solve();
  // }

	this->solveOneWeekWithPenalties();

}

void StochasticSolver::solveOneWeekWithPenalties() {

	Solver* pSolver = setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand());
	pSolver->computeWeightsTotalShiftsForStochastic();
	pSolver->solve();
	solution_ = pSolver->getSolution();

}
