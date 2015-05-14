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

StochasticSolver::StochasticSolver(Scenario* pScenario, Algorithm generationAlgorithm, Algorithm evaluationAlgorithm,
		int nExtraDaysGenerationDemands, int nEvaluationDemands, int nDaysEvaluation, int nMaxGenerationDemands,
		vector<Demand*> demandHistory):
		Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
		generationAlgorithm_(generationAlgorithm), evaluationAlgorithm_(evaluationAlgorithm), nExtraDaysGenerationDemands_(nExtraDaysGenerationDemands),
		nEvaluationDemands_(nEvaluationDemands), nDaysEvaluation_(nDaysEvaluation), nGenerationDemandsMax_(nMaxGenerationDemands), demandHistory_(demandHistory)
{
	bestScore_ = LARGE_SCORE;
	bestSchedule_ = -1;
	nGenerationDemands_ = 0;
	nSchedules_ = 0;
	pEmptyPreferencesForEvaluation_ = new Preferences(pScenario_->nbNurses(), nDaysEvaluation_, pScenario_->nbShifts());
}

StochasticSolver::~StochasticSolver(){
	// delete the demands used for generation
	while (!pGenerationDemands_.empty()) {
		if (pGenerationDemands_.back()) delete pGenerationDemands_.back();
		pGenerationDemands_.pop_back();
	}
	// delete the demands used for evaluation
	while (!pEvaluationDemands_.empty()) {
		if (pEvaluationDemands_.back()) delete pEvaluationDemands_.back();
		pEvaluationDemands_.pop_back();
	}
	// delete the solvers used for generation
	while (!pGenerationSolvers_.empty()) {
		if (pGenerationSolvers_.back()) delete pGenerationSolvers_.back();
		pGenerationSolvers_.pop_back();
	}
	// delete the solvers used for evaluation
	for(int n=0; n<pEvaluationSolvers_.size(); n++){
		while (!pEvaluationSolvers_[n].empty()) {
			if (pEvaluationSolvers_[n].back()) delete pEvaluationSolvers_[n].back();
			pEvaluationSolvers_[n].pop_back();
		}
	}
	// delete the empty preference list
	delete pEmptyPreferencesForEvaluation_;
}



//----------------------------------------------------------------------------
//
// SOLVE FUNCTIONS
// The one is general for the whole process
//
//----------------------------------------------------------------------------

// Main function
double StochasticSolver::solve(vector<Roster> initialSolution){

	// A. Special case of the last week
	//

	if(pScenario_->nbWeeks() == pScenario_->thisWeek()){
		solveOneWeekWithoutPenalties();
	}

	// B. Regular week when using the perturbation + generate single schedule
	//
	else if(evaluationAlgorithm_ == NONE){
		solveOneWeekWithPenalties();
	}

	// C. Regular week when using generation-evaluation
	else {
		solveOneWeekGenerationEvaluation();
	}

	return DBL_MAX;
}

// Does everything for the one week and only keeps the best schedule for it
void StochasticSolver::solveOneWeekWithPenalties() {
	Solver* pSolver = setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(), generationAlgorithm_);
	pSolver->computeWeightsTotalShiftsForStochastic();
//	pSolver->computeWeightsTotalShiftsForPrimalDual();
	pSolver->solve();
	solution_ = pSolver->getSolution();
   /* update nurse States */
	for(int n=0; n<pScenario_->nbNurses_; ++n){
	   theLiveNurses_[n]->roster_ = solution_[n];
	   theLiveNurses_[n]->buildStates();
	}
	status_ = pSolver->getStatus();
}

// Special case of the last week
void StochasticSolver::solveOneWeekWithoutPenalties(){
	Solver* pSolver = setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(), GENCOL);
	pSolver->solve();
	solution_ = pSolver->getSolution();
	/* update nurse States */
   for(int n=0; n<pScenario_->nbNurses_; ++n){
      theLiveNurses_[n]->roster_ = solution_[n];
      theLiveNurses_[n]->buildStates();
   }
	status_ = pSolver->getStatus();
}

// Solves the problem by generation + evaluation of scenarios
void StochasticSolver::solveOneWeekGenerationEvaluation(){
	int prevBestSchedule = -1;
	// TODO: Add a time constraint in the while condition
	while(nSchedules_<nGenerationDemandsMax_){
		addAndSolveNewSchedule();

		// TODO: write the output NOW so that it is not lost
		if(prevBestSchedule != bestSchedule_){

		}

		prevBestSchedule = bestSchedule_;
	}
}

//----------------------------------------------------------------------------
//
// GENERIC FUNCTION TO DO EVERYTHING FOR ONE SCHEDULE
// Includes generation, evaluation of the score, and update of the rankings
// and data.
//
//----------------------------------------------------------------------------

// Do everything for the new schedule (incl. generation, score, ranking)
void StochasticSolver::addAndSolveNewSchedule(){
	generateNewSchedule();
	evaluateSchedule(nSchedules_-1);
	updateRankingsAndScores();
}



//----------------------------------------------------------------------------
//
// GENERATION OF DEMANDS FOR THE CURRENT WEEK (=FOR SCHEDULE GENERATION)
// Note that these demands share a common first week which is the week we currently try to solve.
//
//----------------------------------------------------------------------------

// Generate a new demand for generation
void StochasticSolver::generateSingleGenerationDemand(){
	DemandGenerator dg (1, pScenario_->pWeekDemand()->nbDays_ + nExtraDaysGenerationDemands_, demandHistory_ , pScenario_);
	Demand * pSingleDemand = dg.generateSinglePerturbatedDemand();
	pGenerationDemands_.push_back( pScenario_->pWeekDemand()->append(pSingleDemand) );
	delete pSingleDemand;
}



//----------------------------------------------------------------------------
//
// GENERATION OF SCENARIOS FOR THE FUTURE (=FOR SCHEDULE EVALUATION)
//
//----------------------------------------------------------------------------

// Generate the schedules that are used for evaluation
void StochasticSolver::generateAllEvaluationDemands(){
	DemandGenerator dg (nEvaluationDemands_, nDaysEvaluation_, demandHistory_, pScenario_);
	pEvaluationDemands_ = dg.generatePerturbedDemands();
	// Initialize structures for scores
	for(int j=0; j<nEvaluationDemands_; j++){
		map<double, set<int> > m;
		schedulesFromObjectiveByEvaluationDemand_.push_back(m);
	}
}



//----------------------------------------------------------------------------
//
// GENERATION OF SCHEDULES
// A solution is a potential candidate to be the chosen schedule for the week we are solving.
// A result is, given a solution and a potential future, the value obtained for that couple (solution,demand) [i.e. LP bound for instance]
// A score is, given a solution, the average score it obtains, compared to the other solutions (the precise meaning of "score" should be better defined)
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule GENERATION
Solver* StochasticSolver::setGenerationSolverWithInputAlgorithm(Demand* pDemand){
	Solver* pSolver;
	switch(generationAlgorithm_){
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

// Generate a new schedule
void StochasticSolver::generateNewSchedule(){
	// A. Generate a demand that will be the origin of the scenario generation
	//
	generateSingleGenerationDemand();
	Demand * newDemand = pGenerationDemands_[nGenerationDemands_-1];

	// B. Solve this schedule (in a way that should be defined) so as to have a schedule
	//
	Solver* pGenSolver = setGenerationSolverWithInputAlgorithm( newDemand );
	pGenSolver->computeWeightsTotalShiftsForStochastic();
	pGenSolver->solve();

	// C. Update the data
	//
	pGenerationSolvers_.push_back(pGenSolver);
	nSchedules_ ++;
}



//----------------------------------------------------------------------------
//
// EVALUATION OF SCHEDULES
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule EVALUATION
Solver* StochasticSolver::setEvaluationWithInputAlgorithm(Demand* pDemand, vector<State>* stateEndOfSchedule){
	Solver* pSolver;
	switch(evaluationAlgorithm_){
	case GREEDY:
		pSolver = new Greedy(pScenario_, pDemand, pEmptyPreferencesForEvaluation_, stateEndOfSchedule);
		break;
	case GENCOL:
		pSolver = new MasterProblem(pScenario_, pDemand, pEmptyPreferencesForEvaluation_, stateEndOfSchedule, S_BCP);
		break;
	default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}
	return pSolver;
}

// Initialization
void StochasticSolver::initScheduleEvaluation(int sched){
	// Extend theEvaluationSolvers_
	vector<Solver*> v;
	for(int i=0; i<nEvaluationDemands_; i++){
		Solver * s;
		v.push_back(s);
	}
	pEvaluationSolvers_.push_back(v);

}

// Evaluate 1 schedule on all evaluation instances
void StochasticSolver::evaluateSchedule(int sched){
	for(int j=0; j<nEvaluationDemands_; j++){

		delete pEvaluationSolvers_[sched][j];
		pEvaluationSolvers_[sched][j] = setEvaluationWithInputAlgorithm(pEvaluationDemands_[j], new vector<State> (pGenerationSolvers_[sched]->getFinalStates()));
		pEvaluationSolvers_[sched][j]->evaluate();

		// Insert the solution cost and solution
		//
		// TODO : ici, arondi a l'entier -> peut etre modifie si besoin
		double currentCost = ((int) pEvaluationSolvers_[sched][j]->solutionCost());
		// If already in the costs -> add it to the set of schedules that found that cost
		if(schedulesFromObjectiveByEvaluationDemand_[j].find(currentCost) != schedulesFromObjectiveByEvaluationDemand_[j].end()){
			schedulesFromObjectiveByEvaluationDemand_[j].at(currentCost).insert(sched);
		}
		// Otherwise, add a new pair
		else{
			set<int> s; s.insert(sched);
			schedulesFromObjectiveByEvaluationDemand_[j].insert(pair<double, set<int> >( currentCost, s));
		}
	}
}

// Recompute all scores after one schedule evaluation
void StochasticSolver::updateRankingsAndScores(){

	vector<double> theNewScores;

	for(int j=0; j<nEvaluationDemands_; j++){
		int localRank = 1;
		map<double, set<int> > localCosts = schedulesFromObjectiveByEvaluationDemand_[j];
		for(map<double, set<int> >::iterator it = localCosts.begin(); it != localCosts.end(); ++it){
			for(int sched : it->second){
				theNewScores[sched] += ((double)(localRank + it->second.size() - 1)) / ((double) it->second.size());
				localRank += it->second.size();
			}
		}
	}

	theScores_ = theNewScores;

	bestScore_ = LARGE_SCORE;
	bestSchedule_ = -1;
	for(int i=0; i<nSchedules_; i++){
		if(theScores_[i] < bestScore_){
			bestScore_ = theScores_[i];
			bestSchedule_ = i;
		}
	}





}

















Solver* StochasticSolver::setSubSolverWithInputAlgorithm(Demand* pDemand, Algorithm algorithm) {

	Solver* pSolver;
	switch(algorithm){
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
