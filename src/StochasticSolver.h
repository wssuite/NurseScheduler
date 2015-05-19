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


class StochasticSolverOptions{

public:

	StochasticSolverOptions(){

		SolverParam gp;
		generationParameters_ = gp;

		SolverParam ep;
		evaluationParameters_ = ep;

	};
	~StochasticSolverOptions(){};

	bool withEvaluation_ = false;

	bool generationCostPerturbation_ = true;
	bool evaluationCostPerturbation_ = false;

	Algorithm generationAlgorithm_ = GENCOL;
	Algorithm evaluationAlgorithm_ = NONE;

	int totalTimeLimitSeconds_ = LARGE_TIME;

	int nExtraDaysGenerationDemands_ = 7;

	int nEvaluationDemands_ = 0;
	int nDaysEvaluation_ = 0;
	int nGenerationDemandsMax_ = 1;

	SolverParam generationParameters_;
	SolverParam evaluationParameters_;

};



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

	StochasticSolver(Scenario* pScenario, StochasticSolverOptions options, vector<Demand*> demandHistory);

	~StochasticSolver();






	//----------------------------------------------------------------------------
	//
	// SOLVE FUNCTIONS
	// The one is general for the whole process
	//
	//----------------------------------------------------------------------------

	// Main function
	double solve(vector<Roster> initialSolution = {});

protected:

	void init();
	StochasticSolverOptions options_;



	//----------------------------------------------------------------------------
	//
	// SUBSOLVE FUNCTIONS
	//
	//----------------------------------------------------------------------------
	// Solves the problem by generation + evaluation of scenarios
	void solveOneWeekGenerationEvaluation();
	// Does everything for one schedule (for one week): Includes generation,
	// evaluation of the score, and update of the rankings and data.
	void addAndSolveNewSchedule();
	// Solves the problem by generating a schedule + using cost penalties
	void solveOneWeekNoGenerationEvaluation();
	// Special case of the last week
	void solveOneWeekWithoutPenalties();


	//----------------------------------------------------------------------------
	//
	// GENERATION OF DEMANDS FOR THE CURRENT WEEK (=FOR SCHEDULE GENERATION)
	// Note that these demands share a common first week which is the week we currently try to solve.
	//
	//----------------------------------------------------------------------------

	// History
	vector<Demand *> demandHistory_;
	// Number of demands generated
	int nGenerationDemands_;
	// Vector of random demands that are used to GENERATE the schedules
	vector<Demand*> pGenerationDemands_;
	// Generate a new demand for generation
	void generateSingleGenerationDemand();



	//----------------------------------------------------------------------------
	//
	// GENERATION OF SCENARIOS FOR THE FUTURE (=FOR SCHEDULE EVALUATION)
	//
	//----------------------------------------------------------------------------

	// Vector of random demands that are used to EVAULATE the generated schedules
	vector<Demand*> pEvaluationDemands_;
	// Generate the schedules that are used for evaluation
	void generateAllEvaluationDemands();



	//----------------------------------------------------------------------------
	//
	// GENERATION OF SCHEDULES
	// A solution is a potential candidate to be the chosen schedule for the week we are solving.
	// A result is, given a solution and a potential future, the value obtained for that couple (solution,demand) [i.e. LP bound for instance]
	// A score is, given a solution, the average score it obtains, compared to the other solutions (the precise meaning of "score" should be better defined)
	//
	//----------------------------------------------------------------------------

	// Schedules
	int nSchedules_;
	vector<Solver*> pGenerationSolvers_;

	// Return a solver with the algorithm specified for schedule GENERATION
	Solver * setGenerationSolverWithInputAlgorithm(Demand* pDemand);
	// Generate a new schedule
	void generateNewSchedule();



	//----------------------------------------------------------------------------
	//
	// EVALUATION OF SCHEDULES
	//
	//----------------------------------------------------------------------------

	// Empty preferences -> only 1 to avoid multiplying them
	Preferences * pEmptyPreferencesForEvaluation_;
	// Evaluation
	vector<vector<Solver*> > pEvaluationSolvers_;
	vector<map<double, set<int> > > schedulesFromObjectiveByEvaluationDemand_;
	// Scores
	vector<double> theScores_;
	int bestSchedule_;
	double bestScore_;

	// Return a solver with the algorithm specified for schedule EVALUATION
	Solver * setEvaluationWithInputAlgorithm(Demand* pDemand, vector<State> * stateEndOfSchedule);
	// Initialization
	void initScheduleEvaluation(int sched);
	// Evaluate 1 schedule and store the corresponding detailed results
	void evaluateSchedule(int sched);
	// Recompute all scores after one schedule evaluation
	void updateRankingsAndScores();
	// Getter
	double valueOfEvaluation(int sched, int evalDemand){return pEvaluationSolvers_[sched][evalDemand]->solutionCost();}




protected:

	// Update the weights
	// For now, the update depends only on the initial states and on the contract
	// of the nurses, on the number of days on the demand, on the number of weeks
	// already treated and on the number of weeks left
	//
	void computeWeightsTotalShifts();

	Solver * setSubSolverWithInputAlgorithm(Demand* pDemand, Algorithm algo);










};

#endif
