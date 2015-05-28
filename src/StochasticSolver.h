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

enum RankingStrategy {RK_MEAN, RK_SCORE, RK_NONE};

class StochasticSolverOptions{

public:

	StochasticSolverOptions(){

		SolverParam gp;
		generationParameters_ = gp;
		generationParameters_.weightStrategy_ =  RANDOMMEANMAX;


		SolverParam ep;
		evaluationParameters_ = ep;
		evaluationParameters_.stopAfterXSolution_ = 0;
		evaluationParameters_.weightStrategy_ =  BOUNDRATIO;

	};
	~StochasticSolverOptions(){};

	bool withEvaluation_ = true;
	bool withIterativeDemandIncrease_ = false;

	bool generationCostPerturbation_ = true;
	bool evaluationCostPerturbation_ = true;

	bool withResolveForGeneration_ = false;
	Algorithm generationAlgorithm_ = GENCOL;
	bool withResolveForEvaluation_ = true;
	Algorithm evaluationAlgorithm_ = GENCOL;
	RankingStrategy rankingStrategy_ = RK_SCORE;

	int totalTimeLimitSeconds_ = LARGE_TIME;

	int nExtraDaysGenerationDemands_ = 7;
	int nEvaluationDemands_ = 1;
	int nDaysEvaluation_ = 14;
	int nGenerationDemandsMax_ = 100;

	string logfile_ = "";

	SolverParam generationParameters_;
	SolverParam evaluationParameters_;


	int verbose = 1;

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

	StochasticSolver(Scenario* pScenario, StochasticSolverOptions options, vector<Demand*> demandHistory, double costPreviousWeeks=0);

	~StochasticSolver();






	//----------------------------------------------------------------------------
	//
	// SOLVE FUNCTIONS
	// The one is general for the whole process
	//
	//----------------------------------------------------------------------------

	// Main function
	double solve(vector<Roster> initialSolution = {});

	//get the number of generated schedules
	//
	int getNbSchedules() { return schedules_.size(); }

protected:

	void init();

	// Options that characterize the execution of the stochastic solver
	StochasticSolverOptions options_;

	// Timer started at the creation of the solver and stopped at destruction
	Tools::Timer* timerTotal_;

	// Log file that can be useful when calling the solver through simulator
	Tools::LogOutput* pLogStream_;



	//----------------------------------------------------------------------------
	//
	// SUBSOLVE FUNCTIONS
	//
	//----------------------------------------------------------------------------
	// Solves the problem by generation + evaluation of scenarios
	void solveOneWeekGenerationEvaluation();
	// Does everything for one schedule (for one week): Includes generation,
	// evaluation of the score, and update of the rankings and data.
	// Returns false if time has run out
	bool addAndSolveNewSchedule();
	// Iterative solution process in which the week is first solved by itsef,
	// before adding one perturbebd week demand and solving the new extended
	// demand demand until no time is left
	void solveIterativelyWithIncreasingDemand();
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
	// For reusable solvers
	Solver * pReusableGenerationSolver_;
	vector<vector<Roster> > schedules_;
	vector<vector<State> > finalStates_;

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
	vector<Solver*> pReusableEvaluationSolvers_;
	vector<map<double, set<int> > > schedulesFromObjectiveByEvaluationDemand_;
	vector<map<double, set<int> > > schedulesFromObjectiveByEvaluationDemandGreedy_;
	// Scores
	vector<double> theScores_;
	vector<double> theScoresGreedy_;


	int bestSchedule_;
	double bestScore_;
	double costPreviousWeeks_;

	// Return a solver with the algorithm specified for schedule EVALUATION
	Solver * setEvaluationWithInputAlgorithm(Demand* pDemand, vector<State> * stateEndOfSchedule);
	// Initialization
	void initScheduleEvaluation(int sched);
	// Evaluate 1 schedule and store the corresponding detailed results (returns false if time has run out)
	bool evaluateSchedule(int sched);
	// Recompute all scores after one schedule evaluation
	void updateRankingsAndScores(RankingStrategy strategy);
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
