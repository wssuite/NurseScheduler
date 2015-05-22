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

// #define COMPARE_EVALUATIONS

//-----------------------------------------------------------------------------
//
//  C l a s s   S t o c h a s t i c S o l v e r
//
//  Solves the problem with uncertainty on the demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

StochasticSolver::StochasticSolver(Scenario * pScenario, StochasticSolverOptions options, vector<Demand*> demandHistory):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
options_(options), demandHistory_(demandHistory){
	std::cout << "# New stochastic solver created!" << endl;

	bestScore_ = LARGE_SCORE;
	bestSchedule_ = -1;
	nGenerationDemands_ = 0;
	nSchedules_ = 0;

	pEmptyPreferencesForEvaluation_ = new Preferences(pScenario_->nbNurses(), options_.nDaysEvaluation_, pScenario_->nbShifts());

	// create the timer that records the life time of the solver and start it
	timerTotal_ = new Tools::Timer();
    timerTotal_->init();
    timerTotal_->start();

    // initialize the log output
    pLogStream_ = new Tools::LogOutput(options_.logfile_);

    if (!options_.generationParameters_.logfile_.empty()) {
		FILE * pFile;
		pFile = fopen (options_.generationParameters_.logfile_.c_str(),"w");
		fprintf(pFile,"HERE COME THE LOGS OF THE MASTER PROBLEM\n\n");
		fclose(pFile);
   }
}

StochasticSolver::~StochasticSolver(){
	// kill the timer and display the total time spent in the algorithm
	timerTotal_->stop();
	(*pLogStream_) << "Total time spent in the algorithm : " << timerTotal_->dSinceInit() << endl;
	delete timerTotal_;

	// Delete the output
    delete pLogStream_;

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

	// Special case of the last week -> always to optimality with no time limit
	//
	if(pScenario_->nbWeeks()-1 == pScenario_->thisWeek()){
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solving week no. " << pScenario_->thisWeek() << " as the LAST WEEK (hence, to optimality !)" << std::endl;
		// General options
		options_.nExtraDaysGenerationDemands_ = 0;
		options_.withEvaluation_ = false;
		options_.generationCostPerturbation_ = false;
		// Options for the generation algo (-> optimality, no time limit, write every solution)
		options_.generationParameters_.solveToOptimality_ = true;
		//options_.generationParameters_.maxSolvingTimeSeconds_ = LARGE_TIME;
		options_.generationParameters_.printEverySolution_ = true;
	}

	// A. No generation-evaluation
	//
	if(! options_.withEvaluation_){
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solving week no. " << pScenario_->thisWeek() << " with PERTURBATIONS." << std::endl;
		solveOneWeekNoGenerationEvaluation();
		while(status_ == INFEASIBLE or status_ == UNSOLVED){
			(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Status is INFEASIBLE or UNSOLVED..." << std::endl;
			(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solving week no. " << pScenario_->thisWeek() << " with PERTURBATIONS -> trying again." << std::endl;
			solveOneWeekNoGenerationEvaluation();
		}
	}

	// B. Generation-evaluation
	else {
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solving week no. " << pScenario_->thisWeek() << " with GENERATION-EVALUATION." << std::endl;
		solveOneWeekGenerationEvaluation();
	}

	/* update nurse States */
	for(int n=0; n<pScenario_->nbNurses_; ++n){
		theLiveNurses_[n]->roster_ = solution_[n];
		theLiveNurses_[n]->buildStates();
	}

	return solutionCost();
}

// Does everything for the one week and only keeps the best schedule for it

void StochasticSolver::solveOneWeekNoGenerationEvaluation() {


	Solver * pSolver;
	// Need to extend the current demand?
	//
	if(options_.nExtraDaysGenerationDemands_ > 0){
		generateSingleGenerationDemand();
		pSolver = setSubSolverWithInputAlgorithm(pGenerationDemands_[0], options_.generationAlgorithm_);
	} else {
		pSolver = setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(), options_.generationAlgorithm_);
	}

	// Need to perturb the costs?
	//
	if(options_.generationCostPerturbation_){
//		pSolver->computeWeightsTotalShiftsForStochastic();
		pSolver->computeWeightsTotalShiftsForPrimalDual();
	}

	// Solve
	//
	(*pLogStream_) << "# Solve without evaluation\n";
	pSolver->solve(options_.generationParameters_);
	solution_ = pSolver->getSolution();
	status_ = pSolver->getStatus();
}

// Special case of the last week
void StochasticSolver::solveOneWeekWithoutPenalties(){
	Solver* pSolver = setSubSolverWithInputAlgorithm(pScenario_->pWeekDemand(), GENCOL);
	pSolver->solve(options_.generationParameters_);
	solution_ = pSolver->getSolution();
	status_ = pSolver->getStatus();
}

// Solves the problem by generation + evaluation of scenarios
void StochasticSolver::solveOneWeekGenerationEvaluation(){

	while(nSchedules_<options_.nGenerationDemandsMax_){

		// get the time left to solve another schedule
		double timeLeft = options_.totalTimeLimitSeconds_;
		if (nSchedules_ > 0) {
			timeLeft = options_.totalTimeLimitSeconds_-timerTotal_->dSinceInit();
			if (timeLeft < 1.0) break;
			options_.generationParameters_.maxSolvingTimeSeconds_  = (timeLeft-1.0)/2.0;
			options_.evaluationParameters_.maxSolvingTimeSeconds_  = (timeLeft-1.0)/(2.0*options_.nEvaluationDemands_);
		}
		(*pLogStream_) << "# Time left: " << timeLeft << std::endl;

		// This the main function that finds a new schedule and evaluates it
		addAndSolveNewSchedule();

		// Get the new best schedule
		//
		int newBestSchedule = -1;
		double newBestScore = LARGE_SCORE;
		for(int i=0; i<nSchedules_; i++){
			if(theScores_[i] < newBestScore){
				newBestScore = theScores_[i];
				newBestSchedule = i;
			}
		}

		// write the output NOW so that it is not lost
		//
		bestScore_ = newBestScore;
		if(newBestSchedule != bestSchedule_){
			bestSchedule_ = newBestSchedule;

			solution_ = pGenerationSolvers_[bestSchedule_]->getSolutionAtDay(6);
			loadSolution(solution_);
			Tools::LogOutput outStream(options_.generationParameters_.outfile_);
			outStream << solutionToString();

			(*pLogStream_) << "# New best is schedule n°" << bestSchedule_ << " (score: " << bestScore_ << ")" << std::endl;
			(*pLogStream_) << "# The new best solution was written in " << options_.generationParameters_.outfile_ << std::endl;
		}
		else {
			(*pLogStream_) << "# Best schedule did not change and is no. " << bestSchedule_ << " (score: " << bestScore_ << ")" << std::endl;

		}

		// Stop if the average time per schedule is smaller than the time left
		// no need to start building a schedule if there a risk that we won't
		// have any time left to evaluate it
		timeLeft = options_.totalTimeLimitSeconds_-timerTotal_->dSinceInit();
		double avgTimePerSchedule = timerTotal_->dSinceInit()/nSchedules_;
		if (timeLeft < avgTimePerSchedule) break;
	}
	#ifdef COMPARE_EVALUATIONS
	for(int i=0; i<nSchedules_; i++){
		(*pLogStream_) << " The score of schedule " << i << ". GENCOL : " << theScores_[i] << " ; GREEDY : " << theScoresGreedy_[i] << std::endl;
	}
	#endif
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
	if(nSchedules_ == 1){
		generateAllEvaluationDemands();
	}
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

	int nDaysInDemand = options_.nExtraDaysGenerationDemands_;
	bool isFeasible = false;
	Demand * pCompleteDemand;
	Demand * pSingleDemand;

	(*pLogStream_) << "# Generating new schedule" << std::endl;

	while(!isFeasible){
		DemandGenerator dg (1, nDaysInDemand, demandHistory_ , pScenario_);
		pSingleDemand = dg.generateSinglePerturbatedDemand();
		pCompleteDemand = pScenario_->pWeekDemand()->append(pSingleDemand);
		isFeasible = dg.checkDemandFeasibility(pCompleteDemand);
		if(!isFeasible){
			delete pCompleteDemand;
			delete pSingleDemand;
		}
	}

	pGenerationDemands_.push_back( pScenario_->pWeekDemand()->append(pSingleDemand) );
	nGenerationDemands_ ++;
	delete pSingleDemand;
	(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Generation demand no. " << nGenerationDemands_ << " created (over " << pGenerationDemands_[nGenerationDemands_-1]->nbDays_ << " days)." << std::endl;
}



//----------------------------------------------------------------------------
//
// GENERATION OF SCENARIOS FOR THE FUTURE (=FOR SCHEDULE EVALUATION)
//
//----------------------------------------------------------------------------

// Generate the schedules that are used for evaluation
void StochasticSolver::generateAllEvaluationDemands(){
	DemandGenerator dg (options_.nEvaluationDemands_, options_.nDaysEvaluation_, demandHistory_, pScenario_);
	pEvaluationDemands_ = dg.generatePerturbedDemands();
	// Initialize structures for scores
	for(int j=0; j<options_.nEvaluationDemands_; j++){
		map<double, set<int> > m;
		schedulesFromObjectiveByEvaluationDemand_.push_back(m);

		#ifdef COMPARE_EVALUATIONS
		schedulesFromObjectiveByEvaluationDemandGreedy_.push_back(m);
		#endif

		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Evaluation demand no. " << j << " created (over " << options_.nDaysEvaluation_ << " days)." << std::endl;
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
	switch(options_.generationAlgorithm_){
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
	pGenSolver->solve(options_.generationParameters_);

	// C. Update the data
	//
	pGenerationSolvers_.push_back(pGenSolver);
	nSchedules_ ++;

	// D. Display
	//
	(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Candidate schedule no. " << (nSchedules_-1) << " generated: (length: " << pGenerationSolvers_[nSchedules_-1]->getNbDays() << " days)" << std::endl;
}



//----------------------------------------------------------------------------
//
// EVALUATION OF SCHEDULES
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for schedule EVALUATION
Solver* StochasticSolver::setEvaluationWithInputAlgorithm(Demand* pDemand, vector<State> * stateEndOfSchedule){
	Solver* pSolver;
	Scenario * pScen = new Scenario (*pScenario_);
	pScen->linkWithDemand(new Demand ());

	// update the scenario to treat next week
	pScen->updateNewWeek(pDemand, *pEmptyPreferencesForEvaluation_, *stateEndOfSchedule);

	switch(options_.evaluationAlgorithm_){
	case GREEDY:
		pSolver = new Greedy(pScen, pDemand, pEmptyPreferencesForEvaluation_, stateEndOfSchedule);
		break;
	case GENCOL:
		pSolver = new MasterProblem(pScen, pDemand, pEmptyPreferencesForEvaluation_, stateEndOfSchedule, S_BCP);
		break;
	default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}
	return pSolver;
}

// Initialization
void StochasticSolver::initScheduleEvaluation(int sched){
	// Extend pEvaluationSolvers_
	vector<Solver*> v;
	for(int j=0; j<options_.nEvaluationDemands_; j++){
		Solver * s;
		v.push_back(s);
	}
	pEvaluationSolvers_.push_back(v);
}

// Evaluate 1 schedule on all evaluation instances
void StochasticSolver::evaluateSchedule(int sched){

	(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Evaluation of the schedule no. " << sched << std::endl;

	#ifdef COMPARE_EVALUATIONS
	vector<Solver*> pGreedyEvaluators;
	for(int j=0; j<options_.nEvaluationDemands_; j++){
		Solver * s;
		pGreedyEvaluators.push_back(s);
	}
	#endif

	initScheduleEvaluation(sched);
	vector<State> initialStates = pGenerationSolvers_[sched]->getStatesOfDay(6);
	for (int i = 0; i < pScenario_->nbNurses_; i++) {
		initialStates[i].dayId_ = 0;
	}

	// set the time per evaluation to the ratio of the time left over the number of evaluations
	double timeLeft = options_.totalTimeLimitSeconds_-timerTotal_->dSinceInit();
	options_.evaluationParameters_.maxSolvingTimeSeconds_ = (timeLeft-1.0)/(double)options_.nEvaluationDemands_;


	for(int j=0; j<options_.nEvaluationDemands_; j++){

		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Starting evaluation of schedule no. " << sched << " over evaluation demand no. " << j << std::endl;

		pEvaluationSolvers_[sched][j] = setEvaluationWithInputAlgorithm(pEvaluationDemands_[j], & initialStates);

		#ifdef COMPARE_EVALUATIONS
		options_.evaluationAlgorithm_ = GREEDY;
		pGreedyEvaluators[j] = setEvaluationWithInputAlgorithm(pEvaluationDemands_[j], & initialStates);
		options_.evaluationAlgorithm_ = GENCOL;
		#endif

		if(pEvaluationSolvers_[sched][j]->getNbDays() + (7*pScenario_->thisWeek()+1) < 7* pScenario_->nbWeeks_){
			pEvaluationSolvers_[sched][j]->computeWeightsTotalShiftsForStochastic();

			#ifdef COMPARE_EVALUATIONS
			pGreedyEvaluators[j]->computeWeightsTotalShiftsForStochastic();
			#endif
		}

		// If the first schedule took more than half the available time to be solved, there is
		// no need to go through evaluation since there will not be any time left for 
		// a second schedule
		bool isTimeForMoreThanOneSchedule = true;
		if ( (nSchedules_==0) && (timerTotal_->dSinceInit() > options_.totalTimeLimitSeconds_/2.0) ) {
			isTimeForMoreThanOneSchedule = false;
		}

		// Only perform the evaluation if the schedule is feasible and 
		// there is time for more than one schedule
		double currentCost, currentCostGreedy;
		if (pGenerationSolvers_[sched]->getStatus() == INFEASIBLE && isTimeForMoreThanOneSchedule) {
			currentCost = 1.0e6;
			currentCostGreedy = 1.0e6;
		}
		else {
			// Perform the actual evaluation on demand j by running the chosen algorithm
			// TODO : ici, arondi a l'entier -> peut etre modifie si besoin	
			currentCost = (int) pEvaluationSolvers_[sched][j]->evaluate(options_.evaluationParameters_);
			#ifdef COMPARE_EVALUATIONS
			pGreedyEvaluators[j]->solve();
			currentCostGreedy = pGreedyEvaluators[j]->solutionCost();
			#endif
		}

		// Display
		//
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Schedule no. " << sched << " evaluated over evaluation demand no. " << j << " (solution cost: " << currentCost << ")." << std::endl;

		// Insert the solution cost and solution
		//
		// If already in the costs -> add it to the set of schedules that found that cost
		if(schedulesFromObjectiveByEvaluationDemand_[j].find(currentCost) != schedulesFromObjectiveByEvaluationDemand_[j].end()){
			schedulesFromObjectiveByEvaluationDemand_[j].at(currentCost).insert(sched);
		}
		// Otherwise, add a new pair
		else{
			set<int> s; s.insert(sched);
			schedulesFromObjectiveByEvaluationDemand_[j].insert(pair<double, set<int> >( currentCost, s));
		}
		#ifdef COMPARE_EVALUATIONS
		if(schedulesFromObjectiveByEvaluationDemandGreedy_[j].find(currentCostGreedy) != schedulesFromObjectiveByEvaluationDemandGreedy_[j].end()){
			schedulesFromObjectiveByEvaluationDemandGreedy_[j].at(currentCostGreedy).insert(sched);
		}
		else {
			set<int> s; s.insert(sched);
			schedulesFromObjectiveByEvaluationDemandGreedy_[j].insert(pair<double, set<int> >( currentCostGreedy, s));
		}
		#endif
	}

	(*pLogStream_) << "# Evaluation of schedule no. " << sched << " done!" << std::endl;
}

// Recompute all scores after one schedule evaluation
void StochasticSolver::updateRankingsAndScores(){

	(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Starting the update of the scores and ranking." << std::endl;

	vector<double> theNewScores;
	Tools::initDoubleVector(&theNewScores, nSchedules_, 0);

	for(int j=0; j<options_.nEvaluationDemands_; j++){
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solution costs for demand no. " << j << endl;
		int localRank = 1;
		map<double, set<int> > localCosts = schedulesFromObjectiveByEvaluationDemand_[j];
		for(map<double, set<int> >::iterator it = localCosts.begin(); it != localCosts.end(); ++it){
			for(int sched : it->second){
				theNewScores[sched] += (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size());
				// If is infeasible -> double that amount to get more robust
				if((pEvaluationSolvers_[sched][j])->getStatus() == INFEASIBLE){
					theNewScores[sched] += (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size());
				}
				(*pLogStream_) << "#     | sched " << sched << " -> " << it->first << " (score += " << (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size()) << ")" << endl;
			}
			localRank += it->second.size();
		}
	}

	#ifdef COMPARE_EVALUATIONS
	vector<double> theNewScoresGreedy;
	Tools::initDoubleVector(&theNewScoresGreedy, nSchedules_, 0);
	for(int j=0; j<options_.nEvaluationDemands_; j++){
		(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Solution costs for demand no. " << j << endl;
		int localRank = 1;
		map<double, set<int> > localCosts = schedulesFromObjectiveByEvaluationDemandGreedy_[j];
		for(map<double, set<int> >::iterator it = localCosts.begin(); it != localCosts.end(); ++it){
			for(int sched : it->second){
				theNewScoresGreedy[sched] += (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size());
				// If is infeasible -> double that amount to get more robust
				if((pEvaluationSolvers_[sched][j])->getStatus() == INFEASIBLE){
					theNewScoresGreedy[sched] += (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size());
				}
				(*pLogStream_) << "#     | sched " << sched << " -> " << it->first << " (score += " << (double)localRank + ((double)(it->second.size() - 1)) / ((double) it->second.size()) << ")" << endl;
			}
			localRank += it->second.size();
		}
	}
	theScoresGreedy_ = theNewScoresGreedy;
	#endif

	theScores_ = theNewScores;
	(*pLogStream_) << "# [week=" << pScenario_->thisWeek() << "] Update of the scores and ranking done!" << std::endl;

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
