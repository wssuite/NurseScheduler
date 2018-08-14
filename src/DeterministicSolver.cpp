/*
* DeterministicSolver.cpp
*
*  Created on: April 29, 2015
*      Author: jeremy omer
*/

#include "DeterministicSolver.h"
//#include "Greedy.h"
#include "MasterProblem.h"
#include "InitializeSolver.h"

// #define COMPARE_EVALUATIONS

//-----------------------------------------------------------------------------
//
//  C l a s s   D e t e r m i n i s t i c S o l v e r
//
//  Solves the problem with deterministic demand
//  From a given problem (number of weeks, nurses, etc.), can compute a solution.
//
//-----------------------------------------------------------------------------

DeterministicSolver::DeterministicSolver(Scenario* pScenario,InputPaths inputPaths):
Solver(pScenario,pScenario->pWeekDemand(),pScenario->pWeekPreferences(), pScenario->pInitialState()),
pCompleteSolver_(0), pRollingSolver_(0), pLNSSolver_(0) {

	std::cout << "# New deterministic solver created!" << endl;

	// The nurses must be preprocessed to use their positions
	if (!isPreprocessedNurses_) this->preprocessTheNurses();

	// set the options of the deterministic solver
	// (the corresponding method needs to be changed manually for tests)
	//
	std::cout << "# Set the options" << std::endl;
	if (inputPaths.paramFile().empty()) {
		this->setOptionsToDefault(inputPaths);
	}
	else {
		this->readOptionsFromFile(inputPaths);
	}
	std::cout << std::endl;

	if (!options_.logfile_.empty()) {
		FILE * pFile;
		pFile = fopen (options_.logfile_.c_str(),"w");
		fprintf(pFile,"HERE COME THE LOGS OF THE MASTER PROBLEM\n\n");
		fclose(pFile);
	}
}

DeterministicSolver::~DeterministicSolver(){
	// delete also the reusable solvers
	if(pCompleteSolver_) delete pCompleteSolver_;
	if(pRollingSolver_) delete pRollingSolver_;
	// DBG if (pLNSSolver_) delete pLNSSolver_;
}


//----------------------------------------------------------------------------
//
// PARAMETERS FUNCTIONS
// Set the parameters of every solver
//
//----------------------------------------------------------------------------

// Initialize deterministic options with default values
//
void DeterministicSolver::setOptionsToDefault(InputPaths& inputPaths) {

	// initialize the log files when no path is specified
	string logDeterministic = inputPaths.logPath().empty() ? "" : inputPaths.logPath()+"LogDeterministic.txt";
	string logSolver = inputPaths.logPath().empty() ? "": inputPaths.logPath()+"LogSolver.txt";

	// global options are set to default values in the .h

	// random seed is initialized
	options_.randomSeed_ = inputPaths.randSeed();

	// parameters of the solution to optimality
	SolverParam completeParameters(options_.verbose_,TWO_DIVES,logDeterministic,logSolver);
	completeParameters_ = completeParameters;

	// parameters of the solution in the rolling horizon process
	SolverParam rollingParameters(options_.verbose_,TWO_DIVES,logDeterministic,logSolver);
	rollingParameters_ = rollingParameters;

	// parameters of the solution in the large neighborhood search
	SolverParam lnsParameters(options_.verbose_,TWO_DIVES,logDeterministic,logSolver);
	lnsParameters_ = lnsParameters;
}

// Read deterministic options from a file
//
void DeterministicSolver::readOptionsFromFile(InputPaths& inputPaths) {

	// initialize the log files when no path is specified
	string logDeterministic = inputPaths.logPath().empty() ? "" : inputPaths.logPath()+"LogDeterministic.txt";
	string logSolver = inputPaths.logPath().empty() ? "": inputPaths.logPath()+"LogSolver.txt";

	// random seed is initialized
	options_.randomSeed_ = inputPaths.randSeed();

	// open the file
	//
	std::fstream file;
	std::string fileName = inputPaths.paramFile();
	std::cout << "Reading " << fileName << std::endl;
	file.open(fileName.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read the file " << fileName << std::endl;
		throw Tools::myException("readOptionsFromFile: The input file was not opened properly!",__LINE__);
	}

	// go through all the lines of the parameter file
	std::string title;
	SolverParam param;
	OptimalityLevel completeOptimalityLevel=UNTIL_FEASIBLE;
	OptimalityLevel lnsOptimalityLevel=UNTIL_FEASIBLE;
	OptimalityLevel rollingOptimalityLevel=UNTIL_FEASIBLE;
	while(file.good()){
		Tools::readUntilChar(&file, '=', &title);

		// Read the name of the scenario
		//
		if(Tools::strEndsWith(title, "divideIntoConnexPositions")){
			file >> options_.divideIntoConnexPositions_;
		}
		else if (Tools::strEndsWith(title, "withRollingHorizon")) {
			file >> options_.withRollingHorizon_;
		}
		else if (Tools::strEndsWith(title, "rollingSamplePeriod")) {
			file >> options_.rollingSamplePeriod_;
		}
		else if (Tools::strEndsWith(title, "rollingControlHorizon")) {
			file >> options_.rollingControlHorizon_;
		}
		else if (Tools::strEndsWith(title, "rollingPredictionHorizon")) {
			file >> options_.rollingPredictionHorizon_;
		}
		else if (Tools::strEndsWith(title, "withLNS")) {
			file >> options_.withLNS_;
		}
		else if (Tools::strEndsWith(title, "lnsMaxItWithoutImprovement")) {
			file >> options_.lnsMaxItWithoutImprovement_;
		}
		else if (Tools::strEndsWith(title, "lnsNursesRandomDestroy")) {
			file >> options_.lnsNursesRandomDestroy_;
		}
		else if (Tools::strEndsWith(title, "lnsNursesPositionDestroy")) {
			file >> options_.lnsNursesPositionDestroy_;
		}
		else if (Tools::strEndsWith(title, "lnsNursesContractDestroy")) {
			file >> options_.lnsNursesContractDestroy_;
		}
		else if (Tools::strEndsWith(title, "lnsDestroyOverTwoWeeks")) {
			file >> options_.lnsDestroyOverTwoWeeks_;
		}
		else if (Tools::strEndsWith(title, "lnsDestroyOverFourWeeks")) {
			file >> options_.lnsDestroyOverFourWeeks_;
		}
		else if (Tools::strEndsWith(title, "lnsDestroyOverAllWeeks")) {
			file >> options_.lnsDestroyOverAllWeeks_;
		}
		else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverTwoWeeks")) {
			file >> options_.lnsNbNursesDestroyOverTwoWeeks_;
		}
		else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverFourWeeks")) {
			file >> options_.lnsNbNursesDestroyOverFourWeeks_;
		}
		else if (Tools::strEndsWith(title, "lnsNbNursesDestroyOverAllWeeks")) {
			file >> options_.lnsNbNursesDestroyOverAllWeeks_;
		}
		else if (Tools::strEndsWith(title, "solutionAlgorithm")) {
			std::string algoName;
			file >> algoName;
			options_.solutionAlgorithm_ = AlgorithmsByName[algoName];
		}
		else if (Tools::strEndsWith(title, "solverType")) {
			std::string solverName;
			file >> solverName;
			options_.MySolverType_ = MySolverTypesByName[solverName];
			// DBG
			std::cout << "LP solver :" << solverName << std::endl;
		}
		else if (Tools::strEndsWith(title, "verbose")) {
			file >> options_.verbose_;
		}
		// Below, read the parameters that refer to branch and price solution
		// They are not options of the deterministic solver, but they need to be
		// set at this stage
		else if (Tools::strEndsWith(title, "isStabilization")) {
			file >> param.isStabilization_;
		}
		else if (Tools::strEndsWith(title, "isStabUpdateCost")) {
			file >> param.isStabUpdateCost_;
		}
		else if (Tools::strEndsWith(title, "isStabUpdateBounds")) {
			file >> param.isStabUpdateBounds_;
		}
		else if (Tools::strEndsWith(title, "branchColumnDisjoint")) {
			file >> param.branchColumnDisjoint_;
		}
		else if (Tools::strEndsWith(title, "branchColumnUntilValue")) {
			file >> param.branchColumnUntilValue_;
		}
		else if (Tools::strEndsWith(title, "stopAfterXDegenerateIt")) {
			file >> param.stopAfterXDegenerateIt_;
		}
		else if (Tools::strEndsWith(title, "heuristicMinIntegerPercent")) {
			file >> param.heuristicMinIntegerPercent_;
		}
		else if (Tools::strEndsWith(title, "performHeuristicAfterXNode")) {
			file >> param.performHeuristicAfterXNode_;
		}
		else if (Tools::strEndsWith(title, "rollingOptimalityLevel")) {
			std::string strOpt;
			file >> strOpt;
			rollingOptimalityLevel=stringToOptimalityLevel[strOpt];
		}
		else if (Tools::strEndsWith(title, "lnsOptimalityLevel")) {
			std::string strOpt;
			file >> strOpt;
			lnsOptimalityLevel=stringToOptimalityLevel[strOpt];
		}
		else if (Tools::strEndsWith(title, "completeOptimalityLevel")) {
			std::string strOpt;
			file >> strOpt;
			completeOptimalityLevel=stringToOptimalityLevel[strOpt];		}
		// Subproblem options
		//
		else if (Tools::strEndsWith(title, "spDefaultStrategy")) {
			file >> param.sp_default_strategy_;
		}
		else if (Tools::strEndsWith(title, "spSecondChanceStrategy")) {
			file >> param.sp_secondchance_strategy_;
		}
		else if (Tools::strEndsWith(title, "spNbRotationsPerNurse")) {
			file >> param.sp_nbrotationspernurse_;
		}
		else if (Tools::strEndsWith(title, "spNbNursesToPrice")) {
			file >> param.sp_nbnursestoprice_;
		}
		else if (Tools::strEndsWith(title, "spWithSecondChance")) {
			file >> param.sp_withsecondchance_;
		}
		else if (Tools::strEndsWith(title, "spMaxReducedCostBound")) {
			file >> param.sp_max_reduced_cost_bound_;
		}
	}
	options_.totalTimeLimitSeconds_ = inputPaths.timeOut();
	param.maxSolvingTimeSeconds_ = options_.totalTimeLimitSeconds_;

	// parameters of the solution to optimality
	completeParameters_ = param;
	completeParameters_.initialize(options_.verbose_,completeOptimalityLevel);

	// parameters of the solution in the rolling horizon process
	rollingParameters_ = param;
	rollingParameters_.initialize(options_.verbose_,rollingOptimalityLevel);
	rollingParameters_.performHeuristicAfterXNode_ = -1;

	// parameters of the solution in the large neighborhood search
	lnsParameters_ = param;
	lnsParameters_.initialize(options_.verbose_,lnsOptimalityLevel);
	lnsParameters_.performHeuristicAfterXNode_ = -1;

}

//----------------------------------------------------------------------------
//
// STATISTICS OF THE OVERALL SOLUTION PROCESS
//
//----------------------------------------------------------------------------

void DeterministicSolver::updateInitialStats(MasterProblem* pMaster) {

	BcpModeler* pModel = dynamic_cast<BcpModeler*>(pMaster->getModel());

	stats_.bestUBInitial_= pMaster->computeSolutionCost();
	stats_.bestUB_= pMaster->computeSolutionCost();
	stats_.rootLB_ = pModel->getRootLB();
	stats_.bestLB_ = std::min(stats_.bestUB_, 5.0*ceil( (pModel->get_best_lb()-EPSILON)/5.0));
	stats_.timeInitialSol_= pMaster->getTimerTotal()->dSinceStart();
	pMaster->getTimerTotal()->reset();

	// Details on Branch and price related runtimes
	//
	stats_.timeGenColRoot_= pModel->gettimeFirstRoot();
	stats_.timeGenColMaster_=pModel->getTimeStats().time_lp_solving;
	stats_.timeGenSubProblems_=pModel->getTimeStats().time_var_generation;

	// details on Branch and price iterations
	//
	stats_.itGenColInitial_=pModel->getNbLpIterations();;
	stats_.nodesBBInitial_=pModel->getNbNodes();
	if (!stats_.nodesBBInitial_) stats_.timeGenColRoot_ = stats_.timeInitialSol_;
}
void DeterministicSolver::updateImproveStats(MasterProblem* pMaster) {

	BcpModeler* pModel = dynamic_cast<BcpModeler*>(pMaster->getModel());

	stats_.bestUB_=pMaster->computeSolutionCost();
	stats_.timeImproveSol_=pMaster->getTimerTotal()->dSinceStart();

	// Details on Branch and price related runtimes
	//
	stats_.timeGenColMaster_+=pModel->getTimeStats().time_lp_solving;
	stats_.timeGenSubProblems_+=pModel->getTimeStats().time_var_generation;

	// details on Branch and price iterations
	//
	stats_.itGenColImprove_=pModel->getNbLpIterations();
	stats_.nodesBBImprove_=pModel->getNbNodes();

	// number of problems solved in each phase
	//
	if (options_.withLNS_) {
		stats_.itImproveSol_=0;
	}
}


//----------------------------------------------------------------------------
//
// SOLVE FUNCTIONS
// The one is general for the whole process
//
//----------------------------------------------------------------------------

// Main function
double DeterministicSolver::solve(vector<Roster> initialSolution){

	objValue_ = 0.0;

	// set the random seed to the same fixed value for every solution of a new
	// scenario
	// this is absolutely necessary for reproductibility of the results
	//
	Tools::initializeRandomGenerator(this->options_.randomSeed_);
	srand(this->options_.randomSeed_);

	// DBG
	std::cout << "Next random : " << Tools::randomInt(0, RAND_MAX) << std::endl;
	std::cout << "Next random : " << Tools::randomInt(0, RAND_MAX) << std::endl;

	// Always solve small problems to optimality
	// This can actually save time
	//
	if ( (pScenario_->nbDays() <= 28 && pScenario_->nbNurses() <= 8)
		|| (pScenario_->nbDays() <= 56 && pScenario_->nbNurses() <= 5) ) {
		completeParameters_.setOptimalityLevel(OPTIMALITY);
		objValue_ = this->solveCompleteHorizon();
		if (MasterProblem* pMaster = dynamic_cast<MasterProblem*> (pCompleteSolver_)) {
			this->updateInitialStats(pMaster);
		}
		return objValue_;
	}

	// First divide the scenarion into connex positions if requested
	//
	if (options_.divideIntoConnexPositions_) {
		options_.divideIntoConnexPositions_ = false;
		objValue_ = this->solveByConnexPositions();
	}
	// If the the scenarion is divided into connex positions, the solution
	// of the subproblems goes in the "else" below
	//
	else {
		// Find a good feasible solution using a rolling horizon planning or
		/// solving directly the complete horizon
		//
		if (options_.withRollingHorizon_) {
			objValue_ = this->solveWithRollingHorizon();
			if (MasterProblem* pMaster = dynamic_cast<MasterProblem*> (pRollingSolver_)) {
				this->updateInitialStats(pMaster);
			}
		}
		else {
			objValue_ = this->solveCompleteHorizon();
			if (MasterProblem* pMaster = dynamic_cast<MasterProblem*> (pCompleteSolver_)) {
				this->updateInitialStats(pMaster);
			}
			// do not bother improving the solution if it is already optimal
			if (status_ == OPTIMAL) {
				return objValue_;
			}
		}

		// Improve the solution with an LNS
		//
		if (options_.withLNS_) {
			objValue_ = this->solveWithLNS();

			if (MasterProblem* pMaster = dynamic_cast<MasterProblem*> (pLNSSolver_)) {
				this->updateImproveStats(pMaster);
			}
		}
	}
	return objValue_;
}


//------------------------------------------------------------------------
//
// Solve the complete horizon using the input algorithm
//
//------------------------------------------------------------------------

double DeterministicSolver::solveCompleteHorizon() {

	// Initialize solver and solve
	//
	pCompleteSolver_ = setSolverWithInputAlgorithm(pDemand_);
	pCompleteSolver_->solve(completeParameters_);

	return this->treatResults(pCompleteSolver_);
}


//----------------------------------------------------------------------------
// After the end of a solution process: retrieve status, solution, etc.
//----------------------------------------------------------------------------

double DeterministicSolver::treatResults(Solver* pSolver) {

	// Change status back to feasible a feasible solution was found
	status_ = pSolver->getStatus();
	double timeSinceStart = pTimerTotal_->dSinceStart();
	std::cout << "Time spent until then: " << timeSinceStart << " s ";
	std::cout << "(time limit is "<< options_.totalTimeLimitSeconds_ << " s)" << std::endl;
	if (status_ == TIME_LIMIT) {
		std::cout << "Stop solution process: time limit is reached" << std::endl;
		if (MasterProblem* pMaster = dynamic_cast<MasterProblem*>(pSolver)) {
			if (pMaster->getModel()->nbSolutions() >= 1) {
				status_ = FEASIBLE;
			}
		}
	}

	// Print the solution if required
	if (completeParameters_.printIntermediarySol_) {
		pSolver->printCurrentSol();
	}

	// Update nurses' states when solution is feasible or optimal
	if (status_ == FEASIBLE || status_ == OPTIMAL ) {
		solution_ = pSolver->getSolution();
		for(int n=0; n<pScenario_->nbNurses_; ++n){
			theLiveNurses_[n]->roster_ = solution_[n];
			theLiveNurses_[n]->buildStates();
		}
	}

	objValue_ = this->computeSolutionCost();

	return objValue_;
}


//------------------------------------------------------------------------
//
// Solve the problem using a decomposition of the set nurses by connex
// components of the graph of positions
//
//------------------------------------------------------------------------

double DeterministicSolver::solveByConnexPositions() {

	// DIVIDE THE SCENARIO INTO CONNEX COMPONENTS AND PRINT THE RESULT
	vector<Scenario*> scenariosPerComponent;
	scenariosPerComponent = divideScenarioIntoConnexPositions(pScenario_);

	// SOLVE THE PROBLEM COMPONENT-WISE
	vector<DeterministicSolver*> solverPerComponent;
	for (Scenario* pScenario: scenariosPerComponent) {
		std::cout << "COMPONENT-WISE SCENARIO" << std::endl;
		std::cout << pScenario->toString() << std::endl;

		// SET THE SOLVER AND SOLVE THE SUBPROBLEM
		InputPaths inputPaths;
		solverPerComponent.push_back(new DeterministicSolver(pScenario,inputPaths));
		solverPerComponent.back()->copyParameters(this);

		// set allowed time proportionnally to the number of nurses in each
		// component
		double allowedTime = options_.totalTimeLimitSeconds_*(double)pScenario->nbNurses()/(double)pScenario_->nbNurses();
		// if solving the last component, leave it all the time left
		if (solverPerComponent.size() == scenariosPerComponent.size()) {
			allowedTime = std::max(allowedTime,options_.totalTimeLimitSeconds_-pTimerTotal_->dSinceStart());
		}
		solverPerComponent.back()->setTotalTimeLimit(allowedTime);

		// solve the component
		solverPerComponent.back()->solve();

		// STORE THE SOLUTION
		// Be particularly cautious that the nurse indices are not the same in the
		// initial scenario and in the solvers per component
		for (int n=0; n<pScenario->nbNurses_; ++n) {
			int idNurse = pScenario->theNurses_[n].id_;
			theLiveNurses_[idNurse]->roster_ = solverPerComponent.back()->getSolution()[n];
		}

		// Consolidate the global state of the solver
		stats_.add(solverPerComponent.back()->getGlobalStat());
		Status lastStatus = solverPerComponent.back()->getStatus();

		// in several cases, the new status is the status of the last solver solved
		if (status_ == UNSOLVED || status_ == OPTIMAL ||
			lastStatus == INFEASIBLE || lastStatus == TIME_LIMIT || lastStatus == UNSOLVED) {
			status_ = lastStatus;
		}
		// in all other cases status is unchanged
		else {	}

		// break if the status is not that of a normally finished solution process
		if (status_ == UNSOLVED || status_ == TIME_LIMIT || status_ == INFEASIBLE) {
			std::cout << "Solution process did not terminate normally" << std::endl;
			return -1;
		}

		// the solver of the component can be deleted at this stage
		delete solverPerComponent.back();
		solverPerComponent.back()=0;
	}

	// update nurses' states
	for(int n=0; n<pScenario_->nbNurses_; ++n){
		solution_.push_back(theLiveNurses_[n]->roster_);
		theLiveNurses_[n]->buildStates();
	}

	//  release memory
	for (Scenario* pScenario: scenariosPerComponent) {
		if (pScenario->pWeekDemand()) delete pScenario->pWeekDemand();
		pScenario->setWeekDemand(0);
	}

	return computeSolutionCost();
}

//------------------------------------------------------------------------
// Solve the problem with a receeding horizon algorithm
// The sample period is the number of days for which column variables
// must be integer
//------------------------------------------------------------------------

double DeterministicSolver::solveWithRollingHorizon() {

	std::cout << "SOLVE WITH ROLLING HORIZON" << std::endl;

	int samplePeriod = options_.rollingSamplePeriod_;
	int controlPeriod = options_.rollingControlHorizon_;

	// Initialize the solver that will handle the iterative solution of the
	// receeding horizon
	//
	pRollingSolver_ = setSolverWithInputAlgorithm(pDemand_);

	// Solve the instance iteratively with a rolling horizon
	//
	int firstDay = 0; //first day of the current horizon
	while (firstDay < pDemand_->nbDays_) {
		std::cout << "FIRST DAY = " << firstDay <<  std::endl << std::endl;

		// last days of the sample horizon and of the control horizon
		//
		int lastDaySample = std::min(firstDay+samplePeriod-1,pDemand_->nbDays_-1);
		int lastDayControl = std::min(firstDay+controlPeriod-1,pDemand_->nbDays_-1);

		// Relax the integrality constraints on the variables outside the horizon control
		// (and inside the prediction horizon, but at this stage the prediction horizon
		//   includes the whole horizon)
		//
		vector<bool> isRelaxDay(pDemand_->nbDays_,false);
		for (int day=lastDayControl+1; day < pDemand_->nbDays_; day++) isRelaxDay[day] = true;
		pRollingSolver_->relaxDays(isRelaxDay);

		// Solve the problem with a method that allows for a warm start
		//
		this->rollingSetOptimalityLevel(firstDay);
		pRollingSolver_->rollingSolve(rollingParameters_,firstDay);

		if (rollingParameters_.printIntermediarySol_) {
			pRollingSolver_->printCurrentSol();
		}

		// fix the days of the sample horizon to the values of the solution
		//
		vector<bool> isFixDay(pDemand_->nbDays_,false);
		for (int day=0; day<=lastDaySample; day++) isFixDay[day] = true;
		pRollingSolver_->fixDays(isFixDay);

		// update the first and last day of the sample period
		//
		firstDay = firstDay+samplePeriod;
		lastDayControl = std::min(firstDay+controlPeriod-1,pDemand_->nbDays_-1);

		// Set back the integrality constraints for next sampling period
		vector<bool> isUnrelaxDay(pDemand_->nbDays_,false);
		for (int day=0; day <=lastDayControl; day++) isUnrelaxDay[day] = true;
		pRollingSolver_->unrelaxDays(isUnrelaxDay);

		// stop lns if runtime is exceeded
		//
		double timeSinceStart = pTimerTotal_->dSinceStart();
		std::cout << "Time spent until then: " << timeSinceStart << " s" << std::endl;
		if (timeSinceStart > options_.totalTimeLimitSeconds_) {
			std::cout << "Stop the rolling horizon: time limit is reached!" << std::endl;
			break;
		}
	}
	vector<bool>isUnfixDay(pDemand_->nbDays_,true);
	pRollingSolver_->unfixDays(isUnfixDay);
	pRollingSolver_->storeSolution();
	pRollingSolver_->costsConstrainstsToString();

	std::cout << "END OF ROLLING HORIZON" << std::endl << std::endl;

	return treatResults(pRollingSolver_);
}

// Set the optimality level of the rolling horizon solver
// This function needs to be called before each new solution, and the behavior
// depends on the first day of the horizon
//
void DeterministicSolver::rollingSetOptimalityLevel(int firstDay) {
	rollingParameters_.setOptimalityLevel(TWO_DIVES);

 	if (pDemand_->nbDays_-firstDay <= 21) {
		rollingParameters_.setOptimalityLevel(REPEATED_DIVES);
	}
}


//------------------------------------------------------------------------
//
// Solve the problem with a large neighborhood search
// Iteratively fix the schedule of complete weeks or of a set of nurses and
// solve the rest to improve a initial feasible solution
//
//------------------------------------------------------------------------


// Perform the LNS
//
double DeterministicSolver::solveWithLNS() {

	std::cout << "SOLVE WITH LNS" << std::endl << std::endl;

	// Initialize data structures for LNS
	//
	this->initializeLNS();
	std::vector<double> nursesSelectionWeights(nursesSelectionOperators_.size(),3.0);
	std::vector<double> daysSelectionWeights(daysSelectionOperators_.size(),3.0);
	std::vector<double> repairWeights(repairOperators_.size(),0.1);

	// Initialize the solver that will handle the repair problems
	//
	if (options_.withRollingHorizon_) {
		pLNSSolver_ = pRollingSolver_;
	}
	else {
		pLNSSolver_ = pCompleteSolver_;
	}

	// set the computational time left for the lns after the initialization
	double timeSinceStart = pTimerTotal_->dSinceStart();
	lnsParameters_.maxSolvingTimeSeconds_ = options_.totalTimeLimitSeconds_ - timeSinceStart;

	// pLNSSolver_ = setSolverWithInputAlgorithm(pDemand_);
	// pLNSSolver_->initialize(options_.lnsParameters_,this->solution_);

	// Perform destroy/repair iterations until a given number of iterations
	// without improvement is reached
	//
	int nbItWithoutImprovement=0;
	double bestObjVal=this->computeSolutionCost();
	while (true) { //nbItWithoutImprovement < options_.lnsMaxItWithoutImprovement_) {

		// draw the next destroy operators randomly according to the weights
		int nurseIndex = Tools::drawRandomWithWeights(nursesSelectionWeights);
		int dayIndex = Tools::drawRandomWithWeights(daysSelectionWeights);
		int repairIndex = Tools::drawRandomWithWeights(repairWeights);
		NursesSelectionOperator nurseOperator = nursesSelectionOperators_[nurseIndex];
		DaysSelectionOperator dayOperator = daysSelectionOperators_[dayIndex];

		// apply the destroy operator
		this->adaptiveDestroy(nurseOperator, dayOperator);

		// run the repair operator
		//
		double currentObjVal = pLNSSolver_->LNSSolve(lnsParameters_);

		// stop lns if runtime is exceeded
		//
		timeSinceStart = pTimerTotal_->dSinceStart();
		std::cout << "Time spent until then: " << timeSinceStart << " s" ;
		std::cout << "(time limit is "<< options_.totalTimeLimitSeconds_ << " s)" << std::endl;
		if (pLNSSolver_->getStatus()==TIME_LIMIT && timeSinceStart <= options_.totalTimeLimitSeconds_ - 5.0) {
			Tools::throwError("Error with the timers in LNS!");
		}
		if (pLNSSolver_->getStatus()==TIME_LIMIT || timeSinceStart > options_.totalTimeLimitSeconds_) {
			std::cout << "Stop the lns: time limit is reached" << std::endl;
			break;
		}

		// store the solution
		//
		solution_ = pLNSSolver_->getSolution();
		status_ = pLNSSolver_->getStatus();
		if (lnsParameters_.printIntermediarySol_) {
			pLNSSolver_->printCurrentSol();
		}

		// update the weight of adaptive lns
		//
		if (currentObjVal < bestObjVal-EPSILON) {
			// update stats
			stats_.lnsImprovementValueTotal_+=bestObjVal-currentObjVal;
			stats_.lnsNbIterationsWithImprovement_++;

			// update weights of the destroy operators
			bestObjVal = currentObjVal;
			nbItWithoutImprovement = 0;
			lnsParameters_.setOptimalityLevel(TWO_DIVES);
			nursesSelectionWeights[nurseIndex] = nursesSelectionWeights[nurseIndex]+1.0;
			daysSelectionWeights[nurseIndex] = daysSelectionWeights[nurseIndex]+1.0;

			// update weights of the repair operators
			double timeIteration = pTimerTotal_->dSinceInit()-timeSinceStart;
			repairWeights[repairIndex] += 10.0/timeIteration;


			// update the counters of iterations with improvement
			stats_.nbImprovementsWithNursesSelection_[nurseIndex]++;
			stats_.nbImprovementsWithDaysSelection_[dayIndex]++;
			stats_.nbImprovementsWithRepair_[repairIndex]++;
		}
		else {
			nbItWithoutImprovement++;

			if (nbItWithoutImprovement > std::min(30,options_.lnsMaxItWithoutImprovement_/2) ) {
				lnsParameters_.setOptimalityLevel(OPTIMALITY);
			}
			else if (nbItWithoutImprovement > std::min(10,options_.lnsMaxItWithoutImprovement_/4) ) {
				lnsParameters_.setOptimalityLevel(REPEATED_DIVES);
			}
		}

		std::cout << "**********************************************" << std::endl
		          << "LNS iteration: " << stats_.lnsNbIterations_
		          << "\t" << "Best solution: " << bestObjVal << std::endl
							<< "**********************************************" << std::endl;

		// unfix every nurse and/or days for next iteration
		//
		std::vector<bool> isUnfixNurse(pScenario_->nbNurses_,true);
		std::vector<bool> isUnfixDay(getNbDays(),true);
		pLNSSolver_->unfixNurses(isUnfixNurse);
		pLNSSolver_->unfixDays(isUnfixDay);

		stats_.lnsNbIterations_++;
	}

	std::cout << "END OF LNS" << std::endl << std::endl;

	return treatResults(pLNSSolver_);
}


// Prepare data structures for LNS
//
void DeterministicSolver::initializeLNS() {
	// initialize the set of destroy operators of the solver
	if (options_.lnsNursesRandomDestroy_) {
		nursesSelectionOperators_.push_back(NURSES_RANDOM);
	}
	if (options_.lnsNursesPositionDestroy_) {
		nursesSelectionOperators_.push_back(NURSES_POSITION);
		this->organizeTheLiveNursesByPosition();
	}
	if (options_.lnsNursesContractDestroy_) {
		nursesSelectionOperators_.push_back(NURSES_CONTRACT);
		this->organizeTheLiveNursesByContract();
	}
	if (options_.lnsDestroyOverTwoWeeks_) {
		daysSelectionOperators_.push_back(TWO_WEEKS);
	}
	if ( (options_.lnsDestroyOverFourWeeks_ && getNbDays() > 28) ||
		(options_.lnsDestroyOverAllWeeks_ && getNbDays() <= 28) )  {
		daysSelectionOperators_.push_back(FOUR_WEEKS);
	}
	if (options_.lnsDestroyOverAllWeeks_ && getNbDays() > 28) {
		daysSelectionOperators_.push_back(ALL_WEEKS);
	}

	// initialize the set of repair operators of the lns
	repairOperators_.push_back(REPAIR_TWO_DIVES);
	repairOperators_.push_back(REPAIR_REPEATED_DIVES);
	repairOperators_.push_back(REPAIR_OPTIMALITY);

	// initialize the counters of improvements
	stats_.nbImprovementsWithRepair_.insert(stats_.nbImprovementsWithRepair_.begin(),repairOperators_.size(),0);
	stats_.nbImprovementsWithNursesSelection_.insert(stats_.nbImprovementsWithNursesSelection_.begin(),nursesSelectionOperators_.size(),0);
	stats_.nbImprovementsWithDaysSelection_.insert(stats_.nbImprovementsWithDaysSelection_.begin(),daysSelectionOperators_.size(),0);
}


// Application of the destroy operator
//
void DeterministicSolver::adaptiveDestroy(NursesSelectionOperator nurseOp, DaysSelectionOperator dayOp) {
	// apply the destroy operator
	std::vector<bool> isFixNurse(pScenario_->nbNurses_,true);
	std::vector<bool> isFixDay(pScenario_->nbDays(),true);
	std::vector<int> randIndVector;

	// FIRST SET THE NUMBER OF NURSES AND DAYS THAT MUST BE FIXED
	int nbNursesDestroy = 0;
	int nbDaysDestroy = 0;
	switch (dayOp) {
		case TWO_WEEKS:
			nbNursesDestroy=options_.lnsNbNursesDestroyOverTwoWeeks_;
			nbDaysDestroy=14;
			break;
		case FOUR_WEEKS:
			nbNursesDestroy=options_.lnsNbNursesDestroyOverFourWeeks_;
			nbDaysDestroy=28;
			break;
		case ALL_WEEKS:
			nbNursesDestroy=options_.lnsNbNursesDestroyOverAllWeeks_;
			nbDaysDestroy=getNbDays();
			break;
	}

	// SECOND GENERATE THE NURSES WHOSE PLANNING WILL BE FIXED
	//
	switch (nurseOp) {
		case NURSES_RANDOM: {
			randIndVector = Tools::drawRandomIndices(nbNursesDestroy,0,pScenario_->nbNurses_-1);
			for (int ind:randIndVector) {
				isFixNurse[ind]=false;
			}
			break;
		}
		case NURSES_POSITION: {
			int randPos = Tools::drawRandomWithWeights(positionWeights_);
			int nbNursesWithPos = theLiveNursesByPosition_[randPos].size();
			randIndVector = Tools::drawRandomIndices(nbNursesDestroy,0,nbNursesWithPos-1);
			for (int ind:randIndVector) {
				isFixNurse[theLiveNursesByPosition_[randPos][ind]->id_]=false;
			}
			break;
		}
		case NURSES_CONTRACT: {
			int randContract = Tools::drawRandomWithWeights(contractWeights_);
			int nbNursesWithContract = theLiveNursesByContract_[randContract].size();
			randIndVector = Tools::drawRandomIndices(nbNursesDestroy,0,nbNursesWithContract-1);
			for (int ind:randIndVector) {
				isFixNurse[theLiveNursesByContract_[randContract][ind]->id_]=false;
			}
			break;
		}
	}
	// Fix the nurses that are not destroyed
	pLNSSolver_->fixNurses(isFixNurse);

	// GENERATE THE DAYS THAT WILL BE DESTROYED AND FIX THE OTHERS
	// fix no day if the number of days in the scenario is small
	if (nbDaysDestroy < this->getNbDays()) {
		// draw the first day of the relaxed interval
		std::vector<double> weightDays(pScenario_->nbDays()-nbDaysDestroy-1,1.0);
		weightDays[0] = 7;
		for (int i=1; i< std::max(getNbDays()-nbDaysDestroy-1,6); i++) {
			weightDays[i] = 0.1;
		}
		int firstDay = Tools::drawRandomWithWeights(weightDays); // Tools::randomInt(0, getNbDays()-nbDaysDestroy-1);
		for (int day=0; day <nbDaysDestroy; day++) {
			isFixDay[firstDay+day] = false;
		}
		pLNSSolver_->fixDays(isFixDay);
	}

	// DBG
	std::cout << "REPAIR NURSES: ";
	for (int i=0; i < isFixNurse.size(); i++) {
		if (!isFixNurse[i]) std::cout << i << "\t";
	}
	std::cout << std::endl;
	std::cout << "REPAIR DAYS: ";
	for (int i=0; i < isFixDay.size(); i++) {
		if (!isFixDay[i]) std::cout << i << "\t";
	}
	std::cout << std::endl;

}


// Initialize the organized vectors of live nurses
//
void DeterministicSolver::organizeTheLiveNursesByPosition() {
	nbPositions_=0;
	std::vector<LiveNurse*> copyTheLiveNurses = theLiveNurses_;

	// Transfer the nurses with same positions from the copy vector of live nurse
	// to the organized vector of live nurses
	while (!copyTheLiveNurses.empty()) {
		nbPositions_++;
		std::vector<LiveNurse*> theLiveNursesAtPosition;

		// initialize with the first nurse
		theLiveNursesAtPosition.push_back(copyTheLiveNurses[0]);
		const Position* thisPosition = copyTheLiveNurses[0]->pPosition_;
		copyTheLiveNurses.erase(copyTheLiveNurses.begin());

		// search for the nurses with same position in the remaining nurses
		int nbNursesLeft=copyTheLiveNurses.size();
		for (int n=nbNursesLeft-1; n >= 0; n--) {
			if (copyTheLiveNurses[n]->pPosition_->id_ == thisPosition->id_) {
				theLiveNursesAtPosition.push_back(copyTheLiveNurses[n]);
				copyTheLiveNurses.erase(copyTheLiveNurses.begin()+n);
			}
		}
		theLiveNursesByPosition_.push_back(theLiveNursesAtPosition);
	}

	// Initialize the position weights according to the number of nurses with
	// each position
	//
	for (int p = 0; p < nbPositions_; p++) {
		positionWeights_.push_back(theLiveNursesByPosition_[p].size());
	}
}

void DeterministicSolver::organizeTheLiveNursesByContract() {
	nbContracts_=0;
	std::vector<LiveNurse*> copyTheLiveNurses = theLiveNurses_;

	// Transfer the nurses with same contract from the copy vector of live nurse
	// to the organized vector of live nurses
	while (!copyTheLiveNurses.empty()) {
		nbContracts_++;
		std::vector<LiveNurse*> theLiveNursesWithContract;

		// initialize with the first nurse
		theLiveNursesWithContract.push_back(copyTheLiveNurses[0]);
		const Contract* thisContract = copyTheLiveNurses[0]->pContract_;
		copyTheLiveNurses.erase(copyTheLiveNurses.begin());

		// search for the nurses with same contract in the remaining nurses
		int nbNursesLeft=copyTheLiveNurses.size();
		for (int n=nbNursesLeft-1; n >= 0; n--) {
			if (copyTheLiveNurses[n]->pContract_->id_ ==thisContract->id_) {
				theLiveNursesWithContract.push_back(copyTheLiveNurses[n]);
				copyTheLiveNurses.erase(copyTheLiveNurses.begin()+n);
			}
		}
		theLiveNursesByContract_.push_back(theLiveNursesWithContract);
	}

	// Initialize the contract weights according to the number of nurses with
	// each contract
	//
	for (int c = 0; c < nbContracts_; c++) {
		contractWeights_.push_back(theLiveNursesByContract_[c].size());
	}
}


//----------------------------------------------------------------------------
//
// Construct solvers that will be used to really make the job
//
//----------------------------------------------------------------------------

// Return a solver with the algorithm specified for resolution
Solver * DeterministicSolver::setSolverWithInputAlgorithm(Demand* pDemand) {
	Solver* pSolver=NULL;
	switch(options_.solutionAlgorithm_){
		//case GREEDY:
		//pSolver = new Greedy(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState());
		//break;
		case GENCOL:
		pSolver = new MasterProblem(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState(), options_.MySolverType_);
		break;
		default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}
	return pSolver;
}

// Return a solver with the input algorithm
Solver* DeterministicSolver::setSubSolverWithInputAlgorithm(Demand* pDemand, Algorithm algorithm) {

	Solver* pSolver=NULL;
	switch(algorithm){
		//case GREEDY:
		//pSolver = new Greedy(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState());
		//break;
		case GENCOL:
		pSolver = new MasterProblem(pScenario_, pDemand, pScenario_->pWeekPreferences(), pScenario_->pInitialState(), options_.MySolverType_);
		break;
		default:
		Tools::throwError("The algorithm is not handled yet");
		break;
	}
	return pSolver;
}
