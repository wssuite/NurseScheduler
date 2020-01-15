/*
 * MasterProblem.cpp
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#include "MasterProblem.h"
#include "BcpModeler.h"
#include "RotationPricer.h"
#include "TreeManager.h"
#include "OsiClpSolverInterface.hpp"

#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"
#endif

#ifdef USE_CBC
#include "CbcModeler.h"
#endif

#ifdef USE_SCIP
#include "ScipModeler.h"
#endif

/* namespace usage */
using namespace std;



//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

void Rotation::computeCost(Scenario* pScenario, Preferences* pPreferences, const vector<LiveNurse*>& liveNurses, int horizon){
	//check if pNurse points to a nurse
	if(nurseId_ == -1)
		Tools::throwError("LiveNurse = NULL");

	LiveNurse* pNurse = liveNurses[nurseId_];

	/************************************************
	 * Compute all the costs of a rotation:
	 ************************************************/
	//   double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_ ;

	//if first day of the planning, check on the past, otherwise 0 (rest)
	int lastShift = (firstDay_==0) ? pNurse->pStateIni_->shift_ : 0;
	//nbConsShift = number of consecutive shift
	//if first day of the planning, check on the past, otherwise 0
	int nbConsShifts = (firstDay_==0) ? pNurse->pStateIni_->consShifts_ : 0;
	//consShiftCost = cost of be outside of the interval [min,max] of the consecutives shifts
	consShiftsCost_ = 0;

	//nbConsWorked = number of consecutive worked days
	//if first day of the planning, check on the past , otherwise 0
	int nbConsDaysWorked = (firstDay_==0) ? pNurse->pStateIni_->consDaysWorked_ : 0;
	//consWorkedCost = cost of be outside of the interval [min,max] of the consecutives worked days
	consDaysWorkedCost_ = 0 ;

	//cost of not doing the whole weekend
	completeWeekendCost_ = 0;

	//preferencesCost = cost of not respecting preferences
	preferenceCost_ = 0;

	//initial resting cost
	initRestCost_ =  0;

	/*
	 * Compute consShiftCost
	 */

	// if the initial shift has already exceeded the max, substract now the cost that will be readd later
	if( (firstDay_==0) && (lastShift>0) &&
	    (nbConsShifts > pScenario->maxConsShiftsOfTypeOf(lastShift))){
	  consShiftsCost_ -= (nbConsShifts-pScenario->maxConsShiftsOfTypeOf(lastShift))*WEIGHT_CONS_SHIFTS;
	}

	for(int k=firstDay_; k<firstDay_+length_; ++k){
		if(lastShift == shifts_[k]){
			nbConsShifts ++;
			continue;
		}
		if(lastShift > 0){
		  int diff = max(pScenario->minConsShiftsOfTypeOf(lastShift) - nbConsShifts,
				 nbConsShifts-pScenario->maxConsShiftsOfTypeOf(lastShift));
			if(diff>0) {
				consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
			}
		}
		//initialize nbConsShifts and lastShift
		nbConsShifts = 1;
		lastShift = shifts_[k];
	}

	//compute consShiftsCost for the last shift
	int diff = max((firstDay_+length_ == horizon) ? 0 : pScenario->minConsShiftsOfTypeOf(lastShift) - nbConsShifts,
		       nbConsShifts-pScenario->maxConsShiftsOfTypeOf(lastShift));
	if(diff>0) {
		consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
	}


	/*
	 * Compute consDaysWorkedCost
	 */

	// if already worked too much
	double diffDays = nbConsDaysWorked - pNurse->pContract_->maxConsDaysWork_;
	consDaysWorkedCost_ = (diffDays > 0) ? - diffDays * WEIGHT_CONS_DAYS_WORK : 0 ;

	nbConsDaysWorked += length_;
	//check if nbConsDaysWorked < min, if finishes on last day, does not count
	if(nbConsDaysWorked < pNurse->minConsDaysWork() && firstDay_+length_ < horizon)
		consDaysWorkedCost_ += (pNurse->minConsDaysWork() - nbConsDaysWorked) * WEIGHT_CONS_DAYS_WORK;
	//check if nbConsDaysWorked > max
	else if(nbConsDaysWorked > pNurse->maxConsDaysWork())
		consDaysWorkedCost_ += (nbConsDaysWorked - pNurse->maxConsDaysWork()) * WEIGHT_CONS_DAYS_WORK;

	/*
	 * Compute completeWeekendCost
	 */
	if(pNurse->needCompleteWeekends()){
		//if first day is a Sunday, the saturday is not worked
		if(Tools::isSunday(firstDay_))
			completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
		//if last day + 1 is a Sunday, the sunday is not worked
		if(Tools::isSunday(firstDay_+length_))
			completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
	}

	/*
	 * Compute preferencesCost
	 */

	for(int k=firstDay_; k<firstDay_+length_; ++k)
		if(pPreferences->wantsTheShiftOff(nurseId_, k, shifts_[k]))
			preferenceCost_ += WEIGHT_PREFERENCES;

	/*
	 * Compute initial resting cost
	 */

	if(firstDay_==0 && pNurse->pStateIni_->shift_==0){
		int diff = pNurse->minConsDaysOff() - pNurse->pStateIni_->consDaysOff_;
		initRestCost_ = (diff > 0) ? diff*WEIGHT_CONS_DAYS_OFF : 0;
	}

	/*
	 * Compute the sum of the cost and stores it in cost_
	 */
	if(false){
		cout << "# Calcul du cout:" << endl;
		cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
		cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
		cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
		cout << "#       | Preferences       : " << preferenceCost_ << endl;
		cout << "#       | Initial rest      : " << initRestCost_ << endl;
		cout << "# " << endl;
	}

	cost_ = consShiftsCost_ + consDaysWorkedCost_ + completeWeekendCost_ + preferenceCost_ +  initRestCost_;
}


void Rotation::checkDualCost(DualCosts& costs){
	//check if pNurse points to a nurse
	if(nurseId_ == -1)
		Tools::throwError("LiveNurse = NULL");

	/************************************************
	 * Compute all the dual costs of a rotation:
	 ************************************************/

	double dualCost(cost_);

	/* Working dual cost */
	for(int k=firstDay_; k<firstDay_+length_; ++k)
		dualCost -= costs.dayShiftWorkCost(k, shifts_[k]-1);
	/* Start working dual cost */
	dualCost -= costs.startWorkCost(firstDay_);
	/* Stop working dual cost */
	dualCost -= costs.endWorkCost(firstDay_+length_-1);
	/* Working on weekend */
	if(Tools::isSunday(firstDay_))
		dualCost -= costs.workedWeekendCost();
	for(int k=firstDay_; k<firstDay_+length_; ++k)
		if(Tools::isSaturday(k))
			dualCost -= costs.workedWeekendCost();


	// Display: set to true if you want to display the details of the cost

	if(abs(dualCost_ - dualCost) > EPSILON ){
		cout << "# " << endl;
		cout << "# " << endl;
		cout << "Bad dual cost: " << dualCost_ << " != " << dualCost << endl;
		cout << "# " << endl;
		cout << "#   | Base cost     : + " << cost_ << endl;

		cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
		cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
		cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
		cout << "#       | Preferences       : " << preferenceCost_ << endl;
		cout << "#       | Initial rest      : " << initRestCost_ << endl;

		for(int k=firstDay_; k<firstDay_+length_; ++k)
			cout << "#   | Work day-shift: - " << costs.dayShiftWorkCost(k, shifts_[k]-1) << endl;
		cout << "#   | Start work    : - " << costs.startWorkCost(firstDay_) << endl;
		cout << "#   | Finish Work   : - " << costs.endWorkCost(firstDay_+length_-1) << endl;
		if(Tools::isSunday(firstDay_))
			cout << "#   | Weekends      : - " << costs.workedWeekendCost() << endl;
		for(int k=firstDay_; k<firstDay_+length_; ++k)
			if(Tools::isSaturday(k))
				cout << "#   | Weekends      : - " << costs.workedWeekendCost() << endl;
		std::cout << "#   | ROTATION:" << "  cost=" << cost_ << "  dualCost=" << dualCost_ << "  firstDay=" << firstDay_ << "  length=" << length_ << std::endl;
		std::cout << "#               |";
		toString(56);
		cout << "# " << endl;
		cout << "# " << endl;

	}
}

//Compare rotations on index
//
bool Rotation::compareId(const Rotation& rot1, const Rotation& rot2){
	return ( rot1.id_ < rot2.id_ );
}

//Compare rotations on cost
//
bool Rotation::compareCost(const Rotation& rot1, const Rotation& rot2){
	if(rot1.cost_ == DBL_MAX || rot2.cost_ == DBL_MAX)
		Tools::throwError("Rotation cost not computed.");
	return ( rot1.cost_ < rot2.cost_ );
}

//Compare rotations on dual cost
//
bool Rotation::compareDualCost(const Rotation& rot1, const Rotation& rot2){
	if(rot1.dualCost_ == DBL_MAX || rot2.dualCost_ == DBL_MAX)
		Tools::throwError("Rotation cost not computed.");
	return ( rot1.dualCost_ < rot2.dualCost_ );
}

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

// Default constructor
MasterProblem::MasterProblem(Scenario* pScenario, Demand* pDemand,
		Preferences* pPreferences, vector<State>* pInitState, MySolverType solverType):

								   Solver(pScenario, pDemand, pPreferences, pInitState), PrintSolution(),
								   solverType_(solverType), pModel_(0), pPricer_(0), pRule_(0), pTree_(0),
								   positionsPerSkill_(pScenario->nbSkills_), skillsPerPosition_(pScenario->nbPositions()), restsPerDay_(pScenario->nbNurses_),

								   columnVars_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
								   minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
								   minWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedWeekendAvgVars_(pScenario_->nbNurses_),
								   minWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedWeekendContractAvgVars_(pScenario_->nbContracts_),
								   optDemandVars_(pDemand_->nbDays_),numberOfNursesByPositionVars_(pDemand_->nbDays_), skillsAllocVars_(pDemand_->nbDays_),

								   restFlowCons_(pScenario->nbNurses_), workFlowCons_(pScenario->nbNurses_),
								   minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
								   minWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedWeekendAvgCons_(pScenario_->nbNurses_),
								   minWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedWeekendContractAvgCons_(pScenario_->nbContracts_),
								   minDemandCons_(pDemand_->nbDays_), optDemandCons_(pDemand_->nbDays_),numberOfNursesByPositionCons_(pDemand_->nbDays_), feasibleSkillsAllocCons_(pDemand_->nbDays_),
									// STAB
									stabRestFlowPlus_(pScenario->nbNurses_), stabRestFlowMinus_(pScenario->nbNurses_), stabWorkFlowPlus_(pScenario->nbNurses_), stabWorkFlowMinus_(pScenario->nbNurses_),
									stabMinWorkedDaysPlus_(pScenario->nbNurses_), stabMaxWorkedDaysMinus_(pScenario->nbNurses_), stabMaxWorkedWeekendMinus_(pScenario->nbNurses_),
									stabMinDemandPlus_(pDemand_->nbDays_), stabOptDemandPlus_(pDemand_->nbDays_)
{
	// build the model
	this->initializeSolver(solverType);
}

MasterProblem::~MasterProblem(){
	if (pPricer_) delete pPricer_;
	if (pRule_) delete pRule_;
	if (pTree_) delete pTree_;
	if (pModel_) delete pModel_;
}

// Initialize the solver at construction
//
void MasterProblem::initializeSolver(MySolverType solverType) {

	// This OSI interface is created only to retrieve the proper value of
	// infinity for the solver
	OsiSolverInterface* solver=NULL;

	switch(solverType){
	case S_CLP:
		pModel_ = new BcpModeler(this,PB_NAME, CLP);
		solver = new OsiClpSolverInterface();
		break;
	case S_Gurobi:
#ifdef USE_GUROBI
		pModel_ = new BcpModeler(this,PB_NAME, Gurobi);
		solver = new OsiGrbSolverInterface();
#else
	  Tools::throwError("BCP has not been built with Gurobi.");
#endif
		break;
	case S_Cplex:
#ifdef USE_CPLEX
		pModel_ = new BcpModeler(this,PB_NAME, Cplex);
		solver = new OsiCpxSolverInterface();
#else
	  Tools::throwError("BCP has not been built with Cplex.");
#endif
		break;
	case S_CBC:
		pModel_ = new BcpModeler(this,PB_NAME);
		solver = new OsiClpSolverInterface();
		break;
	default:
		Tools::throwError("MasterProblem::initializeSolver: the requested solver is not supported presently.");
		break;
	}

	pModel_->setInfinity(solver->getInfinity());
	// DBG : QUESTION
	// DBG : REPONSE, cet Osi n'est pas celui qui sera utilise pour resoudre les
	// LPs, il ne sert qu'a recuperer la valeur infinity. Il est donc necessaire
	// de le supprimer  pour eviter les fuites
	delete solver;

	this->preprocessData();

	/*
	 * Build the two vectors linking positions and skills
	 */
	for(int p=0; p<skillsPerPosition_.size(); p++){
		vector<int> skills(pScenario_->pPositions()[p]->skills_.size());
		for(int sk=0; sk<pScenario_->pPositions()[p]->skills_.size(); ++sk)
			skills[sk]=pScenario_->pPositions()[p]->skills_[sk];
		skillsPerPosition_[p] = skills;
	}
	for(int sk=0; sk<positionsPerSkill_.size(); sk++){
		vector<int> positions(pScenario_->nbPositions());
		int i(0);
		for(int p=0; p<positions.size(); p++)
			if(find(skillsPerPosition_[p].begin(), skillsPerPosition_[p].end(), sk) != skillsPerPosition_[p].end()){
				positions[i]=p;
				++i;
			}
		positions.resize(i);
		positionsPerSkill_[sk] = positions;
	}

	// initialize the vectors indicating whether the min/max total constraints
	// with averaged bounds are considered
	for (int i=0; i < pScenario_->nbNurses_; i++) {
		isMinWorkedDaysAvgCons_.push_back(false);
		isMaxWorkedDaysAvgCons_.push_back(false);
		isMaxWorkedWeekendAvgCons_.push_back(false);
	}
	for(int p=0; p<pScenario_->nbContracts_; ++p){
		isMinWorkedDaysContractAvgCons_.push_back(false);
		isMaxWorkedDaysContractAvgCons_.push_back(false);
		isMaxWorkedWeekendContractAvgCons_.push_back(false);
	}
}

//build the rostering problem
void MasterProblem::build(SolverParam param){
	/* Rotation constraints */
	buildRotationCons(param);

	/* Min/Max constraints */
	buildMinMaxCons(param);

	/* Skills coverage constraints */
	buildSkillsCoverageCons(param);

	/* We add initial rotations to be always feasible */
	std::string baseName("feasibilityRotation");
	//We add a column with 1 everywhere for each nurse to be always feasible
	//build a map of shift -1 everywhere
	map<int,int> shifts;
	for(int k=0; k<pDemand_->nbDays_; ++k)
		shifts.insert(pair<int,int>( k , -1 ));

	for(int i=0; i<pScenario_->nbNurses_; ++i){
		// DBG: Compute the cost of artificial variables in accordance to the soft
		// constraints
		double artificialCost = WEIGHT_TOTAL_SHIFTS*pScenario_->nbShifts_*(pScenario_->nbDays()-pScenario_->maxTotalShiftsOf(i));
		artificialCost += WEIGHT_CONS_DAYS_WORK*pScenario_->nbShifts_*(pScenario_->nbDays()-pScenario_->maxConsDaysWorkOf(i));
		for (int i = 1; i < pScenario_->nbShifts_; i++) {
		  artificialCost += WEIGHT_CONS_SHIFTS*(pScenario_->nbDays()-pScenario_->maxConsShiftsOfTypeOf(i));
		}
		Rotation rotation(shifts, i, LARGE_SCORE);// artificialCost);//
		addRotation(rotation, baseName.c_str(), true);
	}

	/* Initialize the objects used in the branch and price unless the CBC is used
      to solve the problem
	 */
	if (solverType_ != S_CBC) {
		/* Rotation pricer */
		pPricer_ = new RotationPricer(this, "pricer", param);
		pModel_->addObjPricer(pPricer_);

		/* Tree */
		RestTree* pTree = new RestTree(pScenario_, pDemand_);
		pTree_ = pTree;
		pModel_->addTree(pTree_);

		/* Branching rule */
		pRule_ = new DiveBranchingRule(this, pTree, "branching rule");
		pModel_->addBranchingRule(pRule_);
	}
}

//solve the rostering problem
double MasterProblem::solve(vector<Roster> solution){
	return  solve(solution, true);
}

// Solve the rostering problem with parameters

double MasterProblem::solve(SolverParam param, vector<Roster> solution){
	param.saveFunction_ = this;
	pModel_->setParameters(param);
	return  solve(solution, true);
}

//solve the rostering problem
double MasterProblem::solve(vector<Roster> solution, bool rebuild){

	// build the model first
	if(rebuild)
		this->build(pModel_->getParameters());
	else
		pModel_->reset();

	// input an initial solution
	this->initializeSolution(solution);

	// DBG
	// pModel_->writeProblem("outfiles/model.lp");

	// RqJO: warning, it would be better to define an enumerate type of verbosity
	// levels and create the matching in the Modeler subclasses
	if (solverType_ != S_CBC ) {
		pModel_->setVerbosity(1);
	}
	this->solveWithCatch();

	if (pModel_->getParameters().printBranchStats_ ) {
		pModel_->printStats();
	}

	if(!pModel_->printBestSol()) {
		return pModel_->getRelaxedObjective();
	}

	storeSolution();
	costsConstrainstsToString();

	return pModel_->getObjective();
}

void MasterProblem::solveWithCatch(){
	pModel_->solve();
}

//Resolve the problem with another demand and keep the same preferences
//
double MasterProblem::resolve(Demand* pDemand, SolverParam param, vector<Roster> solution){
	updateDemand(pDemand);
	param.saveFunction_ = this;
	pModel_->setParameters(param);
	return solve(solution, false);
}

// Initialization of the master problem with/without solution
void MasterProblem::initialize(SolverParam param, vector<Roster> solution) {
	this->build(param);
	this->initializeSolution(solution);
	param.saveFunction_ = this;
	pModel_->setParameters(param);

	// in case the initial solution is not empty, fix the corresponding rotations
	// to one and solve the problem to get the solution properly
	if (!solution.empty()) {
		pModel_->fixEveryRotation();
	}
}

//initialize the rostering problem with one column to be feasible if there is no initial solution
//otherwise build the columns corresponding to the initial solution
void MasterProblem::initializeSolution(vector<Roster> solution){
	string baseName("initialRotation");
	//rotations are added for each nurse of the initial solution
	if(solution.size() != 0){
		//build the rotations of each nurse
		for(int i=0; i<pScenario_->nbNurses_; ++i){
			//load the roster of nurse i
			Roster roster = solution[i];

			bool workedLastDay = false;
			int lastShift = 0;
			map<int,int> shifts;
			//build all the successive rotation of this nurse
			for(int k=0; k<pDemand_->nbDays_; ++k){
				//shift=0 => rest
				int shift = roster.shift(k);
				//if work, insert the shift in the map
				if(shift>0){
					shifts.insert(pair<int,int>(k, shift));
					lastShift = shift;
					workedLastDay = true;
				}
				else if(shift<0 && lastShift>0){
					shifts.insert(pair<int,int>(k, lastShift));
					workedLastDay = true;
				}
				//if stop to work, build the rotation
				else if(workedLastDay){
					Rotation rotation(shifts, i);
					rotation.computeCost(pScenario_, pPreferences_, theLiveNurses_, pDemand_->nbDays_);
					rotation.computeTimeDuration(pScenario_);
					pModel_->addActiveColumn(addRotation(rotation, baseName.c_str()));
					shifts.clear();
					lastShift = shift;
					workedLastDay = false;
				}
			}
			//if work on the last day, build the rotation
			if(workedLastDay){
				Rotation rotation(shifts, i);
				rotation.computeCost(pScenario_, pPreferences_, theLiveNurses_,pDemand_->nbDays_);
				rotation.computeTimeDuration(pScenario_);
				pModel_->addActiveColumn(addRotation(rotation, baseName.c_str()));
				shifts.clear();
			}
		}
	}
}

//------------------------------------------------
// Solution with rolling horizon process
//------------------------------------------------

// relax/unrelax the integrality constraints of the variables corresponding to input days
//
void MasterProblem::relaxDays(vector<bool> isRelax) {
	if (isRelaxDay_.empty()) {
		isRelaxDay_.insert(isRelaxDay_.begin(),pDemand_->nbDays_,false);
		isPartialRelaxDays_ = true;
	}

	for (int day=0; day < pDemand_->nbDays_; day++) {
		isRelaxDay_[day]= isRelaxDay_[day]? true:isRelax[day];
	}

	pModel_->relaxRotationsStartingFromDays(isRelax);
}
void MasterProblem::unrelaxDays(vector<bool> isUnrelax) {
	if (isRelaxDay_.empty()) {
		isPartialRelaxDays_= false;
	}
	else {
		for (int day=0; day < pDemand_->nbDays_; day++)
			isRelaxDay_[day]= isRelaxDay_[day]? !isUnrelax[day]:false;

		pModel_->unrelaxRotationsStartingFromDays(isUnrelax);
	}

}

// fix/unfix all the variables corresponding to the input vector of days
void MasterProblem::fixDays(vector<bool> isFix) {
	// initialize the list of fixed days if empty
	if (isFixDay_.empty()) {
		isFixDay_.insert(isFixDay_.begin(),pDemand_->nbDays_,false);
		isPartialFixDays_ =true;
	}

	// set the list of fixed day
	// + forbid the generation of rotations with these starting days
	for (int day=0; day < pDemand_->nbDays_; day++) {
		isFixDay_[day]= isFixDay_[day]?true:isFix[day];
		if (isFixDay_[day]) pPricer_->forbidStartingDay(day);
	}

	// actually fix the active columns starting with the input days with their current values
	pModel_->fixRotationsStartingFromDays(isFix);
}
void MasterProblem::unfixDays(vector<bool> isUnfix) {
	// easy to treat the case where no day is fixed yet
	if (isFixDay_.empty()) {
		isPartialFixDays_ = false;
		pPricer_->clearForbiddenStartingDays();
	}
	else {
		// set the list of unfixed day
		// + authorize the generation of rotations with these starting days
		for (int day=0; day < pDemand_->nbDays_; day++) {
			isFixDay_[day]= isFixDay_[day]?!isUnfix[day]:false;
			if (!isFixDay_[day]) pPricer_->authorizeStartingDay(day);
		}

		// actually unfix the active columns starting with the input days
		pModel_->unfixRotationsStartingFromDays(isUnfix);
	}
}

// fix/unfix all the variables corresponding to the input vector of nurse ids
void MasterProblem::fixNurses(vector<bool> isFix){
	// initialize the list of fixed nurses if empty
	if (isFixNurse_.empty()) {
		isFixNurse_.insert(isFixNurse_.begin(),getNbNurses(),false);
		isPartialFixNurses_ = true;
	}
	// set the list of fixed day
	// + forbid the generation of rotations of the input nurses
	for (LiveNurse* pNurse: theLiveNurses_){
		int n = pNurse->id_;
		isFixNurse_[n] = isFixNurse_[n]?true:isFix[n];
		if (isFixNurse_[n]) pPricer_->forbidNurse(n);
	}

	// actually fix the active columns of the input nurses with their current values
	pModel_->fixRotationsOfNurses(isFix);
}
void MasterProblem::unfixNurses(vector<bool> isUnfix){
	// easy to treat the case where no nurse is fixed yet
	if (isFixNurse_.empty()) {
		isPartialFixNurses_ = false;
		pPricer_->clearForbiddenNurses();
	}
	else {
		// set the list of unfixed nurses
		// + authorize the generation of rotations for the input nurses
	for (LiveNurse* pNurse: theLiveNurses_){
			int n = pNurse->id_;
			isFixNurse_[n]= isFixNurse_[n]?!isUnfix[n]:false;
			if (!isFixNurse_[n]) pPricer_->authorizeNurse(n);
		}

		// actually unfix the active columns of the input nurses
		pModel_->unfixRotationsOfNurses(isUnfix);
	}
}

//------------------------------------------------------------------------------
// Solve the problem with a method that allows for a warm start
//------------------------------------------------------------------------------
double MasterProblem::rollingSolve(SolverParam param, int firstDay) {

	// build the model and initialize with artificial columns at the first iteration
	if(firstDay == 0) {
		initialize(param);
	}
	else {
		pModel_->reset(true);
		param.saveFunction_ = this;
		pModel_->setParameters(param);
	}

	// solve the problem
	if (solverType_ != S_CBC ) {
		pModel_->setVerbosity(1);
	}

	solveWithCatch();
	pModel_->loadBestSol();

	// output information and save the solution
	if (pModel_->getParameters().printBranchStats_) {
		pModel_->printStats();
	}

	return pModel_->getObjective();
}


//------------------------------------------------------------------------------
// Solve the problem with a method that can be specific to our implementation
// of LNS
//------------------------------------------------------------------------------
double MasterProblem::LNSSolve(SolverParam param) {

	// in lns, we always re-optimize, and we assume that the best solution until
	// there is already loaded
	pModel_->reset(true);
	param.saveFunction_ = this;
	pModel_->setParameters(param);

	// solve the problem
	pModel_->setVerbosity(1);

	solveWithCatch();
	pModel_->loadBestSol();

	// output information and save the solution
	if (pModel_->getParameters().printBranchStats_) {
		pModel_->printStats();
	}

	storeSolution();
	costsConstrainstsToString();

	return pModel_->getObjective();
}


//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void MasterProblem::storeSolution(){
	//retrieve a feasible allocation of skills
	vector< vector< vector< vector<double> > > > skillsAllocation(pDemand_->nbDays_);

	for(int k=0; k<pDemand_->nbDays_; ++k){
		vector< vector< vector<double> > > skillsAllocation2(pScenario_->nbShifts_-1);

		for(int s=0; s<pScenario_->nbShifts_-1; ++s){
			vector< vector<double> > skillsAllocation3(pScenario_->nbSkills_);

			for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
				skillsAllocation3[sk] = pModel_->getVarValues(skillsAllocVars_[k][s][sk]);

			skillsAllocation2[s] = skillsAllocation3;
		}
		skillsAllocation[k] = skillsAllocation2;
	}

	//build the rosters
	for(LiveNurse* pNurse: theLiveNurses_)
		pNurse->roster_.reset();

	for(MyVar* var: pModel_->getActiveColumns()){
		if(pModel_->getVarValue(var) > EPSILON){
			Rotation rot(var->getPattern());
			LiveNurse* pNurse = theLiveNurses_[rot.nurseId_];
			for(int k=rot.firstDay_; k<rot.firstDay_+rot.length_; ++k){
				bool assigned = false;
				for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
					if(skillsAllocation[k][rot.shifts_[k]-1][sk][pNurse->pPosition_->id_] > EPSILON){
						pNurse->roster_.assignTask(k,rot.shifts_[k],sk);
						skillsAllocation[k][rot.shifts_[k]-1][sk][pNurse->pPosition_->id_] --;
						assigned = true;
						break;
					}
				if(!assigned){
					char error[255];
					sprintf(error, "No skill found for Nurse %d on day %d on shift %d", pNurse->id_, k, rot.shifts_[k]);
					Tools::throwError((const char*) error);
				}
			}
		}
	}

	//build the states of each nurse
	solution_.clear();
	for(LiveNurse* pNurse: theLiveNurses_){
		pNurse->buildStates();
		solution_.push_back(pNurse->roster_);
	}
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void MasterProblem::save(vector<int>& weekIndices, string outfile){
	storeSolution();
	// initialize the log stream
	// first, concatenate the week numbers
	int nbWeeks = weekIndices.size();

	// write separately the solutions of each week in the required output format
	int firstDay = pDemand_->firstDay_;
	for(int w=0; w<nbWeeks; ++w){
		//		string solutionFile = outdir+std::to_string(weekIndices[w])+".txt";
		Tools::LogOutput solutionStream(outfile);
		solutionStream << solutionToString(firstDay, 7, pScenario_->thisWeek()+w);
		firstDay += 7;
	}
}

//------------------------------------------------------------------------------
// Build the, possibly fractional, roster corresponding to the solution
// currently stored in the model
//------------------------------------------------------------------------------
vector<vector<vector<double>>> MasterProblem::getFractionalRoster() {
	vector<vector<vector<double>>> fractionalRoster(getNbNurses());
	for(vector<vector<double>>& fractionalRoster2: fractionalRoster) {
		Tools::initDoubleVector2D(&fractionalRoster2,getNbDays(),pDemand_->nbShifts_-1,0);
	}

	// Retrieve current fractional roster for each nurse
	// Warning, the working shifts are numbered from 0 to nbShifts_-1 instead of
	// 1 to nbShifts_ in this vector
	double value = 0.0;
	for(MyVar* var : pModel_->getActiveColumns()){
		if (var->getPattern().empty()) continue;
		Rotation rot(var->getPattern());
		vector<vector<double>>& fractionalRoster2 = fractionalRoster[rot.nurseId_];
		value = pModel_->getVarValue(var);
		for(pair<int,int> p : rot.shifts_)
			fractionalRoster2[p.first][p.second-1] += value;
	}

	return fractionalRoster;
}

void MasterProblem::printCurrentSol(){
	allocationToString();
		coverageToString();
}


//------------------------------------------------------------------------------
// Build the variable of the rotation as well as all the affected constraints
// with their coefficients. if s=-1, the nurse i works on all shifts
//------------------------------------------------------------------------------
MyVar* MasterProblem::addRotation(Rotation& rotation, const char* baseName, bool coreVar){
	//nurse index
	int nurseId = rotation.nurseId_;

	//Column var, its name, and affected constraints with their coefficients
	MyVar* var;
	char name[255];
	vector<MyCons*> cons;
	vector<double> coeffs;

	/* Min/Max constraints */
	int nbWeekends = Tools::containsWeekend(rotation.firstDay_, rotation.firstDay_+rotation.length_-1);
	//addMinMaxConsToCol(cons, coeffs, nurseId, rotation.length_, nbWeekends);
	addMinMaxConsToCol(cons, coeffs, nurseId, rotation.timeDuration_, nbWeekends);   // pour prendre en compte les heures plutôt que les jours

	/* Skills coverage constraints */
	for(int k=rotation.firstDay_; k<rotation.firstDay_+rotation.length_; ++k)
		addSkillsCoverageConsToCol(cons, coeffs, nurseId, k, rotation.shifts_[k]);

	sprintf(name, "%s_N%d_%ld",baseName , nurseId, rotation.id_);
	if(coreVar){
		// DBG
		// The artificial variables are taken out of the flow constraints to
		// allow them to be used even after branching
		// Otherwise, the problem might mistakenly appear infeasible after
		// branching
		// addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_, true, false);
		// addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_+rotation.length_-1, false, true);

		pModel_->createPositiveVar(&var, name, rotation.cost_, rotation.getCompactPattern());
		for(int i=0; i<cons.size(); i++)
			pModel_->addCoefLinear(cons[i], var, coeffs[i]);
	}
	else {
		/* Rotation constraints
			They are added only for real rotations to be sure that the artificial variables
			can always be used to create a feasible solution
		 */
		addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_, true, false);
		addRotationConsToCol(cons, coeffs, nurseId, rotation.firstDay_+rotation.length_-1, false, true);

		if (this->isRelaxDay(rotation.firstDay_)) {
			pModel_->createPositiveColumn(&var, name, rotation.cost_, rotation.getCompactPattern(), rotation.dualCost_, cons, coeffs);
		}
		else {
			pModel_->createIntColumn(&var, name, rotation.cost_, rotation.getCompactPattern(), rotation.dualCost_, cons, coeffs);
		}
	}
	return var;
}

/*
 * Rotation constraints
 */
void MasterProblem::buildRotationCons(SolverParam param){
	char name[255];
	//build the rotation network for each nurse
	for(int i=0; i<pScenario_->nbNurses_; i++){
		int minConsDaysOff(theLiveNurses_[i]->minConsDaysOff()),
				maxConsDaysOff(theLiveNurses_[i]->maxConsDaysOff()),
				initConsDaysOff(theLiveNurses_[i]->pStateIni_->consDaysOff_);
		//=true if we have to compute a cost for resting days exceeding the maximum allowed
		//=false otherwise
		bool const maxRest = (maxConsDaysOff < pDemand_->nbDays_ + initConsDaysOff);
		//number of long resting arcs as function of maxRest
		int const nbLongRestingArcs((maxRest) ? maxConsDaysOff : minConsDaysOff);
		//first day when a rest arc exists =
		//nbLongRestingArcs - number of consecutive worked days in the past
		int const firstRestArc( min( max( 0, nbLongRestingArcs - initConsDaysOff ), pDemand_->nbDays_-1 ) );
		//first day when a restingVar exists: at minimun 1
		//if firstRestArc=0, the first resting arc is a longRestingVar
		int const indexStartRestArc = max(1, firstRestArc);
		//number of resting arcs
		int const nbRestingArcs( pDemand_->nbDays_- indexStartRestArc );

		//initialize vectors
		vector< vector< MyVar* > > restsPerDay2(pDemand_->nbDays_);
		vector< MyVar* > restingVars2(nbRestingArcs);
		vector< vector<MyVar*> > longRestingVars2(pDemand_->nbDays_);
		vector<MyCons*> restFlowCons2(pDemand_->nbDays_);
		vector<MyCons*> workFlowCons2(pDemand_->nbDays_);
		vector<MyVar*> stabRestFlowPlus2(pDemand_->nbDays_);
		vector<MyVar*> stabRestFlowMinus2(pDemand_->nbDays_);
		vector<MyVar*> stabWorkFlowPlus2(pDemand_->nbDays_);
		vector<MyVar*> stabWorkFlowMinus2(pDemand_->nbDays_);

		/*****************************************
		 * Creating arcs
		 *****************************************/
		for(int k=0; k<pDemand_->nbDays_; ++k){
			/*****************************************
			 * first long resting arcs
			 *****************************************/
			if(k==0){
				//number of min long resting arcs
				int nbMinRestArcs( max(0, minConsDaysOff - initConsDaysOff) );
				//initialize cost
				int cost (nbMinRestArcs * WEIGHT_CONS_DAYS_OFF);
				Rotation rot = computeInitStateRotation(theLiveNurses_[i]);

				//initialize vectors
				//Must have a minimum of one long resting arcs
				vector<MyVar*> longRestingVars3_0(indexStartRestArc);

				//create minRest arcs
				for(int l=1; l<=nbMinRestArcs; ++l){
					cost -= WEIGHT_CONS_DAYS_OFF;
					sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
					pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, cost+rot.cost_, rot.getCompactPattern());
					initialStateVars_.push_back(longRestingVars3_0[l-1]);
					//add this resting arc for each day of rest
					for(int k1=0; k1<l; ++k1)
						restsPerDay2[k1].push_back(longRestingVars3_0[l-1]);
				}

				//create maxRest arcs, if maxRest=true
				if(maxRest){
					for(int l=1+nbMinRestArcs; l<=firstRestArc; ++l){
						sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
						pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, rot.cost_, rot.getCompactPattern());
						initialStateVars_.push_back(longRestingVars3_0[l-1]);
						//add this resting arc for each day of rest
						for(int k1=0; k1<l; ++k1)
							restsPerDay2[k1].push_back(longRestingVars3_0[l-1]);
					}
				}

				//create the only resting arc (same as a short resting arcs)
				if(firstRestArc == 0){
					sprintf(name, "restingVars_N%d_%d_%d", i, 0, 1);
					pModel_->createPositiveVar(&longRestingVars3_0[0], name, (maxRest) ? WEIGHT_CONS_DAYS_OFF+rot.cost_ : rot.cost_, rot.getCompactPattern());
					initialStateVars_.push_back(longRestingVars3_0[0]);
					//add this resting arc for the first day of rest
					restsPerDay2[0].push_back(longRestingVars3_0[0]);
				}
				//store vectors
				longRestingVars2[0] = longRestingVars3_0;
			}
			/*****************************************
			 * long resting arcs without the first ones
			 *****************************************/
			else{
				//number of long resting arcs = min(nbLongRestingArcs, number of possible long resting arcs)
				int nbLongRestingArcs2( min(nbLongRestingArcs, pDemand_->nbDays_-k) );
				//initialize cost
				//if the arc finishes the last day, the cost is 0. Indeed it will be computed on the next planning
				int cost = minConsDaysOff * WEIGHT_CONS_DAYS_OFF;

				//initialize vectors
				vector<MyVar*> longRestingVars3(nbLongRestingArcs2);

				//create minRest arcs
				for(int l=1; l<=minConsDaysOff; ++l){
					bool doBreak = false;
					cost -= WEIGHT_CONS_DAYS_OFF;
					sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
					//if arc ends before the last day: normal cost
					if(l < pDemand_->nbDays_-k)
						pModel_->createPositiveVar(&longRestingVars3[l-1], name, cost);
					//otherwise, arc finishes on last day
					//so: cost=0 and we break the loop
					else{
						pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
						doBreak = true;
					}
					//add this resting arc for each day of rest
					for(int k1=k; k1<k+l; ++k1)
						restsPerDay2[k1].push_back(longRestingVars3[l-1]);
					if(doBreak)
						break;
				}
				//create maxRest arcs, if maxRest=true
				if(maxRest)
					for(int l=1+minConsDaysOff; l<=maxConsDaysOff; ++l){
						//if exceed last days, break
						if(l > pDemand_->nbDays_-k)
							break;
						sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
						pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
						//add this resting arc for each day of rest
						for(int k1=k; k1<k+l; ++k1)
							restsPerDay2[k1].push_back(longRestingVars3[l-1]);
					}
				//store vectors
				longRestingVars2[k] = longRestingVars3;
			}
			/*****************************************
			 * short resting arcs
			 *****************************************/
			if(k>=indexStartRestArc){
				sprintf(name, "restingVars_N%d_%d_%d", i, k, k+1);
				pModel_->createPositiveVar(&restingVars2[k-indexStartRestArc], name, (maxRest) ? WEIGHT_CONS_DAYS_OFF : 0);
				//add this resting arc for this day of rest
				restsPerDay2[k].push_back(restingVars2[k-indexStartRestArc]);
			}
		}

		/*****************************************
		 * Resting nodes constraints
		 *****************************************/
		for(int k=0; k<pDemand_->nbDays_; ++k){
			vector<double> coeffs(longRestingVars2[k].size());
			for(int l=0; l<longRestingVars2[k].size(); ++l)
				coeffs[l] = 1;
			sprintf(name, "restingNodes_N%d_%d", i, k);
			//Create flow constraints. out flow = 1 if source node (k=0)
			pModel_->createEQConsLinear(&restFlowCons2[k], name, (k==0) ? 1 : 0,
					longRestingVars2[k], coeffs);

			// STAB:Add stabilization variables
			//
			if (param.isStabilization_)	{
				sprintf(name,"stabRestFlowPlus2_%i",k);
				pModel_->createPositiveVar(&stabRestFlowPlus2[k],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
				sprintf(name,"stabRestFlowMinus2_%i",k);
				pModel_->createPositiveVar(&stabRestFlowMinus2[k],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
				pModel_->addCoefLinear(restFlowCons2[k],stabRestFlowPlus2[k],1.0);
				pModel_->addCoefLinear(restFlowCons2[k],stabRestFlowMinus2[k],-1.0);
			}

		}

		/*****************************************
		 * Working nodes constraints
		 *****************************************/
		for(int k=1; k<=pDemand_->nbDays_; ++k){
			//take the min between the number of long resting arcs and the number of possible in arcs
			int nbLongRestingArcs2 = min(nbLongRestingArcs,k);

			vector<MyVar*> vars;
			vector<double> coeffs;
			//add long resting arcs
			for(int l=0; l<nbLongRestingArcs2; ++l){
				//if the long resting arc starts on the source node,
				//check if there exists such an arc
				if( (l > k-1) || (l >= longRestingVars2[k-1-l].size()) )
					break;
				vars.push_back(longRestingVars2[k-1-l][l]);
				//compute in-flow for the sink
				if(k==pDemand_->nbDays_)
					coeffs.push_back(1);
				//compute out-flow
				else
					coeffs.push_back(-1);
			}
			//add resting arcs
			//just 1 out, if first restingVar
			//compute out-flow
			if(k==indexStartRestArc){
				vars.push_back(restingVars2[0]);
				coeffs.push_back(1);
			}
			//just 1 in, if last resting arcs
			//compute in-flow for the sink
			else if (k==pDemand_->nbDays_){
				vars.push_back(restingVars2[restingVars2.size()-1]);
				coeffs.push_back(1);
			}
			//2 otherwise: 1 in and 1 out
			//compute out-flow
			else if(k>indexStartRestArc){
				vars.push_back(restingVars2[k-1-indexStartRestArc]);
				coeffs.push_back(-1);
				vars.push_back(restingVars2[k-indexStartRestArc]);
				coeffs.push_back(1);
			}
			sprintf(name, "workingNodes_N%d_%d", i, k);
			//Create flow constraints. in flow = 1 if sink node (k==pDemand_->nbDays_)
			pModel_->createEQConsLinear(&workFlowCons2[k-1], name, (k==pDemand_->nbDays_) ? 1 : 0,
					vars, coeffs);

			// STAB:Add stabilization variables
			//
			if (param.isStabilization_) {
				sprintf(name,"stabWorkFlowPlus2_%i",k-1);
				pModel_->createPositiveVar(&stabWorkFlowPlus2[k-1],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
				sprintf(name,"stabWorkFlowMinus2_%i",k-1);
				pModel_->createPositiveVar(&stabWorkFlowMinus2[k-1],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
				pModel_->addCoefLinear(workFlowCons2[k-1],stabWorkFlowPlus2[k-1],1.0);
				pModel_->addCoefLinear(workFlowCons2[k-1],stabWorkFlowMinus2[k-1],-1.0);
			}
		}

		//store vectors
		restsPerDay_[i] = restsPerDay2;
		restingVars_[i] = restingVars2;
		longRestingVars_[i] = longRestingVars2;
		restFlowCons_[i] = restFlowCons2;
		workFlowCons_[i] = workFlowCons2;

		// STAB
		stabRestFlowPlus_[i] = stabRestFlowPlus2;
		stabRestFlowMinus_[i] = stabRestFlowMinus2;
		stabWorkFlowPlus_[i] = stabWorkFlowPlus2;
		stabWorkFlowMinus_[i] = stabWorkFlowMinus2;
	}
}

int MasterProblem::addRotationConsToCol(vector<MyCons*>& cons, vector<double>& coeffs, int i, int k, bool firstDay, bool lastDay){
	//check if the rotation starts on day k
	if(firstDay){
		//compute out-flow
		coeffs.push_back(1.0);
		//add to source constraint
		if(k==0)
			cons.push_back(restFlowCons_[i][0]);
		//add to work node constraint
		else
			cons.push_back(workFlowCons_[i][k-1]);

		return 1;
	}

	//check if the rotation finishes on day k
	else if(lastDay){
		//add to sink constraint
		//compute in-flow
		if(k==pDemand_->nbDays_-1){
			coeffs.push_back(1.0);
			cons.push_back(workFlowCons_[i][pDemand_->nbDays_-1]);
		}
		//add to rest node constraint
		//compute out-flow
		else{
			coeffs.push_back(-1.0);
			cons.push_back(restFlowCons_[i][k+1]);
		}

		return 1;
	}

	return 0;
}

/*
 * Min/Max constraints
 */
void MasterProblem::buildMinMaxCons(SolverParam param){
	char name[255];
	for(int i=0; i<pScenario_->nbNurses_; i++){
		sprintf(name, "minWorkedDaysVar_N%d", i);
		pModel_->createPositiveVar(&minWorkedDaysVars_[i], name, weightTotalShiftsMin_[i]);
		sprintf(name, "maxWorkedDaysVar_N%d", i);
		pModel_->createPositiveVar(&maxWorkedDaysVars_[i], name, weightTotalShiftsMax_[i]);

		sprintf(name, "minWorkedDaysCons_N%d", i);
		vector<MyVar*> vars1 = {minWorkedDaysVars_[i]};
		vector<double> coeffs1 = {1};
		pModel_->createGEConsLinear(&minWorkedDaysCons_[i], name, minTotalShifts_[i], vars1, coeffs1);

		// STAB:Add stabilization variable
		//
		if (param.isStabilization_) {
			sprintf(name,"stabMinWorkedDaysPlus_%i",i);
			pModel_->createPositiveVar(&stabMinWorkedDaysPlus_[i],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
			pModel_->addCoefLinear(minWorkedDaysCons_[i],stabMinWorkedDaysPlus_[i],1.0);
		}

		sprintf(name, "maxWorkedDaysCons_N%d", i);
		vector<MyVar*> vars2 = {maxWorkedDaysVars_[i]};
		vector<double> coeffs2 = {-1};
		pModel_->createLEConsLinear(&maxWorkedDaysCons_[i], name, maxTotalShifts_[i], vars2, coeffs2);

		// STAB:Add stabilization variable
		//
		if (param.isStabilization_) {
			sprintf(name,"stabMaxWorkedDaysMinus_%i",i);
			pModel_->createPositiveVar(&stabMaxWorkedDaysMinus_[i],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
			pModel_->addCoefLinear(maxWorkedDaysCons_[i],stabMaxWorkedDaysMinus_[i],-1.0);
		}

		// add constraints on the total number of shifts to satisfy bounds that
		// correspond to the global bounds averaged over the weeks
		//
		// STAB: not implemented there yet
		if (!minTotalShiftsAvg_.empty() && !maxTotalShiftsAvg_.empty() && !weightTotalShiftsAvg_.empty()) {

			// only add the constraint if is tighter than the already added constraint
			if (minTotalShiftsAvg_[i] > minTotalShifts_[i]) {
				sprintf(name, "minWorkedDaysAvgVar_N%d", i);
				pModel_->createPositiveVar(&minWorkedDaysAvgVars_[i], name, weightTotalShiftsAvg_[i]);

				sprintf(name, "minWorkedDaysAvgCons_N%d", i);
				vector<MyVar*> varsAvg1 = {minWorkedDaysVars_[i], minWorkedDaysAvgVars_[i]};
				vector<double> coeffsAvg1 = {1,1};
				pModel_->createGEConsLinear(&minWorkedDaysAvgCons_[i], name, minTotalShiftsAvg_[i], varsAvg1, coeffsAvg1);

				isMinWorkedDaysAvgCons_[i] = true;
			}

			if (maxTotalShiftsAvg_[i] < maxTotalShifts_[i]) {
				sprintf(name, "maxWorkedDaysAvgVar_N%d", i);
				pModel_->createPositiveVar(&maxWorkedDaysAvgVars_[i], name, weightTotalShiftsAvg_[i]);

				sprintf(name, "maxWorkedDaysAvgCons_N%d", i);
				vector<MyVar*> varsAvg2 = {maxWorkedDaysVars_[i],maxWorkedDaysAvgVars_[i]};
				vector<double> coeffsAvg2 = {-1,-1};
				pModel_->createLEConsLinear(&maxWorkedDaysAvgCons_[i], name, maxTotalShiftsAvg_[i], varsAvg2, coeffsAvg2);

				isMaxWorkedDaysAvgCons_[i] = true;
			}
		}

		sprintf(name, "maxWorkedWeekendVar_N%d", i);
		pModel_->createPositiveVar(&maxWorkedWeekendVars_[i], name, weightTotalWeekendsMax_[i]);

		sprintf(name, "maxWorkedWeekendCons_N%d", i);
		vector<MyVar*> vars3 = {maxWorkedWeekendVars_[i]};
		vector<double> coeffs3 = {-1};
		pModel_->createLEConsLinear(&maxWorkedWeekendCons_[i], name, maxTotalWeekends_[i],
				vars3, coeffs3);

		// STAB:Add stabilization variable
		//
		if (param.isStabilization_) {
			sprintf(name,"stabMaxWorkedWeekendMinus_%i",i);
			pModel_->createPositiveVar(&stabMaxWorkedWeekendMinus_[i],name,-param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
			pModel_->addCoefLinear(maxWorkedWeekendCons_[i],stabMaxWorkedWeekendMinus_[i],-1.0);
		}

		// STAB: not implemented there yet
		if ( !maxTotalWeekendsAvg_.empty()  && !weightTotalWeekendsAvg_.empty()
				&& maxTotalWeekendsAvg_[i] < theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) {

			sprintf(name, "maxWorkedWeekendAvgVar_N%d", i);
			pModel_->createPositiveVar(&maxWorkedWeekendAvgVars_[i], name, weightTotalWeekendsAvg_[i]);

			sprintf(name, "maxWorkedWeekendAvgCons_N%d", i);
			vector<MyVar*> varsAvg3 = {maxWorkedWeekendVars_[i],maxWorkedWeekendAvgVars_[i]};
			vector<double> coeffsAvg3 = {-1,-1};
			pModel_->createLEConsLinear(&maxWorkedWeekendAvgCons_[i], name, maxTotalWeekendsAvg_[i]- theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
					varsAvg3, coeffsAvg3);

			isMaxWorkedWeekendAvgCons_[i] = true;

		}

	}

	// WEEKEND CUTS
	// sprintf(name, "sumMaxWorkedWeekendCons");
	// std::vector<MyVar*> varsSum3;
	// std::vector<double> coeffsSum3;
	// for(int i=0; i<pScenario_->nbNurses_; i++) {
	// 	varsSum3.push_back(maxWorkedWeekendVars_[i]);
	// 	coeffsSum3.push_back(1);
	// }
	// pModel_->createGEConsLinear(&sumMaxWorkedWeekendCons_, name, 4, varsSum3, coeffsSum3);



	for(int p=0; p<pScenario_->nbContracts_; ++p){

		if(!minTotalShiftsContractAvg_.empty() && !maxTotalShiftsContractAvg_.empty()  && !weightTotalShiftsContractAvg_.empty()){
			sprintf(name, "minWorkedDaysContractAvgVar_P%d", p);
			pModel_->createPositiveVar(&minWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);
			sprintf(name, "maxWorkedDaysContractAvgVar_P%d", p);
			pModel_->createPositiveVar(&maxWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);

			sprintf(name, "minWorkedDaysContractAvgCons_P%d", p);
			vector<MyVar*> vars1 = {minWorkedDaysContractAvgVars_[p]};
			vector<double> coeffs1 = {1};
			pModel_->createGEConsLinear(&minWorkedDaysContractAvgCons_[p], name, minTotalShiftsContractAvg_[p], vars1, coeffs1);

			sprintf(name, "maxWorkedDaysContractAvgCons_P%d", p);
			vector<MyVar*> vars2 = {maxWorkedDaysContractAvgVars_[p]};
			vector<double> coeffs2 = {-1};
			pModel_->createLEConsLinear(&maxWorkedDaysContractAvgCons_[p], name, maxTotalShiftsContractAvg_[p], vars2, coeffs2);

			isMinWorkedDaysContractAvgCons_[p] = true;
			isMaxWorkedDaysContractAvgCons_[p] = true;
		}

		if(!maxTotalWeekendsContractAvg_.empty()  && !weightTotalWeekendsContractAvg_.empty()){
			sprintf(name, "maxWorkedWeekendContractAvgVar_P%d", p);
			pModel_->createPositiveVar(&maxWorkedWeekendContractAvgVars_[p], name, weightTotalWeekendsContractAvg_[p]);

			sprintf(name, "maxWorkedWeekendContractAvgCons_C%d", p);
			vector<MyVar*> varsAvg3 = {maxWorkedWeekendContractAvgVars_[p]};
			vector<double> coeffsAvg3 = {-1 };
			pModel_->createLEConsLinear(&maxWorkedWeekendContractAvgCons_[p], name, maxTotalWeekendsContractAvg_[p],
					varsAvg3, coeffsAvg3);

			isMaxWorkedWeekendContractAvgCons_[p] = true;
		}
	}
}

int MasterProblem::addMinMaxConsToCol(vector<MyCons*>& cons, vector<double>& coeffs, int i, int nbDays, int nbWeekends){
	int nbCons(0);
	int p = theLiveNurses_[i]->pContract_->id_;
	++nbCons;
	cons.push_back(minWorkedDaysCons_[i]);
	coeffs.push_back(nbDays);
	++nbCons;
	cons.push_back(maxWorkedDaysCons_[i]);
	coeffs.push_back(nbDays);
	if (isMinWorkedDaysAvgCons_[i]) {
		++nbCons;
		cons.push_back(minWorkedDaysAvgCons_[i]);
		coeffs.push_back(nbDays);
	}
	if (isMaxWorkedDaysAvgCons_[i]) {
		++nbCons;
		cons.push_back(maxWorkedDaysAvgCons_[i]);
		coeffs.push_back(nbDays);
	}
	if (isMinWorkedDaysContractAvgCons_[p]) {
		++nbCons;
		cons.push_back(minWorkedDaysContractAvgCons_[p]);
		coeffs.push_back(nbDays);
	}
	if (isMaxWorkedDaysContractAvgCons_[p]) {
		++nbCons;
		cons.push_back(maxWorkedDaysContractAvgCons_[p]);
		coeffs.push_back(nbDays);
	}


	if(nbWeekends){
		++nbCons;
		cons.push_back(maxWorkedWeekendCons_[i]);
		coeffs.push_back(nbWeekends);

		if (isMaxWorkedWeekendAvgCons_[i]) {
			++nbCons;
			cons.push_back(maxWorkedWeekendAvgCons_[i]);
			coeffs.push_back(nbWeekends);
		}

		if (isMaxWorkedWeekendContractAvgCons_[p]) {
			++nbCons;
			cons.push_back(maxWorkedWeekendContractAvgCons_[p]);
			coeffs.push_back(nbWeekends);
		}
	}



	return nbCons;
}

/*
 * Skills coverage constraints
 */
void MasterProblem::buildSkillsCoverageCons(SolverParam param){
	char name[255];
	for(int k=0; k<pDemand_->nbDays_; k++){
		//initialize vectors
		vector< vector<MyVar*> > optDemandVars1(pScenario_->nbShifts_-1);
		vector< vector<MyVar*> > numberOfNursesByPositionVars1(pScenario_->nbShifts_-1);
		vector< vector< vector<MyVar*> > > skillsAllocVars1(pScenario_->nbShifts_-1);
		vector< vector<MyCons*> > minDemandCons1(pScenario_->nbShifts_-1);
		vector< vector<MyCons*> > optDemandCons1(pScenario_->nbShifts_-1);
		vector< vector<MyCons*> > numberOfNursesByPositionCons1(pScenario_->nbShifts_-1);
		vector< vector<MyCons*> > feasibleSkillsAllocCons1(pScenario_->nbShifts_-1);

		// STAB
		vector< vector<MyVar*> > stabMinDemandPlus1(pScenario_->nbShifts_-1);
		vector< vector<MyVar*> > stabOptDemandPlus1(pScenario_->nbShifts_-1);


		//forget s=0, it's a resting shift
		for(int s=1; s<pScenario_->nbShifts_; s++){
			//initialize vectors
			vector<MyVar*> optDemandVars2(pScenario_->nbSkills_);
			vector<MyVar*> numberOfNursesByPositionVars2(pScenario_->nbPositions());
			vector< vector<MyVar*> > skillsAllocVars2(pScenario_->nbSkills_);
			vector<MyCons*> minDemandCons2(pScenario_->nbSkills_);
			vector<MyCons*> optDemandCons2(pScenario_->nbSkills_);
			vector<MyCons*> numberOfNursesByPositionCons2(pScenario_->nbPositions());
			vector<MyCons*> feasibleSkillsAllocCons2(pScenario_->nbPositions());

			// STAB
			vector<MyVar*> stabMinDemandPlus2(pScenario_->nbSkills_);
			vector<MyVar*> stabOptDemandPlus2(pScenario_->nbSkills_);

			for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				//initialize vectors
				vector<MyVar*> skillsAllocVars3(pScenario_->nbPositions());

				//create variables
				sprintf(name, "optDemandVar_%d_%d_%d", k, s, sk);
				pModel_->createPositiveVar(&optDemandVars2[sk], name, WEIGHT_OPTIMAL_DEMAND);
				for(int p=0; p<pScenario_->nbPositions(); p++){
					sprintf(name, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
					// DBG
					pModel_->createPositiveVar(&skillsAllocVars3[p], name, 0);
					//pModel_->createIntVar(&skillsAllocVars3[p], name, 0);
				}
				//store vectors
				skillsAllocVars2[sk] = skillsAllocVars3;

				//adding variables and building minimum demand constraints
				vector<MyVar*> vars1(positionsPerSkill_[sk].size());
				vector<double> coeffs1(positionsPerSkill_[sk].size());
				for(int p=0; p<positionsPerSkill_[sk].size(); ++p){
					vars1[p] = skillsAllocVars3[positionsPerSkill_[sk][p]];
					coeffs1[p] = 1;
				}
				sprintf(name, "minDemandCons_%d_%d_%d", k, s, sk);
				pModel_->createFinalGEConsLinear(&minDemandCons2[sk], name, pDemand_->minDemand_[k][s][sk],
						vars1, coeffs1);

				// STAB:Add stabilization variable
				if (param.isStabilization_) {
					sprintf(name,"stabMinDemandPlus_%d_%d_%d",k,s,sk);
					pModel_->createPositiveVar(&stabMinDemandPlus2[sk],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
					pModel_->addCoefLinear(minDemandCons2[sk],stabMinDemandPlus2[sk],1.0);
				}

				//adding variables and building optimal demand constraints
				vars1.push_back(optDemandVars2[sk]);
				coeffs1.push_back(1);
				sprintf(name, "optDemandCons_%d_%d_%d", k, s, sk);
				pModel_->createFinalGEConsLinear(&optDemandCons2[sk], name, pDemand_->optDemand_[k][s][sk],
						vars1, coeffs1);

				// STAB:Add stabilization variable
				if (param.isStabilization_) {
					sprintf(name,"stabOptDemandPlus_%d_%d_%d",k,s,sk);
					pModel_->createPositiveVar(&stabOptDemandPlus2[sk],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
					pModel_->addCoefLinear(optDemandCons2[sk],stabOptDemandPlus2[sk],1.0);
				}
			}

			for(int p=0; p<pScenario_->nbPositions(); p++){
				//creating variables
				sprintf(name, "nursesNumber_%d_%d_%d", k, s, p);
				// DBG
				// pModel_->createIntVar(&numberOfNursesByPositionVars2[p], name, 0);
				pModel_->createPositiveVar(&numberOfNursesByPositionVars2[p], name, 0);
				//adding variables and building number of nurses constraints
				vector<MyVar*> vars3;
				vector<double> coeff3;
				vars3.push_back(numberOfNursesByPositionVars2[p]);
				coeff3.push_back(-1);
				sprintf(name, "nursesNumberCons_%d_%d_%d", k, s, p);
				pModel_->createEQConsLinear(&numberOfNursesByPositionCons2[p], name, 0,
						vars3, coeff3);

				//adding variables and building skills allocation constraints
				int const nonZeroVars4(1+skillsPerPosition_[p].size());
				vector<MyVar*> vars4(nonZeroVars4);
				vector<double> coeff4(nonZeroVars4);
				vars4[0] = numberOfNursesByPositionVars2[p];
				coeff4[0] = 1;
				for(int sk=1; sk<nonZeroVars4; ++sk){
					vars4[sk] = skillsAllocVars2[skillsPerPosition_[p][sk-1]][p];
					coeff4[sk] =-1;
				}
				sprintf(name, "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
				pModel_->createEQConsLinear(&feasibleSkillsAllocCons2[p], name, 0,
						vars4, coeff4);
			}

			//store vectors
			optDemandVars1[s-1] = optDemandVars2;
			numberOfNursesByPositionVars1 [s-1] = numberOfNursesByPositionVars2;
			skillsAllocVars1[s-1] = skillsAllocVars2;
			minDemandCons1[s-1] = minDemandCons2;
			optDemandCons1[s-1] = optDemandCons2;
			numberOfNursesByPositionCons1 [s-1] = numberOfNursesByPositionCons2;
			feasibleSkillsAllocCons1[s-1] = feasibleSkillsAllocCons2;

			// STAB
			stabMinDemandPlus1[s-1] = stabMinDemandPlus2;
			stabOptDemandPlus1[s-1] = stabOptDemandPlus2;
		}

		//store vectors
		optDemandVars_[k] = optDemandVars1;
		numberOfNursesByPositionVars_[k] = numberOfNursesByPositionVars1;
		skillsAllocVars_[k] = skillsAllocVars1;
		minDemandCons_[k] = minDemandCons1;
		optDemandCons_[k] = optDemandCons1;
		numberOfNursesByPositionCons_[k] = numberOfNursesByPositionCons1;
		feasibleSkillsAllocCons_[k] = feasibleSkillsAllocCons1;

		// STAB
		stabMinDemandPlus_[k] = stabMinDemandPlus1;
		stabOptDemandPlus_[k] = stabOptDemandPlus1;
	}
}

int MasterProblem::addSkillsCoverageConsToCol(vector<MyCons*>& cons, vector<double>& coeffs, int i, int k, int s){
	int nbCons(0);

	int p(theLiveNurses_[i]->pPosition_->id_);
	if(s==-1){
		for(int s0=1; s0<pScenario_->nbShifts_; ++s0){
			++nbCons;
			cons.push_back(numberOfNursesByPositionCons_[k][s0-1][p]);
			coeffs.push_back(1.0);
		}
	}
	else{
		++nbCons;
		cons.push_back(numberOfNursesByPositionCons_[k][s-1][p]);
		coeffs.push_back(1.0);
	}

	return nbCons;
}

void MasterProblem::updateDemand(Demand* pDemand){
	if(pDemand->nbDays_ != pDemand_->nbDays_)
		Tools::throwError("The new demand must have the same size than the old one, so that's ");

	//set the pointer
	pDemand_ = pDemand;

	//modify the associated constraints
	for(int k=0; k<pDemand_->nbDays_; k++)
		for(int s=1; s<pScenario_->nbShifts_; s++)
			for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				minDemandCons_[k][s-1][sk]->setLhs(pDemand_->minDemand_[k][s][sk]);
				optDemandCons_[k][s-1][sk]->setLhs(pDemand_->optDemand_[k][s][sk]);
			}
}

string MasterProblem::costsConstrainstsToString(){
	stringstream rep;

	double initStateRestCost = getRotationCosts(INIT_REST_COST);
	char buffer[100];
	sprintf(buffer, "%-30s %10.0f \n", "Column costs", getRotationCosts() - initStateRestCost);
	rep << buffer;
	rep << "-----------------------------------------\n";
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. shifts costs", getRotationCosts(CONS_SHIFTS_COST));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. worked days costs", getRotationCosts(CONS_WORKED_DAYS_COST));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Complete weekend costs", getRotationCosts(COMPLETE_WEEKEND_COST));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Preferences costs", getRotationCosts(PREFERENCE_COST));
	rep << buffer;
	rep << "-----------------------------------------\n";
	double initStateCost = getRotationCosts(TOTAL_COST, true);
	sprintf(buffer, "%-30s %10.0f \n", "History costs (counted)", initStateCost + initStateRestCost);
	rep << buffer;
	sprintf(buffer, "%-30s %10.0f \n", "Resting costs", pModel_->getTotalCost(restingVars_)+pModel_->getTotalCost(longRestingVars_) + initStateRestCost - initStateCost);
	rep << buffer;
	sprintf(buffer, "%-30s %10.0f \n", "Min worked days costs", pModel_->getTotalCost(minWorkedDaysVars_));
	rep << buffer;
	sprintf(buffer, "%-30s %10.0f \n", "Max worked days costs", pModel_->getTotalCost(maxWorkedDaysVars_));
	rep << buffer;
	sprintf(buffer, "%-30s %10.0f \n", "Max worked weekend costs", pModel_->getTotalCost(maxWorkedWeekendVars_));
	rep << buffer;
	sprintf(buffer, "%-30s %10.0f \n", "Coverage costs", pModel_->getTotalCost(optDemandVars_));//, true));
	rep << buffer;
	rep << "-----------------------------------------\n";
	rep << "\n";

	cout << rep.str();

	return rep.str();
}

double MasterProblem::getRotationCosts(CostType costType, bool initStateRotation){
	//if(initStateRotation), search for empty rotation (=rotation for initial state)
	if(initStateRotation)  return getRotationCosts(costType, initialStateVars_);
	return getRotationCosts(costType, pModel_->getActiveColumns());
}

double MasterProblem::getRotationCosts(CostType costType, const vector<MyVar*>& vars){
	double cost = 0;
	for(MyVar* var: vars){
		double value = pModel_->getVarValue(var);
		if(value > EPSILON){
			Rotation rot(var->getPattern());
			rot.computeCost(pScenario_, pPreferences_, theLiveNurses_, pDemand_->nbDays_);
			switch(costType){
			case CONS_SHIFTS_COST: cost += rot.consShiftsCost_*value;
			break;
			case CONS_WORKED_DAYS_COST: cost += rot.consDaysWorkedCost_*value;
			break;
			case COMPLETE_WEEKEND_COST: cost += rot.completeWeekendCost_*value;
			break;
			case PREFERENCE_COST: cost += rot.preferenceCost_*value;
			break;
			case INIT_REST_COST: cost += rot.initRestCost_*value;
			break;
			default: cost += rot.cost_*value;
			//            if(!initStateRotation && rot.second.length_>0){
			//               rot.second.toString(pDemand_->nbDays_);
			//               pModel_->toString(rot.first);
			//            }
			break;
			}
		}
	}
	return cost;
}

Rotation MasterProblem::computeInitStateRotation(LiveNurse* pNurse){
	//initialize rotation
	Rotation rot(map<int,int>(), pNurse->id_);

	//compute cost for previous cons worked shifts and days
	int lastShift = pNurse->pStateIni_->shift_;
	if(lastShift>0){
		int nbConsWorkedDays = pNurse->pStateIni_->consDaysWorked_;
		int diff = pNurse->minConsDaysWork() - nbConsWorkedDays;
		rot.consDaysWorkedCost_ += (diff>0) ? diff*WEIGHT_CONS_DAYS_WORK : 0;

		int nbConsShifts = pNurse->pStateIni_->consShifts_;
		int diff2 = pScenario_->minConsShiftsOfTypeOf(lastShift) - nbConsShifts;
		rot.consShiftsCost_ += (diff2>0) ? diff2*WEIGHT_CONS_SHIFTS : 0;
	}
	rot.cost_ = rot.consDaysWorkedCost_ + rot.consShiftsCost_;

	return rot;
}



string MasterProblem::allocationToString(bool printInteger){
	stringstream rep;

	int nbNurses = pScenario_->nbNurses_;
	int nbShifts = pScenario_->nbShifts_;
	int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nbDays_;

	rep << std::endl;
	rep << "Allocations of the (potentially fractional) current solution:" << std::endl;
	rep << "\t\t  ";
	for (int day = firstDay; day < firstDay+nbDays; day++) {
		rep << "| " << Tools::intToDay(day).at(0) << " ";
	}
	rep << "|" << std::endl;
	rep << "-------------------------------------"<< std::endl;

	for (int n = 0; n < nbNurses; n ++) {
		LiveNurse* pNurse = theLiveNurses_[n];
		vector<vector<double>> fractionalRoster; Tools::initDoubleVector2D(&fractionalRoster,nbDays,nbShifts-1,0);
		for(MyVar* var : pModel_->getActiveColumns()){
			if(var->getPattern()[0] != pNurse->id_)
				continue;
			Rotation rot(var->getPattern());
			for(pair<int,int> p : rot.shifts_)
				fractionalRoster[p.first][p.second-1] += pModel_->getVarValue(var);
		}
		rep << pNurse->name_ << "\t";
		for(int s=1; s<nbShifts; ++s){
			if (s>1) rep << "\t";
			rep << pScenario_->intToShift_[s].at(0) << "\t";
			for (int day = firstDay; day < firstDay+nbDays; day++){
				double shiftValue = fractionalRoster[day][s-1];
				if(shiftValue > 1-EPSILON){
					rep << "|  1 ";
				} else if(shiftValue > EPSILON){
					char buffer[100];
					sprintf(buffer, "|%1.2f", shiftValue);
					rep << buffer;
				} else {
					rep << "| -- ";
				}
				if(Tools::isSunday(day)) rep << "| ";
			}
			rep << "|" << std::endl;
		}
		rep << std::endl;
	}
	rep << std::endl;

	cout << rep.str();

	return rep.str();
}

string MasterProblem::coverageToString(bool printInteger){
	stringstream rep;

	int nbShifts = pScenario_->nbShifts_;
	int nbSkills = pScenario_->nbSkills_;
	int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nbDays_;

	rep << std::endl;
	rep << "Coverage of the (potentially fractional) current solution:" << std::endl;
	rep << "\t\t   ";
	for (int day = firstDay; day < firstDay+nbDays; day++) {
		rep << "| " << Tools::intToDay(day).at(0) << " ";
	}
	rep << "|" << std::endl;
	rep << "-------------------------------------"<< std::endl;

	string tab = "\t";

	for(int s=1; s<nbShifts; ++s){
		rep << pScenario_->intToShift_[s].at(0) << "\t";
		for(int sk=0; sk<nbSkills; sk++){
			if(sk!=0) rep << "\t";

			char buffer0[20];
			string skill = pScenario_->intToSkill_[sk];
			sprintf(buffer0, "%-12s", skill.c_str());
			rep << buffer0;

			for (int day = firstDay; day < firstDay+nbDays; day++){
				double shiftValue = pModel_->getVarValue(skillsAllocVars_[day][s-1][sk]);
				char buffer[100];
				if(abs(shiftValue - round(shiftValue)) < EPSILON)
					sprintf(buffer, "|%4d", (int) round(shiftValue));
				else sprintf(buffer, "|%2.2f", shiftValue);
				rep << buffer;
				if(Tools::isSunday(day)) rep << "| ";
			}

			rep << "|" << std::endl;
		}
	}
	rep << std::endl;

	cout << rep.str();

	return rep.str();
}

string MasterProblem::workedWeekendsToString(bool printInteger){
	return "";
}

//---------------------------------------------------------------------------
//
// STAB: Methods required to implement stabilization in the column generation
//
//---------------------------------------------------------------------------

// STAB
// Multiply the upper bound of the input variable by the input factor
void MasterProblem::multiplyUbInSolver(MyVar* pVar, OsiSolverInterface* solver, double factor) {
	int varind = pVar->getIndex();
	double ub = pVar->getUB();

	if (ub != solver->getColUpper()[varind]) {
		Tools::throwError("multiplyUbInSolver: the upper bound stored in the variable is not the same as that in the solver!");
	}

	solver->setColUpper(varind,factor*ub);
	pVar->setUB(factor*ub);
}

// STAB
// Set the bound of the input variable to the input value
void MasterProblem::updateVarUbInSolver(MyVar* pVar, OsiSolverInterface* solver, double value) {
	int varind = pVar->getIndex();
	double ub = pVar->getUB();

	if (ub != solver->getColUpper()[varind]) {
		Tools::throwError("updateVarUbInSolver: the upper bound stored in the variable is not the same as that in the solver!");
	}

	solver->setColUpper(varind,value);
	pVar->setUB(value);
}

// STAB
// Set the cost of the input variable to the input value
void MasterProblem::updateVarCostInSolver(MyVar* pVar, OsiSolverInterface* solver, double value) {
	int varind = pVar->getIndex();
	double cost = pVar->getCost();

	if (cost != solver->getObjCoefficients()[varind]) {
		Tools::throwError("updateVarCostInSolver: the cost stored in the variable is not the same as that in the solver!");
	}

	solver->setObjCoeff(varind,value);
	pVar->setCost(value);
}

// STAB
// Update all the upper bounds of the stabilization variables by multiplying
// them by an input factor
void MasterProblem::stabUpdateBound(OsiSolverInterface* solver, double factor) {
	for(int i=0; i<pScenario_->nbNurses_; i++){
		// stabilization variables corresponding to the global constraints of
		// of the nurses
		multiplyUbInSolver(stabMinWorkedDaysPlus_[i], solver, factor);
		multiplyUbInSolver(stabMaxWorkedDaysMinus_[i], solver, factor);
		multiplyUbInSolver(stabMaxWorkedWeekendMinus_[i], solver, factor);

		// stabilization variables corresponding to the flow constraints
		for(int k=0; k<pDemand_->nbDays_; ++k) {
			multiplyUbInSolver(stabRestFlowMinus_[i][k], solver, factor);
			multiplyUbInSolver(stabRestFlowPlus_[i][k], solver, factor);
			multiplyUbInSolver(stabWorkFlowMinus_[i][k], solver, factor);
			multiplyUbInSolver(stabWorkFlowPlus_[i][k], solver, factor);
		}
	}

	// stabilization variables corresponding to the cover constraints
	for(int k=0; k<pDemand_->nbDays_; k++){
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				multiplyUbInSolver(stabMinDemandPlus_[k][s-1][sk], solver, factor);
				multiplyUbInSolver(stabOptDemandPlus_[k][s-1][sk], solver, factor);
			}
		}
	}
}

// STAB
// Update all the costs of the stabilization variables to the values
// corresponding dual variables with a small margin in input
void MasterProblem::stabUpdateCost(OsiSolverInterface* solver, double margin) {
	for(int i=0; i<pScenario_->nbNurses_; i++){
		// stabilization variables corresponding to the global constraints of
		// of the nurses
		double minWorkedDaysDual = pModel_->getDual(minWorkedDaysCons_[i], true);
		double maxWorkedDaysDual = pModel_->getDual(maxWorkedDaysCons_[i], true);
		double maxWorkedWeekendDual = pModel_->getDual(maxWorkedWeekendCons_[i], true);
		updateVarCostInSolver(stabMinWorkedDaysPlus_[i], solver, minWorkedDaysDual+margin);
		updateVarCostInSolver(stabMaxWorkedDaysMinus_[i], solver, -maxWorkedDaysDual+margin);
		updateVarCostInSolver(stabMaxWorkedWeekendMinus_[i], solver, -maxWorkedWeekendDual+margin);

		// stabilization variables corresponding to the flow constraints
		for(int k=0; k<pDemand_->nbDays_; ++k) {
			double restFlowDual = pModel_->getDual(restFlowCons_[i][k], true);
			double workFlowDual = pModel_->getDual(workFlowCons_[i][k], true);
			updateVarCostInSolver(stabRestFlowMinus_[i][k], solver, -restFlowDual+margin);
			updateVarCostInSolver(stabRestFlowPlus_[i][k], solver, restFlowDual+margin);
			updateVarCostInSolver(stabWorkFlowMinus_[i][k], solver, -workFlowDual+margin);
			updateVarCostInSolver(stabWorkFlowPlus_[i][k], solver, workFlowDual+margin);
		}
	}

	// stabilization variables corresponding to the cover constraints
	for(int k=0; k<pDemand_->nbDays_; k++){
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				double minDemandDual = pModel_->getDual(minDemandCons_[k][s-1][sk], true);
				double optDemandDual = pModel_->getDual(optDemandCons_[k][s-1][sk], true);
				updateVarCostInSolver(stabMinDemandPlus_[k][s-1][sk], solver, minDemandDual+margin);
				updateVarCostInSolver(stabOptDemandPlus_[k][s-1][sk], solver, optDemandDual+margin);
			}
		}
	}
}

// STAB
// Check the stopping criterion of the relaxation solution specific to the
// the stabilization
// The point is that current solution can be infeasible if  stabilization
// variables are non zero
bool MasterProblem::stabCheckStoppingCriterion() {
	if (!pModel_->getParameters().isStabilization_) {
		return true;
	}
	else {
		for(int i=0; i<pScenario_->nbNurses_; i++){
			// stabilization variables corresponding to the global constraints of
			// of the nurses
			if (pModel_->getVarValue(stabMinWorkedDaysPlus_[i]) > EPSILON ||
				pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]) >EPSILON ||
				pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]) > EPSILON) {
				return false;
			}

			// stabilization variables corresponding to the flow constraints
			for(int k=0; k<pDemand_->nbDays_; ++k) {
				if (pModel_->getVarValue(stabRestFlowMinus_[i][k]) > EPSILON ||
					pModel_->getVarValue(stabRestFlowPlus_[i][k]) > EPSILON ||
					pModel_->getVarValue(stabWorkFlowMinus_[i][k]) > EPSILON ||
					pModel_->getVarValue(stabWorkFlowPlus_[i][k]) > EPSILON ) {
					return false;
				}
			}
		}

		// stabilization variables corresponding to the cover constraints
		for(int k=0; k<pDemand_->nbDays_; k++){
			for(int s=1; s<pScenario_->nbShifts_; s++){
				for(int sk=0; sk<pScenario_->nbSkills_; sk++){
					if (pModel_->getVarValue(stabMinDemandPlus_[k][s-1][sk]) > EPSILON ||
						pModel_->getVarValue(stabOptDemandPlus_[k][s-1][sk]) > EPSILON) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

// STAB: compute the lagrangian bound
//
double MasterProblem::computeLagrangianBound(double objVal,double sumRedCost) {
	double stabSumCostValue = 0.0;
	if (!pModel_->getParameters().isStabilization_) {
		return objVal+sumRedCost;
	}
	else {
		for(int i=0; i<pScenario_->nbNurses_; i++){
			// stabilization variables corresponding to the global constraints of
			// of the nurses
			stabSumCostValue += stabMinWorkedDaysPlus_[i]->getCost()*pModel_->getVarValue(stabMinWorkedDaysPlus_[i]);
			stabSumCostValue += stabMaxWorkedDaysMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]);
			stabSumCostValue += stabMaxWorkedWeekendMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]);

			// stabilization variables corresponding to the flow constraints
			for(int k=0; k<pDemand_->nbDays_; ++k) {
				stabSumCostValue += stabRestFlowPlus_[i][k]->getCost()*pModel_->getVarValue(stabRestFlowPlus_[i][k]);
				stabSumCostValue += stabRestFlowMinus_[i][k]->getCost()*pModel_->getVarValue(stabRestFlowMinus_[i][k]);
				stabSumCostValue += stabWorkFlowPlus_[i][k]->getCost()*pModel_->getVarValue(stabWorkFlowPlus_[i][k]);
				stabSumCostValue += stabWorkFlowMinus_[i][k]->getCost()*pModel_->getVarValue(stabWorkFlowMinus_[i][k]);
			}
		}

		// stabilization variables corresponding to the cover constraints
		for(int k=0; k<pDemand_->nbDays_; k++){
			for(int s=1; s<pScenario_->nbShifts_; s++){
				for(int sk=0; sk<pScenario_->nbSkills_; sk++){
					stabSumCostValue += stabMinDemandPlus_[k][s-1][sk]->getCost()*pModel_->getVarValue(stabMinDemandPlus_[k][s-1][sk]);
					stabSumCostValue += stabOptDemandPlus_[k][s-1][sk]->getCost()*pModel_->getVarValue(stabOptDemandPlus_[k][s-1][sk]);
				}
			}
		}
	}
	return objVal+sumRedCost-stabSumCostValue;
}

// STAB: reset the costs and bounds of the stabilization variables
//
void MasterProblem::stabResetBoundAndCost(OsiSolverInterface* solver, SolverParam param) {
	for(int i=0; i<pScenario_->nbNurses_; i++){
		// stabilization variables corresponding to the global constraints of
		// of the nurses
		updateVarCostInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabCostIni_+param.stabCostMargin_);
		updateVarCostInSolver(stabMaxWorkedDaysMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
		updateVarCostInSolver(stabMaxWorkedWeekendMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
		updateVarUbInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabBoundIni_);
		updateVarUbInSolver(stabMaxWorkedDaysMinus_[i], solver, param.stabBoundIni_);
		updateVarUbInSolver(stabMaxWorkedWeekendMinus_[i], solver, param.stabBoundIni_);

		// stabilization variables corresponding to the flow constraints
		for(int k=0; k<pDemand_->nbDays_; ++k) {
			updateVarCostInSolver(stabRestFlowMinus_[i][k], solver, -param.stabCostIni_+param.stabCostMargin_);
			updateVarCostInSolver(stabRestFlowPlus_[i][k], solver, param.stabCostIni_+param.stabCostMargin_);
			updateVarCostInSolver(stabWorkFlowMinus_[i][k], solver, -param.stabCostIni_+param.stabCostMargin_);
			updateVarCostInSolver(stabWorkFlowPlus_[i][k], solver, param.stabCostIni_+param.stabCostMargin_);
			updateVarUbInSolver(stabRestFlowMinus_[i][k], solver, param.stabBoundIni_);
			updateVarUbInSolver(stabRestFlowPlus_[i][k], solver, param.stabBoundIni_);
			updateVarUbInSolver(stabWorkFlowMinus_[i][k], solver, param.stabBoundIni_);
			updateVarUbInSolver(stabWorkFlowPlus_[i][k], solver, param.stabBoundIni_);
		}
	}

	// stabilization variables corresponding to the cover constraints
	for(int k=0; k<pDemand_->nbDays_; k++){
		for(int s=1; s<pScenario_->nbShifts_; s++){
			for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				updateVarCostInSolver(stabMinDemandPlus_[k][s-1][sk], solver, param.stabCostIni_+param.stabCostMargin_);
				updateVarCostInSolver(stabOptDemandPlus_[k][s-1][sk], solver, param.stabCostIni_+param.stabCostMargin_);
				updateVarUbInSolver(stabMinDemandPlus_[k][s-1][sk], solver, param.stabBoundIni_);
				updateVarUbInSolver(stabOptDemandPlus_[k][s-1][sk], solver, param.stabBoundIni_);
			}
		}
	}
}
