/*
 * MasterProblem.cpp
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#include "MasterProblem.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"

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
using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;


//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

// Default constructor
MasterProblem::MasterProblem(PScenario pScenario, PDemand pDemand,
		PPreferences pPreferences, vector<State>* pInitState, MySolverType solverType):

								   Solver(pScenario, pDemand, pPreferences, pInitState), PrintSolution(),
								   pModel_(0), positionsPerSkill_(pScenario->nbSkills_),
								   skillsPerPosition_(pScenario->nbPositions()), pPricer_(0), pTree_(0), pRule_(0),
								   solverType_(solverType), columnVars_(pScenario->nbNurses_),
								   minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
								   minWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedDaysAvgVars_(pScenario->nbNurses_), maxWorkedWeekendAvgVars_(pScenario_->nbNurses_),
								   minWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedDaysContractAvgVars_(pScenario->nbContracts_), maxWorkedWeekendContractAvgVars_(pScenario_->nbContracts_),
								   optDemandVars_(pDemand_->nbDays_),numberOfNursesByPositionVars_(pDemand_->nbDays_), skillsAllocVars_(pDemand_->nbDays_),

								   minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
								   minWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedDaysAvgCons_(pScenario->nbNurses_), maxWorkedWeekendAvgCons_(pScenario_->nbNurses_),
								   minWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedDaysContractAvgCons_(pScenario->nbContracts_), maxWorkedWeekendContractAvgCons_(pScenario_->nbContracts_),
								   minDemandCons_(pDemand_->nbDays_), optDemandCons_(pDemand_->nbDays_),numberOfNursesByPositionCons_(pDemand_->nbDays_), feasibleSkillsAllocCons_(pDemand_->nbDays_),
									// STAB
									stabMinWorkedDaysPlus_(pScenario->nbNurses_), stabMaxWorkedDaysMinus_(pScenario->nbNurses_), stabMaxWorkedWeekendMinus_(pScenario->nbNurses_),
									stabMinDemandPlus_(pDemand_->nbDays_), stabOptDemandPlus_(pDemand_->nbDays_)
{
	// build the model
	this->initializeSolver(solverType);
}

MasterProblem::~MasterProblem(){
  delete pPricer_;
  delete pRule_;
  delete pTree_;
  delete pModel_;
}

// Initialize the solver at construction
//
void MasterProblem::initializeSolver(MySolverType solverType) {

	// This OSI interface is created only to retrieve the proper value of
	// infinity for the solver
  double inf = -1;
	switch(solverType){
	case S_CLP:
    pModel_ = new BcpModeler(this, PB_NAME, CLP);
    inf = OsiClpSolverInterface().getInfinity();
		break;
	case S_Gurobi:
#ifdef USE_GUROBI
		pModel_ = new BcpModeler(this,PB_NAME, Gurobi);
		inf = OsiGrbSolverInterface().getInfinity();
#else
	  Tools::throwError("BCP has not been built with Gurobi.");
#endif
		break;
	case S_Cplex:
#ifdef USE_CPLEX
		pModel_ = new BcpModeler(this,PB_NAME, Cplex);\
		inf = OsiCpxSolverInterface().getInfinity();
#else
	  Tools::throwError("BCP has not been built with Cplex.");
#endif
		break;
	case S_CBC:
    pModel_ = new BcpModeler(this,PB_NAME);
    inf = OsiClpSolverInterface().getInfinity();
		break;
	default:
		Tools::throwError("MasterProblem::initializeSolver: the requested solver is not supported presently.");
		break;
	}

//  pModel_->setInfinity(inf);

	this->preprocessData();

	/*
	 * Build the two vectors linking positions and skills
	 */
	for(unsigned int p=0; p<skillsPerPosition_.size(); p++){
		vector<int> skills(pScenario_->pPositions()[p]->skills_.size());
		for(unsigned int sk=0; sk<pScenario_->pPositions()[p]->skills_.size(); ++sk)
			skills[sk]=pScenario_->pPositions()[p]->skills_[sk];
		skillsPerPosition_[p] = skills;
	}
	for(unsigned int sk=0; sk<positionsPerSkill_.size(); sk++){
		vector<int> positions(pScenario_->nbPositions());
		int i(0);
		for(unsigned int p=0; p<positions.size(); p++)
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

//solve the rostering problem
double MasterProblem::solve(vector<Roster> solution){
	return  solve(solution, true);
}

// Solve the rostering problem with parameters

double MasterProblem::solve(const SolverParam& param, vector<Roster> solution){
	pModel_->setParameters(param, this);
	return  solve(solution, true);
}

//solve the rostering problem
double MasterProblem::solve(vector<Roster> solution, bool rebuild) {

  // build the model first
  if (rebuild) {
    pModel_->clear();
    this->build(pModel_->getParameters());
  } else
    pModel_->reset();

  // input an initial solution
  this->initializeSolution(solution);

  // DBG
#ifdef  DBG
  pModel_->writeProblem("master_model");
#endif

  // RqJO: warning, it would be better to define an enumerate type of verbosity
  // levels and create the matching in the Modeler subclasses
  if (solverType_ != S_CBC) {
    pModel_->setVerbosity(1);
  }
  this->solveWithCatch();

  if (pModel_->getParameters().printBranchStats_) {
    pModel_->printStats();
  }

  if (!pModel_->printBestSol()) {
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
double MasterProblem::resolve(PDemand pDemand, const SolverParam& param, vector<Roster> solution){
	updateDemand(pDemand);
	pModel_->setParameters(param, this);
	return solve(solution, false);
}

// Initialization of the master problem with/without solution
void MasterProblem::initialize(const SolverParam& param, vector<Roster> solution) {
	this->build(param);
	this->initializeSolution(solution);
	pModel_->setParameters(param, this);

	// in case the initial solution is not empty, fix the corresponding rotations
	// to one and solve the problem to get the solution properly
	if (!solution.empty()) {
		pModel_->fixEveryRotation();
	}
}

//build the rostering problem
void MasterProblem::build(const SolverParam& param){
  /* Min/Max constraints */
  buildMinMaxCons(param);

  /* Skills coverage constraints */
  buildSkillsCoverageCons(param);

  /* Initialize the objects used in the branch and price unless the CBC is used
      to solve the problem
   */
  if (solverType_ != S_CBC) {
    /* Rotation pricer */
    pPricer_ = new RCPricer(this, "pricer", param);
    pModel_->addObjPricer(pPricer_);

    /* Tree */
    RestTree* pTree =new RestTree(pScenario_, pDemand_);
    pTree_ = pTree;
    pModel_->addTree(pTree_);

    /* Branching rule */
    pRule_ = new DiveBranchingRule(this, pTree, "branching rule");
    pModel_->addBranchingRule(pRule_);
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
    for (int day = 0; day < pDemand_->nbDays_; day++)
      isRelaxDay_[day] = isRelaxDay_[day] ? !isUnrelax[day] : false;
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
	for (PLiveNurse pNurse: theLiveNurses_){
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
	for (PLiveNurse pNurse: theLiveNurses_){
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
double MasterProblem::rollingSolve(const SolverParam& param, int firstDay) {

	// build the model and initialize with artificial columns at the first iteration
	if(firstDay == 0) {
		initialize(param);
	}
	// otherwise reset the solution but keep the best solution
	else {
	  // load and store the best solution
	  pModel_->loadBestSol();
	  storeSolution();
	  // reset the model
		pModel_->reset();
		pModel_->setParameters(param, this);
    // add the best solution  back in the model
		initializeSolution(solution_);
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
double MasterProblem::LNSSolve(const SolverParam& param) {
  // load and store the best solution
  pModel_->loadBestSol();
  storeSolution();
  // reset the model
  pModel_->reset();
  pModel_->setParameters(param, this);
  // add the best solution  back in the model
  initializeSolution(solution_);

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
	vector4D<double> skillsAllocation;
	Tools::initVector4D(skillsAllocation, pDemand_->nbDays_, pScenario_->nbShifts_-1,
	    pScenario_->nbSkills_, 0, .0);

	for(int k=0; k<pDemand_->nbDays_; ++k)
	  for(int s=0; s<pScenario_->nbShifts_-1; ++s)
	    for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
				skillsAllocation[k][s][sk] = pModel_->getVarValues(skillsAllocVars_[k][s][sk]);

	//build the rosters
	for(PLiveNurse pNurse: theLiveNurses_)
		pNurse->roster_.reset();

	for(MyVar* var: pModel_->getActiveColumns()){
		if(pModel_->getVarValue(var) > EPSILON){
		  PPattern pat = getPattern(var->getPattern());
			PLiveNurse pNurse = theLiveNurses_[pat->nurseId_];
			for(int k=pat->firstDay_; k<pat->firstDay_+pat->length_; ++k){
        int s = pat->getShift(k);
        if(s == 0) continue; // nothing to do
        // assign a skill to the nurse for the shift
				bool assigned = false;
				for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
					if(skillsAllocation[k][s-1][sk][pNurse->pPosition_->id_] > EPSILON){
						pNurse->roster_.assignTask(k,s,sk);
						skillsAllocation[k][s-1][sk][pNurse->pPosition_->id_] --;
						assigned = true;
						break;
					}
				if(!assigned){
					char error[255];
					sprintf(error, "No skill found for Nurse %d on day %d on shift %d", pNurse->id_, k, s);
					Tools::throwError((const char*) error);
				}
			}
		}
	}

	//build the states of each nurse
	solution_.clear();
	for(PLiveNurse pNurse: theLiveNurses_){
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
vector3D<double> MasterProblem::getFractionalRoster() {
	vector3D<double> fractionalRoster;
  Tools::initVector3D(fractionalRoster, getNbNurses(), getNbDays(), pDemand_->nbShifts_, .0);

	// Retrieve current fractional roster for each nurse
	for(MyVar* var : pModel_->getActiveColumns()){
		if (var->getPattern().empty()) continue;
    double value = pModel_->getVarValue(var);
		if(value < EPSILON) continue;
		PPattern pat = getPattern(var->getPattern());
    vector2D<double>& fractionalRoster2 = fractionalRoster[pat->nurseId_];
		for(int k=pat->firstDay_; k<pat->firstDay_+pat->length_; ++k)
			fractionalRoster2[k][pat->getShift(k)] += value;
	}

	return fractionalRoster;
}

void MasterProblem::checkIfPatternAlreadyPresent(const std::vector<double>& pattern) const {
  for(MyVar* var: pModel_->getActiveColumns()) {
    bool equal = true;
    for (int j = 0; j < pattern.size(); ++j)
      if (abs(pattern[j] - var->getPattern()[j]) > EPSILON) {
        equal = false;
        break;
      }
    if(equal) {
      string name = var->name_;
      Tools::throwError("Pattern already present as column: " + name);
    }
  }
}

void MasterProblem::printCurrentSol(){
	allocationToString();
		coverageToString();
}

/******************************************************
 * Get the duals values per day for a nurse
 ******************************************************/
// build a DualCosts structure
DualCosts MasterProblem::buildDualCosts(PLiveNurse  pNurse) const {
  return DualCosts(getShiftsDualValues(pNurse), getStartWorkDualValues(pNurse),
      getEndWorkDualValues(pNurse), getWorkedWeekendDualValue(pNurse));
}

vector2D<double> MasterProblem::getShiftsDualValues(PLiveNurse  pNurse) const {
  vector2D<double> dualValues(pDemand_->nbDays_);
  int i = pNurse->id_;
  int p = pNurse->pContract_->id_;

  /* Min/Max constraints */
  double minWorkedDays = pModel_->getDual(minWorkedDaysCons_[i], true);
  double maxWorkedDays = pModel_->getDual(maxWorkedDaysCons_[i], true);

  double minWorkedDaysAvg = isMinWorkedDaysAvgCons_[i] ? pModel_->getDual(minWorkedDaysAvgCons_[i], true):0.0;
  double maxWorkedDaysAvg = isMaxWorkedDaysAvgCons_[i] ? pModel_->getDual(maxWorkedDaysAvgCons_[i], true):0.0;

  double minWorkedDaysContractAvg = isMinWorkedDaysContractAvgCons_[p] ?
                                    pModel_->getDual(minWorkedDaysContractAvgCons_[p], true):0.0;
  double maxWorkedDaysContractAvg = isMaxWorkedDaysContractAvgCons_[p] ?
                                    pModel_->getDual(maxWorkedDaysContractAvgCons_[p], true):0.0;

  for(int k=0; k<pDemand_->nbDays_; ++k){
    //initialize vector
    vector<double> dualValues2(pScenario_->nbShifts_-1);

    for(int s=1; s<pScenario_->nbShifts_; ++s){
      /* Min/Max constraints */
      dualValues2[s-1] = minWorkedDays + minWorkedDaysAvg + minWorkedDaysContractAvg;
      dualValues2[s-1] += maxWorkedDays + maxWorkedDaysAvg + maxWorkedDaysContractAvg;

      // pour ajuster les valeurs duales en fonction des heures travaillees

      dualValues2[s-1] *= pScenario_->timeDurationToWork_[s];

      /* Skills coverage */
      dualValues2[s-1] += pModel_->getDual(
          numberOfNursesByPositionCons_[k][s-1][pNurse->pPosition_->id_], true);
    }

    //store vector
    dualValues[k] = dualValues2;
  }

  return dualValues;
}


vector<double> MasterProblem::getStartWorkDualValues(PLiveNurse pNurse) const {
  return vector<double>(pDemand_->nbDays_);
}

vector<double> MasterProblem::getEndWorkDualValues(PLiveNurse pNurse) const {
  return vector<double>(pDemand_->nbDays_);
}

double MasterProblem::getWorkedWeekendDualValue(PLiveNurse pNurse) const{
  int id = pNurse->id_;
  double dualVal = pModel_->getDual(maxWorkedWeekendCons_[id], true);
  if (isMaxWorkedWeekendAvgCons_[id]) {
    dualVal += pModel_->getDual(maxWorkedWeekendAvgCons_[id], true);
  }
  if (isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_]) {
    dualVal += pModel_->getDual(
        maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
  }

  return dualVal;
}

/*
 * Min/Max constraints
 */
void MasterProblem::buildMinMaxCons(const SolverParam& param){
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
      pModel_->createLEConsLinear(&maxWorkedWeekendAvgCons_[i], name,
                                  maxTotalWeekendsAvg_[i] - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
                                  varsAvg3, coeffsAvg3);

			isMaxWorkedWeekendAvgCons_[i] = true;

		}

	}

	for(int p=0; p<pScenario_->nbContracts_; ++p){

		if(!minTotalShiftsContractAvg_.empty() && !maxTotalShiftsContractAvg_.empty()  && !weightTotalShiftsContractAvg_.empty()){
			sprintf(name, "minWorkedDaysContractAvgVar_P%d", p);
			pModel_->createPositiveVar(&minWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);
			sprintf(name, "maxWorkedDaysContractAvgVar_P%d", p);
			pModel_->createPositiveVar(&maxWorkedDaysContractAvgVars_[p], name, weightTotalShiftsContractAvg_[p]);

			sprintf(name, "minWorkedDaysContractAvgCons_P%d", p);
			vector<MyVar*> vars1 = {minWorkedDaysContractAvgVars_[p]};
			vector<double> coeffs1 = {1};
      pModel_->createGEConsLinear(&minWorkedDaysContractAvgCons_[p], name, minTotalShiftsContractAvg_[p], vars1,
                                  coeffs1);

			sprintf(name, "maxWorkedDaysContractAvgCons_P%d", p);
			vector<MyVar*> vars2 = {maxWorkedDaysContractAvgVars_[p]};
			vector<double> coeffs2 = {-1};
      pModel_->createLEConsLinear(&maxWorkedDaysContractAvgCons_[p], name, maxTotalShiftsContractAvg_[p], vars2,
                                  coeffs2);

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
void MasterProblem::buildSkillsCoverageCons(const SolverParam& param){
	char name[255];
  //initialize vectors
  Tools::initVector3D(optDemandVars_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
      pScenario_->nbSkills_, (MyVar*)nullptr);
  Tools::initVector3D(stabMinDemandPlus_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
                      pScenario_->nbSkills_, (MyVar*)nullptr);
  Tools::initVector3D(stabOptDemandPlus_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
                      pScenario_->nbSkills_, (MyVar*)nullptr);
  Tools::initVector3D(numberOfNursesByPositionVars_, pDemand_->nbDays_,
      pScenario_->nbShifts_-1, pScenario_->nbPositions(), (MyVar*)nullptr);
  Tools::initVector4D(skillsAllocVars_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
      pScenario_->nbSkills_, pScenario_->nbPositions(), (MyVar*)nullptr);
  Tools::initVector3D(minDemandCons_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
                      pScenario_->nbSkills_, (MyCons*)nullptr);
  Tools::initVector3D(optDemandCons_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
                      pScenario_->nbSkills_, (MyCons*)nullptr);
  Tools::initVector3D(numberOfNursesByPositionCons_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
      pScenario_->nbPositions(), (MyCons*)nullptr);
  Tools::initVector3D(feasibleSkillsAllocCons_, pDemand_->nbDays_, pScenario_->nbShifts_-1,
      pScenario_->nbPositions(), (MyCons*)nullptr);

  for(int k=0; k<pDemand_->nbDays_; k++){
		//forget s=0, it's a resting shift
		for(int s=1; s<pScenario_->nbShifts_; s++){
		  for(int sk=0; sk<pScenario_->nbSkills_; sk++){
				//create variables
				sprintf(name, "optDemandVar_%d_%d_%d", k, s, sk);
				pModel_->createPositiveVar(&optDemandVars_[k][s-1][sk], name, WEIGHT_OPTIMAL_DEMAND);
				for(int p=0; p<pScenario_->nbPositions(); p++){
					sprintf(name, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
					// DBG
					pModel_->createPositiveVar(&skillsAllocVars_[k][s-1][sk][p], name, 0);
					//pModel_->createIntVar(&skillsAllocVars3[p], name, 0);
				}
				//adding variables and building minimum demand constraints
				vector<MyVar*> vars1(positionsPerSkill_[sk].size());
				vector<double> coeffs1(positionsPerSkill_[sk].size());
				for(unsigned int p=0; p<positionsPerSkill_[sk].size(); ++p){
					vars1[p] = skillsAllocVars_[k][s-1][sk][positionsPerSkill_[sk][p]];
					coeffs1[p] = 1;
				}
				sprintf(name, "minDemandCons_%d_%d_%d", k, s, sk);
        pModel_->createGEConsLinear(&minDemandCons_[k][s - 1][sk], name, pDemand_->minDemand_[k][s][sk],
                                    vars1, coeffs1);

				// STAB:Add stabilization variable
				if (param.isStabilization_) {
					sprintf(name,"stabMinDemandPlus_%d_%d_%d",k,s,sk);
					pModel_->createPositiveVar(&stabMinDemandPlus_[k][s-1][sk],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
					pModel_->addCoefLinear(minDemandCons_[k][s-1][sk],stabMinDemandPlus_[k][s-1][sk],1.0);
				}

				//adding variables and building optimal demand constraints
				vars1.push_back(optDemandVars_[k][s-1][sk]);
				coeffs1.push_back(1);
				sprintf(name, "optDemandCons_%d_%d_%d", k, s, sk);
        pModel_->createGEConsLinear(&optDemandCons_[k][s - 1][sk], name, pDemand_->optDemand_[k][s][sk],
                                    vars1, coeffs1);

				// STAB:Add stabilization variable
				if (param.isStabilization_) {
					sprintf(name,"stabOptDemandPlus_%d_%d_%d",k,s,sk);
					pModel_->createPositiveVar(&stabOptDemandPlus_[k][s-1][sk],name,param.stabCostIni_+param.stabCostMargin_,DEFAULT_PATTERN,0,param.stabBoundIni_);
					pModel_->addCoefLinear(optDemandCons_[k][s-1][sk],stabOptDemandPlus_[k][s-1][sk],1.0);
				}
			}

			for(int p=0; p<pScenario_->nbPositions(); p++){
				//creating variables
				sprintf(name, "nursesNumber_%d_%d_%d", k, s, p);
				// DBG
				// pModel_->createIntVar(&numberOfNursesByPositionVars2[p], name, 0);
				pModel_->createPositiveVar(&numberOfNursesByPositionVars_[k][s-1][p], name, 0);
				//adding variables and building number of nurses constraints
				vector<MyVar*> vars3;
				vector<double> coeff3;
				vars3.push_back(numberOfNursesByPositionVars_[k][s-1][p]);
				coeff3.push_back(-1);
				sprintf(name, "nursesNumberCons_%d_%d_%d", k, s, p);
        pModel_->createEQConsLinear(&numberOfNursesByPositionCons_[k][s - 1][p], name, 0,
                                    vars3, coeff3);

				//adding variables and building skills allocation constraints
				int const nonZeroVars4(1+skillsPerPosition_[p].size());
				vector<MyVar*> vars4(nonZeroVars4);
				vector<double> coeff4(nonZeroVars4);
				vars4[0] = numberOfNursesByPositionVars_[k][s-1][p];
				coeff4[0] = 1;
				for(int sk=1; sk<nonZeroVars4; ++sk){
					vars4[sk] = skillsAllocVars_[k][s-1][skillsPerPosition_[p][sk-1]][p];
					coeff4[sk] =-1;
				}
				sprintf(name, "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
        pModel_->createEQConsLinear(&feasibleSkillsAllocCons_[k][s - 1][p], name, 0,
                                    vars4, coeff4);
			}
		}
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

void MasterProblem::updateDemand(PDemand pDemand){
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
	std::stringstream rep;

	char buffer[100];
	sprintf(buffer, "%-30s %10.0f \n", "Column costs", getColumnsCost(TOTAL_COST, false));
	rep << buffer;
	rep << "-----------------------------------------\n";
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. shifts costs", getColumnsCost(CONS_SHIFTS_COST, false));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. worked days costs", getColumnsCost(CONS_WORKED_DAYS_COST, false));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Complete weekend costs", getColumnsCost(COMPLETE_WEEKEND_COST, false));
	rep << buffer;
	sprintf(buffer, "%5s%-25s %10.0f \n", "", "Preferences costs", getColumnsCost(PREFERENCE_COST, false));
	rep << buffer;
  sprintf(buffer, "%5s%-25s %10.0f \n", "", "Resting costs", getColumnsCost(REST_COST, false));
  rep << buffer;
  sprintf(buffer, "%5s%-25s %10.0f \n", "", "History work costs (counted)", getColumnsCost(TOTAL_COST, true));
  rep << buffer;
  sprintf(buffer, "%5s%-25s %10.0f \n", "", "History rest costs (counted)", getColumnsCost(REST_COST,  true));
  rep << buffer;
	rep << "-----------------------------------------\n";
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


string MasterProblem::allocationToString(bool printInteger){
	std::stringstream rep;

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

  vector3D<double> fractionalRoster = getFractionalRoster();
	for (int n = 0; n < nbNurses; n ++) {
		PLiveNurse pNurse = theLiveNurses_[n];
		const vector2D<double>& fnurseFractionalRoster = fractionalRoster[n];
		rep << pNurse->name_ << "\t";
		for(int s=1; s<nbShifts; ++s){
			if (s>1) rep << "\t";
			rep << pScenario_->intToShift_[s] << "\t";
			for (int day = firstDay; day < firstDay+nbDays; day++){
				double shiftValue = fnurseFractionalRoster[day][s];
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
	std::stringstream rep;

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
		rep << pScenario_->intToShift_[s] << "\t";
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
bool MasterProblem::stabCheckStoppingCriterion() const {
	if (!pModel_->getParameters().isStabilization_)
		return true;

  for(int i=0; i<pScenario_->nbNurses_; i++)
    // stabilization variables corresponding to the global constraints of
    // of the nurses
    if (pModel_->getVarValue(stabMinWorkedDaysPlus_[i]) > EPSILON ||
      pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]) >EPSILON ||
      pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]) > EPSILON)
      return false;

  // stabilization variables corresponding to the cover constraints
  for(int k=0; k<pDemand_->nbDays_; k++)
    for(int s=1; s<pScenario_->nbShifts_; s++)
      for(int sk=0; sk<pScenario_->nbSkills_; sk++)
        if (pModel_->getVarValue(stabMinDemandPlus_[k][s-1][sk]) > EPSILON ||
          pModel_->getVarValue(stabOptDemandPlus_[k][s-1][sk]) > EPSILON)
          return false;

	return true;
}

// STAB: compute the lagrangian bound
//
double MasterProblem::computeLagrangianBound(double objVal,double sumRedCost) const {
  return -LARGE_SCORE;
//	double stabSumCostValue = 0.0;
//	if (!pModel_->getParameters().isStabilization_) {
//		return objVal+sumRedCost;
//	}
//	else {
//		for(int i=0; i<pScenario_->nbNurses_; i++){
//			// stabilization variables corresponding to the global constraints of
//			// of the nurses
//			stabSumCostValue += stabMinWorkedDaysPlus_[i]->getCost()*pModel_->getVarValue(stabMinWorkedDaysPlus_[i]);
//			stabSumCostValue += stabMaxWorkedDaysMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedDaysMinus_[i]);
//			stabSumCostValue += stabMaxWorkedWeekendMinus_[i]->getCost()*pModel_->getVarValue(stabMaxWorkedWeekendMinus_[i]);
//		}
//
//		// stabilization variables corresponding to the cover constraints
//		for(int k=0; k<pDemand_->nbDays_; k++)
//			for(int s=1; s<pScenario_->nbShifts_; s++)
//				for(int sk=0; sk<pScenario_->nbSkills_; sk++){
//					stabSumCostValue += stabMinDemandPlus_[k][s-1][sk]->getCost()*pModel_->getVarValue(stabMinDemandPlus_[k][s-1][sk]);
//					stabSumCostValue += stabOptDemandPlus_[k][s-1][sk]->getCost()*pModel_->getVarValue(stabOptDemandPlus_[k][s-1][sk]);
//				}
//	}
//	return objVal+sumRedCost-stabSumCostValue;
}

// STAB: reset the costs and bounds of the stabilization variables
//
void MasterProblem::stabResetBoundAndCost(OsiSolverInterface* solver, const SolverParam& param) {
	for(int i=0; i<pScenario_->nbNurses_; i++){
		// stabilization variables corresponding to the global constraints of
		// of the nurses
		updateVarCostInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabCostIni_+param.stabCostMargin_);
		updateVarCostInSolver(stabMaxWorkedDaysMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
		updateVarCostInSolver(stabMaxWorkedWeekendMinus_[i], solver, -param.stabCostIni_+param.stabCostMargin_);
		updateVarUbInSolver(stabMinWorkedDaysPlus_[i], solver, param.stabBoundIni_);
		updateVarUbInSolver(stabMaxWorkedDaysMinus_[i], solver, param.stabBoundIni_);
		updateVarUbInSolver(stabMaxWorkedWeekendMinus_[i], solver, param.stabBoundIni_);
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
