/*
 * RotationPricer.cpp
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#include "RotationPricer.h"
#include "BcpModeler.h"

/* namespace usage */
using namespace std;

//////////////////////////////////////////////////////////////
//
// R O T A T I O N   P R I C E R
//
//////////////////////////////////////////////////////////////


static char const * baseName = "rotation";

/* Constructs the pricer object. */
RotationPricer::RotationPricer(MasterProblem* master, const char* name, SolverParam param):
  MyPricer(name), pMaster_(master), pScenario_(master->pScenario_), nbDays_(master->pDemand_->nbDays_), pModel_(master->getModel()), nursesToSolve_(master->theNursesSorted_), nbMaxRotationsToAdd_(20), nbSubProblemsToSolve_(15), nb_int_solutions_(0)
{
	// Initialize the parameters
	initPricerParameters(param);

	/* sort the nurses */
	//   random_shuffle( nursesToSolve_.begin(), nursesToSolve_.end());
}

/* Destructs the pricer object. */
RotationPricer::~RotationPricer() {
	for(pair<const Contract*, SubProblem*> p: subProblems_)
		delete p.second;
}

void RotationPricer::initPricerParameters(SolverParam param){
	// Here: doing a little check to be sure that secondchance is activated only when the parameters are different...
	withSecondchance_ = param.sp_withsecondchance_ && (param.sp_default_strategy_ != param.sp_secondchance_strategy_);
	nbMaxRotationsToAdd_ = param.sp_nbrotationspernurse_;
	nbSubProblemsToSolve_ = param.sp_nbnursestoprice_;
	defaultSubprobemStrategy_ = param.sp_default_strategy_;
	secondchanceSubproblemStrategy_ = param.sp_secondchance_strategy_;
    shortSubproblem_ = param.sp_short_;
	currentSubproblemStrategy_ = defaultSubprobemStrategy_;
}

/******************************************************
 * Perform pricing
 ******************************************************/
vector<MyVar*> RotationPricer::pricing(double bound, bool before_fathom){

	// Reset all rotations, columns, counters, etc.
	resetSolutions();


	// count and store the nurses whose subproblems produced rotations.
	// DBG: why minDualCost? Isn't it more a reduced cost?
	double minDualCost = 0;
	vector<LiveNurse*> nursesSolved;

	for(vector<LiveNurse*>::iterator it0 = nursesToSolve_.begin(); it0 != nursesToSolve_.end();){

		// RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
		LiveNurse* pNurse = *it0;
		bool nurseForbidden = isNurseForbidden(pNurse->id_);

		// cout << "NURSE # " << pNurse->id_ << "    " << pNurse->name_ << endl;

		// IF THE NURSE IS NOT FORBIDDEN, SOLVE THE SUBPROBLEM
		if(!nurseForbidden){

			// BUILD OR RE-USE THE SUBPROBLEM
			SubProblem* subProblem = retriveSubproblem(pNurse);

			// RETRIEVE DUAL VALUES
			vector< vector<double> > workDualCosts(getWorkDualValues(pNurse));
			vector<double> startWorkDualCosts(getStartWorkDualValues(pNurse));
			vector<double> endWorkDualCosts(getEndWorkDualValues(pNurse));
			double workedWeekendDualCost = getWorkedWeekendDualValue(pNurse);
			DualCosts dualCosts (workDualCosts, startWorkDualCosts, endWorkDualCosts, workedWeekendDualCost, true);

			// UPDATE FORBIDDEN SHIFTS
			if (pModel_->getParameters().isColumnDisjoint_) {
				addForbiddenShifts();
			}
			set<pair<int,int> > nurseForbiddenShifts(forbiddenShifts_);
			pModel_->addForbiddenShifts(pNurse, nurseForbiddenShifts);

			// SET SOLVING OPTIONS
			SubproblemParam sp_param (currentSubproblemStrategy_,pNurse);

			// DBG ***
			// generateRandomForbiddenStartingDays();

			// SOLVE THE PROBLEM
			++ nbSPTried_;
			subProblem->solve(pNurse, &dualCosts, sp_param, nurseForbiddenShifts, forbiddenStartingDays_, true ,
					bound);

			// RETRIEVE THE GENERATED ROTATIONS
			newRotationsForNurse_ = subProblem->getRotations();

			// DBG ***
			// checkForbiddenStartingDays();
			// recordSPStats(subProblem);
			// for(Rotation& rot: newRotationsForNurse_){
			// 	rot.checkDualCost(dualCosts);
			// }

			// ADD THE ROTATIONS TO THE MASTER PROBLEM
			addRotationsToMaster();



		}

		// CHECK IF THE SUBPROBLEM GENERATED NEW ROTATIONS
		// If yes, store the nures
		if(newRotationsForNurse_.size() > 0 && !nurseForbidden){
			++nbSPSolvedWithSuccess_;
			if(newRotationsForNurse_[0].dualCost_ < minDualCost)
				minDualCost = newRotationsForNurse_[0].dualCost_;

			nursesToSolve_.erase(it0);
			nursesSolved.push_back(pNurse);
		}
		// Otherwise (no rotation generated or nurse is forbidden), try the next nurse
		else {
			++it0;

			// If it was the last nurse to search AND no improving column was found AND we may want to solve SP with
			// different parameters -> change these parameters and go for another loop of solving
			if( it0 == nursesToSolve_.end() && allNewColumns_.empty() && withSecondchance_){
				if(currentSubproblemStrategy_ == defaultSubprobemStrategy_){
					nursesToSolve_.insert(nursesToSolve_.end(), nursesSolved.begin(), nursesSolved.end());
					it0 = nursesToSolve_.begin();
					currentSubproblemStrategy_ = secondchanceSubproblemStrategy_;
				} else if (currentSubproblemStrategy_ == secondchanceSubproblemStrategy_) {
					currentSubproblemStrategy_ = defaultSubprobemStrategy_;
				}
			}
		}

		//if the maximum number of subproblem solved is reached, break.
		if(nbSPSolvedWithSuccess_ == nbSubProblemsToSolve_)
			break;
	}

	//Add the nurse in nursesSolved at the end
	nursesToSolve_.insert(nursesToSolve_.end(), nursesSolved.begin(), nursesSolved.end());

	//set statistics
	BcpModeler* model = dynamic_cast<BcpModeler*>(pModel_);
	if(model){
		model->setLastNbSubProblemsSolved(nbSPTried_);
		model->setLastMinDualCost(minDualCost);
	}

	if(allNewColumns_.empty())
		print_current_solution_();

	//   std::cout << "# -------  END  ------- Subproblems!" << std::endl;

	//   return optimal;
	return allNewColumns_;
}

/******************************************************
 * Get the duals values per day for a nurse
 ******************************************************/
vector< vector<double> > RotationPricer::getWorkDualValues(LiveNurse* pNurse){
	vector< vector<double> > dualValues(nbDays_);
	int i = pNurse->id_;
	int p = pNurse->pContract_->id_;

	/* Min/Max constraints */
	double minWorkedDays = pModel_->getDual(pMaster_->minWorkedDaysCons_[i], true);
	double maxWorkedDays = pModel_->getDual(pMaster_->maxWorkedDaysCons_[i], true);

	double minWorkedDaysAvg = pMaster_->isMinWorkedDaysAvgCons_[i] ? pModel_->getDual(pMaster_->minWorkedDaysAvgCons_[i], true):0.0;
	double maxWorkedDaysAvg = pMaster_->isMaxWorkedDaysAvgCons_[i] ? pModel_->getDual(pMaster_->maxWorkedDaysAvgCons_[i], true):0.0;

	double minWorkedDaysContractAvg = pMaster_->isMinWorkedDaysContractAvgCons_[p] ?
			pModel_->getDual(pMaster_->minWorkedDaysContractAvgCons_[p], true):0.0;
	double maxWorkedDaysContractAvg = pMaster_->isMaxWorkedDaysContractAvgCons_[p] ?
			pModel_->getDual(pMaster_->maxWorkedDaysContractAvgCons_[p], true):0.0;

	for(int k=0; k<nbDays_; ++k){
		//initialize vector
		vector<double> dualValues2(pScenario_->nbShifts_-1);

		for(int s=1; s<pScenario_->nbShifts_; ++s){
			/* Min/Max constraints */
			dualValues2[s-1] = minWorkedDays + minWorkedDaysAvg + minWorkedDaysContractAvg;
			dualValues2[s-1] += maxWorkedDays + maxWorkedDaysAvg + maxWorkedDaysContractAvg;

			// pour ajuster les valeurs duales en fonction des heures travaillees
			
			dualValues2[s-1] *= pScenario_->hoursToWork_[s];
			
			/* Skills coverage */
			dualValues2[s-1] += pModel_->getDual(
					pMaster_->numberOfNursesByPositionCons_[k][s-1][pNurse->pPosition_->id_], true);
		}

		//store vector
		dualValues[k] = dualValues2;
	}

	return dualValues;
}


vector<double> RotationPricer::getStartWorkDualValues(LiveNurse* pNurse){
	int i = pNurse->id_;
	vector<double> dualValues(nbDays_);

	//get dual value associated to the source
	dualValues[0] =  pModel_->getDual(pMaster_->restFlowCons_[i][0], true);
	//get dual values associated to the work flow constraints
	//don't take into account the last which is the sink
	for(int k=1; k<nbDays_; ++k)
		dualValues[k] = pModel_->getDual(pMaster_->workFlowCons_[i][k-1], true);

	return dualValues;
}

vector<double> RotationPricer::getEndWorkDualValues(LiveNurse* pNurse){
	int i = pNurse->id_;
	vector<double> dualValues(nbDays_);

	//get dual values associated to the work flow constraints
	//don't take into account the first which is the source
	//take into account the cost, if the last day worked is k
	for(int k=0; k<nbDays_-1; ++k)
		dualValues[k] = -pModel_->getDual(pMaster_->restFlowCons_[i][k+1], true);

	//get dual value associated to the sink
	dualValues[nbDays_-1] =  pModel_->getDual(
			pMaster_->workFlowCons_[i][nbDays_-1], true);

	return dualValues;
}

double RotationPricer::getWorkedWeekendDualValue(LiveNurse* pNurse){
	int id = pNurse->id_;
	double dualVal = pModel_->getDual(pMaster_->maxWorkedWeekendCons_[id], true);
	if (pMaster_->isMaxWorkedWeekendAvgCons_[id]) {
		dualVal += pModel_->getDual(pMaster_->maxWorkedWeekendAvgCons_[id], true);
	}
	if (pMaster_->isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_]) {
		dualVal += pModel_->getDual(pMaster_->maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
	}

	return dualVal;
}

/******************************************************
 * add some forbidden shifts
 ******************************************************/
void RotationPricer::addForbiddenShifts(){
	//search best rotation
	vector<Rotation>::iterator bestRotation;
	double bestDualcost = DBL_MAX;
	for(vector<Rotation>::iterator it = newRotationsForNurse_.begin(); it != newRotationsForNurse_.end(); ++it)
		if(it->dualCost_ < bestDualcost){
			bestDualcost = it->dualCost_;
			bestRotation = it;
		}

	//forbid shifts of the best rotation
	if(bestDualcost != DBL_MAX)
		for(pair<int,int> pair: bestRotation->shifts_)
			forbiddenShifts_.insert(pair);
}

// Returns a pointer to the right subproblem
SubProblem* RotationPricer::retriveSubproblem(LiveNurse* pNurse){
	SubProblem* subProblem;
	map<const Contract*, SubProblem*>::iterator it =  subProblems_.find(pNurse->pContract_);
	// Each contract has one subproblem. If it has not already been created, create it.
	if( it == subProblems_.end() ){
	  if (shortSubproblem_)
	    subProblem = new SubProblemShort(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitState_);
	  else
	    subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitState_);
	  // then build the graph
	  subProblem->build();
	    
	  //	  subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitState_, noShort);
		subProblems_.insert(it, pair<const Contract*, SubProblem*>(pNurse->pContract_, subProblem));
	} else {
		subProblem = it->second;
	}
	return subProblem;
}

// Add the rotations to the master problem
void RotationPricer::addRotationsToMaster(){

    // COMPUTE THE COST OF THE ROTATIONS
	for(Rotation& rot: newRotationsForNurse_){
		rot.computeCost(pScenario_, pMaster_->pPreferences_, pMaster_->theLiveNurses_,nbDays_);
		rot.computeTimeDuration(pScenario_);
		rot.treeLevel_ = pModel_->getCurrentTreeLevel();
	}

	// SORT THE ROTATIONS
	sortNewlyGeneratedRotations();

	// SECOND, ADD THE ROTATIONS TO THE MASTER PROBLEM (in the previously computed order)
	int nbRotationsAdded = 0;
	for(Rotation& rot: newRotationsForNurse_){
		allNewColumns_.push_back(pMaster_->addRotation(rot, baseName));
		++nbRotationsAdded;
		// DBG
//		cout << rot.toString(nbDays_, pScenario_->shiftIDToShiftTypeID_) << endl;
		if(nbRotationsAdded >= nbMaxRotationsToAdd_)
			break;
	}
}

// Sort the rotations that just were generated for a nurse. Default option is sort by increasing reduced cost but we
// could try something else (involving disjoint columns for ex.)
void RotationPricer::sortNewlyGeneratedRotations(){
	std::stable_sort(newRotationsForNurse_.begin(), newRotationsForNurse_.end(), Rotation::compareDualCost);

}

// Set the subproblem options depending on the parameters
//
//void RotationPricer::setSubproblemOptions(vector<SolveOption>& options, int& maxRotationLengthForSubproblem,
//		LiveNurse* pNurse){
//
//	// Default option that should be used
//	options.push_back(SOLVE_ONE_SINK_PER_LAST_DAY);
//
//	// Options are added depending on the chosen parameters
//	//
//	if (currentPricerParam_.isExhaustiveSearch()) {
//		options.push_back(SOLVE_SHORT_ALL);
//		maxRotationLengthForSubproblem = LARGE_TIME;
//	}
//	else if (currentPricerParam_.nonPenalizedRotationsOnly()) {
//		options.push_back(SOLVE_SHORT_DAY_0_ONLY);
//		maxRotationLengthForSubproblem = pNurse->maxConsDaysWork();
//	}
//	else {
//		if(currentPricerParam_.withShortRotations()) options.push_back(SOLVE_SHORT_ALL);
//		maxRotationLengthForSubproblem = currentPricerParam_.maxRotationLength();
//	}
//
//}





// ------------------------------------------
//
// PRINT functions
//
// ------------------------------------------

void RotationPricer::print_current_solution_(){
	//pMaster_->costsConstrainstsToString();
	//pMaster_->allocationToString();
	if(nb_int_solutions_ < pMaster_->getModel()->nbSolutions()){
		nb_int_solutions_ = pMaster_->getModel()->nbSolutions();
		//pMaster_->getModel()->printBestSol();
		//pMaster_->costsConstrainstsToString();
	}
}
void RotationPricer::printStatSPSolutions(){
	double tMeanSubproblems = timeInExSubproblems_ / ((double)nbExSubproblems_);
	double tMeanS = timeForS_ / ((double)nbS_);
	double tMeanNL = timeForNL_ / ((double)nbNL_);
	double tMeanN = timeForN_ / ((double)nbN_);
	string sepLine = "+-----------------+------------------------+-----------+\n";
	printf("\n");
	printf("%s", sepLine.c_str());
	printf("| %-15s |%10s %12s |%10s |\n", "type", "time", "number", "mean time");
	printf("%s", sepLine.c_str());
	printf("| %-15s |%10.2f %12d |%10.4f |\n",	"Ex. Subproblems", timeInExSubproblems_, nbExSubproblems_, tMeanSubproblems);
	printf("| %-15s |%10.2f %12d |%10.4f |\n",	"Short rotations", timeForS_, nbS_, tMeanS);
	printf("| %-15s |%10.2f %12d |%10.4f |\n",	"NL rotations", timeForNL_, nbNL_, tMeanNL);
	printf("%s", sepLine.c_str());
	printf("| %-15s |%10.2f %12d |%10.4f |\n", "N rotations", timeForN_, nbN_, tMeanN);
	printf("%s", sepLine.c_str());
	printf("\n");
}





// ------------------------------------------
//
// DBG functions - may be useful
//
// ------------------------------------------

//void RotationPricer::recordSPStats(SubProblem* sp){
//	if (currentPricerParam_.isExhaustiveSearch()) {
//		nbExSubproblems_++; nbSubproblems_++; nbS_++; nbNL_++;
//		timeForS_ += sp->timeInS_->dSinceStart();
//		timeForNL_ += sp->timeInNL_->dSinceStart();
//		timeInExSubproblems_ += sp->timeInS_->dSinceStart() + sp->timeInNL_->dSinceStart();
//	}
//	else if (currentPricerParam_.nonPenalizedRotationsOnly()) {
//		nbSubproblems_++; nbN_++;
//		timeForN_ += sp->timeInNL_->dSinceStart();
//	}
//	else {
//	}
//}

void RotationPricer::generateRandomForbiddenStartingDays(){
	set<int> randomForbiddenStartingDays;
	for(int m=0; m<5; m++){
		int k = Tools::randomInt(0, nbDays_-1);
		randomForbiddenStartingDays.insert(k);
	}
	forbiddenStartingDays_ = randomForbiddenStartingDays;
}
void RotationPricer::checkForbiddenStartingDays(){
	for(Rotation& rot: newRotationsForNurse_){
		int startingDay = rot.firstDay_;
		if(forbiddenStartingDays_.find(startingDay) != forbiddenStartingDays_.end()){
			cout << "# On a généré une rotation qui commence un jour interdit !" << endl;
			getchar();
		}
	}
}
