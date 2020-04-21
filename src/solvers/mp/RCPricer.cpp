/*
 * RCPricer.cpp
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#include "solvers/mp/RCPricer.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/rcspp/ShortSP.h"


/* namespace usage */
using namespace std;

//////////////////////////////////////////////////////////////
//
// R C   P R I C E R
//
//////////////////////////////////////////////////////////////

/* Constructs the pricer object. */
RCPricer::RCPricer(MasterProblem* master, const char* name, const SolverParam& param):
  MyPricer(name), pMaster_(master), pScenario_(master->getScenario()), nbDays_(master->getNbDays()),
  pModel_(master->getModel()), nursesToSolve_(master->getSortedLiveNurses()),
  nbMaxColumnsToAdd_(20), nbSubProblemsToSolve_(15), nb_int_solutions_(0)
{
	// Initialize the parameters
	initPricerParameters(param);

	/* sort the nurses */
	//   random_shuffle( nursesToSolve_.begin(), nursesToSolve_.end());
}

/* Destructs the pricer object. */
RCPricer::~RCPricer() {
	for(pair<const Contract*, SubProblem*> p: subProblems_)
		delete p.second;
}

void RCPricer::initPricerParameters(const SolverParam& param){
	// Here: doing a little check to be sure that secondchance is activated only when the parameters are different...
	withSecondchance_ = param.sp_withsecondchance_ && (param.sp_default_strategy_ != param.sp_secondchance_strategy_);
	nbMaxColumnsToAdd_ = param.sp_nbrotationspernurse_;
	nbSubProblemsToSolve_ = param.sp_nbnursestoprice_;
	defaultSubprobemStrategy_ = param.sp_default_strategy_;
	secondchanceSubproblemStrategy_ = param.sp_secondchance_strategy_;
    shortSubproblem_ = param.sp_short_;
	currentSubproblemStrategy_ = defaultSubprobemStrategy_;
}

/******************************************************
 * Perform pricing
 ******************************************************/
vector<MyVar*> RCPricer::pricing(double bound, bool before_fathom) {

	// Reset all rotations, columns, counters, etc.
	resetSolutions();


	// count and store the nurses whose subproblems produced rotations.
	// DBG: why minDualCost? Isn't it more a reduced cost?
	double minDualCost = 0;
	vector<LiveNurse*> nursesSolved;

	for(auto it0 = nursesToSolve_.begin(); it0 != nursesToSolve_.end();){

		// RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
		LiveNurse* pNurse = *it0;
		bool nurseForbidden = isNurseForbidden(pNurse->id_);

		// cout << "NURSE # " << pNurse->id_ << "    " << pNurse->name_ << endl;

		// IF THE NURSE IS NOT FORBIDDEN, SOLVE THE SUBPROBLEM
		if(!nurseForbidden){

			// BUILD OR RE-USE THE SUBPROBLEM
			SubProblem* subProblem = retriveSubproblem(pNurse);

			// RETRIEVE DUAL VALUES
			DualCosts dualCosts = pMaster_->buildDualCosts(pNurse);

			// UPDATE FORBIDDEN SHIFTS
			if (pModel_->getParameters().isColumnDisjoint_) {
				addForbiddenShifts();
			}
			set<pair<int,int> > nurseForbiddenShifts(forbiddenShifts_);
			pModel_->addForbiddenShifts(pNurse, nurseForbiddenShifts);

			// SET SOLVING OPTIONS
			SubproblemParam sp_param (currentSubproblemStrategy_,pNurse);

			// SOLVE THE PROBLEM
			++ nbSPTried_;
			subProblem->solve(pNurse, &dualCosts, sp_param, nurseForbiddenShifts, forbiddenStartingDays_, true ,
					bound);

			// RETRIEVE THE GENERATED ROTATIONS
			newSolutionsForNurse_ = subProblem->getSolutions();

			// ADD THE ROTATIONS TO THE MASTER PROBLEM
      addColumnsToMaster(pNurse->id_);
		}

		// CHECK IF THE SUBPROBLEM GENERATED NEW ROTATIONS
		// If yes, store the nures
		if(!newSolutionsForNurse_.empty() && !nurseForbidden){
			++nbSPSolvedWithSuccess_;
			if(newSolutionsForNurse_.front().cost < minDualCost)
				minDualCost = newSolutionsForNurse_.front().cost;

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
	BcpModeler* model = static_cast<BcpModeler*>(pModel_);
	model->setLastNbSubProblemsSolved(nbSPTried_);
	model->setLastMinDualCost(minDualCost);

	if(allNewColumns_.empty())
		print_current_solution_();

	//   std::cout << "# -------  END  ------- Subproblems!" << std::endl;

	//   return optimal;
	return allNewColumns_;
}

/******************************************************
 * add some forbidden shifts
 ******************************************************/
void RCPricer::addForbiddenShifts(){
	//search best rotation
	RCSolution bestSolution;
	double bestDualcost = DBL_MAX;
	for(const RCSolution &sol: newSolutionsForNurse_)
		if(sol.cost < bestDualcost){
			bestDualcost = sol.cost;
      bestSolution = sol;
		}

	//forbid shifts of the best rotation
	if(bestDualcost != DBL_MAX) {
	  int k = bestSolution.firstDay;
    for (int s: bestSolution.shifts)
      forbiddenShifts_.insert(pair<int,int>(k++,s));
  }
}

// Returns a pointer to the right subproblem
SubProblem* RCPricer::retriveSubproblem(LiveNurse* pNurse){
	SubProblem* subProblem;
	auto it = subProblems_.find(pNurse->pContract_);
	// Each contract has one subproblem. If it has not already been created, create it.
	if( it == subProblems_.end() ){
	  if (shortSubproblem_)
	    subProblem = new ShortSP(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitialStates());
	  else
	    subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitialStates());
	  // then build the rcspp
	  subProblem->build();
	    
	  //	  subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitState_, noShort);
		subProblems_.insert(it, pair<const Contract*, SubProblem*>(pNurse->pContract_, subProblem));
	} else {
		subProblem = it->second;
	}
	return subProblem;
}

// Add the rotations to the master problem
int RCPricer::addColumnsToMaster(int nurseId){
	// SORT THE SOLUTIONS
  sortNewlyGeneratedSolutions();

	// SECOND, ADD THE ROTATIONS TO THE MASTER PROBLEM (in the previously computed order)
	int nbcolumnsAdded = 0;
	for(const RCSolution& sol: newSolutionsForNurse_){
		allNewColumns_.emplace_back(pMaster_->addColumn(nurseId, sol));
		++nbcolumnsAdded;
		if(nbcolumnsAdded >= nbMaxColumnsToAdd_)
			break;
	}

  return nbcolumnsAdded;
}

// Sort the rotations that just were generated for a nurse. Default option is sort by increasing reduced cost but we
// could try something else (involving disjoint columns for ex.)
void RCPricer::sortNewlyGeneratedSolutions(){
	std::stable_sort(newSolutionsForNurse_.begin(), newSolutionsForNurse_.end(),
	    [](const RCSolution& sol1, const RCSolution& sol2) { return sol1.cost < sol2.cost; });

}

// Set the subproblem options depending on the parameters
//
//void RCPricer::setSubproblemOptions(vector<SolveOption>& options, int& maxRotationLengthForSubproblem,
//		pLiveNurse pNurse){
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

void RCPricer::print_current_solution_(){
	//pMaster_->costsConstrainstsToString();
	//pMaster_->allocationToString();
	if(nb_int_solutions_ < pMaster_->getModel()->nbSolutions()){
		nb_int_solutions_ = pMaster_->getModel()->nbSolutions();
		//pMaster_->getModel()->printBestSol();
		//pMaster_->costsConstrainstsToString();
	}
}
void RCPricer::printStatSPSolutions(){
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

//void RCPricer::recordSPStats(SubProblem* sp){
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

void RCPricer::generateRandomForbiddenStartingDays(){
	set<int> randomForbiddenStartingDays;
	for(int m=0; m<5; m++){
		int k = Tools::randomInt(0, nbDays_-1);
		randomForbiddenStartingDays.insert(k);
	}
	forbiddenStartingDays_ = randomForbiddenStartingDays;
}
