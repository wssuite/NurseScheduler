/*
 * RCPricer.cpp
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#include "solvers/mp/RCPricer.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/rcspp/SubProblemShort.h"


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
	for(auto& p: subProblems_)
		delete p.second;
}

void RCPricer::initPricerParameters(const SolverParam& param){
	nbMaxColumnsToAdd_ = param.sp_nbrotationspernurse_;
	nbSubProblemsToSolve_ = param.sp_nbnursestoprice_;
	defaultSubprobemStrategy_ = param.sp_default_strategy_;
	shortSubproblem_ = param.sp_short_;
	Tools::initVector(currentSubproblemStrategy_, pMaster_->getNbNurses(), defaultSubprobemStrategy_);
}

/******************************************************
 * Perform pricing
 ******************************************************/
vector<MyVar*> RCPricer::pricing(double bound, bool before_fathom, bool after_fathom, bool backtracked) {
  // reset the current strategies at the beginning of a node
  if(after_fathom) // first pricing for a new node
    Tools::initVector(currentSubproblemStrategy_, pMaster_->getNbNurses(), defaultSubprobemStrategy_);

	// Reset all rotations, columns, counters, etc.
	resetSolutions();


	// count and store the nurses whose subproblems produced rotations.
	// DBG: why minDualCost? Isn't it more a reduced cost?
	double minDualCost = 0;
	vector<PLiveNurse> nursesSolved, nursesIncreasedStrategy;

	for(auto it0 = nursesToSolve_.begin(); it0 != nursesToSolve_.end();){

		// RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
		PLiveNurse pNurse = *it0;

    // try next nurse if forbidden
		if(isNurseForbidden(pNurse->id_)) {
      ++it0;
      continue;
		}

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
    SubproblemParam sp_param(currentSubproblemStrategy_[pNurse->id_], pNurse);

    // SOLVE THE PROBLEM
    ++ nbSPTried_;
    subProblem->solve(pNurse, &dualCosts, sp_param, nurseForbiddenShifts, forbiddenStartingDays_, true ,
        bound);

    // RETRIEVE THE GENERATED ROTATIONS
    newSolutionsForNurse_ = subProblem->getSolutions();

    // ADD THE ROTATIONS TO THE MASTER PROBLEM
    addColumnsToMaster(pNurse->id_);

		// CHECK IF THE SUBPROBLEM GENERATED NEW ROTATIONS
		// If yes, store the nures
		if(!newSolutionsForNurse_.empty()){
			++nbSPSolvedWithSuccess_;
			if(newSolutionsForNurse_.front().cost < minDualCost)
				minDualCost = newSolutionsForNurse_.front().cost;

			// erase the  current nurse (and thus try next one in the loop)
			nursesToSolve_.erase(it0);
			nursesSolved.push_back(pNurse);

      //if the maximum number of subproblem solved is reached, break.
      if(nbSPSolvedWithSuccess_ == nbSubProblemsToSolve_)
        break;
      // otherwise continue
      else
        continue;
		}
		// Otherwise (no solution generated),
    // try next nurse
    // if not the most exaustive, increase the strategy for next round
    else if (currentSubproblemStrategy_[pNurse->id_]  <  SubproblemParam::maxSubproblemStrategyLevel_) {
      currentSubproblemStrategy_[pNurse->id_]++;
      nursesIncreasedStrategy.push_back(pNurse);
      // try next nurse
      nursesToSolve_.erase(it0);
    }
    // just try next nurse
    else {
      ++it0;
    }

    // If it was the last nurse to search AND no improving column was found AND we have increased strategy level
    // try to solve for these nurses
    if( it0 == nursesToSolve_.end() && allNewColumns_.empty() && !nursesIncreasedStrategy.empty()){
      // add the nurses left at the beginning of the nursesSolved vector for next loop
      nursesSolved.insert(nursesSolved.begin(), nursesToSolve_.begin(), nursesToSolve_.end());
      // then update the nurse to solve and restart loop
      nursesToSolve_ = nursesIncreasedStrategy;
      it0 = nursesToSolve_.begin();
      // remove all the nurses that have just been added back
      nursesIncreasedStrategy.clear();
    }
	}

  //Add the nurse in nursesIncreasedStrategy at the beginning (first to be tried next round)
  nursesToSolve_.insert(nursesToSolve_.begin(), nursesIncreasedStrategy.begin(), nursesIncreasedStrategy.end());

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
SubProblem* RCPricer::retriveSubproblem(PLiveNurse pNurse){
	SubProblem* subProblem;
	auto it = subProblems_.find(pNurse->pContract_);
	// Each contract has one subproblem. If it has not already been created, create it.
	if( it == subProblems_.end() ){
	  if (shortSubproblem_)
	    subProblem = new SubProblemShort(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitialStates());
	  else
	    subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitialStates());
	  // then build the rcspp
	  subProblem->build();
	    
	  //	  subProblem = new SubProblem(pScenario_, nbDays_, pNurse->pContract_, pMaster_->pInitState_, noShort);
		subProblems_[pNurse->pContract_] = subProblem;
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
