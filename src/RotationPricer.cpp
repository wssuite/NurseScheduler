/*
 * RotationPricer.cpp
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#include "RotationPricer.h"

/* namespace usage */
using namespace std;

//////////////////////////////////////////////////////////////
//
// R O T A T I O N   P R I C E R
//
//////////////////////////////////////////////////////////////


static char* baseName = "rotation";

/* Constructs the pricer object. */
RotationPricer::RotationPricer(MasterProblem* master, const char* name):
                        MyPricer(name), nbMaxRotationsToAdd_(20), nbSubProblemsToSolve_(15), nursesToSolve_(master->theLiveNurses_),
                        master_(master), pScenario_(master->pScenario_), pDemand_(master->pDemand_), pModel_(master->getModel())
{
   /* sort the nurses */
//   random_shuffle( nursesToSolve_.begin(), nursesToSolve_.end());
}

/* Destructs the pricer object. */
RotationPricer::~RotationPricer() {
   for(pair<const Contract*, SubProblem*> p: subProblems_)
      delete p.second;
}

/******************************************************
 * Perform pricing
 ******************************************************/
bool RotationPricer::pricing(double bound){
   //=false if once optimality hasn't be proven
   bool optimal = true;
   //forbidden shifts
   set<pair<int,int> > forbiddenShifts;
   //computed new rotations
   vector<Rotation> rotations;

   std::cout << "# ------- BEGIN ------- Subproblems..." << std::endl;

   //count and store the nurses for whom their subproblem has generated rotations.
   int nbSubProblemSolved = 0;
   vector<LiveNurse*> nursesSolved;
   for(vector<LiveNurse*>::iterator it0 = nursesToSolve_.begin(); it0 != nursesToSolve_.end();){
      LiveNurse* pNurse = *it0;

      /* Build or re-use a subproblem */
      SubProblem* subProblem;
      //search the contract
      map<const Contract*, SubProblem*>::iterator it =  subProblems_.find(pNurse->pContract_);
      //if doesn't find => create new subproblem
      if( it == subProblems_.end() ){
         subProblem = new SubProblem(pScenario_, pDemand_, pNurse->pContract_, master_->pInitState_);
         subProblems_.insert(it, pair<const Contract*, SubProblem*>(pNurse->pContract_, subProblem));
      }

      //otherwise retrieve the subproblem associated to the contract
      else
         subProblem = it->second;

      /* Retrieves dual values */
      vector< vector<double> > workDualCosts; workDualCosts = getWorkDualValues(pNurse);
      vector<double> startWorkDualCosts; startWorkDualCosts = getStartWorkDualValues(pNurse);
      vector<double> endWorkDualCosts; endWorkDualCosts = getEndWorkDualValues(pNurse);
      double workedWeekendDualCost = getWorkedWeekendDualValue(pNurse);

      Costs costs (&workDualCosts, &startWorkDualCosts, &endWorkDualCosts, workedWeekendDualCost);

      /* Compute forbidden */
      computeForbiddenShifts(forbiddenShifts, rotations);

	   /* Solve options */
	   vector<SolveOption> options;
	   //options.push_back(SOLVE_ONE_SINK_PER_FIRST_DAY);
	   options.push_back(SOLVE_ONE_SINK_PER_LAST_DAY);
	   options.push_back(SOLVE_FORBIDDEN_RESET);

//	   cout << "#  SP " << pNurse->name_ << " begins" << endl;

	   /* Solve subproblems */
	   if( subProblem->solve(pNurse, &costs, options, forbiddenShifts, false, 10) )
		   optimal = false;
	   else
		   subProblem->solve(pNurse, &costs, options, forbiddenShifts, true);

//	   cout << "#  SP " << pNurse->name_ << " solved" << endl;

		/* Retrieve rotations */
		rotations = subProblem->getRotations();
		std::sort(rotations.begin(), rotations.end(), Rotation::compareDualCost);
		// add them to the master problem
		int nbRotationsAdded = 0;
		for(Rotation rot: rotations){
			double c = rot.cost_;
			rot.computeCost(pScenario_, master_->pPreferences_, master_->pDemand_->nbDays_);
			rot.computeDualCost(workDualCosts, startWorkDualCosts, endWorkDualCosts, workedWeekendDualCost);
			master_->addRotation(rot, baseName);
			++nbRotationsAdded;
			if(nbRotationsAdded > nbMaxRotationsToAdd_)
			   break;
		}

//		cout << "#  SP " << pNurse->name_ << " added columns" << endl;

      //count if the subproblem has generated some rotations and then store the nurse
      if(rotations.size() > 0){
         ++nbSubProblemSolved;
         nursesToSolve_.erase(it0);
         nursesSolved.push_back(pNurse);
      }
      //try the next nurse
      else {
         ++it0;
      }
		//if the maximum number of subproblem solved is reached, break.
		if(nbSubProblemSolved == nbSubProblemsToSolve_)
		   break;
   }
   
   //Add the nurse in nursesSolved at the end
   nursesToSolve_.insert(nursesToSolve_.end(), nursesSolved.begin(), nursesSolved.end());

   std::cout << "# -------  END  ------- Subproblems!" << std::endl;

   return optimal;
}

/******************************************************
 * Get the duals values per day for a nurse
 ******************************************************/
vector< vector<double> > RotationPricer::getWorkDualValues(LiveNurse* pNurse){
   vector< vector<double> > dualValues(pDemand_->nbDays_);
   int i = pNurse->id_;

   /* Min/Max constraints */
   double minWorkedDays = pModel_->getDual(master_->minWorkedDaysCons_[i], true);
   double maxWorkedDays = pModel_->getDual(master_->maxWorkedDaysCons_[i], true);

   for(int k=0; k<pDemand_->nbDays_; ++k){
      //initialize vector
      vector<double> dualValues2(pScenario_->nbShifts_-1);

      for(int s=1; s<pScenario_->nbShifts_; ++s){
         /* Min/Max constraints */
         dualValues2[s-1] = minWorkedDays;
         dualValues2[s-1] += maxWorkedDays;

         /* Skills coverage */
         dualValues2[s-1] += pModel_->getDual(
            master_->feasibleSkillsAllocCons_[k][s-1][pNurse->pPosition_->id_], true);
      }

      //store vector
      dualValues[k] = dualValues2;
   }

   return dualValues;
}


vector<double> RotationPricer::getStartWorkDualValues(LiveNurse* pNurse){
   int i = pNurse->id_;
   vector<double> dualValues(pDemand_->nbDays_);

   //get dual value associated to the source
   dualValues[0] =  pModel_->getDual(master_->restFlowCons_[i][0], true);
   //get dual values associated to the work flow constraints
   //don't take into account the last which is the sink
   for(int k=1; k<pDemand_->nbDays_; ++k)
      dualValues[k] = pModel_->getDual(master_->workFlowCons_[i][k-1], true);

   return dualValues;
}

vector<double> RotationPricer::getEndWorkDualValues(LiveNurse* pNurse){
   int i = pNurse->id_;
   vector<double> dualValues(pDemand_->nbDays_);

   //get dual values associated to the work flow constraints
   //don't take into account the first which is the source
   //take into account the cost, if the last day worked is k
   for(int k=0; k<pDemand_->nbDays_-1; ++k)
      dualValues[k] = -pModel_->getDual(master_->restFlowCons_[i][k+1], true);

   //get dual value associated to the sink
   dualValues[pDemand_->nbDays_-1] =  pModel_->getDual(
      master_->workFlowCons_[i][pDemand_->nbDays_-1], true);

   return dualValues;
}

double RotationPricer::getWorkedWeekendDualValue(LiveNurse* pNurse){
   return pModel_->getDual(master_->maxWorkedWeekendCons_[pNurse->id_], true);
}

/******************************************************
 * add some forbidden shifts
 ******************************************************/
void RotationPricer::computeForbiddenShifts(
   set<pair<int,int>>& forbiddenShifts, vector<Rotation> rotations){
   //search best rotation
   vector<Rotation>::iterator bestRotation(0);
   double bestDualcost = DBL_MAX;
   for(vector<Rotation>::iterator it = rotations.begin(); it != rotations.end(); ++it)
      if(it->dualCost_ < bestDualcost){
         bestDualcost = it->dualCost_;
         bestRotation = it;
      }

   //forbid shifts of the best rotation
   if(bestDualcost != DBL_MAX)
      for(pair<int,int> pair: bestRotation->shifts_)
         forbiddenShifts.insert(pair);

}


//////////////////////////////////////////////////////////////
//
// B R A N C H I N G  R U L E
//
//////////////////////////////////////////////////////////////

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem* master, const char* name):
                        MyBranchingRule(name), master_(master), pModel_(master->getModel())
{ }

//remove all bad candidates from fixingCandidates
void DiveBranchingRule::logical_fixing(vector<MyObject*>& fixingCandidates){
   if(searchStrategy_ == DepthFirstSearch){
      //look for fractional columns
      //Fix all column above BRANCH_LB
      vector<MyObject*> candidatesToFix;
      //set lb to 1 for the colcumn which is the closest to 1
      MyObject* bestFixingCandidate(0);
      double bestValue = 0;

      //search the good candidates
      for(MyObject* var: fixingCandidates){
         double value = pModel_->getVarValue(var);
         //if var not fractional, continue
         if( pModel_->isInteger(var) )
            continue;
         //if value > BRANCH_LB, add this candidate to candidatesToFix
         if( value > BRANCH_LB)
            candidatesToFix.push_back(var);
         //else if value > bestValue, choose this candidate for the moment
         else if( value > bestValue){
            bestFixingCandidate = var;
            bestValue = value;
         }
      }

      //Clear all the candidates and add the chosen one
      fixingCandidates.clear();
      if(candidatesToFix.size()>0)
         fixingCandidates = candidatesToFix;
      else if(bestFixingCandidate)
         fixingCandidates.push_back(bestFixingCandidate);
   }
   //otherwise clear all the candidates
   else
      fixingCandidates.clear();
}

void DiveBranchingRule::branching_candidates(vector<MyObject*>& branchingCandidates){
   MyObject* bestVar(0);
   double bestValue = DBL_MAX;

   //manage integrality on the skill allocation variables
   switch(searchStrategy_){
   case DepthFirstSearch:
      //variable closest to upper integer
      for(MyObject* var: pModel_->getIntegerCoreVars()){
         if(pModel_->isInteger(var))
            continue;

         double value = pModel_->getVarValue(var);
         double frac = value - floor(value);
         double closeToInt = 1-frac;
         if(closeToInt < bestValue){
            bestVar = var;
            bestValue = closeToInt;
            if(closeToInt<EPSILON)
               break;
         }
      }
      break;
   default:
      //variable closest to .5
      for(MyObject* var: pModel_->getIntegerCoreVars()){
         if(pModel_->isInteger(var))
            continue;

         double value = pModel_->getVarValue(var);
         double frac = value - floor(value);
         double closeTo5 = abs(0.5-frac);
         if(closeTo5 < bestValue){
            bestVar = var;
            bestValue = closeTo5;
            if(closeTo5<EPSILON)
               break;
         }
      }
      break;
   }

   if(bestVar != 0)
      branchingCandidates.push_back(bestVar);
}

