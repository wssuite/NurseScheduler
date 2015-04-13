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
                     MyPricer(name),
                     master_(master), pScenario_(master->pScenario_), pDemand_(master->pDemand_), pModel_(master->getModel())
{ }

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

   /* sort the nurses */
   vector<LiveNurse*> sortedNurses = sortNurses();

   for(LiveNurse* pNurse: sortedNurses){
      /* Build or re-use a subproblem */
      SubProblem* subProblem;
      //search the contract
      map<const Contract*, SubProblem*>::iterator it =  subProblems_.find(pNurse->pContract_);
      //if doesn't find => create new subproblem
      if( it == subProblems_.end() ){
         subProblem = new SubProblem(pScenario_, pDemand_, pNurse->pContract_);
         subProblems_.insert(it, pair<const Contract*, SubProblem*>(pNurse->pContract_, subProblem));
      }

		// SR: Ici, je modifie pour en creer un nouveau a chaque fois car pour l'instant, j'ai des problemes si je relance l'algo de plus courts
		//     chemins sur un SP qui existait deja avant (a changer)
		// TODO : commenter ces lignes pour que le SP ne soit pas recree a chaque fois, à implementer plus tard
		else if(true) {
			subProblem = new SubProblem(pScenario_, pDemand_, pNurse->pContract_);
			it->second = subProblem;
		}

      //otherwise retrieve the subproblem associated to the contract
      else
         subProblem = it->second;

      /* Retrieves dual values */
      vector< vector<double> > workDualCosts = getWorkDualValues(pNurse);
      vector<double> startWorkDualCosts = getStartWorkDualValues(pNurse);
      vector<double> endWorkDualCosts = getEndWorkDualValues(pNurse);
      double workedWeekendDualCost = getWorkedWeekendDualValue(pNurse);
      Costs costs (&workDualCosts, &startWorkDualCosts, &endWorkDualCosts, workedWeekendDualCost);

      /* Compute forbidden */
      computeForbiddenShifts(&forbiddenShifts, rotations);

      /* Solve subproblems */
      if( subProblem->solve(pNurse, &costs, forbiddenShifts, false) )
         optimal = false;
      else
         subProblem->solve(pNurse, &costs, forbiddenShifts, true);


      // SR - TODO : calcul du cout a chaque fois, car pas fait dans le SP
		/* Retrieve rotations and add them to the master problem*/
		rotations = subProblem->getRotations();
		for(Rotation rot: rotations){
			std::cout << "# Cost update check : " << rot.cost_;
			rot.computeCost(pScenario_, master_->pPreferences_, master_->pDemand_->nbDays_);
			std::cout << "  ->  " << rot.cost_ << std::endl;
			master_->addRotation(rot, baseName);

		}

      /* Retrieve rotations and add them to the master problem*/
      /*
      rotations = subProblem->getRotations();
      for(Rotation rot: rotations)
         master_->addRotation(rot, baseName);
      */

   }

   return optimal;
}

/******************************************************
 * Sort the nurses to give an order to solve subproblems
 ******************************************************/
vector<LiveNurse*> RotationPricer::sortNurses(){
   return master_->theLiveNurses_;
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
   set<pair<int,int> >* forbiddenShifts, vector<Rotation> rotations){

}


//////////////////////////////////////////////////////////////
//
// R O T A T I O N   P R I C E R
//
//////////////////////////////////////////////////////////////

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem* master, const char* name):
                     MyBranchingRule(name), master_(master), pModel_(master->getModel()),
                     checkChildren(false)
{ }

void DiveBranchingRule::logical_fixing(vector<MyObject*>& fixingCandidates){
   //look for fractional columns
   //if value > branchUB, then set lb to 1
      for(int i=0; i<master_->getRotations().size(); ++i){
         for(pair<MyObject*, Rotation> pair: master_->getRotations()[i]){
            //if var not fractional, continue
            if( pModel_->isInteger(pair.first) )
               continue;

            //if value > branchLB, branch on it
            if( BRANCH_LB < pModel_->getVarValue(pair.first) )
               fixingCandidates.push_back(pair.first);
      }
   }//end search fractional variables
}

void DiveBranchingRule::branching_candidates(vector<MyObject*>& branchingCandidates){
   //if we have checked children, this is the end, no more branching
   //otherwise look for fractional columns
   if(!checkChildren){
      for(int i=0; i<master_->getRotations().size(); ++i){
         MyObject* bestVariable(0); //store best variable for branching for nurse i
         for(pair<MyObject*, Rotation> pair: master_->getRotations()[i]){// check all the columns
            //if var not fractional, continue
            if( pModel_->isInteger(pair.first) )
               continue;

            /* else, branch for each nurse on the column with the highest value */
            //check if already a variable in bestVariable
            //update the variable for the nurse if the variable is higher
            if( (bestVariable == 0)  ||
               (pModel_->getVarValue(pair.first) > pModel_->getVarValue(bestVariable)) )
               bestVariable = pair.first;
         }// end rotations

         //if branching candidates -> store it
         if(bestVariable != 0)
            branchingCandidates.push_back(bestVariable);
      }

      //if branching candidates -> branching just once
      checkChildren = (branchingCandidates.size() > 0);
   }//end search fractional variables
}

