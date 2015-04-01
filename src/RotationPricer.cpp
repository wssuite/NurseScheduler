/*
 * RotationPricer.cpp
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#include "RotationPricer.h"

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

      /* Retrieve rotations and add them to the master problem*/
      rotations = subProblem->getRotations();
      for(Rotation rot: rotations)
         master_->addRotation(rot, baseName);
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

