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

   //   std::cout << "# ------- BEGIN ------- Subproblems..." << std::endl;

   for(LiveNurse* pNurse: sortedNurses){

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
      computeForbiddenShifts(&forbiddenShifts, rotations);

<<<<<<< HEAD
      /* Solve options */
      vector<SolveOption> options;
      //options.push_back(SOLVE_ONE_SINK_PER_FIRST_DAY);
      //options.push_back(SOLVE_ONE_SINK_PER_LAST_DAY);
      options.push_back(SOLVE_FORBIDDEN_RESET);
=======
	   /* Solve options */
	   vector<SolveOption> options;
	   options.push_back(SOLVE_ONE_SINK_PER_FIRST_DAY);
	   //options.push_back(SOLVE_ONE_SINK_PER_LAST_DAY);
	   options.push_back(SOLVE_FORBIDDEN_RESET);
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses

<<<<<<< HEAD
      /* Solve subproblems */
      if( subProblem->solve(pNurse, &costs, options, forbiddenShifts, false, 4) )
         optimal = false;
      else
         subProblem->solve(pNurse, &costs, options, forbiddenShifts, true);
=======
	   cout << "#  SP " << pNurse->name_ << " begins" << endl;

	   /* Solve subproblems */
	   if( subProblem->solve(pNurse, &costs, options, forbiddenShifts, false, 10) )
		   optimal = false;
	   else
		   subProblem->solve(pNurse, &costs, options, forbiddenShifts, true);

	   cout << "#  SP " << pNurse->name_ << " solved" << endl;
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses

<<<<<<< HEAD
      /* Retrieve rotations and add them to the master problem*/
      rotations = subProblem->getRotations();
      for(Rotation rot: rotations){
         double c = rot.cost_;
         rot.computeCost(pScenario_, master_->pPreferences_, master_->pDemand_->nbDays_);
         rot.computeDualCost(workDualCosts, startWorkDualCosts, endWorkDualCosts, workedWeekendDualCost);
         master_->addRotation(rot, baseName);
      }
=======
		/* Retrieve rotations and add them to the master problem*/
		rotations = subProblem->getRotations();
		for(Rotation rot: rotations){
			double c = rot.cost_;
			rot.computeCost(pScenario_, master_->pPreferences_, master_->pDemand_->nbDays_);
			rot.computeDualCost(workDualCosts, startWorkDualCosts, endWorkDualCosts, workedWeekendDualCost);
			master_->addRotation(rot, baseName);
		}

		cout << "#  SP " << pNurse->name_ << " added columns" << endl;
>>>>>>> branch 'master' of https://github.com/jeremyomer/RosterDesNurses
   }

   //   std::cout << "# -------  END  ------- Subproblems!" << std::endl;

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
   set<pair<int,int> >* forbiddenShifts, vector<Rotation> rotations){

}


//////////////////////////////////////////////////////////////
//
// B R A N C H I N G  R U L E
//
//////////////////////////////////////////////////////////////

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem* master, const char* name):
                        MyBranchingRule(name), master_(master), pModel_(master->getModel()),
                        checkChildren(false)
{ }

//remove all bad candidates from fixingCandidates
void DiveBranchingRule::logical_fixing(vector<MyObject*>& fixingCandidates){
   //look for fractional columns
   //if value > branchUB, then set lb to 1
   vector<MyObject*> fixingCandidates2;
   for(MyObject* var: fixingCandidates){
      double value = pModel_->getVarValue(var);
      //if var not fractional, continue
      if( pModel_->isInteger(var) )
         continue;
      //if value < branchLB, remove it
      else if( value < BRANCH_LB)
         continue;
      fixingCandidates2.push_back(var);
   }
   fixingCandidates = fixingCandidates2;
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

