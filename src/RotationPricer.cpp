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
bool RotationPricer::pricing(double bound, bool before_fathom){
   //=false if once optimality hasn't be proven
   bool optimal = true;
   //forbidden shifts
   set<pair<int,int> > forbiddenShifts;
   //computed new rotations
   vector<Rotation> rotations;

//   std::cout << "# ------- BEGIN ------- Subproblems..." << std::endl;

   //count and store the nurses for whom their subproblem has generated rotations.
   int nbSubProblemSolved = 0, nbIteration = 0;
   double minDualCoast = 0;
   vector<LiveNurse*> nursesSolved;
   for(vector<LiveNurse*>::iterator it0 = nursesToSolve_.begin(); it0 != nursesToSolve_.end();){
      ++nbIteration;
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
	   //options.push_back(SOLVE_VERY_SHORT_ONLY);
	   //options.push_back(SOLVE_ONE_SINK_PER_LAST_DAY);
	   options.push_back(SOLVE_FORBIDDEN_RESET);
	   options.push_back(SOLVE_NEGATIVE_ALLVALUES);

//	   cout << "#  SP " << pNurse->name_ << " begins" << endl;

	   /* Solve subproblems */
	   if( subProblem->solve(pNurse, &costs, options, forbiddenShifts, false, 10) )
		   optimal = false;
	   else
		   subProblem->solve(pNurse, &costs, options, forbiddenShifts, true);

//	   cout << "#  SP " << pNurse->name_ << " solved" << endl;

	   /*
	    * Rotations
	    */

		/* Retrieve rotations */
		rotations = subProblem->getRotations();
		/* sort rotations */
      for(Rotation& rot: rotations){
         rot.computeCost(pScenario_, master_->pPreferences_, master_->pDemand_->nbDays_);
         rot.computeDualCost(workDualCosts, startWorkDualCosts, endWorkDualCosts, workedWeekendDualCost);
      }
		std::sort(rotations.begin(), rotations.end(), Rotation::compareDualCost);
		/* add them to the master problem */
		int nbRotationsAdded = 0;
		for(Rotation& rot: rotations){
			master_->addRotation(rot, baseName);
			++nbRotationsAdded;
			if(nbRotationsAdded > nbMaxRotationsToAdd_)
			   break;
		}

//		cout << "#  SP " << pNurse->name_ << " added columns" << endl;

      //count if the subproblem has generated some rotations and then store the nurse
      if(rotations.size() > 0){
         ++nbSubProblemSolved;
         if(rotations[0].dualCost_ < minDualCoast)
            minDualCoast = rotations[0].dualCost_;
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

   //set statistics
   BcpModeler* model = dynamic_cast<BcpModeler*>(pModel_);
   if(model){
      model->setLastNbSubProblemsSolved(nbIteration);
      model->setLastMinDualCost(minDualCoast);
   }

//   std::cout << "# -------  END  ------- Subproblems!" << std::endl;

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
            master_->numberOfNursesByPositionCons_[k][s-1][pNurse->pPosition_->id_], true);
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

/*************************************************************
 * Diving branching rule: dive then close to .5
 *************************************************************/

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem* master, const char* name):
                        MyBranchingRule(name), master_(master), pModel_(master->getModel()), nbBranchingCandidates_(5)
{ }

//add all good candidates
void DiveBranchingRule::logical_fixing(vector<MyObject*>& fixingCandidates){
   //look for fractional columns
   //Fix all column above BRANCH_LB
   //search the good candidates
   for(int i=0; i<master_->getRotations().size(); ++i)
      for(pair<MyObject*, Rotation> var: master_->getRotations()[i]){
         double value = pModel_->getVarValue(var.first);
         //if var not fractional or the rotation is not a real rotation (length = 0), continue
         if( var.second.length_==0 || pModel_->isInteger(var.first) )
            continue;
         //if value > BRANCH_LB, add this candidate to candidatesToFix
         if( value > BRANCH_LB)
            fixingCandidates.push_back(var.first);
      }
}

void DiveBranchingRule::branching_candidates(vector<MyObject*>& branchingCandidates){
   //search all candidates
   vector<pair<MyObject*, double>> candidates;
   for(int i=0; i<master_->getRotations().size(); ++i)
      for(pair<MyObject*, Rotation> var: master_->getRotations()[i])
         //if var is fractional and the rotation is a real rotation (length > 0)
         if(var.second.length_>0 && !pModel_->isInteger(var.first) )
            candidates.push_back(pair<MyObject*, double>(var.first, pModel_->getVarValue(var.first)));

   switch(searchStrategy_){
   case DepthFirstSearch:
      sort(candidates.begin(), candidates.end(), compareColumnCloseToInt);
      break;
   default:
      sort(candidates.begin(), candidates.end(), compareColumnCloseTo5);
   }

   for(int i=nbBranchingCandidates_-1; i>=0; --i)
      branchingCandidates.push_back(candidates[i].first);
}

bool DiveBranchingRule::compareColumnCloseToInt(pair<MyObject*, double>& obj1, pair<MyObject*, double>& obj2){
   double frac1 = obj1.second - floor(obj1.second), frac2 = obj2.second - floor(obj2.second);
   double closeToInt1 = 1-frac1, closeToInt2 = 1-frac2;
   return (closeToInt1 < closeToInt2);
}

bool DiveBranchingRule::compareColumnCloseTo5(pair<MyObject*, double>& obj1, pair<MyObject*, double>& obj2){
   double frac1 = obj1.second - floor(obj1.second), frac2 = obj2.second - floor(obj2.second);
   double closeTo5_1 = abs(0.5-frac1), closeTo5_2 = abs(0.5-frac2);
   return (closeTo5_1 < closeTo5_2);
}

//if(mediumCandidates_.size() == 0)
//   for(MyObject* var: pModel_->getIntegerCoreVars()){
//      string str2 = "nursesNumber";
//      string str0(var->name_);
//      string str1 = str0.substr(0,str2.size());
//      if(strcmp(str1.c_str(), str2.c_str()) == 0)
//         bestCandidates_.push_back(var);
//      else
//         mediumCandidates_.push_back(var);
//   }
//
//MyObject *bestVar(0);
//double bestValue = DBL_MAX;
//
////manage integrality on the skill allocation variables
//switch(searchStrategy_){
//case DepthFirstSearch:
//   //variable closest to upper integer
//   for(MyObject* var: bestCandidates_){
//      if(pModel_->isInteger(var))
//         continue;
//
//      double value = pModel_->getVarValue(var);
//      double frac = value - floor(value);
//      double closeToInt = 1-frac;
//
//         if(closeToInt < bestValue){
//            bestVar = var;
//            bestValue = closeToInt;
//            if(closeToInt<EPSILON)
//               break;
//         }
//   }
//
//   if(bestVar != 0)
//      break;
//
//   for(MyObject* var: mediumCandidates_){
//      if(pModel_->isInteger(var))
//         continue;
//
//      double value = pModel_->getVarValue(var);
//      double frac = value - floor(value);
//      double closeToInt = 1-frac;
//
//         if(closeToInt < bestValue){
//            bestVar = var;
//            bestValue = closeToInt;
//            if(closeToInt<EPSILON)
//               break;
//         }
//   }
//
//   break;
//default:
//   //variable closest to .5
//   for(MyObject* var: bestCandidates_){
//      if(pModel_->isInteger(var))
//         continue;
//
//      double value = pModel_->getVarValue(var);
//      double frac = value - floor(value);
//      double closeTo5 = abs(0.5-frac);
//
//      if(closeTo5 < bestValue){
//         bestVar = var;
//         bestValue = closeTo5;
//         if(closeTo5<EPSILON)
//            break;
//      }
//   }
//
//   if(bestVar != 0)
//      break;
//
//   for(MyObject* var: mediumCandidates_){
//      if(pModel_->isInteger(var))
//         continue;
//
//      double value = pModel_->getVarValue(var);
//      double frac = value - floor(value);
//      double closeTo5 = abs(0.5-frac);
//
//      if(closeTo5 < bestValue){
//         bestVar = var;
//         bestValue = closeTo5;
//         if(closeTo5<EPSILON)
//            break;
//      }
//   }
//
//   break;
//}
//
//if(bestVar != 0)
//   branchingCandidates.push_back(bestVar);

/*************************************************************
 * CorePriority branching rule: branch on core variables first
 *************************************************************/

/* Constructs the branching rule object. */
CorePriorityBranchingRule::CorePriorityBranchingRule(Modeler* pModel, const char* name):
                        MyBranchingRule(name), pModel_(pModel)
{ }

//remove all bad candidates from fixingCandidates while keeping the order
void CorePriorityBranchingRule::logical_fixing(vector<MyObject*>& fixingCandidates){
   //choose the var nursesNumber
   for(MyObject* var: pModel_->getIntegerCoreVars()){
      string str2 = "nursesNumber";
      string str0(var->name_);
      string str1 = str0.substr(0,str2.size());
      if(strcmp(str1.c_str(), str2.c_str()) == 0)
         fixingCandidates.push_back(var);
   }
}

//remove all worst/best candidates from fixingCandidates while keeping the order
void CorePriorityBranchingRule::branching_candidates(vector<MyObject*>& branchingCandidates){
   //choose the var nursesNumber
   for(MyObject* var: pModel_->getIntegerCoreVars()){
      string str2 = "skillsAlloc";
      string str0(var->name_);
      string str1 = str0.substr(0,str2.size());
      if(strcmp(str1.c_str(), str2.c_str()) == 0)
         branchingCandidates.push_back(var);
   }
}


