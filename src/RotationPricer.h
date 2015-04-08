/*
 * RotationPricer.h
 *
 * Allow to link the sub problems and the master problem trough scip
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#ifndef ROTATIONPRICER_H_
#define ROTATIONPRICER_H_

#include "MyTools.h"
#include "MasterProblem.h"
#include "SubProblem.h"
#include "Modeler.h"

/* namespace usage */
using namespace std;

class RotationPricer: public MyPricer
{
public:
   RotationPricer(MasterProblem* master, const char* name);
   virtual ~RotationPricer();

   /* perform pricing */
   bool pricing(double bound = 0);

private:
   //Pointer to the master problem to link the master and the sub problems
   //
   MasterProblem* master_;

   //pointers to the data
   //
   Scenario* pScenario_;
   Demand* pDemand_;
   Modeler* pModel_;

   //map of the contract and sub problems
   //one subproblem per type of contract
   //
   map<const Contract*, SubProblem*> subProblems_;

   /*
    * Methods
    */

   //sort the nurses to give an order to solve subproblems
   //
   vector<LiveNurse*> sortNurses();

   //get the duals values per day and per shift for a nurse
   //
   vector< vector<double> > getWorkDualValues(LiveNurse* pNurse);
   vector<double> getStartWorkDualValues(LiveNurse* pNurse);
   vector<double> getEndWorkDualValues(LiveNurse* pNurse);
   double getWorkedWeekendDualValue(LiveNurse* pNurse);

   //compute some forbidden shifts from the lasts rotations and forbidden shifts
   //
   void computeForbiddenShifts(set<pair<int,int> >* forbiddenShifts, vector<Rotation> rotations);
};

class DiveBranchingRule: public MyBranchingRule
{
public:
   DiveBranchingRule(MasterProblem* master, const char* name);
   virtual ~DiveBranchingRule() { }

   /* compute branching decisions */
   void branching_candidates(vector<MyObject*>& branchingCandidates);

   /* compute fixing decisions */
   void logical_fixing(vector<MyObject*>& fixingCandidates);

protected:
   //Pointer to the master problem to link the master and the sub problems
   //
   MasterProblem* master_;

   //pointers to the data
   //
   Modeler* pModel_;

   // Settings
   //
   bool checkChildren; //true = branch once
};

#endif /* ROTATIONPRICER_H_ */
