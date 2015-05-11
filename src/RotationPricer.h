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
   bool pricing(double bound=0, bool before_fathom = true);

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
   vector<LiveNurse*> nursesToSolve_;
   map<const Contract*, SubProblem*> subProblems_;

   /*
    * Settings
    */
   int nbMaxRotationsToAdd_, nbSubProblemsToSolve_;

   /*
    * Methods
    */

   //get the duals values per day and per shift for a nurse
   //
   vector< vector<double> > getWorkDualValues(LiveNurse* pNurse);
   vector<double> getStartWorkDualValues(LiveNurse* pNurse);
   vector<double> getEndWorkDualValues(LiveNurse* pNurse);
   double getWorkedWeekendDualValue(LiveNurse* pNurse);

   //compute some forbidden shifts from the lasts rotations and forbidden shifts
   //
   void computeForbiddenShifts(set<pair<int,int> >& forbiddenShifts, vector<Rotation> rotations);
};

static bool compareObject(const pair<MyObject*,double>& p1, const pair<MyObject*,double>& p2);

class DiveBranchingRule: public MyBranchingRule
{
public:
   DiveBranchingRule(MasterProblem* master, const char* name);
   virtual ~DiveBranchingRule() { }

   /* compute branching decisions */
   void branching_candidates(vector<MyObject*>& branchingCandidates);

   /* branch on the number of nurses */
   void branchOnNumberOfNurses(vector<MyObject*>& branchingCandidates);

   /* branch on a set of resting arcs */
   void branchOnRestingArcs(vector<MyObject*>& branchingCandidates);

   /* compute fixing decisions */
   void logical_fixing(vector<MyObject*>& fixingCandidates);

   /* compare columns */
   static bool compareColumnCloseToInt(pair<MyObject*, double> obj1, pair<MyObject*, double> obj2);

   static bool compareColumnCloseTo5(pair<MyObject*, double> obj1, pair<MyObject*, double> obj2);

protected:
   //Pointer to the master problem to link the master and the sub problems
   //
   MasterProblem* master_;

   //pointers to the data
   //
   Modeler* pModel_;

   //number of candidates
   //
   int nbBranchingCandidates_;

   //vectors of the variables on which we can branch
   //
   vector<MyObject*> bestCandidates_, mediumCandidates_;
};

class CorePriorityBranchingRule: public MyBranchingRule
{
public:
   CorePriorityBranchingRule(Modeler* pModel, const char* name);
   virtual ~CorePriorityBranchingRule() { }

   /* compute branching decisions */
   void branching_candidates(vector<MyObject*>& branchingCandidates);

   /* compute fixing decisions */
   void logical_fixing(vector<MyObject*>& fixingCandidates);

protected:
   //pointers to the data
   //
   Modeler* pModel_;
};

#endif /* ROTATIONPRICER_H_ */
