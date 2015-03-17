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
#include "objscip/objscip.h"
#include "scip/pub_var.h"

/* namespace usage */
using namespace std;
using namespace scip;

class RotationPricer : public ObjPricer
{
public:
   RotationPricer(MasterProblem* master, const char* name);
   virtual ~RotationPricer();

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** perform pricing */
   bool pricing();

private:
   /*
    * Variables
    */

   //Pointer to the master problem to link the master and the sub problems
   //
   MasterProblem* master_;

   //name of the pricer handler
   //
   const char* name_;

   //pointers to the data
   //
   Scenario* pScenario_;
   Demand* pDemand_;
   Modeler* pScip_;

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

#endif /* ROTATIONPRICER_H_ */
