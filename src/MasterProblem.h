/*
 * MasterProblem.h
 *
 *  Created on: 3 févr. 2015
 *      Author: samuel
 */

#ifndef MASTERPROBLEM_H_
#define MASTERPROBLEM_H_

/* Inheritance */
#include "Solver.h"

/* standard library includes */
#include<algorithm>
#include <cfloat>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

/* Tools include */
#include "ScipModeler.h"
#include "MyTools.h"
#include "Nurse.h"

class MasterProblem : public Solver{


public:
   // Specific constructor and destructor
   MasterProblem(Scenario* pScenario, Demand* pDemand,
      Preferences* pPreferences, vector<State>* pInitState, vector<Roster> solution = {});
   ~MasterProblem(){}

   // Main method to solve the rostering problem for a given input
   void solve();

   /*
    * Solving parameters
    */
   char* PB_NAME = "GenCol";
   int solvingTime;
   int bigM = 10000;

private:
   ScipModeler scip_;
   vector2D positionsPerSkill_;//link positions to skills
   vector2D skillsPerPosition_;//link skills to positions

   /*
    * Variables
    */
   vector< vector<SCIP_VAR*> > columnVars_; //binary variables for the columns

   vector< vector<SCIP_VAR*> > restingVars_; //binary variables for the resting arcs in the rotation network
   vector< vector< vector<SCIP_VAR*> > > longRestingVars_; //binary variables for the resting arcs in the rotation network

   vector<SCIP_VAR*> minWorkedDaysVars_; //count the number of missing worked days per nurse
   vector<SCIP_VAR*> maxWorkedDaysVars_; //count the number of exceeding worked days per nurse
   vector<SCIP_VAR*> maxWorkedWeekendVars_; //count the number of exceeding worked weekends per nurse

   vector< vector< vector<SCIP_VAR*> > > optDemandVars_; //count the number of missing nurse to reach the optimal
   vector< vector< vector< vector<SCIP_VAR*> > > > skillsAllocVars_; //makes the allocation of the skills

   /*
    * Constraints
    */
   //transmission of the flow on the resting nodes
   //initialization of the flow constraint at the first position of each restFlowCons_[i] (i=nurse)
   vector< vector<SCIP_CONS*> > restFlowCons_;
   //transmission of the flow on the working nodes
   //end of the flow constraint at the last position of each workFlowCons_[i] (i=nurse)
   vector< vector<SCIP_CONS*> > workFlowCons_;

   vector<SCIP_CONS*> minWorkedDaysCons_; //count the number of missing worked days per nurse
   vector<SCIP_CONS*> maxWorkedDaysCons_; //count the number of exceeding worked days per nurse
   vector<SCIP_CONS*> maxWorkedWeekendCons_; //count the number of exceeding worked weekends per nurse

   vector< vector< vector<SCIP_CONS*> > > minDemandCons_; //ensure a minimal coverage per day, per shift, per skill
   vector< vector< vector<SCIP_CONS*> > > optDemandCons_; //count the number of missing nurse to reach the optimal
   vector< vector< vector<SCIP_CONS*> > > feasibleSkillsAllocCons_; // ensures that each nurse works with the good skill

   /*
    * Methods
    */

   // Main method to build the rostering problem for a given input
   void build();

   //Initialization of the rostering problem with/without solution
   void initialize(vector<Roster> solution);

   //add the correct constraints and coefficients for the nurse i working on day k on shift s
   //if s=-1, the nurse works on all shifts
   //if firstDay, add a working arc starting on day k in the rotation network
   //if lastDay, add a working arc ending on day k+1 in the rotation network
   //return the number of constraints added
   int addConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, int s=-1, bool firstDay = false, bool lastDay = false);

   /* Build each set of constraints - Add also the coefficient of a column for each set */
   void buildRotationCons();
   int addRotationConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, bool firstDay, bool lastDay);
   void buildMinMaxCons();
   int addMinMaxConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k);
   void buildSkillsCoverageCons();
   int addSkillsCoverageConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, int s=-1);

};

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------
struct Rotation {

   // Specific constructors and destructors
   //

   Rotation(map<int,int> shift, double cost = 999999, double dualCost = 999999) :
      cost_(cost), shifts_(shift), dualCost_(dualCost), length_(shift.size())
   {
      firstDay_ = 999;
      for(map<int,int>::iterator itS = shift.begin(); itS != shift.end(); ++itS)
         if(itS->first < firstDay_) firstDay_ = itS->first;
   };

   Rotation(int firstDay, vector<int> shiftSuccession, double cost = 999999, double dualCost = 999999) :
      cost_(cost), dualCost_(dualCost), firstDay_(firstDay), length_(shiftSuccession.size())
   {
      map<int,int> m;
      for(int k=0; k<shiftSuccession.size(); k++) m.insert(pair<int,int>( (firstDay+k) , shiftSuccession[k] ));
      shifts_ = m;
   }

   ~Rotation(){};

   // Cost
   //
   double cost_;

   // Dual cost as found in the subproblem
   //
   double dualCost_;

   // Shifts to be performed
   //
   map<int,int> shifts_;

   // First worked day
   //
   int firstDay_;

   // Duration
   //
   int length_;

};

#endif /* MASTERPROBLEM_H_ */
