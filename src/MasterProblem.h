/*
 * MasterProblem.h
 *
 *  Created on: 3 f��vr. 2015
 *      Author: samuel
 */

#ifndef MASTERPROBLEM_H_
#define MASTERPROBLEM_H_

/* Inheritance */
#include "Solver.h"

/* Tools include */
#include "MyTools.h"

/* My includes */
#include "Nurse.h"
#include "ScipModeler.h"
#include "Modeler.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------
enum CostType {TOTAL_COST, CONS_SHIFTS_COST, CONS_WORKED_DAYS_COST, COMPLETE_WEEKEND_COST, PREFERENCE_COST, INIT_REST_COST};

struct Rotation {

   // Specific constructors and destructors
   //
   Rotation(map<int,int> shift, LiveNurse* nurse = NULL, double cost = 999999, double dualCost = 999999) :
      shifts_(shift), pNurse_(nurse), cost_(cost), dualCost_(dualCost), length_(shift.size())
   {
      firstDay_ = 999;
      for(map<int,int>::iterator itS = shift.begin(); itS != shift.end(); ++itS)
         if(itS->first < firstDay_) firstDay_ = itS->first;
   };

   Rotation(int firstDay, vector<int> shiftSuccession, LiveNurse* nurse = NULL, double cost = 999999, double dualCost = 999999) :
      pNurse_(nurse), cost_(cost), dualCost_(dualCost), firstDay_(firstDay), length_(shiftSuccession.size())
   {
      map<int,int> m;
      for(int k=0; k<shiftSuccession.size(); k++) m.insert(pair<int,int>( (firstDay+k) , shiftSuccession[k] ));
      shifts_ = m;
   }

   ~Rotation(){};

   //the nurse
   //
   LiveNurse* pNurse_;

   // Cost
   //
   double cost_;
   double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_, initRestCost_ ;

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

   //Compute the cost of a rotation
   //
   void computeCost(Scenario* pScenario, Preferences* pPreferences, int horizon);

};


//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------
class MasterProblem : public Solver{
   //allows RotationPricer to access all private arguments and methods of MasterProblem
   friend class RotationPricer;
public:
   // Specific constructor and destructor
   MasterProblem(Scenario* pScenario, Demand* pDemand,
      Preferences* pPreferences, vector<State>* pInitState, vector<Roster> solution = {});
   ~MasterProblem();

   //solve the rostering problem
   void solve();

   //get the pointer to scip
   Modeler* getModel(){
      return pModel_;
   }

   /*
    * Solving parameterdoubles
    */
   char* PB_NAME = "GenCol";
   int solvingTime;
   int bigM = 10000;

private:
   Modeler* pModel_;
   vector2D positionsPerSkill_;//link positions to skills
   vector2D skillsPerPosition_;//link skills to positions
   MyPricer* pPricer_;//prices the rotations
   vector< map<MyVar*, Rotation> > rotations_;//stores the scip variables and the rotations for each nurse

   /*
    * Variables
    */
   vector< vector<MyVar*> > columnVars_; //binary variables for the columns

   vector< vector<MyVar*> > restingVars_; //binary variables for the resting arcs in the rotation network
   vector< vector< vector<MyVar*> > > longRestingVars_; //binary variables for the resting arcs in the rotation network

   vector<MyVar*> minWorkedDaysVars_; //count the number of missing worked days per nurse
   vector<MyVar*> maxWorkedDaysVars_; //count the number of exceeding worked days per nurse
   vector<MyVar*> maxWorkedWeekendVars_; //count the number of exceeding worked weekends per nurse

   vector< vector< vector<MyVar*> > > optDemandVars_; //count the number of missing nurse to reach the optimal
   vector< vector< vector< vector<MyVar*> > > > skillsAllocVars_; //makes the allocation of the skills

   /*
    * Constraints
    */
   //transmission of the flow on the resting nodes
   //initialization of the flow constraint at the first position of each restFlowCons_[i] (i=nurse)
   vector< vector<MyCons*> > restFlowCons_;
   //transmission of the flow on the working nodes
   //end of the flow constraint at the last position of each workFlowCons_[i] (i=nurse)
   vector< vector<MyCons*> > workFlowCons_;

   vector<MyCons*> minWorkedDaysCons_; //count the number of missing worked days per nurse
   vector<MyCons*> maxWorkedDaysCons_; //count the number of exceeding worked days per nurse
   vector<MyCons*> maxWorkedWeekendCons_; //count the number of exceeding worked weekends per nurse

   vector< vector< vector<MyCons*> > > minDemandCons_; //ensure a minimal coverage per day, per shift, per skill
   vector< vector< vector<MyCons*> > > optDemandCons_; //count the number of missing nurse to reach the optimal
   vector< vector< vector<MyCons*> > > feasibleSkillsAllocCons_; // ensures that each nurse works with the good skill

   /*
    * Methods
    */

   // Main method to build the rostering problem for a given input
   void build();

   //Initialization of the rostering problem with/without solution
   void initialize(vector<Roster> solution);

   //solve a solution in the output
   void storeSolution();

   //Create a new rotation variable
   //add the correct constraints and coefficients for the nurse i working on a rotation
   //if s=-1, the nurse works on all shifts
   //store the rotation in rotations_
   void addRotation(Rotation rotation, char* baseName);

   //compute and add the last rotation finishing on the day just before the first one
   Rotation computeInitStateRotation(LiveNurse* pNurse);

   //get the cost of all shosen rotations in solution sol for a certain CostType
   double getRotationCosts(CostType costType = TOTAL_COST, bool initStateRotation = false);

   /* Build each set of constraints - Add also the coefficient of a column for each set */
   void buildRotationCons();
   int addRotationConsToCol(vector<MyCons*>* cons, vector<double>* coeffs, int i, int k, bool firstDay, bool lastDay);
   void buildMinMaxCons();
   int addMinMaxConsToCol(vector<MyCons*>* cons, vector<double>* coeffs, int i, int k, bool weekend = false);
   void buildSkillsCoverageCons();
   int addSkillsCoverageConsToCol(vector<MyCons*>* cons, vector<double>* coeffs, int i, int k, int s=-1);

   /* Display functions */
   string costsConstrainstsToString();

};

#endif /* MASTERPROBLEM_H_ */
