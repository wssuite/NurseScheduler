/*
 * MasterProblem.cpp
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#include "MasterProblem.h"
#include "BcpModeler.h"
#include "CbcModeler.h"
#include "ScipModeler.h"
#include "RotationPricer.h"

/* namespace usage */
using namespace std;



//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

void Rotation::computeCost(Scenario* pScenario, Preferences* pPreferences, int horizon){
   //check if pNurse points to a nurse
   if(pNurse_ == NULL)
      Tools::throwError("LiveNurse = NULL");

   /************************************************
    * Compute all the costs of a rotation:
    ************************************************/
   //   double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_ ;

   //if first day of the planning, check on the past, otherwise 0 (rest)
   int lastShift = (firstDay_==0) ? pNurse_->pStateIni_->shift_ : 0;
   //nbConsShift = number of consecutive shift
   //if first day of the planning, check on the past, otherwise 0
   int nbConsShifts = (firstDay_==0) ? pNurse_->pStateIni_->consShifts_ : 0;
   //consShiftCost = cost of be outside of the interval [min,max] of the consecutives shifts
   consShiftsCost_ = 0;

   //nbConsWorked = number of consecutive worked days
   //if first day of the planning, check on the past , otherwise 0
   int nbConsDaysWorked = (firstDay_==0) ? pNurse_->pStateIni_->consDaysWorked_ : 0;
   //consWorkedCost = cost of be outside of the interval [min,max] of the consecutives worked days
   consDaysWorkedCost_ = 0;

   //cost of not doing the whole weekend
   completeWeekendCost_ = 0;

   //preferencesCost = cost of not respecting preferences
   preferenceCost_ = 0;

   //initial resting cost
   initRestCost_ =  0;

   /*
    * Compute consShiftCost
    */

   // if the initial shift has already exceeded the max, substract now the cost that will be readd later
   if( (firstDay_==0) && (lastShift>0) &&
      (nbConsShifts > pScenario->maxConsShifts_[lastShift])){
      consShiftsCost_ -= (nbConsShifts-pScenario->maxConsShifts_[lastShift])*WEIGHT_CONS_SHIFTS;
   }

   for(int k=firstDay_; k<firstDay_+length_; ++k){
      if(lastShift == shifts_[k]){
         nbConsShifts ++;
         continue;
      }
      if(lastShift > 0){
         int diff = max(pScenario->minConsShifts_[lastShift] - nbConsShifts,
            nbConsShifts-pScenario->maxConsShifts_[lastShift]);
         if(diff>0) {
            consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
         }
      }
      //initialize nbConsShifts and lastShift
      nbConsShifts = 1;
      lastShift = shifts_[k];
   }

   //compute consShiftsCost for the last shift
   int diff = max((firstDay_+length_ == horizon) ? 0 : pScenario->minConsShifts_[lastShift] - nbConsShifts,
      nbConsShifts-pScenario->maxConsShifts_[lastShift]);
   if(diff>0) {
      consShiftsCost_ += diff * WEIGHT_CONS_SHIFTS;
   }


   /*
    * Compute consDaysWorkedCost
    */

   nbConsDaysWorked += length_;
   //check if nbConsDaysWorked < min, if finishes on last day, does not count
   if(nbConsDaysWorked < pNurse_->minConsDaysWork() && firstDay_+length_ < horizon)
      consDaysWorkedCost_ += (pNurse_->minConsDaysWork() - nbConsDaysWorked) * WEIGHT_CONS_DAYS_WORK;
   //check if nbConsDaysWorked > max
   else if(nbConsDaysWorked > pNurse_->maxConsDaysWork())
      consDaysWorkedCost_ += (nbConsDaysWorked - pNurse_->maxConsDaysWork()) * WEIGHT_CONS_DAYS_WORK;

   /*
    * Compute completeWeekendCost
    */
   if(pNurse_->needCompleteWeekends()){
      //if first day is a Sunday, the saturday is not worked
      if(Tools::isSunday(firstDay_))
         completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
      //if last day + 1 is a Sunday, the sunday is not worked
      if(Tools::isSunday(firstDay_+length_))
         completeWeekendCost_ += WEIGHT_COMPLETE_WEEKEND;
   }

   /*
    * Compute preferencesCost
    */

   for(int k=firstDay_; k<firstDay_+length_; ++k)
      if(pPreferences->wantsTheShiftOff(pNurse_->id_, k, shifts_[k]))
         preferenceCost_ += WEIGHT_PREFERENCES;

   /*
    * Compute initial resting cost
    */

   if(firstDay_==0 && pNurse_->pStateIni_->shift_==0){
      int diff = pNurse_->minConsDaysOff() - pNurse_->pStateIni_->consDaysOff_;
      initRestCost_ = (diff > 0) ? diff*WEIGHT_CONS_DAYS_OFF : 0;
   }

   /*
    * Compute the sum of the cost and stores it in cost_
    */
   if(false){
	   cout << "# Calcul du cout:" << endl;
	   cout << "#       | Consecutive shifts: " << consShiftsCost_ << endl;
	   cout << "#       | Consecutive days  : " << consDaysWorkedCost_ << endl;
	   cout << "#       | Complete weekends : " << completeWeekendCost_ << endl;
	   cout << "#       | Preferences       : " << preferenceCost_ << endl;
	   cout << "#       | Initial rest      : " << initRestCost_ << endl;
	   cout << "# " << endl;
   }

   cost_ = consShiftsCost_ + consDaysWorkedCost_ + completeWeekendCost_ + preferenceCost_ +  initRestCost_;
}


void Rotation::computeDualCost(vector< vector<double> > workDualCosts, vector<double> startWorkDualCosts,
      vector<double> endWorkDualCosts, double workedWeekendDualCost){
      //check if pNurse points to a nurse
      if(pNurse_ == NULL)
         Tools::throwError("LiveNurse = NULL");

      /************************************************
       * Compute all the dual costs of a rotation:
       ************************************************/

      double dualCost(cost_);

      /* Working dual cost */
      for(int k=firstDay_; k<firstDay_+length_; ++k)
         dualCost -= workDualCosts[k][shifts_[k]-1];
      /* Start working dual cost */
      dualCost -= startWorkDualCosts[firstDay_];
      /* Stop working dual cost */
      dualCost -= endWorkDualCosts[firstDay_+length_-1];
      /* Working on weekend */
      if(Tools::isSunday(firstDay_))
         dualCost -= workedWeekendDualCost;
      for(int k=firstDay_; k<firstDay_+length_; ++k)
         if(Tools::isSaturday(k))
        	 dualCost -= workedWeekendDualCost;


      // Display: set to true if you want to display the details of the cost

//      if(abs(dualCost_ - dualCost) > EPSILON ){
//    	  cout << "# " << endl;
//    	  cout << "# " << endl;
//    	  cout << "Bad dual cost: " << dualCost_ << " != " << dualCost << endl;
//    	  cout << "# " << endl;
//    	  cout << "#   | Base cost     : + " << cost_ << endl;
//    	  for(int k=firstDay_; k<firstDay_+length_; ++k)
//    		  cout << "#   | Work day-shift: - " << workDualCosts[k][shifts_[k]-1] << endl;
//    	  cout << "#   | Start work    : - " << startWorkDualCosts[firstDay_] << endl;
//    	  cout << "#   | Finish Work   : - " << endWorkDualCosts[firstDay_+length_-1] << endl;
//    	  if(Tools::isSunday(firstDay_))
//    		  cout << "#   | Weekends      : - " << workedWeekendDualCost << endl;
//    	  for(int k=firstDay_; k<firstDay_+length_; ++k)
//    		  if(Tools::isSaturday(k))
//    			  cout << "#   | Weekends      : - " << workedWeekendDualCost << endl;
//    	  std::cout << "#   | ROTATION:" << "  cost=" << cost_ << "  dualCost=" << dualCost_ << "  firstDay=" << firstDay_ << "  length=" << length_ << std::endl;
//    	  std::cout << "#               |";
//    	  vector<int> allTasks (56);
//    	  for(map<int,int>::iterator itTask = shifts_.begin(); itTask != shifts_.end(); ++itTask)
//    		  allTasks[itTask->first] = itTask->second;
//    	  for(int i=0; i<allTasks.size(); i++){
//    		  if(allTasks[i] < 1) std::cout << " |";
//    		  else std::cout << allTasks[i] << "|";
//          }
//    	  cout << "# " << endl;
//    	  cout << "# " << endl;
//    	  //getchar();
//      }
}


//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

MasterProblem::MasterProblem(Scenario* pScenario, Demand* pDemand,
   Preferences* pPreferences, vector<State>* pInitState,
   MySolverType solverType, vector<Roster> solution):

   Solver(pScenario, pDemand, pPreferences, pInitState), solverType_(solverType), pPricer_(0), pRule_(0),
   positionsPerSkill_(pScenario->nbSkills_), skillsPerPosition_(pScenario->nbPositions()), rotations_(pScenario->nbNurses_),

   columnVars_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
   minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
   optDemandVars_(pDemand_->nbDays_), skillsAllocVars_(pDemand_->nbDays_),

   restFlowCons_(pScenario->nbNurses_), workFlowCons_(pScenario->nbNurses_),
   minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
   minDemandCons_(pDemand_->nbDays_), optDemandCons_(pDemand_->nbDays_), feasibleSkillsAllocCons_(pDemand_->nbDays_)
{
   switch(solverType){
   case S_SCIP:
      pModel_ = new BcpModeler(PB_NAME);
      break;
   case S_BCP:
      pModel_ = new BcpModeler(PB_NAME);
      break;
   case S_CBC:
      pModel_ = new CbcModeler(PB_NAME);
   }

   this->preprocessTheNurses();

   /*
    * Build the two vectors linking positions and skills
    */
   for(int p=0; p<skillsPerPosition_.size(); p++){
      vector<int> skills(pScenario->pPositions()[p]->skills_.size());
      for(int sk=0; sk<pScenario->pPositions()[p]->skills_.size(); ++sk)
         skills[sk]=pScenario->pPositions()[p]->skills_[sk];
      skillsPerPosition_[p] = skills;
   }
   for(int sk=0; sk<positionsPerSkill_.size(); sk++){
      vector<int> positions(pScenario->nbPositions());
      int i(0);
      for(int p=0; p<positions.size(); p++)
         if(find(skillsPerPosition_[p].begin(), skillsPerPosition_[p].end(), sk) != skillsPerPosition_[p].end()){
            positions[i]=p;
            ++i;
         }
      positions.resize(i);
      positionsPerSkill_[sk] = positions;
   }

   build();
   initialize(solution);
   pModel_->writeProblem("outfiles/model.lp");
}

MasterProblem::~MasterProblem(){
   delete pPricer_;
   delete pRule_;
   delete pModel_;
}

//build the rostering problem
void MasterProblem::build(){
   /* Rotation constraints */
   buildRotationCons();

   /* Min/Max constraints */
   buildMinMaxCons();

   /* Skills coverage constraints */
   buildSkillsCoverageCons();

   /* Initialize the objects used in the branch and price unless the CBC is used
      to solve the problem
   */
   if (solverType_ != S_CBC) {
     /* Rotation pricer */
     pPricer_ = new RotationPricer(this, "pricer");
     pModel_->addObjPricer(pPricer_);

     /* Branching rule */
     pRule_ = new DiveBranchingRule(this, "branching rule");
     pModel_->addBranchingRule(pRule_);
   }
}

//solve the rostering problem
void MasterProblem::solve(){

  // RqJO: warning, it would be better to define an enumerate type of verbosity
  // levels and create the matching in the Modeler subclasses
  if (solverType_ != S_CBC ) {
    pModel_->setVerbosity(1);
  }
   pModel_->solve();
   pModel_->printStats();
   pModel_->printBestSol();
   storeSolution();
   costsConstrainstsToString();
}
//
//initialize the rostering problem with one column to be feasible if there is no initial solution
//otherwise build the columns corresponding to the initial solution
void MasterProblem::initialize(vector<Roster> solution){
   char* baseName = "initialRotation";

   //if there is no initial solution, we add a column with 1 everywhere for each nurse to be feasible
   if(solution.size() == 0){
      //builg a map of shift -1 everywhere
      map<int,int> shifts;
      for(int k=0; k<pDemand_->nbDays_; ++k)
         shifts.insert(pair<int,int>( k , -1 ));

      for(int i=0; i<pScenario_->nbNurses_; ++i){
         Rotation rotation(shifts, theLiveNurses_[i], bigM);
         addRotation(rotation, baseName);
      }
   }
   //otherwise, rotations are added for each nurse of the initial solution
   else{
      //build the rotations of each nurse
      for(int i=0; i<pScenario_->nbNurses_; ++i){
         //load the roster of nurse i
         Roster roster = solution[i];

         bool workedLastDay = false;
         int lastShift = 0;
         map<int,int> shifts;
         //build all the successive rotation of this nurse
         for(int k=0; k<pDemand_->nbDays_; ++k){
            //shift=0 => rest
            int shift = roster.shift(k);
            //if work, insert the shift in the map
            if(shift>0){
               shifts.insert(pair<int,int>(k, shift));
               lastShift = shift;
               workedLastDay = true;
            }
            else if(shift<0 && lastShift>0){
               shifts.insert(pair<int,int>(k, lastShift));
               workedLastDay = true;
            }
            //if stop to work, build the rotation
            else if(workedLastDay){
               Rotation rotation(shifts, theLiveNurses_[i]);
               rotation.computeCost(pScenario_, pPreferences_, pDemand_->nbDays_);
               addRotation(rotation, baseName);
               shifts.clear();
               lastShift = shift;
               workedLastDay = false;
            }
         }
         //if work on the last day, build the rotation
         if(workedLastDay){
            Rotation rotation(shifts, theLiveNurses_[i]);
            rotation.computeCost(pScenario_, pPreferences_, pDemand_->nbDays_);
            addRotation(rotation, baseName);
            shifts.clear();
         }
      }
   }
}

void MasterProblem::storeSolution(){
   //retrieve a feasible allocation of skills
   vector< vector< vector< vector<double> > > > skillsAllocation(pDemand_->nbDays_);

   for(int k=0; k<pDemand_->nbDays_; ++k){
      vector< vector< vector<double> > > skillsAllocation2(pScenario_->nbShifts_-1);

      for(int s=0; s<pScenario_->nbShifts_-1; ++s){
         vector< vector<double> > skillsAllocation3(pScenario_->nbSkills_);

         for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
            skillsAllocation3[sk] = pModel_->getVarValues(skillsAllocVars_[k][s][sk]);

         skillsAllocation2[s] = skillsAllocation3;
      }
      skillsAllocation[k] = skillsAllocation2;
   }

   //build the rosters
   for(LiveNurse* pNurse: theLiveNurses_)
      for(pair<MyObject*, Rotation> p: rotations_[pNurse->id_])
         if(pModel_->getVarValue(p.first) > 0)
            for(int k=p.second.firstDay_; k<p.second.firstDay_+p.second.length_; ++k){
               bool assigned = false;
               for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
                  if(skillsAllocation[k][p.second.shifts_[k]-1][sk][pNurse->pPosition_->id_] > 0){
                     pNurse->roster_.assignTask(k,p.second.shifts_[k],sk);
                     skillsAllocation[k][p.second.shifts_[k]-1][sk][pNurse->pPosition_->id_] --;
                     assigned = true;
                     break;
                  }
               if(!assigned){
                  char error[255];
                  sprintf(error, "No skill found for Nurse %d on day %d on shift %d", pNurse->id_, k, p.second.shifts_[k]);
                  Tools::throwError((const char*) error);
               }
            }

   //build the states of each nurse
   for(LiveNurse* pNurse: theLiveNurses_)
      pNurse->buildStates();
}

//build the variable of the rotation as well as all the affected constraints with their coefficients
//if s=-1, the nurse i works on all shifts
void MasterProblem::addRotation(Rotation rotation, char* baseName){
	//nurse index
	int i = rotation.pNurse_->id_;

	//Column var, its name, and affected constraints with their coefficients
	MyObject* var;
	char name[255];
	vector<MyObject*> cons;
	vector<double> coeffs;

	/* Rotation constraints */
	addRotationConsToCol(&cons, &coeffs, i, rotation.firstDay_, true, false);
	addRotationConsToCol(&cons, &coeffs, i, rotation.firstDay_+rotation.length_-1, false, true);

	/* Min/Max constraints */
	for(int k=rotation.firstDay_; k<rotation.firstDay_+rotation.length_; ++k){
		//check if the nurse works on a saturday and add it in the constraints
		//not yet added in the constraints maxWorkedWeekend
		if(Tools::isSaturday(k))
			addMinMaxConsToCol(&cons, &coeffs, i, k, true);
		//check if the nurse works on a sunday and does not work on a saturday
		//not yet added in the constraints maxWorkedWeekend
		else if( (k==rotation.firstDay_) && Tools::isSunday(k) )
			addMinMaxConsToCol(&cons, &coeffs, i, k, true);
		//otherwise, do not add in the constraints maxWorkedWeekend
		else
			addMinMaxConsToCol(&cons, &coeffs, i, k, false);
	}

	/* Skills coverage constraints */
	for(int k=rotation.firstDay_; k<rotation.firstDay_+rotation.length_; ++k)
		addSkillsCoverageConsToCol(&cons, &coeffs, i, k, rotation.shifts_[k]);

	sprintf(name, "%s_N%d_%d",baseName , i, rotations_[i].size());
	pModel_->createIntColumn(&var, name, rotation.cost_, rotation.dualCost_, cons, coeffs);
	rotations_[i].insert(pair<MyObject*,Rotation>(var, rotation));
}

/*
 * Rotation constraints
 */
void MasterProblem::buildRotationCons(){
   char name[255];
   //build the rotation network for each nurse
   for(int i=0; i<pScenario_->nbNurses_; i++){
      int minConsDaysOff(theLiveNurses_[i]->minConsDaysOff()),
         maxConsDaysOff(theLiveNurses_[i]->maxConsDaysOff()),
         initConsDaysOff(theLiveNurses_[i]->pStateIni_->consDaysOff_);
      //=true if we have to compute a cost for resting days exceeding the maximum allowed
      //=false otherwise
      bool const maxRest = (maxConsDaysOff < pDemand_->nbDays_ + initConsDaysOff);
      //number of long resting arcs as function of maxRest
      int const nbLongRestingArcs((maxRest) ? maxConsDaysOff : minConsDaysOff);
      //first day when a rest arc exists =
      //nbLongRestingArcs - number of consecutive worked days in the past
      int const firstRestArc( min( max( 0, nbLongRestingArcs - initConsDaysOff ), pDemand_->nbDays_-1 ) );
      //first day when a restingVar exists: at minimun 1
      //if firstRestArc=0, the first resting arc is a longRestingVar
      int const indexStartRestArc = max(1, firstRestArc);
      //number of resting arcs
      int const nbRestingArcs( pDemand_->nbDays_- indexStartRestArc );

      //initialize vectors
      vector< MyObject* > restingVars2(nbRestingArcs);
      vector< vector<MyObject*> > longRestingVars2(pDemand_->nbDays_);
      vector<MyObject*> restFlowCons2(pDemand_->nbDays_);
      vector<MyObject*> workFlowCons2(pDemand_->nbDays_);

      /*****************************************
       * Creating arcs
       *****************************************/
      for(int k=0; k<pDemand_->nbDays_; ++k){
         /*****************************************
          * first long resting arcs
          *****************************************/
         if(k==0){
            //number of min long resting arcs
            int nbMinRestArcs( max(0, minConsDaysOff - initConsDaysOff) );
            //initialize cost
            int cost (nbMinRestArcs * WEIGHT_CONS_DAYS_OFF);
            Rotation rot = computeInitStateRotation(theLiveNurses_[i]);

            //initialize vectors
            //Must have a minimum of one long resting arcs
            vector<MyObject*> longRestingVars3_0(indexStartRestArc);

            //create minRest arcs
            for(int l=1; l<=nbMinRestArcs; ++l){
               cost -= WEIGHT_CONS_DAYS_OFF;
               sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
               pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, cost+rot.cost_);
               rotations_[i].insert(pair<MyObject*,Rotation>(longRestingVars3_0[l-1], rot));
            }

            //create maxRest arcs, if maxRest=true
            if(maxRest){
               for(int l=1+nbMinRestArcs; l<=firstRestArc; ++l){
                  sprintf(name, "longRestingVars_N%d_%d_%d", i, 0, l);
                  pModel_->createPositiveVar(&longRestingVars3_0[l-1], name, rot.cost_);
                  rotations_[i].insert(pair<MyObject*,Rotation>(longRestingVars3_0[l-1], rot));
               }
            }

            //create the only resting arc (same as a short resting arcs)
            if(firstRestArc == 0){
               sprintf(name, "restingVars_N%d_%d_%d", i, 0, 1);
               pModel_->createPositiveVar(&longRestingVars3_0[0], name, (maxRest) ? WEIGHT_CONS_DAYS_OFF+rot.cost_ : rot.cost_);
               rotations_[i].insert(pair<MyObject*,Rotation>(longRestingVars3_0[0], rot));
            }
            //store vectors
            longRestingVars2[0] = longRestingVars3_0;
         }
         /*****************************************
          * long resting arcs without the first ones
          *****************************************/
         else{
            //number of long resting arcs = min(nbLongRestingArcs, number of possible long resting arcs)
            int nbLongRestingArcs2( min(nbLongRestingArcs, pDemand_->nbDays_-k) );
            //initialize cost
            //if the arc finishes the last day, the cost is 0. Indeed it will be computed on the next planning
            int cost = minConsDaysOff * WEIGHT_CONS_DAYS_OFF;

            //initialize vectors
            vector<MyObject*> longRestingVars3(nbLongRestingArcs2);

            //create minRest arcs
            for(int l=1; l<=minConsDaysOff; ++l){
               cost -= WEIGHT_CONS_DAYS_OFF;
               sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
               //if arc ends before the last day: normal cost
               if(l < pDemand_->nbDays_-k)
                  pModel_->createPositiveVar(&longRestingVars3[l-1], name, cost);
               //otherwise, arc finishes on last day
               //so: cost=0 and we break the loop
               else{
                  pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
                  break;
               }
            }
            //create maxRest arcs, if maxRest=true
            if(maxRest)
               for(int l=1+minConsDaysOff; l<=maxConsDaysOff; ++l){
                  //if exceed last days, break
                  if(l > pDemand_->nbDays_-k)
                     break;
                  sprintf(name, "longRestingVars_N%d_%d_%d", i, k, k+l);
                  pModel_->createPositiveVar(&longRestingVars3[l-1], name, 0);
               }
            //store vectors
            longRestingVars2[k] = longRestingVars3;
         }
         /*****************************************
          * short resting arcs
          *****************************************/
         if(k>=indexStartRestArc){
            sprintf(name, "restingVars_N%d_%d_%d", i, k, k+1);
            pModel_->createPositiveVar(&restingVars2[k-indexStartRestArc], name, (maxRest) ? WEIGHT_CONS_DAYS_OFF : 0);
         }
      }

      /*****************************************
       * Resting nodes constraints
       *****************************************/
      for(int k=0; k<pDemand_->nbDays_; ++k){
         vector<double> coeffs(longRestingVars2[k].size());
         for(int l=0; l<longRestingVars2[k].size(); ++l)
            coeffs[l] = 1;
         sprintf(name, "restingNodes_N%d_%d", i, k);
         //Create flow constraints. out flow = 1 if source node (k=0)
         pModel_->createEQConsLinear(&restFlowCons2[k], name, (k==0) ? 1 : 0,
            longRestingVars2[k], coeffs);
      }

      /*****************************************
       * Working nodes constraints
       *****************************************/
      for(int k=1; k<=pDemand_->nbDays_; ++k){
         //take the min between the number of long resting arcs and the number of possible in arcs
         int nbLongRestingArcs2 = min(nbLongRestingArcs,k);

         vector<MyObject*> vars;
         vector<double> coeffs;
         //add long resting arcs
         for(int l=0; l<nbLongRestingArcs2; ++l){
            //if the long resting arc starts on the source node,
            //check if there exists such an arc
            if( (l > k-1) || (l >= longRestingVars2[k-1-l].size()) )
               break;
            vars.push_back(longRestingVars2[k-1-l][l]);
            //compute in-flow for the sink
            if(k==pDemand_->nbDays_)
               coeffs.push_back(1);
            //compute out-flow
            else
               coeffs.push_back(-1);
         }
         //add resting arcs
         //just 1 out, if first restingVar
         //compute out-flow
         if(k==indexStartRestArc){
            vars.push_back(restingVars2[0]);
            coeffs.push_back(1);
         }
         //just 1 in, if last resting arcs
         //compute in-flow for the sink
         else if (k==pDemand_->nbDays_){
            vars.push_back(restingVars2[restingVars2.size()-1]);
            coeffs.push_back(1);
         }
         //2 otherwise: 1 in and 1 out
         //compute out-flow
         else if(k>indexStartRestArc){
            vars.push_back(restingVars2[k-1-indexStartRestArc]);
            coeffs.push_back(-1);
            vars.push_back(restingVars2[k-indexStartRestArc]);
            coeffs.push_back(1);
         }
         sprintf(name, "workingNodes_N%d_%d", i, k);
         //Create flow constraints. in flow = 1 if sink node (k==pDemand_->nbDays_)
         pModel_->createEQConsLinear(&workFlowCons2[k-1], name, (k==pDemand_->nbDays_) ? 1 : 0,
            vars, coeffs);
      }

      //store vectors
      restingVars_[i] = restingVars2;
      longRestingVars_[i] = longRestingVars2;
      restFlowCons_[i] = restFlowCons2;
      workFlowCons_[i] = workFlowCons2;
   }
}

int MasterProblem::addRotationConsToCol(vector<MyObject*>* cons, vector<double>* coeffs, int i, int k, bool firstDay, bool lastDay){
   //check if the rotation starts on day k
   if(firstDay){
      //compute out-flow
      coeffs->push_back(1.0);
      //add to source constraint
      if(k==0)
         cons->push_back(restFlowCons_[i][0]);
      //add to work node constraint
      else
         cons->push_back(workFlowCons_[i][k-1]);

      return 1;
   }

   //check if the rotation finishes on day k
   else if(lastDay){
      //add to sink constraint
      //compute in-flow
      if(k==pDemand_->nbDays_-1){
         coeffs->push_back(1.0);
         cons->push_back(workFlowCons_[i][pDemand_->nbDays_-1]);
      }
      //add to rest node constraint
      //compute out-flow
      else{
         coeffs->push_back(-1.0);
         cons->push_back(restFlowCons_[i][k+1]);
      }

      return 1;
   }

   return 0;
}

/*
 * Min/Max constraints
 */
void MasterProblem::buildMinMaxCons(){
   char name[255];
   for(int i=0; i<pScenario_->nbNurses_; i++){
      sprintf(name, "minWorkedDaysVar_N%d", i);
      pModel_->createPositiveVar(&minWorkedDaysVars_[i], name, WEIGHT_TOTAL_SHIFTS);
      sprintf(name, "maxWorkedDaysVar_N%d", i);
      pModel_->createPositiveVar(&maxWorkedDaysVars_[i], name, WEIGHT_TOTAL_SHIFTS);
      sprintf(name, "maxWorkedWeekendVar_N%d", i);
      pModel_->createPositiveVar(&maxWorkedWeekendVars_[i], name, WEIGHT_TOTAL_WEEKENDS);

      sprintf(name, "minWorkedDaysCons_N%d", i);
      vector<MyObject*> vars1 = {minWorkedDaysVars_[i]};
      vector<double> coeffs1 = {1};
      pModel_->createGEConsLinear(&minWorkedDaysCons_[i], name, theLiveNurses_[i]->minTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_,
         vars1, coeffs1);

      sprintf(name, "maxWorkedDaysCons_N%d", i);
      vector<MyObject*> vars2 = {maxWorkedDaysVars_[i]};
      vector<double> coeffs2 = {-1};
      pModel_->createLEConsLinear(&maxWorkedDaysCons_[i], name, theLiveNurses_[i]->maxTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_,
         vars2, coeffs2);

      sprintf(name, "maxWorkedWeekendCons_N%d", i);
      vector<MyObject*> vars3 = {maxWorkedWeekendVars_[i]};
      vector<double> coeffs3 = {-1};
      pModel_->createLEConsLinear(&maxWorkedWeekendCons_[i], name, theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
         vars3, coeffs3);
   }
}

int MasterProblem::addMinMaxConsToCol(vector<MyObject*>* cons, vector<double>* coeffs, int i, int k, bool weekend){
   int nbCons(0);

   ++nbCons;
   cons->push_back(minWorkedDaysCons_[i]);
   coeffs->push_back(1.0);
   ++nbCons;
   cons->push_back(maxWorkedDaysCons_[i]);
   coeffs->push_back(1.0);
   if(weekend){
      ++nbCons;
      cons->push_back(maxWorkedWeekendCons_[i]);
      coeffs->push_back(1.0);
   }

   return nbCons;
}

/*
 * Skills coverage constraints
 */
void MasterProblem::buildSkillsCoverageCons(){
   char name[255];
   for(int k=0; k<pDemand_->nbDays_; k++){
      //initialize vectors
      vector< vector<MyObject*> > optDemandVars1(pScenario_->nbShifts_-1);
      vector< vector< vector<MyObject*> > > skillsAllocVars1(pScenario_->nbShifts_-1);
      vector< vector<MyObject*> > minDemandCons1(pScenario_->nbShifts_-1);
      vector< vector<MyObject*> > optDemandCons1(pScenario_->nbShifts_-1);
      vector< vector<MyObject*> > feasibleSkillsAllocCons1(pScenario_->nbShifts_-1);

      //forget s=0, it's a resting shift
      for(int s=1; s<pScenario_->nbShifts_; s++){
         //initialize vectors
         vector<MyObject*> optDemandVars2(pScenario_->nbSkills_);
         vector< vector<MyObject*> > skillsAllocVars2(pScenario_->nbSkills_);
         vector<MyObject*> minDemandCons2(pScenario_->nbSkills_);
         vector<MyObject*> optDemandCons2(pScenario_->nbSkills_);
         vector<MyObject*> feasibleSkillsAllocCons2(pScenario_->nbPositions());

         for(int sk=0; sk<pScenario_->nbSkills_; sk++){
            //initialize vectors
            vector<MyObject*> skillsAllocVars3(pScenario_->nbPositions());

            //create variables
            sprintf(name, "optDemandVar_%d_%d_%d", k, s, sk);
            pModel_->createPositiveVar(&optDemandVars2[sk], name, WEIGHT_OPTIMAL_DEMAND);
            for(int p=0; p<pScenario_->nbPositions(); p++){
               sprintf(name, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
               pModel_->createIntVar(&skillsAllocVars3[p], name, 0);
            }
            //store vectors
            skillsAllocVars2[sk] = skillsAllocVars3;

            //adding variables and building minimum demand constraints
            vector<MyObject*> vars1(positionsPerSkill_[sk].size());
            vector<double> coeffs1(positionsPerSkill_[sk].size());
            for(int p=0; p<positionsPerSkill_[sk].size(); ++p){
               vars1[p] = skillsAllocVars3[positionsPerSkill_[sk][p]];
               coeffs1[p] = 1;
            }
            sprintf(name, "minDemandCons_%d_%d_%d", k, s, sk);
            pModel_->createFinalGEConsLinear(&minDemandCons2[sk], name, pDemand_->minDemand_[k][s][sk],
               vars1, coeffs1);

            //adding variables and building optimal demand constraints
            vars1.push_back(optDemandVars2[sk]);
            coeffs1.push_back(1);
            sprintf(name, "optDemandCons_%d_%d_%d", k, s, sk);
            pModel_->createFinalGEConsLinear(&optDemandCons2[sk], name, pDemand_->optDemand_[k][s][sk],
               vars1, coeffs1);
         }

         for(int p=0; p<pScenario_->nbPositions(); p++){
            //adding variables and building skills allocation constraints
            int const nonZeroVars3(skillsPerPosition_[p].size());
            vector<MyObject*> vars3(nonZeroVars3);
            vector<double> coeff3(nonZeroVars3);
            for(int sk=0; sk<nonZeroVars3; ++sk){
               vars3[sk] = skillsAllocVars2[skillsPerPosition_[p][sk]][p];
               coeff3[sk] =-1;
            }
            sprintf(name, "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
            pModel_->createEQConsLinear(&feasibleSkillsAllocCons2[p], name, 0,
               vars3, coeff3);
         }

         //store vectors
         optDemandVars1[s-1] = optDemandVars2;
         skillsAllocVars1[s-1] = skillsAllocVars2;
         minDemandCons1[s-1] = minDemandCons2;
         optDemandCons1[s-1] = optDemandCons2;
         feasibleSkillsAllocCons1[s-1] = feasibleSkillsAllocCons2;
      }

      //store vectors
      optDemandVars_[k] = optDemandVars1;
      skillsAllocVars_[k] = skillsAllocVars1;
      minDemandCons_[k] = minDemandCons1;
      optDemandCons_[k] = optDemandCons1;
      feasibleSkillsAllocCons_[k] = feasibleSkillsAllocCons1;
   }
}

int MasterProblem::addSkillsCoverageConsToCol(vector<MyObject*>* cons, vector<double>* coeffs, int i, int k, int s){
   int nbCons(0);

   int p(theLiveNurses_[i]->pPosition_->id_);
   if(s==-1){
      for(int s0=1; s0<pScenario_->nbShifts_; ++s0){
         ++nbCons;
         cons->push_back(feasibleSkillsAllocCons_[k][s0-1][p]);
         coeffs->push_back(1.0);
      }
   }
   else{
      ++nbCons;
      cons->push_back(feasibleSkillsAllocCons_[k][s-1][p]);
      coeffs->push_back(1.0);
   }

   return nbCons;
}

string MasterProblem::costsConstrainstsToString(){
   stringstream rep;

   double initStateRestCost = getRotationCosts(INIT_REST_COST);
   char buffer[100];
   sprintf(buffer, "%-30s %10.0f \n", "Column costs", getRotationCosts() - initStateRestCost);
   rep << buffer;
   rep << "-----------------------------------------\n";
   sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. shifts costs", getRotationCosts(CONS_SHIFTS_COST));
   rep << buffer;
   sprintf(buffer, "%5s%-25s %10.0f \n", "", "Cons. worked days costs", getRotationCosts(CONS_WORKED_DAYS_COST));
   rep << buffer;
   sprintf(buffer, "%5s%-25s %10.0f \n", "", "Complete weekend costs", getRotationCosts(COMPLETE_WEEKEND_COST));
   rep << buffer;
   sprintf(buffer, "%5s%-25s %10.0f \n", "", "Preferences costs", getRotationCosts(PREFERENCE_COST));
   rep << buffer;
   rep << "-----------------------------------------\n";
   double initStateCost = getRotationCosts(TOTAL_COST, true);
   sprintf(buffer, "%-30s %10.0f \n", "History costs (counted)", initStateCost + initStateRestCost);
   rep << buffer;
   sprintf(buffer, "%-30s %10.0f \n", "Resting costs", pModel_->getTotalCost(restingVars_)+pModel_->getTotalCost(longRestingVars_) + initStateRestCost - initStateCost);
   rep << buffer;
   sprintf(buffer, "%-30s %10.0f \n", "Min worked days costs", pModel_->getTotalCost(minWorkedDaysVars_));
   rep << buffer;
   sprintf(buffer, "%-30s %10.0f \n", "Max worked days costs", pModel_->getTotalCost(maxWorkedDaysVars_));
   rep << buffer;
   sprintf(buffer, "%-30s %10.0f \n", "Max worked weekend costs", pModel_->getTotalCost(maxWorkedWeekendVars_));
   rep << buffer;
   sprintf(buffer, "%-30s %10.0f \n", "Coverage costs", pModel_->getTotalCost(optDemandVars_));
   rep << buffer;
   rep << "-----------------------------------------\n";
   rep << "\n";

   cout << rep.str();

   return rep.str();
}

double MasterProblem::getRotationCosts(CostType costType, bool initStateRotation){
   double cost = 0;
   for(map<MyObject*, Rotation> map0: rotations_){
      for(pair<MyObject*, Rotation> rot: map0){
         //if(initStateRotation), search for empty rotation (=rotation for initial state)
         if(initStateRotation && rot.second.shifts_.size()>0)
            continue;
         double value = pModel_->getVarValue(rot.first);
         switch(costType){
         case CONS_SHIFTS_COST: cost += rot.second.consShiftsCost_*value;
         break;
         case CONS_WORKED_DAYS_COST: cost += rot.second.consDaysWorkedCost_*value;
         break;
         case COMPLETE_WEEKEND_COST: cost += rot.second.completeWeekendCost_*value;
         break;
         case PREFERENCE_COST: cost += rot.second.preferenceCost_*value;
         break;
         case INIT_REST_COST: cost += rot.second.initRestCost_*value;
         break;
         default: cost += rot.second.cost_*value;
         break;
         }
      }
   }
   return cost;
}

Rotation MasterProblem::computeInitStateRotation(LiveNurse* pNurse){
   //initialize rotation
   Rotation rot({}, pNurse);

   //compute cost for previous cons worked shifts and days
   int lastShift = pNurse->pStateIni_->shift_;
   if(lastShift>0){
      int nbConsWorkedDays = pNurse->pStateIni_->consDaysWorked_;
      int diff = pNurse->minConsDaysWork() - nbConsWorkedDays;
      rot.consDaysWorkedCost_ += (diff>0) ? diff*WEIGHT_CONS_DAYS_WORK : 0;

      int nbConsShifts = pNurse->pStateIni_->consShifts_;
      int diff2 = pScenario_->minConsShifts_[lastShift] - nbConsShifts;
      rot.consShiftsCost_ += (diff2>0) ? diff2*WEIGHT_CONS_SHIFTS : 0;
   }
   rot.cost_ = rot.consDaysWorkedCost_ + rot.consShiftsCost_;

   return rot;
}
