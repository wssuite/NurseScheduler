/*
 * MasterProblem.cpp
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#include "MasterProblem.h"
#include "RotationPricer.h"

/* namespace usage */
using namespace std;
using namespace scip;

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------


MasterProblem::MasterProblem(Scenario* pScenario, Demand* pDemand,
   Preferences* pPreferences, vector<State>* pInitState, vector<Roster> solution):
                                                                                                                  Solver(pScenario, pDemand, pPreferences, pInitState),
                                                                                                                  scip_(PB_NAME), positionsPerSkill_(pScenario->nbSkills_), skillsPerPosition_(pScenario->nbPositions()), rotations_(pScenario->nbNurses_),

                                                                                                                  columnVars_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
                                                                                                                  minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
                                                                                                                  optDemandVars_(pDemand_->nbDays_), skillsAllocVars_(pDemand_->nbDays_),

                                                                                                                  restFlowCons_(pScenario->nbNurses_), workFlowCons_(pScenario->nbNurses_),
                                                                                                                  minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
                                                                                                                  minDemandCons_(pDemand_->nbDays_), optDemandCons_(pDemand_->nbDays_), feasibleSkillsAllocCons_(pDemand_->nbDays_)
{

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
   scip_.writeProblem("outfiles/model.lp");
}

MasterProblem::~MasterProblem(){ }

//build the rostering problem
void MasterProblem::build(){
   /* Rotation constraints */
   buildRotationCons();

   /* Min/Max constraints */
   buildMinMaxCons();

   /* Skills coverage constraints */
   buildSkillsCoverageCons();

   /* Rotation pricer */
   pPricer_ = new RotationPricer(this, "pricer");
   scip_.addObjPricer(pPricer_);
}

//solve the rostering problem
void MasterProblem::solve(){
   scip_.solve();
   //   scip_.printStats();
   scip_.printBestSol();
   storeSolution();
   scip_.deleteSCIP();
}

//initialize the rostering problem with one column to be feasible if there is no initial solution
//otherwise build the columns corresponding to the initial solution
void MasterProblem::initialize(vector<Roster> solution){
   int nbVars(0);
   vector<SCIP_VAR*> initialVar(pScenario_->nbNurses_);
   const char* baseName = "initialRotation";
   vector<char*> names(pScenario_->nbNurses_);
   vector<double> initialObjCoeffs(pScenario_->nbNurses_);
   vector<int> nbCons(pScenario_->nbNurses_);
   vector< vector<SCIP_CONS*> > initialCons(pScenario_->nbNurses_);
   vector< vector<double> > initialCoeffs(pScenario_->nbNurses_);

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
   //retieve best solution
   SCIP_SOL* sol = scip_.getBestSol();

   //retrieve a feasible allocation of skills
   vector< vector< vector< vector<double> > > > skillsAllocation(pDemand_->nbDays_);

   for(int k=0; k<pDemand_->nbDays_; ++k){
      vector< vector< vector<double> > > skillsAllocation2(pScenario_->nbShifts_-1);

      for(int s=0; s<pScenario_->nbShifts_-1; ++s){
         vector< vector<double> > skillsAllocation3(pScenario_->nbSkills_);

         for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
            skillsAllocation3[sk] = scip_.getVarValues(sol, skillsAllocVars_[k][s][sk]);

         skillsAllocation2[s] = skillsAllocation3;
      }
      skillsAllocation[k] = skillsAllocation2;
   }

   //build the rosters
   for(LiveNurse* pNurse: theLiveNurses_)
      for(pair<SCIP_VAR*, Rotation> p: rotations_[pNurse->id_])
         if(scip_.getVarValue(sol, p.first) > 0)
            for(int k=p.second.firstDay_; k<p.second.firstDay_+p.second.length_; ++k)
               for(int sk=0; sk<pScenario_->nbSkills_; ++sk)
                  if(skillsAllocation[k][p.second.shifts_[k]-1][sk][pNurse->pPosition_->id_] > 0){
                     pNurse->roster_.assignTask(k,p.second.shifts_[k],sk);
                     skillsAllocation[k][p.second.shifts_[k]-1][sk][pNurse->pPosition_->id_] --;
                  }
}

//build the variable of the rotation as well as all the affected constraints with their coefficients
//if s=-1, the nurse i works on all shifts
void MasterProblem::addRotation(Rotation rotation, const char* baseName){
   int nbCons(0);
   //nurse index
   int i = rotation.pNurse_->id_;

   //SCIP column var, its name, and affected constraints with their coefficients
   SCIP_VAR* var;
   char name[255];
   vector<SCIP_CONS*> cons;
   vector<double> coeffs;

   /* Rotation constraints */
   nbCons += addRotationConsToCol(&cons, &coeffs, i, rotation.firstDay_, true, false);
   nbCons += addRotationConsToCol(&cons, &coeffs, i, rotation.firstDay_+rotation.length_-1, false, true);

   /* Min/Max constraints */
   for(int k=rotation.pNurse_->firstDay_; k<rotation.pNurse_->firstDay_+rotation.length_; ++k){
      //check if the nurse works on a saturday and add it in the constraints
      //not yet added in the constraints maxWorkedWeekend
      if(Tools::isSaturday(k))
         nbCons += addMinMaxConsToCol(&cons, &coeffs, i, k, true);
      //check if the nurse works on a sunday and does not work on a saturday
      //not yet added in the constraints maxWorkedWeekend
      else if( (k==rotation.pNurse_->firstDay_) && Tools::isSunday(k) )
         nbCons += addMinMaxConsToCol(&cons, &coeffs, i, k, true);
      //otherwise, do not add in the constraints maxWorkedWeekend
      else
         nbCons += addMinMaxConsToCol(&cons, &coeffs, i, k);
   }

   /* Skills coverage constraints */
   for(int k=rotation.firstDay_; k<rotation.firstDay_+rotation.length_; ++k)
      nbCons += addSkillsCoverageConsToCol(&cons, &coeffs, i, k, rotation.shifts_[k]);

   SCIPsnprintf(name, 255, "%s_N%d_%d",baseName , i, rotations_[i].size());
   scip_.createPositiveColumn(&var, name, rotation.cost_,
      nbCons, &(cons)[0], &(coeffs)[0]);
   rotations_[i].insert(pair<SCIP_VAR*,Rotation>(var, rotation));
}

/*
 * Rotation constraints
 */
void MasterProblem::buildRotationCons(){
   char name[255];
   //build the rotation network for each nurse
   for(int i=0; i<pScenario_->nbNurses_; i++){
      int minConsDaysOff(theLiveNurses_[i]->minConsDaysOff()),
         maxConsDaysOff(theLiveNurses_[i]->maxConsDaysOff());
      //=true if we have to compute a cost for resting days exceeding the maximum allowed
      //=false otherwise
      bool const maxRest = (maxConsDaysOff <= pDemand_->nbDays_);
      //number of long resting arcs as function of maxRest
      int const nbLongRestingArcs((maxRest) ? maxConsDaysOff : minConsDaysOff);
      //first day when a rest arc exists =
      //nbLongRestingArcs - number of consecutive worked days in the past
      int const firstRestArc(max(0, nbLongRestingArcs - theLiveNurses_[i]->pStateIni_->consDaysOff_));
      //number of resting arcs
      int const nbRestingArcs(max(pDemand_->nbDays_-firstRestArc, pDemand_->nbDays_-1));

      //initialize vectors
      vector< SCIP_VAR* > restingVars2(nbRestingArcs);
      vector< vector<SCIP_VAR*> > longRestingVars2(pDemand_->nbDays_);
      vector<SCIP_CONS*> restFlowCons2(pDemand_->nbDays_);
      vector<SCIP_CONS*> workFlowCons2(pDemand_->nbDays_);

      /*****************************************
       * short resting arcs
       *****************************************/
      for(int k=0; k<nbRestingArcs; ++k){
         SCIPsnprintf(name, 255, "restingVars_N%d_%d_%d", i, firstRestArc+k+1, firstRestArc+k+2);
         scip_.createPositiveVar(&(restingVars2[k]), name, (maxRest) ? WEIGHT_CONS_DAYS_OFF : 0);
      }

      /*****************************************
       * long resting arcs without the first ones
       *****************************************/
      for(int k=1; k<pDemand_->nbDays_; ++k){
         //number of long resting arcs = min(nbLongRestingArcs, number of possible long resting arcs)
         int nbLongRestingArcs2( min(nbLongRestingArcs, pDemand_->nbDays_-k) );
         //initialize cost
         //if the arc finishes the last day, the cost is 0. Indeed it will be computed on the next planning
         int cost = minConsDaysOff * WEIGHT_CONS_DAYS_OFF;

         //initialize vectors
         vector<SCIP_VAR*> longRestingVars3(nbLongRestingArcs2);

         //create minRest arcs
         for(int l=1; l<=minConsDaysOff; ++l){
            cost -= WEIGHT_CONS_DAYS_OFF;
            SCIPsnprintf(name, 255, "longRestingVars_N%d_%d_%d", i, k, k+l);
            //if arc ends before the last day: normal cost
            if(l < pDemand_->nbDays_-k)
               scip_.createPositiveVar(&(longRestingVars3[l-1]), name, cost);
            //otherwise, arc finishes on last day
            //so: cost=0 and we break the loop
            else{
               scip_.createPositiveVar(&(longRestingVars3[l-1]), name, 0);
               break;
            }
         }

         //create maxRest arcs, if maxRest=true
         if(maxRest)
            for(int l=1+minConsDaysOff; l<=maxConsDaysOff; ++l){
               //if exceed last days, break
               if(l > pDemand_->nbDays_-k)
                  break;
               SCIPsnprintf(name, 255, "longRestingVars_N%d_%d_%d", i, k, k+l);
               scip_.createPositiveVar(&(longRestingVars3[l-1]), name, 0);
            }
         //store vectors
         longRestingVars2[k] = longRestingVars3;
      }

      /*****************************************
       * first long resting arcs
       *****************************************/
      //number of min resting arcs
      int nbMinRestArcs( max(0, minConsDaysOff - theLiveNurses_[i]->pStateIni_->consDaysOff_) );
      //initialize cost
      int cost (nbMinRestArcs * WEIGHT_CONS_DAYS_OFF);

      //initialize vectors
      //Must have a minimum of one long resting arcs
      vector<SCIP_VAR*> longRestingVars3_0(max(1,firstRestArc));

      //create minRest arcs
      for(int l=1; l<=nbMinRestArcs; ++l){
         cost -= WEIGHT_CONS_DAYS_OFF;
         SCIPsnprintf(name, 255, "longRestingVars_N%d_%d_%d", i, 0, l);
         scip_.createPositiveVar(&(longRestingVars3_0[l-1]), name, cost);
      }

      //create maxRest arcs, if maxRest=true
      if(maxRest){
         for(int l=1+nbMinRestArcs; l<=firstRestArc; ++l){
            SCIPsnprintf(name, 255, "longRestingVars_N%d_%d_%d", i, 0, l);
            scip_.createPositiveVar(&(longRestingVars3_0[l-1]), name, 0);
         }
      }

      //create the only resting arc (same as a short resting arcs)
      if(firstRestArc == 0){
         SCIPsnprintf(name, 255, "restingVars_N%d_%d_%d", i, 0, 1);
         scip_.createPositiveVar(&(longRestingVars3_0[0]), name, (maxRest) ? WEIGHT_CONS_DAYS_OFF : 0);
      }
      //store vectors
      longRestingVars2[0] = longRestingVars3_0;

      /*****************************************
       * Resting nodes constraints
       *****************************************/
      for(int k=0; k<pDemand_->nbDays_; ++k){
         vector<double> coeff(longRestingVars2[k].size());
         for(int l=0; l<longRestingVars2[k].size(); ++l)
            coeff[l] = 1;
         SCIPsnprintf(name, 255, "restingNodes_N%d_%d", i, k);
         //Create flow constraints. out flow = 1 if source node (k=0)
         scip_.createEQConsLinear(&(restFlowCons2[k]), name, (k==0) ? 1 : 0,
            longRestingVars2[k].size(), &(longRestingVars2[k])[0], &coeff[0]);
      }

      /*****************************************
       * Working nodes constraints
       *****************************************/
      for(int k=0; k<pDemand_->nbDays_; ++k){
         int restingArcs(0);
         //add resting arcs, if day k >= first resting arcs
         //just 1 arc if this is the first node or the last
         //2 otherwise
         if(k==firstRestArc || k==pDemand_->nbDays_-1)
            restingArcs = 1;
         else if(k>firstRestArc)
            restingArcs = 2;
         //take the min between the number of long resting arcs and the number of possible in arcs
         int nbLongRestingArcs2 = min(nbLongRestingArcs,k+1);

         vector<SCIP_VAR*> vars;
         vector<double> coeff;
         //add long resting arcs
         for(int l=0; l<nbLongRestingArcs2; ++l){
            //if the long resting arc starts on the source node,
            //check if there exists such an arc
            if( (l == k) && (l >= longRestingVars2[k-l].size()) )
               break;
            vars.push_back(longRestingVars2[k-l][l]);
            //compute in-flow for the sink
            if(k==pDemand_->nbDays_-1)
               coeff.push_back(1);
            //compute out-flow
            else
               coeff.push_back(-1);
         }
         //add resting arcs
         //just 1 out, if first resting arcs
         //compute out-flow
         if(k==firstRestArc){
            vars.push_back(restingVars2[k-firstRestArc]);
            coeff.push_back(1);
         }
         //just 1 in, if last resting arcs
         //compute in-flow for the sink
         else if (k==pDemand_->nbDays_-1){
            vars.push_back(restingVars2[k-1-firstRestArc]);
            coeff.push_back(1);
         }
         //2 otherwise: 1 in and 1 out
         //compute out-flow
         else if(k>firstRestArc){
            vars.push_back(restingVars2[k-1-firstRestArc]);
            coeff.push_back(-1);
            vars.push_back(restingVars2[k-firstRestArc]);
            coeff.push_back(1);
         }
         SCIPsnprintf(name, 255, "workingNodes_N%d_%d", i, k+1);
         //Create flow constraints. in flow = 1 if sink node (k==pDemand_->nbDays_-1)
         scip_.createEQConsLinear(&(workFlowCons2[k]), name, (k==pDemand_->nbDays_-1) ? 1 : 0,
            vars.size(), &vars[0], &coeff[0]);
      }

      //store vectors
      restingVars_[i] = restingVars2;
      longRestingVars_[i] = longRestingVars2;
      restFlowCons_[i] = restFlowCons2;
      workFlowCons_[i] = workFlowCons2;
   }
}

int MasterProblem::addRotationConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, bool firstDay, bool lastDay){
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
      SCIPsnprintf(name, 255, "minWorkedDaysVar_N%d", i);
      scip_.createPositiveVar(&(minWorkedDaysVars_[i]), name, WEIGHT_TOTAL_SHIFTS);
      SCIPsnprintf(name, 255, "maxWorkedDaysVar_N%d", i);
      scip_.createPositiveVar(&(maxWorkedDaysVars_[i]), name, WEIGHT_TOTAL_SHIFTS);
      SCIPsnprintf(name, 255, "maxWorkedWeekendVar_N%d", i);
      scip_.createPositiveVar(&(maxWorkedWeekendVars_[i]), name, WEIGHT_TOTAL_WEEKENDS);

      double coeff1(1);
      SCIPsnprintf(name, 255, "minWorkedDaysCons_N%d", i);
      scip_.createGEConsLinear(&(minWorkedDaysCons_[i]), name, theLiveNurses_[i]->minTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_,
         1, &minWorkedDaysVars_[i], &coeff1);
      double coeff2(-1);
      SCIPsnprintf(name, 255, "maxWorkedDaysCons_N%d", i);
      scip_.createLEConsLinear(&(maxWorkedDaysCons_[i]), name, theLiveNurses_[i]->maxTotalShifts() - theLiveNurses_[i]->pStateIni_->totalDaysWorked_,
         1, &maxWorkedDaysVars_[i], &coeff2);
      double coeff3(-1);
      SCIPsnprintf(name, 255, "maxWorkedWeekendCons_N%d", i);
      scip_.createLEConsLinear(&(maxWorkedWeekendCons_[i]), name, theLiveNurses_[i]->maxTotalWeekends() - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
         1, &maxWorkedWeekendVars_[i], &coeff3);
   }
}

int MasterProblem::addMinMaxConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, bool weekend){
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
      vector< vector<SCIP_VAR*> > optDemandVars1(pScenario_->nbShifts_-1);
      vector< vector< vector<SCIP_VAR*> > > skillsAllocVars1(pScenario_->nbShifts_-1);
      vector< vector<SCIP_CONS*> > minDemandCons1(pScenario_->nbShifts_-1);
      vector< vector<SCIP_CONS*> > optDemandCons1(pScenario_->nbShifts_-1);
      vector< vector<SCIP_CONS*> > feasibleSkillsAllocCons1(pScenario_->nbShifts_-1);

      //forget s=0, it's a resting shift
      for(int s=1; s<pScenario_->nbShifts_; s++){
         //initialize vectors
         vector<SCIP_VAR*> optDemandVars2(pScenario_->nbSkills_);
         vector< vector<SCIP_VAR*> > skillsAllocVars2(pScenario_->nbSkills_);
         vector<SCIP_CONS*> minDemandCons2(pScenario_->nbSkills_);
         vector<SCIP_CONS*> optDemandCons2(pScenario_->nbSkills_);
         vector<SCIP_CONS*> feasibleSkillsAllocCons2(pScenario_->nbPositions());

         for(int sk=0; sk<pScenario_->nbSkills_; sk++){
            //initialize vectors
            vector<SCIP_VAR*> skillsAllocVars3(pScenario_->nbPositions());

            //create variables
            SCIPsnprintf(name, 255, "optDemandVar_%d_%d_%d", k, s, sk);
            scip_.createPositiveVar(&(optDemandVars2[sk]), name, WEIGHT_OPTIMAL_DEMAND);
            for(int p=0; p<pScenario_->nbPositions(); p++){
               SCIPsnprintf(name, 255, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
               scip_.createIntVar(&(skillsAllocVars3[p]), name, 0);
            }
            //store vectors
            skillsAllocVars2[sk] = skillsAllocVars3;

            //adding variables and building minimum demand constraints
            int const nonZeroVars1(positionsPerSkill_[sk].size());
            vector<SCIP_VAR*> vars1(nonZeroVars1);
            vector<double> coeff1(nonZeroVars1);
            for(int p=0; p<positionsPerSkill_[sk].size(); ++p){
               vars1[p] = skillsAllocVars3[positionsPerSkill_[sk][p]];
               coeff1[p] = 1;
            }
            SCIPsnprintf(name, 255, "minDemandCons_%d_%d_%d", k, s, sk);
            scip_.createFinalGEConsLinear(&(minDemandCons2[sk]), name, pDemand_->minDemand_[k][s][sk],
               nonZeroVars1, &vars1[0], &coeff1[0]);

            //adding variables and building optimal demand constraints
            int const nonZeroVars2(nonZeroVars1+1);
            vector<SCIP_VAR*> vars2(vars1);
            vector<double> coeff2(coeff1);
            vars2.push_back(optDemandVars2[sk]);
            coeff2.push_back(1);
            SCIPsnprintf(name, 255, "optDemandCons_%d_%d_%d", k, s, sk);
            scip_.createFinalGEConsLinear(&(optDemandCons2[sk]), name, pDemand_->optDemand_[k][s][sk],
               nonZeroVars2, &vars2[0], &coeff2[0]);
         }

         for(int p=0; p<pScenario_->nbPositions(); p++){
            //adding variables and building skills allocation constraints
            int const nonZeroVars3(skillsPerPosition_[p].size());
            vector<SCIP_VAR*> vars3(nonZeroVars3);
            vector<double> coeff3(nonZeroVars3);
            for(int sk=0; sk<skillsPerPosition_[p].size(); ++sk){
               vars3[sk] = skillsAllocVars2[skillsPerPosition_[p][sk]][p];
               coeff3[sk] =-1;
            }
            SCIPsnprintf(name, 255, "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
            scip_.createEQConsLinear(&(feasibleSkillsAllocCons2[p]), name, 0,
               nonZeroVars3, &vars3[0], &coeff3[0]);
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

int MasterProblem::addSkillsCoverageConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, int s){
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

void Rotation::computeCost(Scenario* pScenario, Preferences* pPreferences, int horizon){
   //check if pNurse points to a nurse
   if(pNurse_ == NULL)
      Tools::throwError("LiveNurse = NULL");

   /************************************************
    * Compute all the costs of a rotation:
    ************************************************/

   //if first day of the planning, check on the past, otherwise 0 (rest)
   int lastShift = (firstDay_==0) ? pNurse_->pStateIni_->shift_ : 0;
   //nbConsShift = number of consecutive shift
   //if first day of the planning, check on the past, otherwise 0
   int nbConsShifts = (firstDay_==0) ? pNurse_->pStateIni_->consShifts_ : 0;
   //consShiftCost = cost of be outside of the interval [min,max] of the consecutives shifts
   double consShiftsCost = 0;

   //nbConsWorked = number of consecutive worked days
   //if first day of the planning, check on the past , otherwise 0
   int nbConsDaysWorked = (firstDay_==0) ? pNurse_->pStateIni_->consDaysWorked_ : 0;
   //consWorkedCost = cost of be outside of the interval [min,max] of the consecutives worked days
   double consDaysWorkedCost = 0;

   //cost of not doing the whole weekend
   double completeWeekendCost = 0;

   //preferencesCost = cost of not respecting preferences
   double preferenceCost = 0;

   /*
    * Compute consShiftCost
    */

   for(int k=firstDay_; k<firstDay_+length_; ++k){
      if(lastShift == shifts_[k]){
         nbConsShifts ++;
         continue;
      }
      if(lastShift > 0){
         int diff = max(pScenario->minConsShifts_[lastShift] - nbConsShifts,
            nbConsShifts-pScenario->maxConsShifts_[lastShift]);
         if(diff>0)
            consShiftsCost += diff * WEIGHT_CONS_SHIFTS;
      }
      //initialize nbConsShifts and lastShift
      nbConsShifts = 1;
      lastShift = shifts_[k];
   }

   //compute consShiftsCost for the last shift
   int diff = max((firstDay_+length_-1 == horizon) ? 0 : pScenario->minConsShifts_[lastShift] - nbConsShifts,
      nbConsShifts-pScenario->maxConsShifts_[lastShift]);
   if(diff>0)
      consShiftsCost += diff * WEIGHT_CONS_SHIFTS;

   /*
    * Compute consDaysWorkedCost
    */

   nbConsDaysWorked += length_;
   //check if nbConsDaysWorked < min
   if(nbConsDaysWorked < pNurse_->minConsDaysWork())
      consDaysWorkedCost += (pNurse_->minConsDaysWork() - nbConsDaysWorked) * WEIGHT_CONS_DAYS_WORK;
   //check if nbConsDaysWorked > max
   else if(nbConsDaysWorked > pNurse_->maxConsDaysWork())
      consDaysWorkedCost += (nbConsDaysWorked - pNurse_->maxConsDaysWork()) * WEIGHT_CONS_DAYS_WORK;

   /*
    * Compute completeWeekendCost
    */
   if(pNurse_->pContract_->needCompleteWeekends_){
      //if first day is a Sunday, the saturday is not worked
      if(Tools::isSunday(firstDay_))
         completeWeekendCost += WEIGHT_COMPLETE_WEEKEND;
      //if last day + 1 is a Sunday, the sunday is not worked
      if(Tools::isSunday(firstDay_+length_))
         completeWeekendCost += WEIGHT_COMPLETE_WEEKEND;
   }

   /*
    * Compute preferencesCost
    */

   for(int k=firstDay_; k<firstDay_+length_; ++k)
      if(pPreferences->wantsTheShiftOff(pNurse_->id_, k, shifts_[k]))
         preferenceCost += WEIGHT_PREFERENCES;

   int nbDaysOff = pPreferences->howManyDaysOff(pNurse_->id_, firstDay_, firstDay_+length_-1);

   /*
    * Compute the sum of the cost and stores it in cost_
    */

   cost_ = consShiftsCost + consDaysWorkedCost + completeWeekendCost + preferenceCost;
}


