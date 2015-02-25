/*
 * MasterProblem.cpp
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#include "MasterProblem.h"

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
                                    scip_(PB_NAME), positionsPerSkill_(pScenario->nbSkills_), skillsPerPosition_(pScenario->nbPositions()),

                                    columnVars_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
                                    minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
                                    optDemandVars_(pDemand_->nbDays_), skillsAllocVars_(pDemand_->nbDays_),

                                    restFlowCons_(pScenario->nbNurses_), workFlowCons_(pScenario->nbNurses_), sourceFlowCons_(pScenario->nbNurses_), sinkFlowCons_(pScenario->nbNurses_),
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

//build the rostering problem
void MasterProblem::build(){
   /* Rotation constraints */
   buildRotationCons();

   /* Min/Max constraints */
   buildMinMaxCons();

   /* Skills coverage constraints */
   buildSkillsCoverageCons();
}

//solve the rostering problem
void MasterProblem::solve(){
   scip_.solve();
   //   scip_.printStats();
   scip_.printBestSol();
}

//initialize the rostering problem with one column to be feasible if there is no initial solution
//otherwise build the columns corresponding to the initial solution
void MasterProblem::initialize(vector<Roster> solution){
   int nbVars(0);
   vector<SCIP_VAR*> initialVar(pScenario_->nbNurses_);
   char name[255];
   vector<char*> names(pScenario_->nbNurses_);
   vector<double> initialObjCoeffs(pScenario_->nbNurses_);
   vector<int> nbCons(pScenario_->nbNurses_);
   vector< vector<SCIP_CONS*> > initialCons(pScenario_->nbNurses_);
   vector< vector<double> > initialCoeffs(pScenario_->nbNurses_);

   //if there is no initial solution, we add a column with 1 everywhere to be feasible
   if(solution.size() == 0){
      ++nbVars;
      SCIPsnprintf(name, 255, "InitialGlobalNurse");
      names[0] = name;
      initialObjCoeffs[0] = bigM;
      vector<SCIP_CONS*> initialCons2(0);
      vector<double> initialCoeffs2(0);
      nbCons[0] = 0;
      for(int i=0; i<pScenario_->nbNurses_; i++)
         for(int k=0; k<pDemand_->nbDays_; k++)
            nbCons[0] += addConsToCol(&initialCons2, &initialCoeffs2, i, k);
      initialCons[0] = initialCons2;
      initialCoeffs[0] = initialCoeffs2;
   }
   //otherwise, a column is added for each nurse of the initial solution
   else{

   }

   //we create the initial column(s)
   for(int i=0; i<nbVars; ++i){
      scip_.createBinaryColumn(&initialVar[i], names[i], initialObjCoeffs[i],
         nbCons[i], &(initialCons[i])[0], &(initialCoeffs[i])[0]);
   }
}

//compute the coefficient for each constraint of a column for the nurse i, the day k and the shift s
//if s=-1, the nurse i works on all shifts
int MasterProblem::addConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k, int s){
   int nbCons(0);

   /* Min/Max constraints */
   nbCons += addRotationConsToCol(cons, coeffs, i, k);

   /* Min/Max constraints */
   nbCons += addMinMaxConsToCol(cons, coeffs, i, k);

   /* Skills coverage constraints */
   nbCons += addSkillsCoverageConsToCol(cons, coeffs, i, k, s);

   return nbCons;
}

void MasterProblem::buildRotationCons(){

}

int MasterProblem::addRotationConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k){
   return 0;
}

/*
 * Min/Max constraints
 */
void MasterProblem::buildMinMaxCons(){
   for(int i=0; i<pScenario_->nbNurses_; i++){
      scip_.createPositiveVar(&minWorkedDaysVars_[i], "minWorkedDaysVar_"+i, WEIGHT_CONS_DAYS_WORK);
      scip_.createPositiveVar(&maxWorkedDaysVars_[i], "maxWorkedDaysVar_"+i, WEIGHT_CONS_DAYS_WORK);
      scip_.createPositiveVar(&maxWorkedWeekendVars_[i], "maxWorkedWeekendVar_"+i, WEIGHT_COMPLETE_WEEKEND);

      double coeff1(1);
      scip_.createGEConsLinear(&minWorkedDaysCons_[i], "minWorkedDaysCons_"+i, theLiveNurses_[i]->minWorkDays_,
         1, &minWorkedDaysVars_[i], &coeff1);
      double coeff2(-1);
      scip_.createLEConsLinear(&maxWorkedDaysCons_[i], "maxWorkedDaysCons_"+i, theLiveNurses_[i]->maxWorkDays_,
         1, &maxWorkedDaysVars_[i], &coeff2);
      double coeff3(-1);
      scip_.createLEConsLinear(&maxWorkedWeekendCons_[i], "maxWorkedWeekendCons_"+i, theLiveNurses_[i]->maxWorkDays_,
         1, &maxWorkedWeekendVars_[i], &coeff3);
   }
}

int MasterProblem::addMinMaxConsToCol(vector<SCIP_CONS*>* cons, vector<double>* coeffs, int i, int k){
   int nbCons(0);

   ++nbCons;
   cons->push_back(minWorkedDaysCons_[i]);
   coeffs->push_back(1.0);
   ++nbCons;
   cons->push_back(maxWorkedDaysCons_[i]);
   coeffs->push_back(1.0);
   string day(Tools::intToDay(i));
   if(day == "Sun" || day == "Sat"){
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
      vector< vector<SCIP_VAR*> > optDemandVars1(pScenario_->nbShifts_);
      vector< vector< vector<SCIP_VAR*> > > skillsAllocVars1(pScenario_->nbShifts_);
      vector< vector<SCIP_CONS*> > minDemandCons1(pScenario_->nbShifts_);
      vector< vector<SCIP_CONS*> > optDemandCons1(pScenario_->nbShifts_);
      vector< vector<SCIP_CONS*> > feasibleSkillsAllocCons1(pScenario_->nbShifts_);

      for(int s=0; s<pScenario_->nbShifts_; s++){
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
            scip_.createPositiveVar(&optDemandVars2[sk], name, WEIGHT_OPTIMAL_DEMAND);
            for(int p=0; p<pScenario_->nbPositions(); p++){
               SCIPsnprintf(name, 255, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
               scip_.createIntVar(&skillsAllocVars3[p], name, 0);
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
            scip_.createFinalGEConsLinear(&minDemandCons2[sk], name, pDemand_->minDemand_[k][s][sk],
               nonZeroVars1, &vars1[0], &coeff1[0]);

            //adding variables and building optimal demand constraints
            int const nonZeroVars2(nonZeroVars1+1);
            vector<SCIP_VAR*> vars2(vars1);
            vector<double> coeff2(coeff1);
            vars2.push_back(optDemandVars2[sk]);
            coeff2.push_back(1);
            SCIPsnprintf(name, 255, "optDemandCons_%d_%d_%d", k, s, sk);
            scip_.createFinalGEConsLinear(&optDemandCons2[sk], name, pDemand_->optDemand_[k][s][sk],
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
            scip_.createEQConsLinear(&feasibleSkillsAllocCons2[p], name, 0,
               nonZeroVars3, &vars3[0], &coeff3[0]);
         }

         //store vectors
         optDemandVars1[s] = optDemandVars2;
         skillsAllocVars1[s] = skillsAllocVars2;
         minDemandCons1[s] = minDemandCons2;
         optDemandCons1[s] = optDemandCons2;
         feasibleSkillsAllocCons1[s] = feasibleSkillsAllocCons2;
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
      for(int s0=0; s0<pScenario_->nbShifts_; ++s0){
         ++nbCons;
         cons->push_back(feasibleSkillsAllocCons_[k][s0][p]);
         coeffs->push_back(1.0);
      }
   }
   else{
      ++nbCons;
      cons->push_back(feasibleSkillsAllocCons_[k][s][p]);
      coeffs->push_back(1.0);
   }

   return nbCons;
}


