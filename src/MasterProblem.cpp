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
   Preferences* pPreferences, vector<State>* pInitState):
                     Solver(pScenario, pDemand, pPreferences, pInitState),
                     scip_("GenCol"), positionsPerSkill_(pScenario->nbSkills_), skillsPerPosition_(pScenario->nbPositions()),

                     columnVars_(pScenario->nbNurses_), restingVars_(pScenario->nbNurses_), longRestingVars_(pScenario->nbNurses_),
                     minWorkedDaysVars_(pScenario->nbNurses_), maxWorkedDaysVars_(pScenario->nbNurses_), maxWorkedWeekendVars_(pScenario->nbNurses_),
                     optDemandVars_(7*pScenario->nbWeeks_), skillsAllocVars_(7*pScenario->nbWeeks_),

                     restFlowCons_(pScenario->nbNurses_), workFlowCons_(pScenario->nbNurses_), sourceFlowCons_(pScenario->nbNurses_), sinkFlowCons_(pScenario->nbNurses_),
                     minWorkedDaysCons_(pScenario->nbNurses_), maxWorkedDaysCons_(pScenario->nbNurses_), maxWorkedWeekendCons_(pScenario->nbNurses_),
                     minDemandCons_(7*pScenario->nbWeeks_), optDemandCons_(7*pScenario->nbWeeks_), feasibleSkillsAllocCons_(7*pScenario->nbWeeks_)
{

   this->preprocessTheNurses();

   /*
    * Build the two vectors linking positions and skills
    */
   for(int p=0; p<skillsPerPosition_.size(); p++){
      vector<int> skills(pScenario->nbSkills_);
      for(auto &sk : pScenario->pPositions()[p]->skills_)
         skills.push_back(sk);
      skillsPerPosition_[p] = skills;
   }
   for(int sk=0; sk<positionsPerSkill_.size(); sk++){
      vector<int> positions(pScenario->nbPositions());
      for(int p=0; p<positions.size(); p++)
         if(find(skillsPerPosition_[p].begin(), skillsPerPosition_[p].end(), sk) != skillsPerPosition_[p].end())
            positions.push_back(p);
      positionsPerSkill_[sk] = positions;
   }

   build();
}

void MasterProblem::build(){
   char name[255];

   /* Min/Max constraints */
   for(int i=0; i<pScenario_->nbNurses_; i++){
      scip_.createPositiveVar(minWorkedDaysVars_[i], "minWorkedDaysVar_"+i, WEIGHT_CONS_DAYS_WORK);
      scip_.createPositiveVar(maxWorkedDaysVars_[i], "maxWorkedDaysVar_"+i, WEIGHT_CONS_DAYS_WORK);
      scip_.createPositiveVar(maxWorkedWeekendVars_[i], "maxWorkedWeekendVar_"+i, WEIGHT_COMPLETE_WEEKEND);

      double coeff1 = 1;
      scip_.createGEConsLinear(minWorkedDaysCons_[i], "minWorkedDaysCons_"+i, theLiveNurses_[i]->minWorkDays_,
         1, &minWorkedDaysVars_[i], &coeff1);
      double coeff2 = -1;
      scip_.createLEConsLinear(minWorkedDaysCons_[i], "minWorkedDaysCons_"+i, theLiveNurses_[i]->maxWorkDays_,
         1, &maxWorkedDaysVars_[i], &coeff2);
      double coeff3 = -1;
      scip_.createLEConsLinear(minWorkedDaysCons_[i], "minWorkedDaysCons_"+i, theLiveNurses_[i]->maxWorkDays_,
         1, &maxWorkedWeekendVars_[i], &coeff3);
   }

   /* Skills coverage constraints */
   for(int k=0; k<7*pScenario_->nbWeeks_; k++){
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
         vector<SCIP_CONS*> feasibleSkillsAllocCons2(pScenario_->nbSkills_);

         for(int sk=0; sk<pScenario_->nbSkills_; sk++){
            //initialize vectors
            vector<SCIP_VAR*> skillsAllocVars3(pScenario_->nbSkills_);

            //create variables
            SCIPsnprintf(name, 255, "optDemandVar_%d_%d_%d", k, s, sk);
            scip_.createPositiveVar(optDemandVars2[sk], name, WEIGHT_OPTIMAL_DEMAND);
            for(int p=0; p<pScenario_->nbPositions(); p++){
               SCIPsnprintf(name, 255, "skillsAllocVar_%d_%d_%d_%d", k, s, sk,p);
               scip_.createIntVar(skillsAllocVars3[p], name, 0);
            }
            //store vectors
            skillsAllocVars2.push_back(skillsAllocVars3);

            //adding variables and building minimum demand constraints
            int const nonZeroVars1 = positionsPerSkill_[sk].size();
            vector<SCIP_VAR*> vars1(nonZeroVars1);
            vector<double> coeff1(nonZeroVars1);
            for(int p: positionsPerSkill_[sk]){
               vars1.push_back(skillsAllocVars3[p]);
               coeff1.push_back(1);
            }
            SCIPsnprintf(name, 255, "minDemandCons_%d_%d_%d", k, s, sk);
            scip_.createGEConsLinear(minDemandCons2[sk], name, pDemand_->minDemand_[k][s][sk],
               nonZeroVars1, &vars1[0], &coeff1[0]);

            //adding variables and building optimal demand constraints
            int const nonZeroVars2 = nonZeroVars1+1;
            vector<SCIP_VAR*> vars2(vars1);
            vector<double> coeff2(coeff1);
            vars2.push_back(optDemandVars2[sk]);
            coeff2.push_back(1);
            SCIPsnprintf(name, 255, "optDemandCons_%d_%d_%d", k, s, sk);
            scip_.createGEConsLinear(optDemandCons2[sk], name, pDemand_->optDemand_[k][s][sk],
               nonZeroVars2, &vars2[0], &coeff2[0]);

            //adding variables and building optimal demand constraints
            vector<SCIP_VAR*> vars3(vars1);
            vector<double> coeff3(nonZeroVars1);
            for(int p: positionsPerSkill_[sk])
               coeff3.push_back(-1);
            SCIPsnprintf(name, 255, "feasibleSkillsAllocCons_%d_%d_%d", k, s, sk);
            scip_.createEQConsLinear(feasibleSkillsAllocCons2[sk], name, 0,
               nonZeroVars1, &vars3[0], &coeff3[0]);
         }

         //store vectors
         optDemandVars1.push_back(optDemandVars2);
         skillsAllocVars1.push_back(skillsAllocVars2);
         minDemandCons1.push_back(minDemandCons2);
         optDemandCons1.push_back(optDemandCons2);
         feasibleSkillsAllocCons1.push_back(feasibleSkillsAllocCons2);
      }

      //store vectors
      optDemandVars_.push_back(optDemandVars1);
      skillsAllocVars_.push_back(skillsAllocVars1);
      minDemandCons_.push_back(minDemandCons1);
      optDemandCons_.push_back(optDemandCons1);
      feasibleSkillsAllocCons_.push_back(feasibleSkillsAllocCons1);
   }

}

void MasterProblem::solve(){

}


