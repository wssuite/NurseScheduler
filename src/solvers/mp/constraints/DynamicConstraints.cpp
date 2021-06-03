/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "DynamicConstraints.h"

#include "solvers/mp/MasterProblem.h"


DynamicConstraints::DynamicConstraints(
    MasterProblem *pMaster,
    TotalShiftDurationConstraint *totalShiftDurationConstraint,
    TotalWeekendConstraint *totalWeekendConstraint) :
    ConstraintMP(pMaster, true),
    totalShiftDurationConstraint_(totalShiftDurationConstraint),
    totalWeekendConstraint_(totalWeekendConstraint),
    minWorkedDaysAvgVars_(pScenario_->nNurses()),
    maxWorkedDaysAvgVars_(pScenario_->nNurses()),
    maxWorkedWeekendAvgVars_(pScenario_->nNurses()),
    minWorkedDaysContractAvgVars_(pScenario_->nContracts()),
    maxWorkedDaysContractAvgVars_(pScenario_->nContracts()),
    maxWorkedWeekendContractAvgVars_(pScenario_->nContracts()),
    minWorkedDaysAvgCons_(pScenario_->nNurses()),
    maxWorkedDaysAvgCons_(pScenario_->nNurses()),
    maxWorkedWeekendAvgCons_(pScenario_->nNurses()),
    minWorkedDaysContractAvgCons_(pScenario_->nContracts()),
    maxWorkedDaysContractAvgCons_(pScenario_->nContracts()),
    maxWorkedWeekendContractAvgCons_(pScenario_->nContracts()) {
  // initialize the vectors indicating whether the min/max total constraints
  // with averaged bounds are considered
  Tools::initVector(&isMinWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedWeekendAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMinWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedWeekendContractAvgCons_,
                    pScenario_->nContracts(),
                    false);

  build();
}

// update the dual values of the constraints based on the current solution
void DynamicConstraints::updateDuals() {
  dualValues_.clear();
  weekenDualValues_.clear();
  for (const PLiveNurse &pNurse : pMaster_->liveNurses()) {
    const int i = pNurse->num_;
    const int p = pNurse->pContract_->id_;

    /* Dynamic constraints */
    double minWorkedDaysAvg =
        isMinWorkedDaysAvgCons_[i] ?
        pModel()->getDual(minWorkedDaysAvgCons_[i], true) : 0.0;
    double maxWorkedDaysAvg =
        isMaxWorkedDaysAvgCons_[i] ?
        pModel()->getDual(maxWorkedDaysAvgCons_[i], true) : 0.0;

    double minWorkedDaysContractAvg =
        isMinWorkedDaysContractAvgCons_[p] ?
        pModel()->getDual(minWorkedDaysContractAvgCons_[p], true) : 0.0;
    double maxWorkedDaysContractAvg =
        isMaxWorkedDaysContractAvgCons_[p] ?
        pModel()->getDual(maxWorkedDaysContractAvgCons_[p], true) : 0.0;

    double d = minWorkedDaysAvg + minWorkedDaysContractAvg
        + maxWorkedDaysAvg + maxWorkedDaysContractAvg;
    dualValues_.push_back(d);

    double wd = 0;
    if (isMaxWorkedWeekendAvgCons_[i])
      wd += pModel()->getDual(maxWorkedWeekendAvgCons_[i], true);
    if (isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_])
      wd += pModel()->getDual(
          maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
    weekenDualValues_.push_back(wd);
  }
}

// return the dual cost of a stretch based on its consumption of
// the constraints
double DynamicConstraints::getDualCost(int nurseNum,
                                       const Stretch &st,
                                       const PAbstractShift &prevS) const {
  double d = st.duration() * dualValues_[nurseNum],
      wd = weekenDualValues_[nurseNum];
  int k = st.firstDay();
  bool weekendWorked =
      Tools::isWeekendDayButNotLastOne(k - 1) && prevS && prevS->isWork();
  for (const PShift &pS : st.pShifts()) {
    if (Tools::isWeekend(k) && !weekendWorked && pS->isWork()) {
      weekendWorked = true;
      d += wd;
    }
    if (Tools::isLastWeekendDay(k++))
      weekendWorked = false;
  }
  return d;
}

// add a given constraint to the column
void DynamicConstraints::addConsToCol(std::vector<MyCons *> *cons,
                                      std::vector<double> *coeffs,
                                      const Pattern &col) const {
  int p = pMaster_->liveNurses()[col.nurseNum()]->pContract_->id_;
  int duration = col.duration();
  if (isMinWorkedDaysAvgCons_[col.nurseNum()]) {
    cons->push_back(minWorkedDaysAvgCons_[col.nurseNum()]);
    coeffs->push_back(duration);
  }
  if (isMaxWorkedDaysAvgCons_[col.nurseNum()]) {
    cons->push_back(maxWorkedDaysAvgCons_[col.nurseNum()]);
    coeffs->push_back(duration);
  }
  if (isMinWorkedDaysContractAvgCons_[p]) {
    cons->push_back(minWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(duration);
  }
  if (isMaxWorkedDaysContractAvgCons_[p]) {
    cons->push_back(maxWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(duration);
  }

  int nWeekends = Tools::nWeekendsInInterval(col.firstDay(), col.lastDay());
  if (nWeekends) {
    if (isMaxWorkedWeekendAvgCons_[col.nurseNum()]) {
      cons->push_back(maxWorkedWeekendAvgCons_[col.nurseNum()]);
      coeffs->push_back(nWeekends);
    }

    if (isMaxWorkedWeekendContractAvgCons_[p]) {
      cons->push_back(maxWorkedWeekendContractAvgCons_[p]);
      coeffs->push_back(nWeekends);
    }
  }
}

void DynamicConstraints::build() {
  char name[255];
  for (int i = 0; i < pScenario_->nNurses(); i++) {
    // add constraints on the total number of shifts to satisfy bounds that
    // correspond to the global bounds averaged over the weeks
    if (!pMaster_->minTotalShiftsAvg_.empty()
        && !pMaster_->maxTotalShiftsAvg_.empty()
        && !pMaster_->weightTotalShiftsAvg_.empty()) {
      // only add the constraint if is tighter than the already added constraint
      if (pMaster_->minTotalShiftsAvg_[i] > pMaster_->minTotalShifts_[i]) {
        snprintf(name, sizeof(name), "minWorkedDaysAvgVar_N%d", i);
        pModel()->createPositiveVar(
            &minWorkedDaysAvgVars_[i],
            name,
            pMaster_->weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "minWorkedDaysAvgCons_N%d", i);
        vector<MyVar*> vars = {minWorkedDaysAvgVars_[i]};
        vector<double> coeffs = {1};
        // if a LB constraint
        if (pMaster_->minTotalShifts_[i] > 0) {
          vars.push_back(totalShiftDurationConstraint_
                             ->getVariables()[i].begin()->second.at(0));
          coeffs.push_back(1);
        }
        pModel()->createGEConsLinear(
            &minWorkedDaysAvgCons_[i],
            name,
            pMaster_->minTotalShiftsAvg_[i],
            vars,
            coeffs);

        isMinWorkedDaysAvgCons_[i] = true;
      }

      if (pMaster_->maxTotalShiftsAvg_[i] < pMaster_->maxTotalShifts_[i]) {
        snprintf(name, sizeof(name), "maxWorkedDaysAvgVar_N%d", i);
        pModel()->createPositiveVar(
            &maxWorkedDaysAvgVars_[i],
            name,
            pMaster_->weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "maxWorkedDaysAvgCons_N%d", i);
        MyVar *maxWorkedDaysVar = totalShiftDurationConstraint_
            ->getVariables()[i].begin()->second.back();
        pModel()->createLEConsLinear(
            &maxWorkedDaysAvgCons_[i],
            name,
            pMaster_->maxTotalShiftsAvg_[i],
            {maxWorkedDaysVar, maxWorkedDaysAvgVars_[i]},
            {-1, -1});

        isMaxWorkedDaysAvgCons_[i] = true;
      }
    }

    // STAB: not implemented there yet
    if (!pMaster_->maxTotalWeekendsAvg_.empty()
        && !pMaster_->weightTotalWeekendsAvg_.empty()
        && pMaster_->maxTotalWeekendsAvg_[i] <
            pMaster_->liveNurses()[i]->maxTotalWeekends()
                - pMaster_->liveNurses()[i]->pStateIni_->totalWeekendsWorked_) {
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgVar_N%d", i);
      pModel()->createPositiveVar(
          &maxWorkedWeekendAvgVars_[i],
          name,
          pMaster_->weightTotalWeekendsAvg_[i]);
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgCons_N%d", i);
      MyVar *maxWorkedWeekendVar = totalWeekendConstraint_
          ->getVariables()[i].begin()->second.at(0);
      pModel()->createLEConsLinear(
          &maxWorkedWeekendAvgCons_[i],
          name,
          pMaster_->maxTotalWeekendsAvg_[i]
              - pMaster_->liveNurses()[i]->pStateIni_->totalWeekendsWorked_,
          {maxWorkedWeekendVar, maxWorkedWeekendAvgVars_[i]},
          {-1, -1});

      isMaxWorkedWeekendAvgCons_[i] = true;
    }
  }

  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    if (!pMaster_->minTotalShiftsContractAvg_.empty()
        && !pMaster_->maxTotalShiftsContractAvg_.empty()
        && !pMaster_->weightTotalShiftsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "minWorkedDaysContractAvgVar_P%d", p);
      pModel()->createPositiveVar(
          &minWorkedDaysContractAvgVars_[p],
          name,
          pMaster_->weightTotalShiftsContractAvg_[p]);
      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgVar_P%d", p);
      pModel()->createPositiveVar(
          &maxWorkedDaysContractAvgVars_[p],
          name,
          pMaster_->weightTotalShiftsContractAvg_[p]);

      snprintf(name, sizeof(name), "minWorkedDaysContractAvgCons_P%d", p);
      pModel()->createGEConsLinear(
          &minWorkedDaysContractAvgCons_[p],
          name,
          pMaster_->minTotalShiftsContractAvg_[p],
          {minWorkedDaysContractAvgVars_[p]},
          {1});

      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgCons_P%d", p);
      pModel()->createLEConsLinear(
          &maxWorkedDaysContractAvgCons_[p],
          name,
          pMaster_->maxTotalShiftsContractAvg_[p],
          {maxWorkedDaysContractAvgVars_[p]},
          {-1});

      isMinWorkedDaysContractAvgCons_[p] = true;
      isMaxWorkedDaysContractAvgCons_[p] = true;
    }

    if (!pMaster_->maxTotalWeekendsContractAvg_.empty()
        && !pMaster_->weightTotalWeekendsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgVar_P%d", p);
      pModel()->createPositiveVar(
          &maxWorkedWeekendContractAvgVars_[p],
          name,
          pMaster_->weightTotalWeekendsContractAvg_[p]);

      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgCons_C%d", p);
      pModel()->createLEConsLinear(
          &maxWorkedWeekendContractAvgCons_[p],
          name,
          pMaster_->maxTotalWeekendsContractAvg_[p],
          {maxWorkedWeekendContractAvgVars_[p]},
          {-1});

      isMaxWorkedWeekendContractAvgCons_[p] = true;
    }
  }
}
