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
    // only add the constraint to the master if there is a version > 0
    ConstraintMP(pMaster, "Dynamic", true,
                 pMaster->getDynamicWeights().version()),
    dynamicWeightsVersion_(0),
    totalShiftDurationConstraint_(totalShiftDurationConstraint),
    totalShiftDurationResources_(pScenario_->nNurses()),
    totalWeekendConstraint_(totalWeekendConstraint),
    totalWeekendResources_(pScenario_->nNurses()),
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
  // initialize resources
  for (int n = 0; n < pScenario_->nNurses(); ++n) {
    if (totalShiftDurationConstraint->getConstraints()[n].size() != 1)
      Tools::throwError("Dynamic constraints work with exactly "
                        "one total duration constraint instead of %d.",
                        totalShiftDurationConstraint->getConstraints().size());
    totalShiftDurationResources_[n] =
        totalShiftDurationConstraint->getConstraints()[n].begin()->first.second;
    if (totalWeekendConstraint->getConstraints()[n].size() != 1)
      Tools::throwError("Dynamic constraints work with exactly "
                        "one total weekend constraint instead of %d.",
                        totalWeekendConstraint->getConstraints().size());
    totalWeekendResources_[n] =
        totalWeekendConstraint->getConstraints()[n].begin()->first.second;
  }

  // update the version
  dynamicWeightsVersion_ = pMaster_->getDynamicWeights().version();

  build();
}

// update the variables and constraints based on the dynamic weights
void DynamicConstraints::update() {
  const DynamicWeights dW = pMaster_->getDynamicWeights();

  // check if dynamic weights has been modified
  if (pMaster_->getDynamicWeights().version() <= dynamicWeightsVersion_)
    return;

  // check if DynamicConstraints has been built
  if (dynamicWeightsVersion_ == 0)
    Tools::throwError("The dynamic constraints have not been created when "
                      "solve() has been executed. These constraints cannot "
                      "be added on resolve(). You need to set the first version"
                      " of dynamic weights before solve().");

  // update the version
  dynamicWeightsVersion_ = pMaster_->getDynamicWeights().version();

  for (int i = 0; i < pScenario_->nNurses(); i++) {
    // Constraints on the number of shifts worked
    minWorkedDaysAvgVars_[i]->setCost(dW.getWeightTotalShiftsAvg()[i]);
    minWorkedDaysAvgCons_[i]->setLhs(dW.getMinTotalShiftsAvg()[i]);

    maxWorkedDaysAvgVars_[i]->setCost(dW.getWeightTotalShiftsAvg()[i]);
    maxWorkedDaysAvgCons_[i]->setRhs(dW.getMaxTotalShiftsAvg()[i]);

    // Constraints on the number of weekends worked
    maxWorkedWeekendAvgVars_[i]->setCost(dW.getWeightTotalWeekendsAvg()[i]);
    maxWorkedWeekendAvgCons_[i]->setRhs(dW.getMaxTotalWeekendsAvg()[i]);
  }

  // add the same bounds for the contracts
  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    // Constraints on the number of shifts worked
    minWorkedDaysContractAvgVars_[p]->setCost(
        dW.getWeightTotalShiftsContractAvg()[p]);
    minWorkedDaysContractAvgCons_[p]->setLhs(
        dW.getMinTotalShiftsContractAvg()[p]);

    maxWorkedDaysContractAvgVars_[p]->setCost(
        dW.getWeightTotalShiftsContractAvg()[p]);
    maxWorkedDaysContractAvgCons_[p]->setRhs(
        dW.getMaxTotalShiftsContractAvg()[p]);

    maxWorkedWeekendContractAvgVars_[p]->setCost(
        dW.getWeightTotalShiftsContractAvg()[p]);
    maxWorkedWeekendContractAvgCons_[p]->setRhs(
        dW.getMaxTotalWeekendsContractAvg()[p]);
  }
}

// update the dual values of the constraints based on the current solution
void DynamicConstraints::updateDuals() {
  dualValues_.clear();
  weekendDualValues_.clear();
  for (const PLiveNurse &pNurse : pMaster_->pLiveNurses()) {
    const int i = pNurse->num_;
    const int p = pNurse->pContract_->id_;

    /* Dynamic constraints */
    double d = 0;
    d += pModel()->getDual(minWorkedDaysAvgCons_[i], true);
    d += pModel()->getDual(maxWorkedDaysAvgCons_[i], true);
    d += pModel()->getDual(minWorkedDaysContractAvgCons_[p], true);
    d += pModel()->getDual(maxWorkedDaysContractAvgCons_[p], true);
    dualValues_.push_back(d);

    double wd = 0;
    wd += pModel()->getDual(maxWorkedWeekendAvgCons_[i], true);
    wd += pModel()->getDual(
        maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
    weekendDualValues_.push_back(wd);
  }
}

// return the consumption for the duration and the weekend resources
std::pair<int, int> DynamicConstraints::getConsumptions(
    int nurseNum, const Stretch &st, const PAbstractShift &prevS) const {
  // consumption for total duration
  auto *pR = totalShiftDurationResources_[nurseNum];
  bool ready = !(prevS && pR->pShift()->includes(*prevS));
  int c = pR->computeConsumption(st, &ready);

  // consumption for total weekend
  auto *pR2 = totalWeekendResources_[nurseNum];
  ready = !(prevS && pR2->pShift()->includes(*prevS));
  int c2 = pR2->computeConsumption(st, &ready);

  return {c, c2};
}

// return the dual cost of a stretch based on its consumption of
// the constraints
double DynamicConstraints::getDualCost(int nurseNum,
                                       const Stretch &st,
                                       const PAbstractShift &prevS) const {
  // consumptions
  auto c = getConsumptions(nurseNum, st, prevS);
  double d = c.first * dualValues_[nurseNum];  // total duration
  d += c.second * weekendDualValues_[nurseNum];  // total weekend

  double dd = st.duration() * dualValues_[nurseNum],
      wd = weekendDualValues_[nurseNum];
  int k = st.firstDayId();
  const PLiveNurse &pN = pMaster_->pLiveNurses()[nurseNum];
  bool weekendWorked = pN->pContract()->isWeekend(k-1) &&
      prevS && prevS->isWork();
  for (const PShift &pS : st.pShifts()) {
    if (pN->pContract()->isWeekend(k) && !weekendWorked && pS->isWork()) {
      weekendWorked = true;
      dd += wd;
    }
    if (pN->pContract()->isLastWeekendDay(k++))
      weekendWorked = false;
  }

  if (fabs(dd - d) > 1e-5)
    int bb = 1;

  return d;
}

// add a given constraint to the column
void DynamicConstraints::addConsToCol(std::vector<MyCons *> *cons,
                                      std::vector<double> *coeffs,
                                      const Column &col) const {
  int n = col.nurseNum();
  auto pContract = pMaster_->pLiveNurses()[n]->pContract_;

  // consumptions
  const PAbstractShift prevS =
      col.firstDayId() > 0 ? pScenario_->pRestShift() :
      pMaster_->pLiveNurses().at(n)->pStateIni_->pShift_;
  auto c = getConsumptions(n, col, prevS);

  int nWeekends = pContract->weekend_.nWeekendsInInterval(
      col.firstDayId(), col.lastDayId());
  if (col.duration() != c.first || nWeekends != c.second)
    int bb = 0;

  // total duration constraints
  cons->push_back(minWorkedDaysAvgCons_[n]);
  coeffs->push_back(c.first);

  cons->push_back(maxWorkedDaysAvgCons_[n]);
  coeffs->push_back(c.first);

  cons->push_back(minWorkedDaysContractAvgCons_[pContract->id_]);
  coeffs->push_back(c.first);

  cons->push_back(maxWorkedDaysContractAvgCons_[pContract->id_]);
  coeffs->push_back(c.first);

  // total weekend constraints
  if (c.second) {
    cons->push_back(maxWorkedWeekendAvgCons_[n]);
    coeffs->push_back(c.second);

    cons->push_back(maxWorkedWeekendContractAvgCons_[pContract->id_]);
    coeffs->push_back(c.second);
  }
}

std::string DynamicConstraints::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  // consumptions
  auto c = getConsumptions(nurseNum, st, nullptr);
  buff << "#   |  Dynamic total duration: "
       << -c.first * dualValues_[nurseNum] << std::endl;
  buff << "#   |  Dynamic total weekend: "
       << -c.second * weekendDualValues_[nurseNum] << std::endl;
  return buff.str();
}


void DynamicConstraints::build() {
  char name[255];
  const DynamicWeights dW = pMaster_->getDynamicWeights();

  for (int i = 0; i < pScenario_->nNurses(); i++) {
    const PNurse &pN = pMaster_->pLiveNurses()[i];
    // add constraints on the total duration of shifts to satisfy bounds that
    // correspond to the global bounds averaged over the weeks

    /** Min total duration */
    // slack variable
    snprintf(name, sizeof(name), "getMinWorkedDaysAvgVar_N%d", i);
    pModel()->createPositiveVar(
        &minWorkedDaysAvgVars_[i],
        name,
        dW.getWeightTotalShiftsAvg()[i]);

    vector<MyVar *> vars = {minWorkedDaysAvgVars_[i]};
    vector<double> coeffs = {1};
    // add the slack of the regular constraint on shift duration to be sure
    // to not count twice a penalty
    MyVar *minWorkedDaysVar =
        totalShiftDurationConstraint_->getVariables()[i].begin()->second.first;
    // create min constraint
    snprintf(name, sizeof(name), "minWorkedDaysAvgCons_N%d", i);
    pModel()->createGEConsLinear(
        &minWorkedDaysAvgCons_[i],
        name,
        dW.getMinTotalShiftsAvg()[i],
        {minWorkedDaysAvgVars_[i], minWorkedDaysVar},
        {1, 1});

    /** Max total duration */
    // slack variable
    snprintf(name, sizeof(name), "maxWorkedDaysAvgVar_N%d", i);
    pModel()->createPositiveVar(
        &maxWorkedDaysAvgVars_[i],
        name,
        dW.getWeightTotalShiftsAvg()[i]);
    // add the slack of the regular constraint on shift duration to be sure
    // to not count twice a penalty
    MyVar *maxWorkedDaysVar =
        totalShiftDurationConstraint_->getVariables()[i].begin()->second.second;
    // create max constraint
    snprintf(name, sizeof(name), "maxWorkedDaysAvgCons_N%d", i);
    pModel()->createLEConsLinear(
        &maxWorkedDaysAvgCons_[i],
        name,
        dW.getMaxTotalShiftsAvg()[i],
        {maxWorkedDaysVar, maxWorkedDaysAvgVars_[i]},
        {-1, -1});

    /** Max total worked weekend */
    // slack variable
    snprintf(name, sizeof(name), "maxWorkedWeekendAvgVar_N%d", i);
    pModel()->createPositiveVar(
        &maxWorkedWeekendAvgVars_[i],
        name,
        dW.getWeightTotalWeekendsAvg()[i]);
    // add the slack of the regular constraint on shift duration to be sure
    // to not count twice a penalty
    MyVar *maxWorkedWeekendVar =
        totalWeekendConstraint_->getVariables()[i].begin()->second.second;
    // create max constraint
    snprintf(name, sizeof(name), "maxWorkedWeekendAvgCons_N%d", i);
    pModel()->createLEConsLinear(
        &maxWorkedWeekendAvgCons_[i],
        name,
        dW.getMaxTotalWeekendsAvg()[i],
        {maxWorkedWeekendVar, maxWorkedWeekendAvgVars_[i]},
        {-1, -1});
  }

  // add the same bounds for the contracts,
  // however the slacks of the regular constraints are not taken into account as
  // we are only bounding the total on a set of nurses
  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    /** Min total duration for all nurses of a given contract */
    snprintf(name, sizeof(name), "minWorkedDaysContractAvgVar_P%d", p);
    pModel()->createPositiveVar(
        &minWorkedDaysContractAvgVars_[p],
        name,
        dW.getWeightTotalShiftsContractAvg()[p]);
    snprintf(name, sizeof(name), "minWorkedDaysContractAvgCons_P%d", p);
    pModel()->createGEConsLinear(
        &minWorkedDaysContractAvgCons_[p],
        name,
        dW.getMinTotalShiftsContractAvg()[p],
        {minWorkedDaysContractAvgVars_[p]},
        {1});

    /** Max total duration for all nurses of a given contract */
    snprintf(name, sizeof(name), "maxWorkedDaysContractAvgVar_P%d", p);
    pModel()->createPositiveVar(
        &maxWorkedDaysContractAvgVars_[p],
        name,
        dW.getWeightTotalShiftsContractAvg()[p]);
    snprintf(name, sizeof(name), "maxWorkedDaysContractAvgCons_P%d", p);
    pModel()->createLEConsLinear(
        &maxWorkedDaysContractAvgCons_[p],
        name,
        dW.getMaxTotalShiftsContractAvg()[p],
        {maxWorkedDaysContractAvgVars_[p]},
        {-1});

    /** Max total worked weekend for all nurses of a given contract */
    snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgVar_P%d", p);
    pModel()->createPositiveVar(
        &maxWorkedWeekendContractAvgVars_[p],
        name,
        dW.getWeightTotalWeekendsContractAvg()[p]);
    snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgCons_C%d", p);
    pModel()->createLEConsLinear(
        &maxWorkedWeekendContractAvgCons_[p],
        name,
        dW.getMaxTotalWeekendsContractAvg()[p],
        {maxWorkedWeekendContractAvgVars_[p]},
        {-1});
  }
}
