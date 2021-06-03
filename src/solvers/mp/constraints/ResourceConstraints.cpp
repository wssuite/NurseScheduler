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

#include "ResourceConstraints.h"

#include <vector>

#include "solvers/mp/MasterProblem.h"


template <class T>
BoundedResourceConstraint<T>::BoundedResourceConstraint(
    MasterProblem *pMaster) :
    ConstraintMP(pMaster, true),
    resourceVars_(pMaster->nNurses()),
    resourceCons_(pMaster->nNurses()) {}

// update the dual values of the constraints based on the current solution
template <class T>
void BoundedResourceConstraint<T>::updateDuals() {
  dualValues_.clear();
  for (const auto &constraints : resourceCons_) {
    std::map<T*, double> duals;
    for (const auto &p : constraints)
      duals[p.first] = pModel()->getDual(p.second);
    dualValues_.push_back(duals);
  }
}

// return the dual cost of a stretch based on its consumption of
// the constraints
template <class T>
double BoundedResourceConstraint<T>::getDualCost(
    int nurseNum,
    const Stretch &st,
    const PAbstractShift &prevS) const {
  double d = 0;
  const auto &duals = dualValues_[nurseNum];
  for (const auto &p : resourceCons_[nurseNum]) {
    // ready if prevS is null or  not included in pShift
    bool ready = !(prevS && p.first->pShift()->includes(*prevS));
    int c = p.first->computeConsumption(st, &ready);
    d += c * duals.at(p.first);
  }
  return d;
}

// add a given constraint to the column
template <class T>
void BoundedResourceConstraint<T>::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Pattern &col) const {
  for (const auto &p : resourceCons_[col.nurseNum()]) {
    bool ready = !p.first->pShift()->includes(*pScenario_->pRestShift());
    int c = p.first->computeConsumption(col, &ready);
    for (MyCons *pC : p.second) {
      cons->push_back(pC);
      coeffs->push_back(c);
    }
  }
}

template <class T>
std::string BoundedResourceConstraint<T>::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  const auto &duals = dualValues_[nurseNum];
  for (const auto &p : resourceCons_[nurseNum]) {
    buff << "#   |  " << p.first->pResource()->name << ": "
         << -getDualCost(nurseNum, st, pScenario_->pRestShift()) << std::endl;
  }
  return buff.str();
}

// update the dual values of the constraints randomly
template <class T>
void BoundedResourceConstraint<T>::createRandomUpdateDuals(double weight) {
  dualValues_.clear();
  for (const auto &constraints : resourceCons_) {
    std::map<T*, double> duals;
    for (const auto &p : constraints)
      duals[p.first] = Tools::randomDouble(0, weight);
    dualValues_.push_back(duals);
  }
}

TotalShiftDurationConstraint::TotalShiftDurationConstraint(
    MasterProblem *pMaster) :
    BoundedResourceConstraint<TotalShiftDuration>(pMaster) {}

// update the dual values of the constraints randomly
void TotalShiftDurationConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  createRandomUpdateDuals(2*pScenario_->weights().WEIGHT_TOTAL_SHIFTS);
}

TotalWeekendConstraint::TotalWeekendConstraint(MasterProblem *pMaster) :
    BoundedResourceConstraint<TotalWeekend>(pMaster) {}

// update the dual values of the constraints randomly
void TotalWeekendConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  createRandomUpdateDuals(2*pScenario_->weights().WEIGHT_TOTAL_WEEKENDS);
}
