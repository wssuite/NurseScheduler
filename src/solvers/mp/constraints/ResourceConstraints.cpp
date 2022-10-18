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

#include <algorithm>
#include <vector>

#include "solvers/mp/MasterProblem.h"


template <class H, class S>
BoundedResourceConstraint<H, S>::BoundedResourceConstraint(
    MasterProblem *pMaster, std::string name) :
    ConstraintMP(pMaster, name, true),
    resourceVars_(pMaster->nNurses()),
    resourceCons_(pMaster->nNurses()) {}

// update the dual values of the constraints based on the current solution
template <class H, class S>
void BoundedResourceConstraint<H, S>::updateDuals() {
  dualValues_.clear();
  for (const auto &constraints : resourceCons_) {
    std::map<std::pair<H*, S*>, double> duals;
    for (const auto &p : constraints)
      duals[p.first] = pModel()->getDual(p.second);
    dualValues_.push_back(duals);
  }
}

// compute consumption
template <class H, class S>
int BoundedResourceConstraint<H, S>::computeConsumption(
    const Stretch &st,
    const std::pair<H*, S*> &p,
    const PAbstractShift &prevS) const {
  auto pS = p.first ? p.first->pShift() : p.second->pShift();
  bool ready = !(prevS && pS->includes(*prevS));
  int c = p.first ? p.first->computeConsumption(st, &ready):
          p.second->computeConsumption(st, &ready);
  return c;
}

// return the dual cost of a stretch based on its consumption of
// the constraints
template <class H, class S>
double BoundedResourceConstraint<H, S>::getDualCost(
    int nurseNum,
    const Stretch &st,
    const PAbstractShift &prevS) const {
  double d = 0;
  const auto &duals = dualValues_[nurseNum];
  for (const auto &p : resourceCons_[nurseNum]) {
    int c = computeConsumption(st, p.first, prevS);
    d += c * duals.at(p.first);
  }
  return d;
}

// add a given constraint to the column
template <class H, class S>
void BoundedResourceConstraint<H, S>::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Column &col) const {
  for (const auto &p : resourceCons_[col.nurseNum()]) {
    int c = computeConsumption(col, p.first, pScenario_->pRestShift());
    cons->push_back(p.second.first);
    cons->push_back(p.second.second);
    coeffs->push_back(c);
    coeffs->push_back(c);
  }
}

template <class H, class S>
std::string BoundedResourceConstraint<H, S>::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  const auto &duals = dualValues_[nurseNum];
  for (const auto &p : resourceCons_[nurseNum]) {
    const std::pair<H*, S*> &p1 = p.first;
    std::string name = p1.first ? p1.first->name : p1.second->name;
    buff << "#   |  " << name << ": "
         << -getDualCost(nurseNum, st, pScenario_->pRestShift()) << std::endl;
  }
  return buff.str();
}

// update the dual values of the constraints randomly
template <class H, class S>
void BoundedResourceConstraint<H, S>::createRandomUpdateDuals(double weight) {
  dualValues_.clear();
  for (const auto &constraints : resourceCons_) {
    std::map<std::pair<H*, S*>, double> duals;
    for (const auto &p : constraints)
      duals[p.first] = Tools::randomDouble(0, weight);
    dualValues_.push_back(duals);
  }
}

// compute bounds for constraint and slack
template <class H, class S>
std::pair<std::pair<int, int>, std::pair<int, int>>
BoundedResourceConstraint<H, S>::computeBounds(H *pHR, S *pSR) const {
  // compute bounds and slacks
  int max = pHR ? pHR->maxConsumptionPerDay() : pSR->maxConsumptionPerDay();
  max *= pMaster_->nDays();
  int lb = 0, ub = max, lb_slack = 0, ub_slack = 0;
  // update bounds
  if (pHR) {
    lb = pHR->getLb();
    ub = pHR->getUb();
  } else {
    lb_slack = std::max(0, pSR->getLb() - lb);
    if (pSR->getLb() > lb)
      lb = pSR->getLb();
    // if no hard constraint on the UB, the slack is not bounded
    ub_slack = ub < max ? std::max(0, ub - pSR->getUb()) : max;
    if (pSR->getUb() < ub)
      ub = pSR->getUb();
  }
  return {{lb, ub}, {lb_slack, ub_slack}};
}


// build a constraint for the LB and the UB of these bounded resources:
// pHR = pointer to a Hard Resource, pSR = pointer to a Soft Resource
// WARNING: The two resources MUST represent the same resource if both defined
// The goal of having both is to be able to define both soft and hard bounds
// on the same constraint
// If there is only a soft constraint (SR), set pHR to nullptr
// If there is only a hard constraint (HR), set pSR to nullptr
template <class H, class S>
void BoundedResourceConstraint<H, S>::build(
    H *pHR, S *pSR, const LiveNurse &pN) {
  // build both constraints
  char name[255];
  auto pBounds = computeBounds(pHR, pSR);
  std::pair<int, int> bounds = pBounds.first, slack_bounds = pBounds.second;

  /* create constraints and slack variables if soft constraint */
  const char *cName = pSR ? pSR->name.c_str() : pHR->name.c_str();
  MyVar * vLB = nullptr, *vUB = nullptr;
  MyCons *cLB = nullptr, *cUB = nullptr;

  // LB slack and constraint
  vector<double> lbCoeffs;
  vector<MyVar*> lbVars;
  if (pSR) {
    snprintf(name, sizeof(name), "%s_lb_N%d_slack", cName, pN.num_);
    pModel()->createPositiveVar(
        &vLB, name, pSR->getLbCost(), {}, 0, slack_bounds.first);
    lbCoeffs = {1};
    lbVars = {vLB};
  }
  snprintf(name, sizeof(name), "%s_lb_N%d", cName, pN.num_);
  pModel()->createGEConsLinear(&cLB, name, bounds.first, lbVars, lbCoeffs);

  // UB slack and constraint
  vector<double> ubCoeffs;
  vector<MyVar*> ubVars;
  if (pSR) {
    snprintf(name, sizeof(name), "%s_ub_N%d_slack", cName, pN.num_);
    pModel()->createPositiveVar(
        &vUB, name, pSR->getUbCost(), {}, 0, slack_bounds.second);
    ubCoeffs = {-1};
    ubVars = {vUB};
  }
  snprintf(name, sizeof(name), "%s_ub_N%d", cName, pN.num_);
  pModel()->createLEConsLinear(&cUB, name, bounds.second, ubVars, ubCoeffs);

  // store the variables and constraints
  resourceVars_[pN.num_][{pHR, pSR}] = {vLB, vUB};
  resourceCons_[pN.num_][{pHR, pSR}] = {cLB, cUB};
}

// update the bounds and costs of slack variables as well as
// bounds of the constraints
template <class H, class S>
void BoundedResourceConstraint<H, S>::update() {
  for (int n = 0; n < resourceCons_.size(); ++n) {
    for (const auto &p : resourceCons_[n]) {
      // compute bounds and slacks
      const std::pair<H*, S*> &p1 = p.first;
      auto pBounds = computeBounds(p1.first, p1.second);
      std::pair<int, int> bounds = pBounds.first, slack_bounds = pBounds.second;

      // update bounds
      const std::pair<MyCons*, MyCons*> &p2 = p.second;
      p2.first->setLhs(bounds.first);
      p2.second->setRhs(bounds.second);

      // update slack
      if (p1.second) {
        std::pair<MyVar *, MyVar *> slacks = resourceVars_[n][p1];
        slacks.first->setUB(slack_bounds.first);
        slacks.second->setUB(slack_bounds.second);
        slacks.first->setCost(p1.second->getLbCost());
        slacks.second->setCost(p1.second->getUbCost());
      }
    }
  }
}

TotalShiftDurationConstraint::TotalShiftDurationConstraint(
    MasterProblem *pMaster) :
    BoundedResourceConstraint<
        HardTotalShiftDurationResource, SoftTotalShiftDurationResource>(
            pMaster, "Worked duration") {}

void TotalShiftDurationConstraint::addConstraintFor(
    const shared_ptr<HardTotalShiftDurationResource> &pHR,
    const shared_ptr<SoftTotalShiftDurationResource> &pSR,
    const LiveNurse &pN) {
  build(pHR.get(), pSR.get(), pN);
}

// update the dual values of the constraints randomly
void TotalShiftDurationConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  createRandomUpdateDuals(2*pScenario_->weights().totalShifts);
}

TotalWeekendConstraint::TotalWeekendConstraint(MasterProblem *pMaster) :
    BoundedResourceConstraint<
        HardTotalWeekendsResource, SoftTotalWeekendsResource>(
                                  pMaster, "Worked weekends") {}

void TotalWeekendConstraint::addConstraintFor(
    const shared_ptr<HardTotalWeekendsResource> &pHR,
    const shared_ptr<SoftTotalWeekendsResource> &pSR,
    const LiveNurse &pN) {
  build(pHR.get(), pSR.get(), pN);
}

// update the dual values of the constraints randomly
void TotalWeekendConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  createRandomUpdateDuals(2*pScenario_->weights().totalWeekends);
}

// ConsWeekendConstraint::ConsWeekendConstraint(MasterProblem *pMaster) :
//    BoundedResourceConstraint<
//        HardConsWeekendShiftResource, SoftConsWeekendShiftResource>(
//        pMaster, "Worked weekends") {}
//
// void ConsWeekendConstraint::addConstraintFor(
//    const shared_ptr<HardConsWeekendShiftResource> &pHR,
//    const shared_ptr<SoftConsWeekendShiftResource> &pSR,
//    const LiveNurse &pN) {
//  build(pHR.get(), pSR.get(), pN);
// }
//
//// update the dual values of the constraints randomly
// void ConsWeekendConstraint::randomUpdateDuals(
//    bool useInputData, int nPerturbations) {
//  createRandomUpdateDuals(2*pScenario_->weights().totalWeekends);
//}
//
//// compute consumption
// int ConsWeekendConstraint::computeConsumption(
//    const Stretch &st,
//    const std::pair<HardConsWeekendShiftResource*,
//                    SoftConsWeekendShiftResource*> &p,
//    const PAbstractShift &prevS) const {
//  BoundedResource *pR = p.first;
//  if (p.first == nullptr) pR = p.second;
//  auto f = [pR](ResourceValues *vChild) {
//    if (vChild->consumption > 0) {
//      // initialize cyclicConsumption on the first reset
//      if (vChild->cyclicConsumption == -1)
//        vChild->cyclicConsumption = vChild->consumption;
//      // check the bounds
//      if (vChild->consumption < pR->getLb() ||
//          vChild->consumption > pR->getUb())
//        return false;
//    }
//    return true;
//  };
//
//  auto pS = p.first ? p.first->pShift() : p.second->pShift();
//  ResourceValues vC;
//  vC.readyToConsume = !(prevS && pS->includes(*prevS));
//  p.first ? p.first->computeConsumption(st, &vC, f) :
//  p.second->computeConsumption(st, &vC, f);
//  return vC.consumption;
//}
