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

#include "RCLabel.h"

#include <list>

#include "RCGraph.h"
#include "solvers/Solver.h"

std::map<CostType, double> initCostPerType() {
  std::map<CostType, double> costPerType;
  for (int t = NO_COST + 1; t < END_INDEX_COST; t++)
    costPerType[(CostType) t] = 0;
  return costPerType;
}

std::map<int, double> initCostPerIntType() {
  std::map<int, double> costPerType;
  for (int t = NO_COST + 1; t < END_INDEX_COST; t++)
    costPerType[t] = 0;
  return costPerType;
}

RCLabel::RCLabel(): num_(-1),
                    pNode_(nullptr),
                    pInArc_(nullptr),
                    pOutArc_(nullptr),
                    pPreviousLabel_(nullptr),
                    pNextLabel_(nullptr),
                    cost_(0),
                    baseCost_(0) {
  baseCost_ = 0;
#ifdef NS_DEBUG
  dualCost_ = 0;
  consShiftCost_ = 0;
  consWeekendShiftCost_ = 0;
  totalShiftCost_ = 0;
  totalWeekendCost_ = 0;
  preferencesCost_ = 0;
#endif
}

RCLabel::RCLabel(int nResources)
    : RCLabel() {
  resourceValues_.resize(nResources);
}

RCLabel::RCLabel(const vector<PResource> &resources):
    RCLabel(static_cast<int>(resources.size())) {
  int ind = 0;
  for (const auto& r : resources) {
    resourceValues_[ind].consumption = 0;
    ind++;
  }
}

RCLabel::RCLabel(const vector<PResource> &resources,
                 const State &initialState)
    : RCLabel(static_cast<int>(resources.size())) {
  int ind = 0;
  for (const auto& r : resources)
    resourceValues_[ind++].consumption = r->getConsumption(initialState);
}

RCLabel::RCLabel(const RCLabel &l): RCLabel() {
  copy(l);
}

void RCLabel::setAsNext(PRCLabel pLPrevious, const PRCArc &pArc) {
  copyValues(*pLPrevious);
  pNode_ = pArc->target;
  pInArc_ = pArc;
  pOutArc_ = nullptr;
  pPreviousLabel_ = std::move(pLPrevious);
  pNextLabel_ = nullptr;
}

void RCLabel::setAsPrevious(PRCLabel pLNext, const PRCArc &pArc) {
  copyValues(*pLNext);
  pNode_ = pArc->origin;
  pInArc_ = nullptr;
  pOutArc_ = pArc;
  pPreviousLabel_ = nullptr;
  pNextLabel_ = std::move(pLNext);
}

void RCLabel::setAsMerged(const PRCLabel &pLForward,
                          const PRCLabel &pLBackward) {
  pPreviousLabel_ = pLForward->getPreviousLabel();
  pNextLabel_ = pLBackward->getNextLabel();
  cost_ = pLForward->cost() + pLBackward->cost();
  baseCost_ = pLForward->baseCost() + pLBackward->baseCost();
#ifdef NS_DEBUG
  dualCost_ = pLForward->dualCost() + pLBackward->dualCost();
#endif
  resourceValues_ = pLForward->allResourceValues();
  for (auto& r : resourceValues_) {
    r.consumption = 0;
  }
  pNode_ = pLForward->getNode();
  pInArc_ = pLForward->getInArc();
  pOutArc_ = pLBackward->getOutArc();
}

void RCLabel::copy(const RCLabel &l) {
  copyValues(l);
  pNode_ = l.pNode_;
  pInArc_ = l.pInArc_;
  pOutArc_ = l.pOutArc_;
  pPreviousLabel_ = l.pPreviousLabel_;
  pNextLabel_ = l.pNextLabel_;
}

void RCLabel::copyValues(const RCLabel &l) {
  cost_ = l.cost_;
  baseCost_ = l.baseCost_;
  resourceValues_ = l.resourceValues_;
  baseCost_ = l.baseCost_;
#ifdef NS_DEBUG
  dualCost_ = l.dualCost_;
  consShiftCost_ = l.consShiftCost_;
  consWeekendShiftCost_ = l.consWeekendShiftCost_;
  totalShiftCost_ = l.totalShiftCost_;
  totalWeekendCost_ = l.totalWeekendCost_;
  preferencesCost_ = l.preferencesCost_;
#endif
}

std::string RCLabel::toString(const vector<PResource> &pResources) const {
  std::stringstream rep;
  rep.setf(std::ios::fixed,  std::ios::floatfield);
  rep.setf(std::ios::showpoint);
  rep << "Label: cost=" << cost() << ", baseCost=" << baseCost();
#ifdef NS_DEBUG
  if (consShiftCost_ > EPSILON)
    rep << ", consShiftCost="
        << std::setprecision(DECIMALS) << consShiftCost_;
  if (consWeekendShiftCost_ > EPSILON)
    rep << ", consWeekendShiftCost="
        << std::setprecision(DECIMALS) << consWeekendShiftCost_;
  if (totalShiftCost_ > EPSILON)
    rep << ", totalShiftCost="
        << std::setprecision(DECIMALS) << totalShiftCost_;
  if (totalWeekendCost_ > EPSILON)
    rep << ", totalWeekendCost="
        << std::setprecision(DECIMALS) << totalWeekendCost_;
#endif
  if (!pResources.empty()) {
    rep << ", Resources:";
    for (const auto &pR : pResources) {
      int c = resourceValues_[pR->id()].consumption;
      if (c)
        rep << " " << pR->name << "=" << c << ",";
    }
  }
  rep << std::endl;
  if (pOutArc_) rep << "Arc out: " << pOutArc_->toString() << std::endl;
  if (pInArc_) rep << "Arc in: " << pInArc_->toString() << std::endl;
  return rep.str();
}

std::string RCLabel::toStringRecursive(
    const vector<PResource> &pResources) const {
  std::stringstream rep;
  // if has any next labels (bidirectionnal has been used)
  if (pOutArc_ != nullptr) {
    const RCLabel *pLabel = this;
    int lvl = 1;
    while (pLabel->pNextLabel_.get() != nullptr) {
      rep << "+" << lvl << " -- " << pLabel->toString(pResources);
      pLabel = pLabel->pNextLabel_.get();
      lvl++;
    }
    rep << std::endl;
  }
  // iterate through previous labels
  if (pInArc_ != nullptr) {
    const RCLabel *pLabel = this;
    int lvl = 1;
    while (pLabel != nullptr) {
      rep << "-" << lvl << " -- " << pLabel->toString(pResources);
      pLabel = pLabel->pPreviousLabel_.get();
      lvl++;
    }
  }
  return rep.str();
}

std::string RCSolution::toString() const {
  std::stringstream rep;
  rep.setf(std::ios_base::fixed, std::ios_base::floatfield);
  rep << "RCSolution: ";
  if (cost_ < DBL_MAX - 1) rep << "cost=" << std::setprecision(0) << cost_;
  if (reducedCost_ < DBL_MAX - 1)
    rep << "  dualCost=" << std::setprecision(2) << reducedCost_;
  rep << "; " << Stretch::toString();
  return rep.str();
}

std::string RCSolution::costsToString() const {
  if (costs_.empty()) {
    std::cout << "WARNING: costs are not computed by cost types." << std::endl;
    return "";
  }
  std::stringstream rep;
  rep.setf(std::ios_base::fixed, std::ios_base::floatfield);
  rep.precision(0);
  for (const auto &p : costs_)
    if (abs(p.second) > 1e-1)
      rep << "#     | " <<  std::setw(25) << std::left
          << prettyNamesByCostType.at(p.first)
          << " : " << p.second << std::endl;
  return rep.str();
}

// Compare rotations on cost
bool RCSolution::compareCost(
    const RCSolution &sol1, const RCSolution &sol2) {
  if (sol1.cost_ == DBL_MAX || sol2.cost_ == DBL_MAX)
    Tools::throwError("Pattern cost not computed.");
  return (sol1.cost_ < sol2.cost_);
}

// Compare rotations on dual cost
bool RCSolution::compareReducedCost(
    const RCSolution &sol1, const RCSolution &sol2) {
  if (sol1.reducedCost_ == DBL_MAX || sol2.reducedCost_ == DBL_MAX)
    Tools::throwError("Pattern dual cost not computed.");
  return (sol1.reducedCost_ < sol2.reducedCost_);
}


PExpander Resource::initialize(const AbstractShift &prevAShift,
                               const Stretch &stretch,
                               const PRCArc &pArc,
                               int indResource) {
  PExpander pE = init(prevAShift, stretch, pArc, indResource);
  if (pE != nullptr)
    pArc->expanders.push_back(pE);
  return pE;
}

DominationStatus Resource::dominates(
    RCLabel *pL1,  RCLabel *pL2, double *cost) const {
  return pL1->getConsumption(id_) <= pL2->getConsumption(id_) ?
         DOMINATED : NOT_DOMINATED;
}

bool Resource::merge(const ResourceValues &vForward,
                     const ResourceValues &vBack,
                     ResourceValues *vMerged,
                     const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  return true;
}

PExpander Resource::init(const AbstractShift &prevAShift,
                         const Stretch &stretch,
                         const PRCArc &pArc,
                         int indResource) {
  return PExpander();
}

DominationStatus BoundedResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  double conso1 = pL1->getConsumption(id_), conso2 = pL2->getConsumption(id_);
  if (conso1 == conso2) return DOMINATED;
  if (conso1 < conso2) return conso1 >= lb_ ? DOMINATED : UB_DOMINATED;
  return NOT_DOMINATED;
}

DominationStatus SoftBoundedResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  if  (useDefaultDomination_ || cost == nullptr)
    return BoundedResource::dominates(pL1, pL2, cost);

  double ubDiff = pL1->getWorstUbCost(id_) - pL2->getWorstUbCost(id_),
      lbDiff = pL1->getWorstLbCost(id_) - pL2->getWorstLbCost(id_);
#ifdef NS_DEBUG
  if (((lbDiff > 0) && (ubDiff > 0)) || ((lbDiff < 0) && (ubDiff < 0)))
    Tools::throwError("ubDiff and lbDiff should never have the same sign");
#endif
  if (lbDiff <= ubDiff) {
    *cost += ubDiff;
    return DOMINATED;
  }
  *cost += lbDiff;
  return UB_DOMINATED;
}

bool SoftBoundedResource::merge(const ResourceValues &vForward,
                                const ResourceValues &vBack,
                                ResourceValues *vMerged,
                                const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  pLMerged->addBaseCost(getCost(vMerged->consumption));
  return true;
}

bool HardBoundedResource::merge(const ResourceValues &vForward,
                                const ResourceValues &vBack,
                                ResourceValues *vMerged,
                                const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  return (vMerged->consumption <= ub_) && (vMerged->consumption >= lb_);
}
