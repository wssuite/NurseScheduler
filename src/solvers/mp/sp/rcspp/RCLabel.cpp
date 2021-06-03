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

#include "RCGraph.h"
#include "solvers/Solver.h"



RCLabel::RCLabel(): num_(-1),
                    pNode_(nullptr),
                    pInArc_(nullptr),
                    pOutArc_(nullptr),
                    pPreviousLabel_(nullptr),
                    pNextLabel_(nullptr),
                    cost_(0),
                    baseCost_(0) {
#ifdef DBG
  baseCost_ = 0;
  dualCost_ = 0;
  consShiftCost_ = 0;
  consWeekendShiftCost_ = 0;
  totalShiftCost_ = 0;
  totalWeekendCost_ = 0;
#endif
}

RCLabel::RCLabel(int nResources)
    : RCLabel() {
  resourceValues_.resize(nResources);
}

RCLabel::RCLabel(const vector<shared_ptr<Resource>> &resources,
                 const State &initialState)
    : RCLabel(resources.size()) {
  for (const auto& r : resources)
    resourceValues_[r->id()].consumption = r->getConsumption(initialState);
}

RCLabel::RCLabel(const RCLabel &l): RCLabel() {
  copy(l);
}

void RCLabel::setAsNext(const PRCLabel &pLPrevious, const PRCArc &pArc) {
  copyValues(*pLPrevious);
  pNode_ = pArc->target;
  pInArc_ = pArc;
  pOutArc_ = nullptr;
  pPreviousLabel_ = pLPrevious;
  pNextLabel_ = nullptr;
}

void RCLabel::setAsPrevious(const shared_ptr<RCLabel> &pLNext,
                            const shared_ptr<RCArc> &pArc) {
  copyValues(*pLNext);
  pNode_ = pArc->origin;
  pInArc_ = nullptr;
  pOutArc_ = pArc;
  pPreviousLabel_ = nullptr;
  pNextLabel_ = pLNext;
}

void RCLabel::setAsMerged(const shared_ptr<RCLabel> &pLForward,
                          const shared_ptr<RCLabel> &pLBackward) {
  pPreviousLabel_ = pLForward->getPreviousLabel();
  pNextLabel_ = pLBackward->getNextLabel();
  cost_ = pLForward->cost() + pLBackward->cost();
  baseCost_ = pLForward->baseCost() + pLBackward->baseCost();
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
#ifdef DBG
  baseCost_ = l.baseCost_;
  dualCost_ = l.dualCost_;
  consShiftCost_ = l.consShiftCost_;
  consWeekendShiftCost_ = l.consWeekendShiftCost_;
  totalShiftCost_ = l.totalShiftCost_;
  totalWeekendCost_ = l.totalWeekendCost_;
#endif
}

std::string RCLabel::toString(const vector<PResource> &pResources) const {
  std::stringstream rep;
  rep << "Label: cost=" << cost() << ", baseCost=" << baseCost();
#ifdef DBG
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
  if (pInArc_) rep << "Arc in: " << pInArc_->toString() << std::endl;
  if (pOutArc_) rep << "Arc out: " << pOutArc_->toString() << std::endl;
  return rep.str();
}

PExpander Resource::initialize(const AbstractShift &prevAShift,
                               const Stretch &stretch,
                               const PRCArc &pArc) {
  PExpander pE = init(prevAShift, stretch, pArc);
  if (pE != nullptr)
    pArc->expanders.push_back(pE);
  return pE;
}

bool Resource::dominates(const PRCLabel &pL1,
                         const PRCLabel &pL2,
                         double *cost) {
  return pL1->getConsumption(id_) <= pL2->getConsumption(id_);
}

bool Resource::merge(const ResourceValues &vForward,
                     const ResourceValues &vBack,
                     ResourceValues *vMerged,
                     const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  return true;
}

bool BoundedResource::dominates(const PRCLabel &pL1,
                                const PRCLabel &pL2,
                                double *cost) {
  double conso1 = pL1->getConsumption(id_), conso2 = pL2->getConsumption(id_);
  if (conso1 == conso2) return true;
  if (conso1 < conso2) return conso1 >= lb_;
  return false;
}

bool SoftBoundedResource::dominates(const PRCLabel &pL1,
                                    const PRCLabel &pL2,
                                    double *cost) {
  if  (useDefaultDomination_)
    return BoundedResource::dominates(pL1, pL2, cost);

  double ubDiff = pL1->getWorstUbCost(id_) - pL2->getWorstUbCost(id_),
      lbDiff = pL1->getWorstLbCost(id_) - pL2->getWorstLbCost(id_);
#ifdef DBG
  if (((lbDiff > 0) && (ubDiff > 0)) || ((lbDiff < 0) && (ubDiff < 0)))
        Tools::throwError("ubDiff and lbDiff should never have the same sign");
#endif
  *cost += std::max(ubDiff, lbDiff);
  return true;
}

bool SoftBoundedResource::merge(const ResourceValues &vForward,
                                const ResourceValues &vBack,
                                ResourceValues *vMerged,
                                const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  pLMerged->addBaseCost(getUbCost(vMerged->consumption));
  pLMerged->addBaseCost(getLbCost(vMerged->consumption));
  return true;
}

bool HardBoundedResource::merge(const ResourceValues &vForward,
                                const ResourceValues &vBack,
                                ResourceValues *vMerged,
                                const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  return (vMerged->consumption <= ub_) && (vMerged->consumption >= lb_);
}