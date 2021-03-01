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

#include "MyRCLabel.h"
#include "MyRCGraph.h"

bool BoundedResource::dominates(int conso1,
                                int conso2) {
  if (conso1 == conso2) return true;
  if (conso1 <= conso2) return conso1 >= lb_;
  return false;
}

RCLabel::RCLabel(): num_(-1),
                    pNode_(nullptr),
                    pInArc_(nullptr),
                    pParentLabel_(nullptr),
                    cost_(0),
                    baseCost_(0),
                    dualCost_(0),
                    consShiftCost_(0),
                    totalShiftCost_(0),
                    totalWeekendCost_(0)  { }

RCLabel::RCLabel(int nResLabels)
    : RCLabel() {
  nResLabels_ = nResLabels;
  resourceValues_.resize(nResLabels_);
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

void RCLabel::copy(const RCLabel &l) {
  cost_ = l.cost_;
  consShiftCost_ = l.consShiftCost_;
  totalShiftCost_ = l.totalShiftCost_;
  totalWeekendCost_ = l.totalWeekendCost_;
  baseCost_ = l.baseCost_;
  dualCost_ = l.dualCost_;
  pNode_ = l.pNode_;
  pInArc_ = l.pInArc_;
  pParentLabel_ = l.pParentLabel_;
  nResLabels_ = l.nResLabels_;
  resourceValues_ = l.resourceValues_;
}

void RCLabel::setInArc(shared_ptr<RCArc> pArc) {
  pNode_ = pArc->target;
  pInArc_ = std::move(pArc);
}

void Resource::initialize(const Shift &prevShift,
                          const Stretch &stretch,
                          const PRCArc &pArc) {
  PExpander pE = init(prevShift, stretch, *pArc);
  if (pE != nullptr)
    pArc->expanders.push_back(std::move(pE));
}
