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

#include "TotalShiftResource.h"

int SoftTotalShiftResource::getConsumption(const State &initialState) const {
  return initialState.totalTimeWorked_;
}

PExpander SoftTotalShiftResource::init(const Shift &prevShift,
                                       const Stretch &stretch,
                                       const RCArc &arc) {
  // we only need to count the number of corresponding shift

  // Number of days left since the day of the target node of the arc
  int nDaysLeft = totalNbDays_ - (stretch.firstDay() + stretch.nDays());
  int consumption = 0;
  for (auto pShift : stretch.pShifts())
    if (shift_->includes(*pShift))
      consumption += 1;

  // if nothing happens
  if (stretch.nDays() == 0)
    return nullptr;

  if (arc.target->type == SINK_NODE)
    return std::make_shared<SoftTotalShiftExpander>(
        *this, consumption, nDaysLeft, true);

  return std::make_shared<SoftTotalShiftExpander>(
      *this, consumption, nDaysLeft);
}

bool SoftTotalShiftExpander::expand(const ResourceValues &vParent,
                                    const PRCLabel &pLChild,
                                    ResourceValues *vChild) {
  if (consumption == 0) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                   nDaysLeft);
    return true;
  }

  // update consumption
  vChild->consumption += consumption;

  // pay for excess of consumption due to this expansion
  if (vChild->consumption  > resource_.getUb()) {
    pLChild->addTotalShiftCost(resource_.getUbCost(vChild->consumption));
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }
  if (arcToSink && vChild->consumption < resource_.getLb())
    pLChild->addTotalShiftCost(resource_.getLbCost(vChild->consumption));

  // Setting 'worst case cost'
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                 nDaysLeft);

  return true;
}


