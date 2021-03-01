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

#include "TotalWeekendsResource.h"

int SoftTotalWeekendsResource::getConsumption(
    const State & initialState) const {
  return initialState.totalWeekendsWorked_;
}

PExpander SoftTotalWeekendsResource::init(const Shift &prevShift,
                                          const Stretch &stretch,
                                          const RCArc &arc) {
  // check if expander is active
  if (!Tools::containsWeekend(
      stretch.firstDay(), stretch.firstDay() + stretch.nDays() - 1))
    return nullptr;

  // Computing the number of Sundays and Saturdays left from the end of the
  // stretch
  int start = stretch.firstDay()+stretch.nDays(), end = totalNbDays()-1;

  int nWeekendsLeft = Tools::containsWeekend(start,  end);
  // remove one saturday if start on a sunday
  int nSaturdaysLeft = nWeekendsLeft - Tools::isSunday(start);
  // remove one sunday if end on a saturday
  int nSundaysLeft = nWeekendsLeft - Tools::isSaturday(end);

  return std::make_shared<SoftTotalWeekendsExpander>(
      *this,
      stretch,
      nSaturdaysLeft,
      nSundaysLeft,
      arc.target->type == SINK_NODE);
}

bool SoftTotalWeekendsExpander::expand(const ResourceValues &vParent,
                                       const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  int nbWeekendsLeft = nSundaysLeft_;

  // The number of non-working weekends remaining is either given by the
  // number of Saturdays or given by the number of Sundays
  auto itShift = stretch_.pShifts().begin();
  for (int i = stretch_.firstDay();
       i < stretch_.firstDay()+stretch_.nDays();
       i++, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if ((*itShift)->isWork() && vChild->readyToConsume) {
        vChild->consumption++;
        vChild->readyToConsume = false;
        nbWeekendsLeft = nSaturdaysLeft_;
      }
      // reset weekend flag
      if (Tools::isLastWeekendDay(i))
        vChild->readyToConsume = true;
    }
  }

  if (vChild->consumption > resource_.getUb()) {
    pLChild->addTotalWeekendCost(resource_.getUbCost(vChild->consumption));
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  if (arcToSink_)
    pLChild->addTotalWeekendCost(resource_.getLbCost(vChild->consumption));

  // Setting 'worst case cost'
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                 nbWeekendsLeft);

  return true;
}

bool SoftTotalWeekendsResource::dominates(int conso1,
                                          int conso2) {
  return BoundedResource::dominates(conso1, conso2);
}



