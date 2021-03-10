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


template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r) {
  // check if expander is active
  if (!Tools::nWeekendsInInterval(stretch.firstDay(), stretch.lastDay()))
    return nullptr;

  // Computing the number of weekends after the last day of the stretch
  int start = stretch.lastDay()+1, end = r.totalNbDays()-1;
  int nWeekendsAfter = Tools::nWeekendsInInterval(start, end);

  // Computing the number of weekends before the first day of the stretch
  start = 0, end = stretch.firstDay() - 1;
  int nWeekendsBefore = Tools::nWeekendsInInterval(start, end);

  return std::make_shared<E>(
      r,
      stretch,
      nWeekendsBefore,
      nWeekendsAfter,
      pArc->target->type == SINK_NODE);
}

int SoftTotalWeekendsResource::getConsumption(
    const State & initialState) const {
  return initialState.totalWeekendsWorked_;
}

PExpander SoftTotalWeekendsResource::init(const AbstractShift &prevAShift,
                                          const Stretch &stretch,
                                          const PRCArc &pArc) {
  return initExpander<SoftTotalWeekendsExpander, SoftTotalWeekendsResource>(
      prevAShift, stretch, pArc, *this);
}

bool SoftTotalWeekendsExpander::expand(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  // go through the days of the stretch and count the weekends that have not
  // been counted yet
  auto itShift = stretch_.pShifts().begin();
  for (int i = stretch_.firstDay(); i <= stretch_.lastDay(); i++, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if ((*itShift)->isWork() && vChild->readyToConsume) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (Tools::isLastWeekendDay(i))
        vChild->readyToConsume = true;
    }
  }
#ifdef DBG
  pLChild->addTotalWeekendCost(resource_.getUbCost(vChild->consumption));
#endif
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addCost(resource_.getUbCost(vChild->consumption));

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  if (arcToSink_) {
    pLChild->addCost(resource_.getLbCost(vChild->consumption));
#ifdef DBG
    pLChild->addTotalWeekendCost(resource_.getLbCost(vChild->consumption));
#endif
  }
  // Setting 'worst case cost'
  // if the resource is not ready to be consumed, it means that the current
  // weekend must not be counted in the potentially remaining weekends
  vChild->worstUbCost =
      resource_.getWorstUbCost(vChild->consumption, nWeekendsAfter_);

  return true;
}


bool SoftTotalWeekendsExpander::expandBack(const PRCLabel &pLChild,
                                           ResourceValues *vChild) {
  // The number of non-working weekends remaining is either given by the
  // number of Saturdays or given by the number of Sundays
  auto itShift = stretch_.pShifts().rbegin();
  int i = stretch_.lastDay();
  for (; i >= stretch_.firstDay(); i--, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if ((*itShift)->isWork() && vChild->readyToConsume) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (Tools::isFirstWeekendDay(i))
        vChild->readyToConsume = true;
    }
  }
  // check UB cost at every node, but LB cost only when merging with the
  // initial label
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addCost(resource_.getUbCost(vChild->consumption));

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  // Setting worst-case costs
  // if the resource is not ready to be consumed, it means that the current
  // weekend must not be counted in the potentially remaining weekends
  vChild->worstUbCost =
      resource_.getWorstUbCost(vChild->consumption, nWeekendsBefore_);

  return true;
}

int HardTotalWeekendsResource::getConsumption(
    const State & initialState) const {
  return initialState.totalWeekendsWorked_;
}

PExpander HardTotalWeekendsResource::init(const AbstractShift &prevAShift,
                                          const Stretch &stretch,
                                          const PRCArc &pArc) {
  return initExpander<HardTotalWeekendsExpander, HardTotalWeekendsResource>(
      prevAShift, stretch, pArc, *this);
}

bool HardTotalWeekendsExpander::expand(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  // go through the days of the stretch and count the weekends that have not
  // been counted yet
  auto itShift = stretch_.pShifts().begin();
  for (int i = stretch_.firstDay(); i <= stretch_.lastDay(); i++, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if ((*itShift)->isWork() && vChild->readyToConsume) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (Tools::isLastWeekendDay(i))
        vChild->readyToConsume = true;
    }
  }

  if (vChild->consumption > resource_.getUb())
    return false;

  // if reach the end while remaining lower than LB -> return false
  // otherwise true
  return vChild->consumption + nWeekendsAfter_ >= resource_.getLb();
}

bool HardTotalWeekendsExpander::expandBack(const PRCLabel &pLChild,
                                           ResourceValues *vChild) {
  // The number of non-working weekends remaining is either given by the
  // number of Saturdays or given by the number of Sundays
  auto itShift = stretch_.pShifts().rbegin();
  int i = stretch_.lastDay();
  for (; i >= stretch_.firstDay(); i--, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if ((*itShift)->isWork() && vChild->readyToConsume) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (Tools::isFirstWeekendDay(i))
        vChild->readyToConsume = true;
    }
  }

  return vChild->consumption <= resource_.getUb();
}
