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

int TotalWeekend::computeConsumption(const Stretch &stretch,
                                     bool *ready) const {
  // go through the days of the stretch and count the weekends that have not
  // been counted yet
  int consumption = 0;
  auto itShift = stretch.pShifts().begin();
  // reset weekend flag
  if (!Tools::isWeekend(stretch.firstDay()-1)) *ready = true;
  if (!Tools::isWeekend(stretch.firstDay())) *ready = true;
  for (int i = stretch.firstDay(); i <= stretch.lastDay(); i++, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if (*ready && pShift_->includes(**itShift)) {
        consumption++;
        *ready = false;
      }
      // reset weekend flag
      if (Tools::isLastWeekendDay(i))
        *ready = true;
    }
  }
  return consumption;
}

int TotalWeekend::computeBackConsumption(const Stretch &stretch,
                                         bool *ready) const {
  // go through the days of the stretch in reverse and
  // count the weekends that have not been counted yet
  int consumption = 0;
  auto itShift = stretch.pShifts().rbegin();
  int i = stretch.lastDay();
  for (; i >= stretch.firstDay(); i--, itShift++) {
    // if a weekend
    if (Tools::isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if (*ready && pShift_->includes(**itShift)) {
        consumption++;
        *ready = false;
      }
      // reset weekend flag
      if (Tools::isFirstWeekendDay(i))
        *ready = true;
    }
  }
  return consumption;
}

template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r,
                           const std::function<double(int)> &ubCost) {
  // check if expander is active
  std::pair<int, int> firstLastDays = r.getFirstLastDays(stretch);
  if (!Tools::nWeekendsInInterval(firstLastDays.first, firstLastDays.second))
    return nullptr;

  // Computing the number of weekends after the last day of the stretch
  int start = firstLastDays.second+1, end = r.firstDay() + r.totalNbDays() - 1;
  int nWeekendsAfter = Tools::nWeekendsInInterval(start, end);

  // Computing the number of weekends before the first day of the stretch
  start = r.firstDay(), end = firstLastDays.first - 1;
  int nWeekendsBefore = Tools::nWeekendsInInterval(start, end);

  int consumption = 0;
//      r.computeConsumption(firstLastDays.first, stretch, prevAShift);
//  if (consumption > r.getUb()) {
//    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
//    pArc->addBaseCost(ubCost(consumption));
//
//   // beware: we never need to store a consumption larger than the upper bound
//    consumption = r.getUb();
//  }

  return std::make_shared<E>(
      r,
      stretch,
      consumption,
      nWeekendsBefore,
      nWeekendsAfter,
      pArc->target->type == SINK_NODE);
}

int SoftTotalWeekendsResource::getConsumption(
    const State & initialState) const {
  return std::min(ub_, initialState.totalWeekendsWorked_);
}

PExpander SoftTotalWeekendsResource::init(const AbstractShift &prevAShift,
                                          const Stretch &stretch,
                                          const PRCArc &pArc) {
  return initExpander<SoftTotalWeekendsExpander, SoftTotalWeekendsResource>(
      prevAShift, stretch, pArc, *this, [this](int c) {
        return this->getUbCost(c);
      });
}

bool SoftTotalWeekendsExpander::expand(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  vChild->consumption +=
      resource_.computeConsumption(stretch_, &vChild->readyToConsume);

  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef DBG
    pLChild->addTotalWeekendCost(resource_.getUbCost(vChild->consumption));
#endif

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  if (arcToSink_) {
    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
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
  vChild->consumption +=
      resource_.computeBackConsumption(stretch_, &vChild->readyToConsume);

  // check UB cost at every node, but LB cost only when merging with the
  // initial label
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef DBG
    pLChild->addTotalWeekendCost(resource_.getUbCost(vChild->consumption));
#endif

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
      prevAShift, stretch, pArc, *this, nullptr);
}

bool HardTotalWeekendsExpander::expand(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  vChild->consumption +=
      resource_.computeConsumption(stretch_, &vChild->readyToConsume);

  if (vChild->consumption > resource_.getUb())
    return false;

  // if reach the end while remaining lower than LB -> return false
  // otherwise true
  return vChild->consumption + nWeekendsAfter_ >= resource_.getLb();
}

bool HardTotalWeekendsExpander::expandBack(const PRCLabel &pLChild,
                                           ResourceValues *vChild) {
  vChild->consumption +=
      resource_.computeBackConsumption(stretch_, &vChild->readyToConsume);

  return vChild->consumption <= resource_.getUb();
}
