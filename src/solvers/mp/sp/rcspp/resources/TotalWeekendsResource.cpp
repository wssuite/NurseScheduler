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
  if (!weekend_.isWeekend(stretch.firstDayId())) *ready = true;
  else if (!weekend_.isWeekend(stretch.firstDayId() - 1)) *ready = true;
  for (const auto& pD : stretch.pDays()) {
    // if a weekend
    if (weekend_.isWeekend(pD)) {
      // check if included on this day and
      // if weekend has not been already counted
      if (*ready && pAShift_->includes(**itShift)) {
        consumption++;
        *ready = false;
      }
      // reset weekend flag
      if (weekend_.isLastWeekendDay(pD))
        *ready = true;
    }
    itShift++;
  }
  return consumption;
}

int TotalWeekend::computeBackConsumption(const Stretch &stretch,
                                         bool *ready) const {
  // go through the days of the stretch in reverse and
  // count the weekends that have not been counted yet
  int consumption = 0;
  auto itShift = stretch.pShifts().rbegin();
  int i = stretch.lastDayId();
  for (; i >= stretch.firstDayId(); i--, itShift++) {
    // if a weekend
    if (weekend_.isWeekend(i)) {
      // check if working on this day and
      // if weekend has not been already counted
      if (*ready && pAShift_->includes(**itShift)) {
        consumption++;
        *ready = false;
      }
      // reset weekend flag
      if (weekend_.isFirstWeekendDay(i))
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
                           const std::function<double(int)> &getCost,
                           const int indResource) {
  // check if expander is active
  std::pair<int, int> firstLastDays = r.getFirstLastDays(stretch);
  if (!stretch.nWeekends(r.weekend()))
    return nullptr;

  // Computing the number of weekends after the last day of the stretch
  int start = firstLastDays.second+1, end =
      r.firstDayId() + r.totalNbDays() - 1;
  int nWeekendsAfter = r.weekend().nWeekendsInInterval(start, end);

  // Computing the number of weekends before the first day of the stretch
  start = r.firstDayId(), end = firstLastDays.first - 1;
  int nWeekendsBefore = r.weekend().nWeekendsInInterval(start, end);
  int consumption = 0;
  // TODO(JO): we should certainly put back at least a part of the code below
  //  in the initialization. The only stretches for which it is not correct
  //  are those that start from a weekend day other than the first, and even
  //  in this case, we never need to count more than the first weekend met in
  //  the stretch.

//      r.computeConsumption(firstLastDays.first, stretch, prevAShift);
//  if (consumption > r.getUb()) {
//    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
//    pArc->addBaseCost(ubCost(consumption));
//
//   // beware: we never need to store a consumption larger than the upper bound
//    consumption = r.getUb();
//  }

  return std::make_shared<E>(
      indResource,
      r,
      stretch,
      consumption,
      nWeekendsBefore,
      nWeekendsAfter,
      pArc->target->type == SINK_NODE);
}

int SoftTotalWeekendsResource::getConsumption(
    const State & initialState) const {
  // always 0, as bounds have been modified according to this value
  return 0;
}

PExpander SoftTotalWeekendsResource::init(const AbstractShift &prevAShift,
                                          const Stretch &stretch,
                                          const shared_ptr<RCArc> &pArc,
                                          int indResource) {
  return initExpander<SoftTotalWeekendsExpander, SoftTotalWeekendsResource>(
      prevAShift, stretch, pArc, *this, [this](int c) {
        return this->getUbCost(c);
      }, indResource);
}

bool SoftTotalWeekendsResource::merge(
    const ResourceValues &vForward,
    const ResourceValues &vBack,
    ResourceValues *vMerged,
    const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  // if both not consuming
  if (vMerged->consumption == 0)
    return true;

  // behave differently if merging on a weekend
  bool mergeOnWeekend =
      weekend().isWeekendDayButNotLastOne(pLMerged->getNode()->dayId);

  // if merging on a weekend, remove one weekend if already counted twice
  if (mergeOnWeekend && !vForward.readyToConsume && !vBack.readyToConsume)
    --vMerged->consumption;

  pLMerged->addBaseCost(getCost(vMerged->consumption));
  return true;
}

bool SoftTotalWeekendsExpander::expand(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  vChild->consumption +=
      resource_.computeConsumption(stretch_, &vChild->readyToConsume);

  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addTotalWeekendCost(resource_.getUbCost(vChild->consumption));
#endif

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  if (arcToSink_) {
    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
#ifdef NS_DEBUG
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
#ifdef NS_DEBUG
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
  // always 0, as bounds have been modified according to this value
  return 0;
}

PExpander HardTotalWeekendsResource::init(const AbstractShift &prevAShift,
                                          const Stretch &stretch,
                                          const shared_ptr<RCArc> &pArc,
                                          int indResource) {
  return initExpander<HardTotalWeekendsExpander, HardTotalWeekendsResource>(
      prevAShift, stretch, pArc, *this, nullptr, indResource);
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
