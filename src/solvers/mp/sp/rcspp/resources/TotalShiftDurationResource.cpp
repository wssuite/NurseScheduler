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

#include "TotalShiftDurationResource.h"

template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r) {
  // we only need to count the number of corresponding shift

  // Number of days left since the day of the target node of the arc
  int nDaysLeft = r.totalNbDays() - (stretch.firstDay() + stretch.nDays());
  int consumption = 0;
  for (const auto &pShift : stretch.pShifts())
    if (r.pShift()->includes(*pShift)) {
      if (r.hasDefaultDuration()) consumption += r.defaultDuration();
      else
        consumption += pShift->duration;
    }

  // if nothing happens
  if (stretch.nDays() == 0)
    return nullptr;

  return std::make_shared<E>(
      r, consumption, nDaysLeft, pArc->target->type == SINK_NODE);
}

int SoftTotalShiftDurationResource::getConsumption(
    const State &initialState) const {
  return initialState.totalTimeWorked_;
}

PExpander SoftTotalShiftDurationResource::init(const AbstractShift &prevAShift,
                                               const Stretch &stretch,
                                               const PRCArc &pArc) {
  return initExpander<SoftTotalShiftDurationExpander,
                      SoftTotalShiftDurationResource>(
      prevAShift, stretch, pArc, *this);
}

bool SoftTotalShiftDurationExpander::expand(const PRCLabel &pLChild,
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
    pLChild->addCost(resource_.getUbCost(vChild->consumption));
#ifdef DBG
    pLChild->addTotalShiftCost(resource_.getUbCost(vChild->consumption));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }
  if (arcToSink_ && vChild->consumption < resource_.getLb()) {
    pLChild->addCost(resource_.getLbCost(vChild->consumption));
#ifdef DBG
    pLChild->addTotalShiftCost(resource_.getLbCost(vChild->consumption));
#endif
  }

  // Setting 'worst case cost'
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                 nDaysLeft);

  return true;
}
bool SoftTotalShiftDurationExpander::expandBack(const PRCLabel &pLChild,
                                                ResourceValues *vChild) {
  if (consumption == 0) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                   resource_.totalNbDays() -
                                                       nDaysLeft);
    return true;
  }

  // update consumption
  vChild->consumption += consumption;

  // pay for excess of consumption due to this expansion
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addCost(resource_.getUbCost(vChild->consumption));

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }
  // check soft LB cost only when merging labels with the source label

  // Setting 'worst case cost'
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                 resource_.totalNbDays() -
                                                     nDaysLeft);

  return true;
}

bool SoftTotalShiftDurationExpander::merge(const ResourceValues &vForward,
                                           const ResourceValues &vBack,
                                           ResourceValues *vMerged,
                                           const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;
  pLMerged->addCost(resource_.getUbCost(vMerged->consumption));
  pLMerged->addCost(resource_.getLbCost(vMerged->consumption));
  return true;
}

int HardTotalShiftDurationResource::getConsumption(
    const State &initialState) const {
  return initialState.totalTimeWorked_;
}

PExpander HardTotalShiftDurationResource::init(const AbstractShift &prevAShift,
                                               const Stretch &stretch,
                                               const PRCArc &pArc) {
  return initExpander<HardTotalShiftDurationExpander,
                      HardTotalShiftDurationResource>(
      prevAShift, stretch, pArc, *this);
}

bool HardTotalShiftDurationExpander::expand(const PRCLabel &pLChild,
                                            ResourceValues *vChild) {
  // update consumption
  vChild->consumption += consumption;

  // if exceed UB
  if (vChild->consumption  > resource_.getUb())
    return false;

  // if reach the end while remaining lower than LB -> return false
  // otherwise true
  return !(arcToSink_ && vChild->consumption < resource_.getLb());
}

bool HardTotalShiftDurationExpander::expandBack(const PRCLabel &pLChild,
                                                ResourceValues *vChild) {
  // update consumption
  vChild->consumption += consumption;

  // if exceed UB return false, but do not check LB at this stage, only when
  // merging labels with the label at source
  return vChild->consumption <= resource_.getUb();
}
bool HardTotalShiftDurationExpander::merge(const ResourceValues &vForward,
                                           const ResourceValues &vBack,
                                           ResourceValues *vMerged,
                                           const PRCLabel &pLMerged) {
  vMerged->consumption = vForward.consumption + vBack.consumption;

  return (vMerged->consumption <= resource_.getUb()) &&
      (vMerged->consumption >= resource_.getLb());
}


