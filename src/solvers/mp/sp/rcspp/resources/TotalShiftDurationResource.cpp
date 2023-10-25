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

int TotalShiftDuration::computeConsumption(
    const Stretch &stretch, bool *ready) const {
  int consumption = 0;
  for (const auto &pShift : stretch.pShifts())
    if (pAShift__->includes(*pShift))
      consumption += duration(pShift);
  return consumption;
}

template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r,
                           const int indResource) {
  // if nothing happens
  if (stretch.nDays() == 0)
    return nullptr;

  // we only need to count the number of corresponding shift
  int consumption = r.computeConsumption(stretch);

  // number of days before the start of the stretch (beware that indices of
  // days start at 0)
  std::pair<int, int> firstLastDays = r.getFirstLastDays(stretch);
  int nDaysBefore = firstLastDays.first - r.firstDayId();
  // Number of days left since the day of the target node of the arc
  int nDaysLeft = r.totalNbDays() + r.firstDayId() - firstLastDays.second - 1;

  return std::make_shared<E>(indResource,
      r, consumption, nDaysBefore * r.maxDuration(),
      nDaysLeft  * r.maxDuration(), pArc->target->type == SINK_NODE);
}

int SoftTotalShiftDurationResource::getConsumption(
    const State &initialState) const {
  // always 0, as bounds have been modified according to this value
  return 0;
}

PExpander SoftTotalShiftDurationResource::init(const AbstractShift &prevAShift,
                                               const Stretch &stretch,
                                               const shared_ptr<RCArc> &pArc,
                                               int indResource) {
  if (isPreprocessed_) return nullptr;
  return initExpander<SoftTotalShiftDurationExpander,
                      SoftTotalShiftDurationResource>(
      prevAShift, stretch, pArc, *this, indResource);
}

void SoftTotalShiftDurationResource::preprocess(const PRCGraph &pRCGraph) {
  if (ub_ > 0) return;

  for (const PRCArc& pA : pRCGraph->pArcs()) {
    double cost = 0;
    preprocess(pA, &cost);
    pA->addBaseCost(cost);
  }
  isPreprocessed_ = true;
}

bool SoftTotalShiftDurationResource::preprocess(
    const PRCArc& pA, double *cost) {
  if (ub_ > 0) return false;
  bool ready = true;
  int c = computeConsumption(pA->stretch, &ready);
  *cost = c * ubCost_;
  return true;
}


bool SoftTotalShiftDurationExpander::expand(const PRCLabel &pLChild,
                                            ResourceValues *vChild) {
  if (consumption_ == 0 && !arcToSink_) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                   maxDurationLeft_);
    return true;
  }

  // update consumption
  vChild->consumption += consumption_;

  // pay for excess of consumption due to this expansion
  if (vChild->consumption  > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addTotalShiftCost(resource_.getUbCost(vChild->consumption));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }
  if (arcToSink_ && vChild->consumption < resource_.getLb()) {
    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addTotalShiftCost(resource_.getLbCost(vChild->consumption));
#endif
  }
  // Setting 'worst case cost'
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption,
                                                 maxDurationLeft_);

  return true;
}

bool SoftTotalShiftDurationExpander::expandBack(const PRCLabel &pLChild,
                                                ResourceValues *vChild) {
  // update consumption
  vChild->consumption += consumption_;

  // pay for excess of consumption due to this expansion
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));

    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }
  // check soft LB cost only when merging labels with the source label

  // Setting 'worst case cost'
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  vChild->worstUbCost =
      resource_.getWorstUbCost(vChild->consumption, maxDurationBefore_);

  return true;
}

int HardTotalShiftDurationResource::getConsumption(
    const State &initialState) const {
  // always 0, as bounds have been modified according to this value
  return 0;
}

DominationStatus HardTotalShiftDurationResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  int c1 = pL1->getConsumption(id_), c2 = pL2->getConsumption(id_);
  if (c1 == c2) return DOMINATED;
  if (c1 < c2) return c1 >= lb_ ? DOMINATED : UB_DOMINATED;
  // do not work if counting rest shifts
  if (pAShift()->isRest()) return NOT_DOMINATED;
  // check if c1 could be exceeded before the end of horizon
  int day = pL1->getNode()->dayId;
  int remainingDays =
      (pL1->getOutArc() == nullptr) ? totalNbDays_ - day: day + 1;
  int remainingDuration = maxDurationForRemainingDays_[remainingDays];
  // if max possible duration worked is lower than the remaining consumption
  // available, the resource ub cannot be violated
  if (remainingDuration <= (ub_ - c1))
    return DOMINATED;
  return NOT_DOMINATED;
}

void HardTotalShiftDurationResource::maxDurationForRemainingDays(
    int maxConsWorkedDays, int minConsRestDays) {
  maxDurationForRemainingDays_.clear();
  int nLength = maxConsWorkedDays + minConsRestDays;
  for (int d=0; d <= totalNbDays_; d++) {
    // max number of remaining work days
    int nSequences = d / nLength;
    int nLeft = d - nSequences * nLength;
    int remainingW = nSequences * maxConsWorkedDays +
        std::min(nLeft, maxConsWorkedDays);
    maxDurationForRemainingDays_.push_back(remainingW * maxDuration_);
  }
}

PExpander HardTotalShiftDurationResource::init(const AbstractShift &prevAShift,
                                               const Stretch &stretch,
                                               const shared_ptr<RCArc> &pArc,
                                               int indResource) {
  return initExpander<HardTotalShiftDurationExpander,
                      HardTotalShiftDurationResource>(
      prevAShift, stretch, pArc, *this, indResource);
}

bool HardTotalShiftDurationExpander::expand(const PRCLabel &pLChild,
                                            ResourceValues *vChild) {
  // update consumption
  vChild->consumption += consumption_;

  // if exceed UB
  if (vChild->consumption  > resource_.getUb())
    return false;

  // it will reach the end while remaining lower than LB -> return false
  // otherwise true
  return vChild->consumption + maxDurationLeft_ >= resource_.getLb();
}

bool HardTotalShiftDurationExpander::expandBack(const PRCLabel &pLChild,
                                                ResourceValues *vChild) {
  // update consumption
  vChild->consumption += consumption_;

  // if exceed UB return false, but do not check LB at this stage, only when
  // merging labels with the label at source
  return vChild->consumption <= resource_.getUb();
}
