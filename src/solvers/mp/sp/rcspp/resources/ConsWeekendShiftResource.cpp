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

#include "ConsWeekendShiftResource.h"

void ConsWeekend::computeConsumption(
    const Stretch &stretch,
    ResourceValues *vChild,
    const std::function<bool(ResourceValues*)>& processReset) const {
  // go through the days of the stretch and count the weekends that have not
  // been counted yet
  auto itShift = stretch.pShifts().begin();
  // reset weekend flag
  if (!weekend_.isWeekend(stretch.firstDayId()) ||
      !weekend_.isWeekend(stretch.firstDayId() - 1))
    vChild->readyToConsume = true;
  for (const auto& pD : stretch.pDays()) {
    // if a weekend
    if (weekend_.isWeekend(pD)) {
      // check if included on this day and
      // if weekend has not been already counted
      if (vChild->readyToConsume && pAShift__->includes(**itShift)) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (weekend_.isLastWeekendDay(pD)) {
        // if there has been a consumption -> do no reset
        if (vChild->readyToConsume) {
          if (!processReset(vChild)) break;
          vChild->consumption = 0;
        }
        vChild->readyToConsume = true;
      }
    }
    itShift++;
  }
}

void ConsWeekend::computeConsumptionBack(
    const Stretch &stretch,
    ResourceValues *vChild,
    const std::function<bool(ResourceValues*)>& processReset) const {
  // go through the days of the stretch and count the weekends that have not
  // been counted yet
  // reset weekend flag
  if (!weekend_.isWeekend(stretch.lastDayId()) ||
      !weekend_.isWeekend(stretch.lastDayId() + 1))
    vChild->readyToConsume = true;
  auto itShift = stretch.pShifts().end();
  auto itDay = stretch.pDays().end();
  while (itDay != stretch.pDays().begin()) {
    itDay--; itShift--;
    // if a weekend
    if (weekend_.isWeekend(*itDay)) {
      // check if included on this day and
      // if weekend has not been already counted
      if (vChild->readyToConsume && pAShift__->includes(**itShift)) {
        vChild->consumption++;
        vChild->readyToConsume = false;
      }
      // reset weekend flag
      if (weekend_.isFirstWeekendDay(*itDay)) {
        // if there has been a consumption -> do no reset
        if (vChild->readyToConsume) {
          if (!processReset(vChild)) break;
          vChild->consumption = 0;
        }
        vChild->readyToConsume = true;
      }
    }
  }
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
  int start = firstLastDays.second + 1,
      end = r.firstDayId() + r.totalNbDays() - 1;
  int nWeekendsAfter = r.weekend().nWeekendsInInterval(start, end);

  // Computing the number of weekends before the first day of the stretch
  start = r.firstDayId(), end = firstLastDays.first - 1;
  int nWeekendsBefore = r.weekend().nWeekendsInInterval(start, end);

  return std::make_shared<E>(
      indResource, r, stretch, r.isCyclic(), nWeekendsBefore, nWeekendsAfter);
}

int SoftConsWeekendShiftResource::getConsumption(
    const State & initialState) const {
  if (pAShift_->isAnyWork())
    return std::min(ub_, initialState.consWeekendWorked_);
  if (pAShift_->isRest())
    return std::min(ub_, initialState.consWeekendOff_);
  return 0;
}

PExpander SoftConsWeekendShiftResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const shared_ptr<RCArc> &pArc,
                                             int indResource) {
  return initExpander<SoftConsWeekendShiftExpander,
                      SoftConsWeekendShiftResource>(
      prevAShift, stretch, pArc, *this, [this](int c) {
        return this->getCost(c);
      }, indResource);
}

bool SoftConsWeekendShiftResource::merge(
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

bool SoftConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  resource_.computeConsumption(
      stretch_, vChild,
      [pLChild, this](ResourceValues *vChild) {
        if (vChild->consumption > 0) {
          // initialize cyclicConsumption on the first reset
          if (vChild->cyclicConsumption == -1)
            vChild->cyclicConsumption = vChild->consumption;
          // pay the cost on the bounds
          pLChild->addBaseCost(resource_.getCost(vChild->consumption));
#ifdef NS_DEBUG
          pLChild->addConsWeekendShiftCost(
              resource_.getCost(vChild->consumption));
#endif
        }
        return true;
      });

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addConsWeekendShiftCost(resource_.getUbCost(vChild->consumption));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  // check if lower than LB at the end
  if (vChild->consumption < resource_.getLb() && vChild->consumption >= 1 &&
      nWeekendsAfter == 0 && resource_.lastDayEndsSequence()) {
    // pay cost on lb
    // pay for violations of soft bounds when at the end of horizon in INRC
    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addConsWeekendShiftCost(
        resource_.getLbCost(vChild->consumption));
#endif
  }

//  // if has consumed, set worst case and return
//  if (vChild->consumption > 0) {
  // Setting 'worst case cost'
  // WARNING: even a consumption of 0 could be dominated in the future
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  // add -1 to the number of weekends if currently in a weekend
  // that still has already been counted
  int nAfter = nWeekendsAfter - !vChild->readyToConsume;
  if (cyclic && vChild->cyclicConsumption > 0)
    nAfter += vChild->cyclicConsumption;
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption, nAfter);
//  if (cyclic && vChild->cyclicConsumption > 0)
//    vChild->worstUbCost += resource_.getUbCost(vChild->cyclicConsumption);
//  } else {
//    // set worst case costs to 0 as no consumption (has reset)
//    vChild->worstLbCost = .0;
//    vChild->worstUbCost = .0;
//  }

  return true;
}

bool SoftConsWeekendShiftExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  resource_.computeConsumptionBack(
      stretch_, vChild,
      [pLChild, this](ResourceValues *vChild) {
        if (vChild->consumption > 0) {
          // initialize cyclicConsumption on the first reset
          if (vChild->cyclicConsumption == -1)
            vChild->cyclicConsumption = vChild->consumption;
          // pay the cost on the bounds
          pLChild->addBaseCost(resource_.getCost(vChild->consumption));
#ifdef NS_DEBUG
          pLChild->addConsWeekendShiftCost(
              resource_.getCost(vChild->consumption));
#endif
        }
        return true;
      });

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsBefore == 0 && vChild->cyclicConsumption > 0)
    Tools::throwError("SoftConsWeekendShiftExpander::expandBack() "
                      "not implemented when cyclic enabled.");
//    vChild->consumption += vChild->cyclicConsumption;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef NS_DEBUG
    pLChild->addConsWeekendShiftCost(resource_.getUbCost(vChild->consumption));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  // check if lower than LB at the end
  // should never happen
//  if (vChild->consumption < resource_.getLb() && vChild->consumption >= 1 &&
//      nWeekendsBefore == 0 && resource_.lastDayEndsSequence()) {
//    // pay cost on lb
//    // pay for violations of soft bounds when at the end of horizon in INRC
//    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
// #ifdef NS_DEBUG
//    pLChild->addConsWeekendShiftCost(
//        resource_.getLbCost(vChild->consumption));
// #endif
//  }

//  // if has consumed, set worst case and return
//  if (vChild->consumption > 0) {
  // Setting 'worst case cost'
  // WARNING: even a consumption of 0 could be dominated in the future
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
  // add -1 to the number of weekends if currently in a weekend
  // that still has already been counted
  int nBefore = nWeekendsBefore - !vChild->readyToConsume;
  if (cyclic && vChild->cyclicConsumption > 0)
    nBefore += vChild->cyclicConsumption;
  vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption, nBefore);
//  if (cyclic && vChild->cyclicConsumption > 0)
//    vChild->worstUbCost += resource_.getUbCost(vChild->cyclicConsumption);
//  } else {
//    // set worst case costs to 0 as no consumption (has reset)
//    vChild->worstLbCost = .0;
//    vChild->worstUbCost = .0;
//  }

  return true;
}

int HardConsWeekendShiftResource::getConsumption(
    const State & initialState) const {
  if (pAShift_->isAnyWork())
    return initialState.consWeekendWorked_;
  if (pAShift_->isRest())
    return initialState.consWeekendOff_;
  return 0;
}

PExpander HardConsWeekendShiftResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const shared_ptr<RCArc> &pArc,
                                             int indResource) {
  return initExpander<HardConsWeekendShiftExpander,
                      HardConsWeekendShiftResource>(
      prevAShift, stretch, pArc, *this, nullptr, indResource);
}

DominationStatus HardConsWeekendShiftResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  if (!cyclic_) return BoundedResource::dominates(pL1, pL2, cost);
  const auto &v1 = pL1->getResourceValues(id_),
      &v2 = pL2->getResourceValues(id_);
  if (v1.consumption > v2.consumption) return NOT_DOMINATED;
  int c1 = v1.consumption + v1.cyclicConsumption,
      c2 = v2.consumption + v2.cyclicConsumption;
  if (v1.consumption == v2.consumption)
    return c1 <= c2 ? DOMINATED : NOT_DOMINATED;
  if (v1.consumption < lb_) return UB_DOMINATED;
  return c1 <= c2 ? DOMINATED : NOT_DOMINATED;
}


bool HardConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  resource_.computeConsumption(
      stretch_, vChild,
      [pLChild, this](ResourceValues *vChild) {
        if (vChild->consumption > 0) {
          // initialize cyclicConsumption on the first reset
          if (vChild->cyclicConsumption == -1)
            vChild->cyclicConsumption = vChild->consumption;
          // check the bounds
          if (vChild->consumption < resource_.getLb() ||
              vChild->consumption > resource_.getUb())
            return false;
        }
        return true;
      });

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // infeasible as it has exceeded the upper bounds
  if (vChild->consumption > resource_.getUb())
    return false;

  // check if a reset or if lower than LB at the end
  if (vChild->consumption < resource_.getLb() && vChild->consumption >= 1)
    if ((nWeekendsAfter == 0 && resource_.lastDayEndsSequence()))
      return false;

  return true;
}

bool HardConsWeekendShiftExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  Tools::throwError(
      "HardConsWeekendShiftExpander::expandBack() not implemented.");
  return false;
}
