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


template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r,
                           const std::function<double(int)> &getCost) {
  // check if expander is active
  std::pair<int, int> firstLastDays = r.getFirstLastDays(stretch);
  if (!Tools::nWeekendsInInterval(firstLastDays.first, firstLastDays.second))
    return nullptr;

  // we need to count the number of times the considered shift appears in
  // the arc's stretch on a weekend
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  int day = firstLastDays.first;
  bool weekendCounted =
      Tools::isWeekend(day-1) && r.pShift()->includes(prevAShift);
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      // on a weekend and same shift not counted -> increment consAfterReset
      if (r.pShift()->includes(**itShift) && !weekendCounted) {
          consBeforeReset++;
          weekendCounted = true;
      }
      if (Tools::isLastWeekendDay(day)) {
        // On last weekend day check if some bounds should be checked
        if (!weekendCounted) {
          // different shifts during whole weekend
          reset = true;
          break;
        } else {
          // reset flag
          weekendCounted = false;
        }
      }
    }
  }

  // if already exceeding UB
  if (consBeforeReset > r.getUb()) {
    // either infeasible
    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
    // if feasible -> pay immediately the current excess
    cost += getCost(consBeforeReset);
    consBeforeReset =  r.getUb();
  }

  // then compute the penalty due to consecutive shifts inside the stretch
  // (should be a monday)
  weekendCounted = false;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      // on a weekend and same shift not counted -> increment consAfterReset
      if (r.pShift()->includes(**itShift) && !weekendCounted) {
          consAfterReset++;
          weekendCounted = true;
      }
      if (Tools::isLastWeekendDay(day)) {
        // On last weekend day check if some bounds should be checked
        // if different shifts during whole weekend
        if (!weekendCounted) {
          if (consAfterReset == 0) continue;
          if (consAfterReset < r.getLb() || consAfterReset > r.getUb()) {
            if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
            cost += getCost(consAfterReset);
          }
          consAfterReset = 0;
        } else {
          // reset flag
          weekendCounted = false;
        }
      }
    }
  }

  // if nothing before and after reset and no reset (price previous consumption)
  // -> initialize nothing
  if (consBeforeReset == 0 && consAfterReset == 0 && !reset)
    return nullptr;

  // Computing the number of weekends after the last day of the stretch
  int start = firstLastDays.second+1, end = r.totalNbDays()-1;
  int nWeekendsAfter = Tools::nWeekendsInInterval(start, end);
  // check if we need to remove current weekend
  if (Tools::isWeekendDayButNotLastOne(firstLastDays.second))
    if (r.pShift()->includes(*stretch.pShifts().back()))
      --nWeekendsAfter;

  if (r.isCyclic() && nWeekendsAfter == 0)
    consAfterReset += consBeforeReset;

  if (consAfterReset > r.getUb()) {
    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
    cost += getCost(consAfterReset);
    consAfterReset = r.getUb();
  }

  // add cost to base cost of the arc (cost != 0 if resource is soft
  pArc->addBaseCost(cost);

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<E>(
      r, false, reset, consBeforeReset, consAfterReset,
      r.isCyclic(), false, nWeekendsAfter);
}

int SoftConsWeekendShiftResource::getConsumption(
    const State & initialState) const {
  if (pShift_->isAnyWork())
    return std::min(ub_, initialState.consWeekendWorked_);
  if (pShift_->isRest())
    return std::min(ub_, initialState.consWeekendOff_);
  return 0;
}

PExpander SoftConsWeekendShiftResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const PRCArc &pArc) {
  return initExpander<SoftConsWeekendShiftExpander,
                      SoftConsWeekendShiftResource>(
      prevAShift, stretch, pArc, *this, [this](int c) {
    return this->getCost(c);
  });
}

bool SoftConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  vChild->consumption += consBeforeReset;

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && !reset && vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(
        resource_.getUbCost() * (vChild->consumption - resource_.getUb()));
#ifdef DBG
    pLChild->addConsWeekendShiftCost(
        resource_.getUbCost() * (vChild->consumption - resource_.getUb()));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  // if not resetting counter, set worst case and return
  if (!reset) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost =
        resource_.getWorstUbCost(vChild->consumption, nWeekendsAfter);
    if (cyclic && vChild->cyclicConsumption > 0)
      vChild->worstUbCost += resource_.getUbCost() * vChild->cyclicConsumption;
    return true;
  }

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // pay for violations of soft bounds when resetting a resource
  // Should check if the resource has been consumed on last weekend
  if (vChild->consumption  > 0) {
    pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
#ifdef DBG
    pLChild->addConsWeekendShiftCost(resource_.getLbCost(vChild->consumption));
#endif
  }

  // initialize cyclicConsumption on the first reset
  if (vChild->cyclicConsumption == -1)
    vChild->cyclicConsumption = vChild->consumption;

  // set new consumption to what is consumed after resetting
  vChild->consumption = consAfterReset;
  // compute worst case costs only if resource is consumed after reset
  // otherwise do nothing (use default value)
  if (consAfterReset) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost =
        resource_.getWorstUbCost(vChild->consumption, nWeekendsAfter);
    if (cyclic && vChild->cyclicConsumption > 0)
      vChild->worstUbCost += resource_.getUbCost() * vChild->cyclicConsumption;
  } else {
    vChild->worstLbCost = .0;
    vChild->worstUbCost = .0;
  }

  return true;
}

bool SoftConsWeekendShiftExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  Tools::throwError("Not implemented.");
  return false;
}

int HardConsWeekendShiftResource::getConsumption(
    const State & initialState) const {
  if (pShift_->isAnyWork())
    return initialState.consWeekendWorked_;
  if (pShift_->isRest())
    return initialState.consWeekendOff_;
  return 0;
}

PExpander HardConsWeekendShiftResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const PRCArc &pArc) {
  return initExpander<HardConsWeekendShiftExpander,
                      HardConsWeekendShiftResource>(
      prevAShift, stretch, pArc, *this, nullptr);
}

bool HardConsWeekendShiftResource::dominates(
    const PRCLabel &pL1, const PRCLabel &pL2, double *cost) {
  if (!cyclic_) return BoundedResource::dominates(pL1, pL2, cost);
  const auto &v1 = pL1->getResourceValues(id_),
      &v2 = pL2->getResourceValues(id_);
  if (v1.consumption > v2.consumption) return false;
  int c1 = v1.consumption + v1.cyclicConsumption,
      c2 = v2.consumption + v2.cyclicConsumption;
  if (v1.consumption == v2.consumption) return c1 <= c2;
  if (v1.consumption < lb_) return false;
  return c1 <= c2;
}


bool HardConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  vChild->consumption += consBeforeReset;

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && !reset && vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // detect infeasibility due to upper bound
  if (vChild->consumption > resource_.getUb())
    return false;

  // if resetting counter
  if (reset) {
    // initialize cyclicConsumption on the first reset
    if (vChild->cyclicConsumption == -1)
      vChild->cyclicConsumption = vChild->consumption;
    // expansion is infeasible if consumption lower than bound at reset
    if (vChild->consumption  > 0 && vChild->consumption < resource_.getLb())
      return false;
    // set new consumption to what is consumed after resetting
    vChild->consumption = consAfterReset;
  }

  // if cyclic, add the initial consumption when close to the end
  if (cyclic && nWeekendsAfter == 0 && vChild->cyclicConsumption > 0) {
    vChild->consumption += vChild->cyclicConsumption;
    // detect infeasibility due to upper bound
    if (vChild->consumption > resource_.getUb())
      return false;
  }

  return true;
}

bool HardConsWeekendShiftExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  Tools::throwError("Not implemented.");
  return false;
}
