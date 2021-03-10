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


int SoftConsWeekendShiftResource::getConsumption(
    const State & initialState) const {
  if (pShift_->isAnyWork())
    return initialState.consWeekendWorked_;
  if (pShift_->isRest())
    return initialState.consWeekendOff_;
  return 0;
}

PExpander SoftConsWeekendShiftResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const PRCArc &pArc) {
  // check if expander is active
  if (!Tools::nWeekendsInInterval(stretch.firstDay(), stretch.lastDay()))
    return nullptr;

  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch

  // number of days before the start of the stretch (beware that indices of
  // days start at 0)
  int nDaysBefore = stretch.firstDay();
  // Number of days left since the day of the target node of the arc
  int nDaysLeft = totalNbDays_ - stretch.lastDay() - 1;
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  int day = stretch.firstDay();
  bool weekendCounted =
      Tools::isWeekend(day-1) && pShift_->includes(prevAShift);
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      if (pShift_->includes(**itShift)) {
        // on a weekend and same shift not counted -> increment consAfterReset
        if (!weekendCounted) {
          consBeforeReset++;
          weekendCounted = true;
        }
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

  // if already exceeding UB -> pay immediately the current excess
  if (consBeforeReset > ub_) {
    cost += getUbCost(consBeforeReset);
    consBeforeReset =  ub_;
  }

  // then compute the penalty due to consecutive shifts inside the stretch
  // -> after reset (should be a monday)
  weekendCounted = false;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      if (pShift_->includes(**itShift)) {
        // on a weekend and same shift not counted -> increment consAfterReset
        if (!weekendCounted) {
          consAfterReset++;
          weekendCounted = true;
        }
      }
      if (Tools::isLastWeekendDay(day)) {
        // On last weekend day check if some bounds should be checked
        if (!weekendCounted) {
          // different shifts during whole weekend
          if (consAfterReset == 0) continue;
          cost += getCost(consAfterReset);
          consAfterReset = 0;
        } else {
          // reset flag
          weekendCounted = false;
        }
      }
    }
  }

  // add cost o base cost of the arc
  pArc->addBaseCost(cost);

  // if nothing before and after reset and no reset (price previous consumption)
  // -> initialize nothing
  if (consBeforeReset == 0 && consAfterReset == 0 && !reset)
    return nullptr;

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<SoftConsWeekendShiftExpander>(
      *this, reset, consBeforeReset, consAfterReset,
      pArc->target->type == SINK_NODE, nDaysBefore, nDaysLeft);
}

bool SoftConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  vChild->consumption += consBeforeReset;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addCost(
        resource_.getUbCost() * (vChild->consumption - resource_.getUb()));
#ifdef DBG
    pLChild->addConsShiftCost(
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
        resource_.getWorstUbCost(vChild->consumption, nDaysLeft);
    return true;
  }

  // pay for violations of soft bounds when resetting a resource
  // Should check if the resource has been consumed on last weekend
  if (vChild->consumption  > 0)
    pLChild->addCost(resource_.getLbCost(vChild->consumption));
#ifdef DBG
  pLChild->addConsShiftCost(resource_.getLbCost(vChild->consumption));
#endif

  // set new consumption to what is consumed after resetting
  vChild->consumption = consAfterReset;
  // compute worst case costs only if resource is consumed after reset
  // otherwise do nothing (use default value)
  if (consAfterReset) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost =
        resource_.getWorstUbCost(vChild->consumption, nDaysLeft);
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
  // we need to count the number of times the considered shift appears in
  // the arc's stretch on a weekend
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  int day = stretch.firstDay();
  bool weekendCounted =
      Tools::isWeekend(day-1) && pShift_->includes(prevAShift);
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      if (pShift_->includes(**itShift)) {
        // on a weekend and same shift not counted -> increment consAfterReset
        if (!weekendCounted) {
          consBeforeReset++;
          weekendCounted = true;
        }
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

  if (consBeforeReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // then compute the penalty due to consecutive shifts inside the stretch
  // (should be a monday)
  weekendCounted = false;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isWeekend(day)) {
      if (pShift_->includes(**itShift)) {
        // on a weekend and same shift not counted -> increment consAfterReset
        if (!weekendCounted) {
          consAfterReset++;
          weekendCounted = true;
        }
      }
      if (Tools::isLastWeekendDay(day)) {
        // On last weekend day check if some bounds should be checked
        if (!weekendCounted) {
          // different shifts during whole weekend
          if (consAfterReset > 0 &&
              (consAfterReset < lb_ || consAfterReset > ub_)) {
            // different shifts
            Tools::throwError("RCSPP arc is infeasible");
          }
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

  if (consAfterReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<HardConsWeekendShiftExpander>(
      *this, false, reset, consBeforeReset, consAfterReset, cost);
}

bool HardConsWeekendShiftExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  Tools::throwError("Not implemented.");
  return false;
}
