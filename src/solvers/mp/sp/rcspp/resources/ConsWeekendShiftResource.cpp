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
  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch

  // number of days before the start of the stretch (beware that indices of
  // days start at 0)
  int nDaysBefore = stretch.firstDay();
  // Number of days left since the day of the target node of the arc
  int nDaysLeft = totalNbDays_ - (stretch.firstDay() + stretch.nDays());
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  int day = stretch.firstDay();
  bool prevDayCounted = Tools::isSunday(day) && pShift_->includes(prevAShift);
  bool saturdayCounted = prevDayCounted;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isSaturday(day)) {
      // same shift ->increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consBeforeReset += 1;
        saturdayCounted = true;
      }
    } else if (Tools::isSunday(day) && !saturdayCounted) {
      // on a sunday and saturday not already counted
      // same shift -> increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consBeforeReset += 1;
      } else if (prevDayCounted || consBeforeReset >= 1) {
        reset = true;
        break;
      }
    } else {
      // reset flag
      saturdayCounted = false;
    }
  }

  // then compute the penalty due to consecutive shifts inside the stretch
  // -> after reset (should be a monday)
  saturdayCounted = false;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isSaturday(day)) {
      // same shift
      if (pShift_->includes(**itShift)) {
        consAfterReset++;
        saturdayCounted = true;
      }
    } else if (Tools::isSunday(day) && !saturdayCounted) {
      // on a sunday and saturday not already counted
      // same shift -> increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consAfterReset += 1;
      } else {
        // different shifts
        if (consAfterReset == 0) continue;
        if (consAfterReset < lb_)
          cost += lbCost_ * (lb_ - consAfterReset);
        else if (consAfterReset > ub_)
          cost += ubCost_ * (consAfterReset - ub_);
        consAfterReset = 0;
      }
    } else {
      // reset flag
      saturdayCounted = false;
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
  return std::make_shared<SoftConsShiftExpander>(
      *this, reset, consBeforeReset, consAfterReset,
      pArc->target->type == SINK_NODE, nDaysBefore, nDaysLeft);
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
  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  int day = stretch.firstDay();
  bool prevDayCounted = Tools::isSunday(day) && pShift_->includes(prevAShift);
  bool saturdayCounted = prevDayCounted;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isSaturday(day)) {
      // same shift ->increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consBeforeReset += 1;
        saturdayCounted = true;
      }
    } else if (Tools::isSunday(day) && !saturdayCounted) {
      // on a sunday and saturday not already counted
      // same shift -> increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consBeforeReset += 1;
      } else if (prevDayCounted || consBeforeReset >= 1) {
        reset = true;
        break;
      }
    } else {
      // reset flag
      saturdayCounted = false;
    }
  }

  // if inactive, initialize nothing
  if (consBeforeReset == 0 && !reset)
    return nullptr;

  if (consBeforeReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // then compute the penalty due to consecutive shifts inside the stretch
  // (should be a monday)
  saturdayCounted = false;
  for (; itShift != stretch.pShifts().end(); itShift++, day++) {
    if (Tools::isSaturday(day)) {
      // same shift
      if (pShift_->includes(**itShift)) {
        consAfterReset++;
        saturdayCounted = true;
      }
    } else if (Tools::isSunday(day) && !saturdayCounted) {
      // on a sunday and saturday not already counted
      // same shift -> increment consBeforeReset
      if (pShift_->includes(**itShift)) {
        consAfterReset += 1;
      } else {
        // different shifts
        if (consAfterReset == 0) continue;
        if (consAfterReset < lb_ || consAfterReset > ub_)
          Tools::throwError("RCSPP arc is infeasible");
      }
    } else {
      // reset flag
      saturdayCounted = false;
    }
  }

  if (consAfterReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<HardConsShiftExpander>(
      *this, false, reset, consBeforeReset, consAfterReset, cost);
}
