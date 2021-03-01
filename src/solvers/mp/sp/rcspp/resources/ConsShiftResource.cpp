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

#include "ConsShiftResource.h"

int SoftConsShiftResource::getConsumption(const State & initialState) const {
  int conso = 0;
  if (pShift_->isAnyWork()) {
    conso = initialState.consDaysWorked_;
  } else if (pShift_->isRest()) {
    conso = initialState.consDaysOff_;
  } else if (pShift_->includes(*initialState.pShift_)) {
    conso = initialState.consShifts_;
  }
  return conso;
}

PExpander SoftConsShiftResource::init(const Shift &prevShift,
                                      const Stretch &stretch,
                                      const RCArc &arc) {
  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch

  // Number of days left since the day of the target node of the arc
  int nDaysLeft = totalNbDays_ - (stretch.firstDay() + stretch.nDays());
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift)
  auto itShift = stretch.pShifts().begin();
  for (; itShift != stretch.pShifts().end(); itShift++) {
    // same shift ->increment consBeforeReset
    if (pShift_->includes(**itShift)) {
      consBeforeReset += 1;
    } else if (pShift_->includes(prevShift) || consBeforeReset >= 1) {
      reset = true;
      break;
    }
  }

  // if inactive, initialize nothing
  if (consBeforeReset == 0 && !reset)
    return nullptr;

  // then compute the penalty due to consecutive shifts inside the stretch
  // -> after reset
  for (; itShift != stretch.pShifts().end(); itShift++) {
    // same shift
    if (pShift_->includes(**itShift)) {
      consAfterReset++;
    } else {
      // different shifts
      if (consAfterReset == 0) continue;
      if (consAfterReset < lb_)
        cost += lbCost_ * (lb_ - consAfterReset);
      else if (consAfterReset > ub_)
        cost += ubCost_ * (consAfterReset - ub_);
      consAfterReset = 0;
    }
  }

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  if (arc.target->type == SINK_NODE)
    return std::make_shared<SoftConsShiftExpander>(
        *this, reset, consBeforeReset, consAfterReset, cost, true, nDaysLeft);
  else
    return std::make_shared<SoftConsShiftExpander>(
        *this, reset, consBeforeReset, consAfterReset, cost, nDaysLeft);
}

bool SoftConsShiftExpander::expand(const ResourceValues &vParent,
                                   const PRCLabel &pLChild,
                                   ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  vChild->consumption += consBeforeReset;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addConsShiftCost(resource_.getUbCost(vChild->consumption));
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  // if not resetting counter, set worst case and return
  if (!reset) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption);
    return true;
  }

  // pay for violations of soft bounds when resetting a resource
  pLChild->addConsShiftCost(resource_.getLbCost(vChild->consumption));

  // set new consumption to what is consumed after resetting
  vChild->consumption = consAfterReset;
  // compute worst case costs only if resource is consumed after reset
  // otherwise do nothing (use default value)
  if (consAfterReset) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    vChild->worstUbCost = resource_.getWorstUbCost(vChild->consumption);
  } else {
    vChild->worstLbCost = .0;
    vChild->worstUbCost = .0;
  }
  // ensure that consumption does not exceed ub_
  if (vChild->consumption > resource_.getUb())
    vChild->consumption = resource_.getUb();

  return true;
}

int HardConsShiftResource::getConsumption(const State &  initialState) const {
  if (pShift_->isAnyWork())
    return initialState.consDaysWorked_;
  if (pShift_->isRest())
    return initialState.consDaysOff_;
  if (pShift_->includes(*initialState.pShift_))
    return initialState.consShifts_;
  return 0;
}

PExpander HardConsShiftResource::init(const Shift &prevShift,
                                      const Stretch &stretch,
                                      const RCArc &arc) {
  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  auto itShift = stretch.pShifts().begin();
  for (; itShift != stretch.pShifts().end(); itShift++) {
    if (pShift_->includes(**itShift)) {
      consBeforeReset += 1;
    } else if (pShift_->includes(prevShift) || consBeforeReset >= 1) {
      reset = true;
      break;
    }
  }

  // if inactive, initialize nothing
  if (consBeforeReset == 0 && !reset)
    return nullptr;

  if (consBeforeReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // then compute the penalty due to consecutive shifts inside the stretch
  for (; itShift != stretch.pShifts().end(); itShift++) {
    if (pShift_->includes(**itShift)) {
      consAfterReset++;
    } else {
      if (!consAfterReset) continue;
      if (consAfterReset < lb_ || consAfterReset > ub_)
        Tools::throwError("RCSPP arc is infeasible");
      consAfterReset = 0;
    }
  }
  if (consAfterReset > ub_)
    Tools::throwError("RCSPP arc is infeasible");

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<HardConsShiftExpander>(
      *this, reset, consBeforeReset, consAfterReset, cost);
}

bool HardConsShiftExpander::expand(const ResourceValues &vParent,
                                   const PRCLabel &pLChild,
                                   ResourceValues *vChild) {
  // consumption before resetting resource if any reset
  vChild->consumption += consBeforeReset;

  // detect infeasibility due to upper bound
  if (vChild->consumption > resource_.getUb())
    return false;

  // if resetting counter
  if (reset) {
    // expansion is infeasible if consumption lower than bound at reset
    if (vChild->consumption < resource_.getLb()) return false;
    // set new consumption to what is consumed after resetting
    vChild->consumption = consAfterReset;
  }
  return true;
}

