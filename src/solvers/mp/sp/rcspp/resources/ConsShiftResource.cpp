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
  if (pShift_->isAnyWork())
    return initialState.consDaysWorked_;
  if (pShift_->isRest())
    return initialState.consDaysOff_;
  if (pShift_->includes(*initialState.pShift_))
    return initialState.consShifts_;
  return 0;
}

PExpander SoftConsShiftResource::init(const AbstractShift &prevAShift,
                                      const Stretch &stretch,
                                      const PRCArc &pArc) {
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
      continue;
    }
    // reset ? either was working on this shift just before or
    // has already worked some same shifts here
    reset = pShift_->includes(prevAShift) || consBeforeReset > 0;
    break;
  }

  // The arc starts the consumption of the resource if the previous shift is
  // not included in the resource abstract shift, and if the resource shift
  // is consumed before reset
  bool start = (!pShift_->includes(prevAShift)) && (consBeforeReset > 0);

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

  // pay for excess of consumption
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (consAfterReset > ub_) {
    cost += ubCost_ * (consAfterReset - ub_);
    consAfterReset = ub_;
  }

  // add cost to base cost of the arc
  pArc->addBaseCost(cost);

  // if nothing before and after reset and no reset (price previous consumption)
  // -> initialize nothing
  if (consBeforeReset == 0 && consAfterReset == 0 && !reset)
    return nullptr;

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<SoftConsShiftExpander>(
      *this, start, reset, consBeforeReset, consAfterReset,
      pArc->target->type == SINK_NODE, nDaysLeft);
}
bool SoftConsShiftResource::merge(const ResourceValues &vForward,
                                  const ResourceValues &vBack,
                                  ResourceValues *vMerged,
                                  const PRCLabel &pLMerged) {
  if (vForward.consumption + vBack.consumption == 0) return true;
  return SoftBoundedResource::merge(vForward, vBack, vMerged, pLMerged);
}

bool SoftConsShiftExpander::expand(const PRCLabel &pLChild,
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

// TODO(JO): forgot to nullify violations due to LB on consecutive shifts
//  when they end at the sink
// TODO(JO): beware of which resource will be active for domination
bool SoftConsShiftExpander::expandBack(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  if (!reset) {
    // consumption is in consBeforeReset if there is no reset
    vChild->consumption += consBeforeReset;
  } else {
    // when going backwards, we first add the consumption after the reset
    vChild->consumption += consAfterReset;
    if (vChild->consumption > 0) {
      // pay for violations of soft bounds when resetting a resource
      // pay lower bound consumption only if the consecutive shifts did not end
      // at a sink
      if (vChild->consumption < nDaysLeft + 1)
        pLChild->addCost(resource_.getLbCost(vChild->consumption));
      // pay for excess of consumption at reset
      if (vChild->consumption > resource_.getUb()) {
        pLChild->addCost(
            resource_.getUbCost() * (vChild->consumption - resource_.getUb()));
      }
    }
    // set new consumption to what is consumed before resetting
    vChild->consumption = consBeforeReset;
  }
  // pay for excess of consumption due to this expansion
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addCost(
        resource_.getUbCost() * (vChild->consumption - resource_.getUb()));
    // beware: we never need to store a consumption larger than the upper
    // bound
    vChild->consumption = resource_.getUb();
  }

  // if it is a start arc, pay for violation of lower bound and reset
  // worst-case costs
  if (start) {
    // pay lower bound consumption only if the consecutive shifts did not end
    // at a sink
    if (vChild->consumption < nDaysLeft + 1)
      pLChild->addCost(resource_.getLbCost(vChild->consumption));
    vChild->consumption = 0;
  }

  // DOC: when expanding back, consumption of resource related to the node we
  // are expanding to is not counted. This means that we know that at least
  // one more unit of resource will be consumed when going out of the node.
  // We can use it to improve the computation of worst-case costs
  vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption + 1);
  vChild->worstUbCost =
      resource_.getWorstUbCost(vChild->consumption + 1,
                               resource_.getTotalNbDays() - nDaysLeft - 1);
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

PExpander HardConsShiftResource::init(const AbstractShift &prevAShift,
                                      const Stretch &stretch,
                                      const PRCArc &pArc) {
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
    } else if (pShift_->includes(prevAShift) || consBeforeReset >= 1) {
      reset = true;
      break;
    }
  }

  // if inactive, initialize nothing
  if (consBeforeReset == 0 && !reset)
    return nullptr;

  // The arc starts the consumption of the resource if the previous shift is
  // not included in the resource abstract shift, and if the resource shift
  // is consumed before reset
  bool start = (!pShift_->includes(prevAShift)) && (consBeforeReset > 0);

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
      *this, start, reset, consBeforeReset, consAfterReset, cost);
}
bool HardConsShiftResource::merge(const ResourceValues &vForward,
                                  const ResourceValues &vBack,
                                  ResourceValues *vMerged,
                                  const PRCLabel &pLMerged) {
  if (vForward.consumption + vBack.consumption == 0) return true;
  return HardBoundedResource::merge(vForward, vBack, vMerged, pLMerged);
}

bool HardConsShiftExpander::expand(const PRCLabel &pLChild,
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
bool HardConsShiftExpander::expandBack(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  // if resetting counter
  if (reset) {
    // first add consumption after reset resource
    vChild->consumption += consAfterReset;
    if (vChild->consumption > 0) {
      // detect infeasibility due to upper bound
      if (vChild->consumption > resource_.getUb()) return false;
      // expansion is infeasible if consumption lower than bound at reset
      if (vChild->consumption < nDaysLeft + 1)
        if (vChild->consumption < resource_.getLb()) return false;
    }
    // set new consumption to what is consumed before resetting
    vChild->consumption = consBeforeReset;
  } else {
    // the consumption is in consBeforeReset if no reset
    vChild->consumption += consBeforeReset;
    // detect infeasibility due to upper bound
    if (vChild->consumption > resource_.getUb()) return false;
  }

  // check lower bound if the arc starts the consumption of the resource
  if (start) {
    if (vChild->consumption < nDaysLeft + 1)
      if (vChild->consumption < resource_.getLb()) return false;
    vChild->consumption = 0;
  }

  return true;
}

