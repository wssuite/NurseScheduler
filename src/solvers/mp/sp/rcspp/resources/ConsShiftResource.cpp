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

#include <set>

template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r,
                           const std::function<double(int)> &getCost,
                           const int indResource) {
  // we need to count the number of times the considered shift appears at the
  // beginning of the arc's stretch
  bool reset = false;
  int consBeforeReset = 0;
  int consAfterReset = 0;
  double cost = 0;

  // Look at consecutive shifts before reset (a different shift); the initial
  // consumption is taken into account when creating the initial label, not
  // at the resource level.
  auto itShift = stretch.pShifts().begin();
  for (; itShift != stretch.pShifts().end(); itShift++) {
    // same shift ->increment consBeforeReset
    if (r.pShift()->includes(**itShift)) {
      consBeforeReset += 1;
      continue;
    }
    // reset ? either was working on this shift just before or
    // has already worked some same shifts here
    reset = r.pShift()->includes(prevAShift) || consBeforeReset > 0;
    if (reset)
      break;
  }

  if (consBeforeReset > r.getUb()) {
    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
    cost += getCost(consBeforeReset);
    consBeforeReset = r.getUb();
  }

  // then compute the penalty due to consecutive shifts inside the stretch
  // -> after reset
  for (; itShift != stretch.pShifts().end(); itShift++) {
    // same shift
    if (r.pShift()->includes(**itShift)) {
      consAfterReset++;
    } else {
      // different shifts
      if (consAfterReset == 0) continue;
      if (consAfterReset < r.getLb() || consAfterReset > r.getUb()) {
        if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
        cost += getCost(consAfterReset);
      }
      consAfterReset = 0;
    }
  }

  // pay for excess of consumption
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  // pay all cost if it's the end of the sequence
  bool arcToSink = pArc->target->type == SINK_NODE;
  if (consAfterReset > r.getUb()) {
    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
    cost += getCost(consAfterReset);
    consAfterReset = r.getUb();
  }

  if (arcToSink && r.lastDayEndsSequence() &&
      consAfterReset > 0 && consAfterReset < r.getLb()) {
    if (r.isHard()) Tools::throwError("RCSPP arc is infeasible");
    // pay for violations of soft bounds when at the end of horizon in INRC
    cost += getCost(consAfterReset);
  }

  // add cost to base cost of the arc
  pArc->addBaseCost(cost);

  // if nothing before and after reset and no reset (price previous consumption)
  // -> initialize nothing
  if (consBeforeReset == 0 && consAfterReset == 0 && !reset)
    return nullptr;

  // The arc starts the consumption of the resource if the previous shift is
  // not included in the resource abstract shift, and if the resource shift
  // is consumed before reset
  bool start = (!r.pShift()->includes(prevAShift)) && (consBeforeReset > 0);
  // number of days before the start of the stretch (beware that indices of
  // days start at 0)
  std::pair<int, int> firstLastDays = r.getFirstLastDays(stretch);
  int nDaysBefore = firstLastDays.first + r.initialConsumption();
  // Number of days left since the day of the target node of the arc
  int nDaysLeft = r.totalNbDays() + r.firstDayId() - firstLastDays.second - 1;

  // if the stretch ends with the considered, we get a non-zero number of
  // consecutive shifts after replenishment
  return std::make_shared<E>(indResource,
      r, start, reset, consBeforeReset, consAfterReset,
      arcToSink, nDaysBefore, nDaysLeft);
}

int SoftConsShiftResource::getConsumption(const State & initialState) const {
  if (pAShift_->isAnyWork())
    return std::min(ub_, initialState.consDaysWorked_);
  if (pAShift_->isRest())
    return std::min(ub_, initialState.consDaysOff_);
  if (pAShift_->includes(*initialState.pShift_))
    return std::min(ub_, initialState.consShifts_);
  return 0;
}

PExpander SoftConsShiftResource::init(const AbstractShift &prevAShift,
                                      const Stretch &stretch,
                                      const shared_ptr<RCArc> &pArc,
                                      int indResource) {
  return initExpander<SoftConsShiftExpander,
                      SoftConsShiftResource>(
      prevAShift, stretch, pArc, *this, [this](int c) {
        return this->getCost(c);
      }, indResource);
}

void SoftConsShiftResource::enumerate(const PRCGraph &pRCGraph,
                                      bool forceEnum) {
  // in general, we do not enumerate these resources, the enumeration is done
  // only if we force it
  if (!forceEnum) return;

  // In case, the abstract shift is not a real shift or a rest
  if (pAShift_->isRest() || pAShift_->isAnyWork()) return;

  // check if ub is active
  bool activeUB = ub_ < totalNbDays_ * maxConsumptionPerDay();

  // enumerate for all possibilities of shifts succession of this type
  std::set<PShift> includedShifts;
  for (const PShift &pS : pRCGraph->pShifts())
    if (pAShift_->includes(*pS))
      includedShifts.insert(pS);

  for (const PRCNode &pOrigin : pRCGraph->pNodes()) {
    // do not treat a sink node
    if (pOrigin->type == SINK_NODE) continue;

    // iterate through the outgoing arcs
    // take a copy as it will be modified when adding new arcs
    vector<PRCArc> outArcs = pOrigin->outArcs;
    for (const PRCArc &pArc : outArcs) {
      // check if the target shift is the same arc, if not continue
      if (!pAShift_->includes(*pArc->target->pAShift)) {
        // add initial cost if any
        if (pAShift_->includes(*pArc->origin->pAShift) && pOrigin->dayId == -1)
          pArc->addBaseCost(getLbCost(initialConsumption_));
        continue;
      }

      // All following arcs added will have pOrigin node as their origin
      // The successor shift is either of the same type as the initial shift
      // if starting from day -1 or of a different type.
      // 1- Compute the initial consumption of this shift type:
      // either 0 or initialConsumption if from source.
      int initialCons(0);
      double initialBaseCost(0);
      if (pAShift_->includes(*pArc->origin->pAShift)) {
        if (pOrigin->dayId == -1) {
          // if starting from day -1, we can work on the same shift type without
          // penalty, but we need to set the initial consumption
          initialCons = initialConsumption_;
        } else {
          // otherwise, we only update the cost of the arc from this node to
          // that with same shift type on next day: it must be taken only if
          // the upper bound is already met, so we add the ubCost to the arc
          pArc->addBaseCost(activeUB ? ubCost_ : 0);
          continue;
        }
      }
      // store initial base cost:
      // it will include the cost associated to the initial state now or
      // it will be added later on when enumerating the other resources
      initialBaseCost = pArc->baseCost;

      // 2 - update the cost of the current arc
      // add lb cost if not the end or if it is always counted
      if (pArc->target->type != SINK_NODE || lastDayEndsSequence_)
        pArc->addBaseCost(getLbCost(pArc->stretch.nDays() + initialCons));

      // 3 - compute all possible stretch length of size > 1 and <= UB or
      // UB - initialCons if starting from day -1 and on the same initial shift
      vector<Stretch> stretches = {pArc->stretch};
      int ub = (activeUB ? ub_ : lb_) - initialCons;
      for (int d = 2; d <= ub && d < totalNbDays_ - pOrigin->dayId; d++) {
        // new stretches
        vector<Stretch> newStretches;

        // add a shift of same type
        for (const PShift &pShift : includedShifts) {
          // Target node of the arc that will be added
          PRCNode pTarget = nullptr;
          for (const PRCNode &pN : pRCGraph->pNodesPerDayId(
              pOrigin->dayId + d)) {
            if (pShift->equals(*pN->pAShift)) {
              pTarget = pN;
              break;
            }
          }
          if (pTarget == nullptr)
            break;

          // Creation of the stretch of the arc that will be added
          for (Stretch st : stretches) {
            if (!st.pShifts().back()->canPrecede(*pShift))
              continue;
            st.pushBack(pShift);  // add one more same shift to the stretch
            newStretches.push_back(st);

            // The new arc is added to the rcGraph (if it is already there, the
            // method only retrieves the arc from the arc list)
            PRCArc pNewArc =
                pRCGraph->addSingleArc(pOrigin, pTarget, st, initialBaseCost);

            // Cost due to penalty of soft lower bound of the corresponding
            // shift (added only if the target node is not a sink node or
            // if always counted).
            // also, add initialCons in case of starting from day -1,
            // generally it's 0.
            if (pTarget->type != SINK_NODE  || lastDayEndsSequence_)
              pNewArc->addBaseCost(getLbCost(st.nDays() + initialCons));
          }
        }

        // set the new stretches
        stretches = newStretches;
      }
    }
  }
  isPreprocessed_ = true;
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

  // if cyclic, add the initial consumption when close to the end
  if (resource_.isCyclic() && nDaysLeft == 0 && !reset &&
      vChild->cyclicConsumption > 0)
    vChild->consumption += vChild->cyclicConsumption;

  // pay for excess of consumption due to this expansion
  // pay attention that consumptions in the label of the source should not
  // exceed the upper bounds
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
#ifdef DBG
    pLChild->addConsShiftCost(resource_.getUbCost(vChild->consumption));
#endif
    // beware: we never need to store a consumption larger than the upper bound
    vChild->consumption = resource_.getUb();
  }

  if (reset) {
    // initialize cyclicConsumption on the first reset
    if (vChild->cyclicConsumption == -1)
      vChild->cyclicConsumption = vChild->consumption;

    // pay for violations of soft bounds when resetting a resource if the
    // consumption is > 0. It could be null when pricing rotations
    if (vChild->consumption > 0) {
      pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
#ifdef DBG
      pLChild->addConsShiftCost(
          resource_.getLbCost(vChild->consumption));
#endif
    }

    // set new consumption to what is consumed after resetting
    vChild->consumption = consAfterReset;
  } else if (arcToSink && resource_.lastDayEndsSequence()
      && (vChild->consumption >= 1)
      && (vChild->consumption < resource_.getLb())) {
    // pay for violations of soft bounds when at the end of horizon in INRC
    // when a reset has happened, this lb cost has been already paid in the init
    pLChild->addBaseCost(
        resource_.getLbCost(vChild->consumption));
#ifdef DBG
    pLChild->addConsShiftCost(
          resource_.getLbCost(vChild->consumption));
#endif
  }

  // compute worst case costs only if resource is consumed
  // otherwise do nothing (use default value)
  if (vChild->consumption > 0) {
    // Setting 'worst case cost'
    vChild->worstLbCost = resource_.getWorstLbCost(vChild->consumption);
    if (resource_.isCyclic() && vChild->cyclicConsumption > 0)
      vChild->worstUbCost = resource_.getWorstUbCost(
          vChild->consumption, nDaysLeft + vChild->cyclicConsumption);
    else
      vChild->worstUbCost =
        resource_.getWorstUbCost(vChild->consumption, nDaysLeft);
  } else {
    vChild->worstLbCost = .0;
    vChild->worstUbCost = .0;
  }

  return true;
}

bool SoftConsShiftExpander::expandBack(const PRCLabel &pLChild,
                                       ResourceValues *vChild) {
  if (!reset) {
    // consumption is in consBeforeReset if there is no reset
    vChild->consumption += consBeforeReset;
    // reset worst-case costs and return if no consumption
    if (vChild->consumption == 0) {
      vChild->worstLbCost = 0;
      vChild->worstUbCost = 0;
      return true;
    }
  } else {
    // when going backwards, we first add the consumption after the reset
    vChild->consumption += consAfterReset;
    if (vChild->consumption > 0) {
      // pay for violations of soft bounds when resetting a resource
      // pay lower bound consumption only if the last sequence always counts or
      // if the consecutive shifts did not end at a sink
      if (resource_.lastDayEndsSequence() ||
          vChild->consumption <= nDaysLeft)
        pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
      // pay for excess of consumption at reset
      if (vChild->consumption > resource_.getUb())
        pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
    }
    // set new consumption to what is consumed before resetting
    vChild->consumption = consBeforeReset;
  }
  // pay for excess of consumption due to this expansion
  if (vChild->consumption > resource_.getUb()) {
    pLChild->addBaseCost(resource_.getUbCost(vChild->consumption));
    // beware: we never need to store a consumption larger than the upper
    // bound
    vChild->consumption = resource_.getUb();
  }

  // if it is a start arc, pay for violation of lower bound and reset
  // worst-case costs
  if (start) {
    // pay lower bound consumption only if always count or
    // if the consecutive shifts did not end at a sink
    if (resource_.lastDayEndsSequence() || vChild->consumption <= nDaysLeft)
      pLChild->addBaseCost(resource_.getLbCost(vChild->consumption));
    vChild->consumption = 0;
    vChild->worstLbCost = 0;
    vChild->worstUbCost = 0;
    return true;
  }

  // DOC: when expanding back, consumption of resource related to the node we
  // are expanding to is not counted. This means that we know that at least
  // one more unit of resource will be consumed when going out of the node.
  // We can use it to improve the computation of worst-case LB cost
  if (!resource_.lastDayEndsSequence() && vChild->consumption == nDaysLeft + 1)
    vChild->worstLbCost = 0;
  else
    vChild->worstLbCost =
        resource_.getWorstLbCost(vChild->consumption+1);
  vChild->worstUbCost =
      resource_.getWorstUbCost(vChild->consumption, nDaysBefore);
  return true;
}

int HardConsShiftResource::getConsumption(const State &  initialState) const {
  if (pAShift_->isAnyWork())
    return initialState.consDaysWorked_;
  if (pAShift_->isRest())
    return initialState.consDaysOff_;
  if (pAShift_->includes(*initialState.pShift_))
    return initialState.consShifts_;
  return 0;
}

PExpander HardConsShiftResource::init(const AbstractShift &prevAShift,
                                      const Stretch &stretch,
                                      const shared_ptr<RCArc> &pArc,
                                      int indResource) {
  return initExpander<HardConsShiftExpander, HardConsShiftResource>(
      prevAShift, stretch, pArc, *this, nullptr, indResource);
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

  // if cyclic, add the initial consumption when close to the end
  if (resource_.isCyclic() && nDaysLeft == 0 &&
      vChild->cyclicConsumption > 0) {
    vChild->consumption += vChild->cyclicConsumption;
    // take into account the reset happening at the end of the cycle
    reset = true;
  }

  // detect infeasibility due to upper bound
  if (vChild->consumption > resource_.getUb())
    return false;

  // if resetting counter
  if (reset) {
#ifdef DBG
    // It could never be equal to 0, as there is at least 1 consumption for the
    // shift of the current node
    if (vChild->consumption  == 0)
      Tools::throwError("The resource value is 0 when it should be > 0.");
#endif
    // initialize cyclicConsumption on the first reset
    if (vChild->cyclicConsumption == -1)
      vChild->cyclicConsumption = vChild->consumption;
    // expansion is infeasible if consumption lower than bound at reset
    if (vChild->consumption < resource_.getLb())
      return false;
    // set new consumption to what is consumed after resetting
    vChild->consumption = consAfterReset;
  } else if (arcToSink && resource_.lastDayEndsSequence()
        && (vChild->consumption >= 1)
        && (vChild->consumption < resource_.getLb())) {
    return false;
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
      // check lower bound only if the consecutive shifts did not end at a sink
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

