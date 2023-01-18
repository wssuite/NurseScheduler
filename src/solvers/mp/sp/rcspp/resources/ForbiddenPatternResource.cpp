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

#include "ForbiddenPatternResource.h"

void ForbiddenPatternResource::preprocessPattern() {
  // for each starting subpattern, detect its repetitions in the pattern
  repeatStartPattern_.clear();
  repeatStartPattern_.resize(pattern_.size_);
  for (int i = 0; i < pattern_.size_ - 1; ++i) {
    for (int j = i + 1; j < pattern_.size_; ++j) {
      bool isRepeated = true;
      for (int k = 0; k <= i; ++k)
        if (!pattern_.pAShifts_[j - i + k]->includes(*pattern_.pAShifts_[k]) ||
            !pattern_.pADays_[j - i + k]->includes(*pattern_.pADays_[k])) {
          isRepeated = false;
          break;
        }
      if (isRepeated) repeatStartPattern_[i].push_back(j);
    }
  }

  // for each ending subpattern, detect its repetitions in the pattern
  repeatEndPattern_.clear();
  repeatEndPattern_.resize(pattern_.size_);
  int n = pattern_.size_ - 1;
  for (int i = 0; i < pattern_.size_ - 1; ++i) {
    for (int j = i + 1; j < pattern_.size_; ++j) {
      bool isRepeated = true;
      for (int k = 0; k <= i; ++k)
        if (!pattern_.pAShifts_[n - j + i - k]->includes(
                *pattern_.pAShifts_[n - k]) ||
            !pattern_.pADays_[n - j + i - k]->includes(
                *pattern_.pADays_[n - k])) {
          isRepeated = false;
          break;
        }
      if (isRepeated) repeatEndPattern_[i].push_back(j);
    }
  }
}

// TODO(JO): we will certainly have a bug here if using history and rolling
//  horizon methods, because we do not consider any initial consumption in
//  the resource
template<typename E, typename R>
shared_ptr<E> initExpander(const AbstractShift &prevAShift,
                           const Stretch &stretch,
                           const PRCArc &pArc,
                           const R &r,
                           const int indResource) {
  /*
   * Initialization for forward label-setting
   */
  // check that the pattern does not appear in the stretch, and get the
  // consumption at the end of the stretch while doing so. The end
  // consumption must start with the first day of the pattern to be
  // considered as such
  int endConsumption = 0;
  auto itS1 = stretch.pShifts().begin();
  auto itD1 = stretch.pDays().begin();
  for (; itS1 != stretch.pShifts().end(); ++itS1, ++itD1) {
    // start checking from here only if it matches the first day-shift of the
    // pattern
    if (!r.getPattern().pAShifts_[0]->includes(**itS1)
        || !r.getPattern().pADays_[0]->includes(**itD1))
      continue;
    auto itSpat = r.getPattern().pAShifts_.begin() + 1;
    auto itDpat = r.getPattern().pADays_.begin() + 1;
    auto itS2 = itS1 + 1;
    auto itD2 = itD1 + 1;
    int conso = 1;
    for (; itS2 < stretch.pShifts().end();
           ++itS2, ++itD2, ++itSpat, ++itDpat) {
      if ((*itSpat)->includes(**itS2) && (*itDpat)->includes(**itD2)) {
        conso += 1;
        if (conso == r.getPattern().size_) {
          // case where the whole forbidden pattern appears in a stretch
          if (r.isHard()) {
            // a. if the resource is hard, we must return, with
            //  forward and backward cons all set to the upper bound
            vector<int> forwardConsumption(r.getUb(), 0);
            vector<int> backwardConsumption(r.getPattern().size_, 0);
            for (int c = 0; c < r.getUb(); ++c) {
              forwardConsumption[c] = r.getUb();
              backwardConsumption[c] = r.getUb();
              return std::make_shared<E>(
                  indResource,
                  r,
                  forwardConsumption,
                  backwardConsumption, 0, 0);
            }
          } else {
            // b. if the resource is soft, we must update arc base cost and
            // continue with the computation of endConsumption
            pArc->addBaseCost(r.getCost());
            conso = 0;
            break;
          }
        }
      } else {
        conso = 0;
        break;
      }
    }
    endConsumption = std::max(endConsumption, conso);
  }

  // compute the final consumption for each possible value of consumption in
  // a label in forward label-setting and store it in a vector with the same
  // size as the pattern
  // one difficulty is that a consumption can also correspond to smaller one
  // if a day-shift appears several times in the pattern
  vector<int> forwardConsumption(r.getPattern().size_, 0);
  for (int conso = 0; conso < r.getPattern().size_; ++conso) {
    auto itSpat = r.getPattern().pAShifts_.begin() + conso;
    auto itDpat = r.getPattern().pADays_.begin() + conso;
    auto itD = stretch.pDays().begin();
    auto itS = stretch.pShifts().begin();
    forwardConsumption[conso] = conso;

    for (; itS != stretch.pShifts().end();
           ++itS, ++itD, ++itSpat, ++itDpat) {
      if (!(*itSpat)->includes(**itS) || !(*itDpat)->includes(**itD)) {
        forwardConsumption[conso] = 0;
        break;
      } else {
        forwardConsumption[conso] += 1;
        if (forwardConsumption[conso] == r.getPattern().size_) {
          // the presence of the complete resource pattern at the beginning
          // of the arc stretch has already been counted so the forward
          // consumption can be set to zero for a zero initial consumption
          if (conso == 0)
            forwardConsumption[conso] = 0;
          break;
        }
      }
    }
    // if the end consumption is larger than what can be achieved from the
    // beginning, just take the end consumption
    forwardConsumption[conso] = std::max(forwardConsumption[conso],
                                         endConsumption);
  }

  // for each possible consumption of the resource, check if a smaller
  // consumption is simultaneously achieved, and in such case, take the
  // larger forward consumption
  // example: if the resource pattern is Day-Day-Late, then
  // repeatStartPattern_[0]=[1] and if the arc stretch is Day-Day, the above
  // loop computed forwardConsumption[0]=2 and forwardConsumption[1]=0, but
  // given that subpattern Day can be seen as a starting subpattern as well
  // as the subpattern from index 1 to 1, then we need to set
  // forwardConsumption[1]=2
  for (int i = 0 ; i < r.getPattern().size_; i++) {
    for (int j : r.repeatStartPattern(i)) {
      forwardConsumption[j] =
          std::max(forwardConsumption[i], forwardConsumption[j]);
    }
  }

  /*
   * Initialization for backward label-setting
   */
  // get the backward consumption at the start of the stretch.
  // The start consumption must end with the last day-shift of the pattern to
  // be considered as such
  int startConsumption = 0;
  auto ritS1 = stretch.pShifts().rbegin();
  auto ritD1 = stretch.pDays().rbegin();
  for (; ritS1 != stretch.pShifts().rend(); ++ritS1, ++ritD1) {
    // start checking from here only if it matches the last day-shift of the
    // pattern
    if (!r.getPattern().pAShifts_.back()->includes(**ritS1)
        || !r.getPattern().pADays_.back()->includes(**ritD1))
      continue;
    auto itSpat = r.getPattern().pAShifts_.rbegin() + 1;
    auto itDpat = r.getPattern().pADays_.rbegin() + 1;
    auto ritS2 = ritS1 + 1;
    auto ritD2 = ritD1 + 1;
    int conso = 1;
    for (; ritS2 != stretch.pShifts().rend();
           ++ritS2, ++ritD2, ++itSpat, ++itDpat) {
      if ((*itSpat)->includes(**ritS2) && (*itDpat)->includes(**ritD2)) {
        conso += 1;
        if (conso == r.getPattern().size_) {
          // case where the whole forbidden pattern appears in a stretch
          if (r.isHard()) {
            // if the resource is hard, this should have been detected while
            // computing end consumptions, and the method should have returned
            Tools::throwError("Forbidden patterns should be detected earlier "
                              "in hard constraints");
          } else {
            // if the resource is soft, the cost of the arc has already
            // been updated, so we just need to go on with the computation of
            // the start consumption
            conso = 0;
            break;
          }
        }
      } else {
        conso = 0;
        break;
      }
    }
    startConsumption = std::max(startConsumption, conso);
  }

  // compute the backward consumption: backward consumption is computed from
  // the end of the pattern
  vector<int> backwardConsumption(r.getPattern().size_, 0);
  for (int conso = 0; conso < r.getPattern().size_; ++conso) {
    auto itSpat = r.getPattern().pAShifts_.rbegin() + conso;
    auto itDpat = r.getPattern().pADays_.rbegin() + conso;
    auto itD = stretch.pDays().rbegin();
    auto itS = stretch.pShifts().rbegin();
    backwardConsumption[conso] = conso;

    for (; itS != stretch.pShifts().rend();
           ++itS, ++itD, ++itSpat, ++itDpat) {
      if (!(*itSpat)->includes(**itS) || !(*itDpat)->includes(**itD)) {
        backwardConsumption[conso] = 0;
        break;
      } else {
        backwardConsumption[conso] += 1;
        if (backwardConsumption[conso] == r.getPattern().size_) {
          // complete patterns have already been counted
          if (conso == 0)
            backwardConsumption[conso] = 0;
          break;
        }
      }
    }
    // if the end consumption is larger than what can be achieved from the
    // beginning, just take the end consumption
    backwardConsumption[conso] = std::max(backwardConsumption[conso],
                                          startConsumption);
  }

  // for each possible consumption of the resource, check if a smaller
  // consumption is simultaneously achieved, and in such case, take the
  // larger forward consumption
  for (int i = 0; i < r.getPattern().size_; i++) {
    for (int j : r.repeatEndPattern(i)) {
      backwardConsumption[j] =
          std::max(backwardConsumption[i], backwardConsumption[j]);
    }
  }

  return std::make_shared<E>(
      indResource,
      r,
      forwardConsumption,
      backwardConsumption,
      endConsumption,
      startConsumption);
}

PExpander SoftForbiddenPatternResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const shared_ptr<RCArc> &pArc,
                                             int indResource) {
  if (isPreprocessed_) return nullptr;
  return initExpander<SoftForbiddenPatternExpander,
                      SoftForbiddenPatternResource>(
      prevAShift, stretch, pArc, *this, indResource);
}

PExpander HardForbiddenPatternResource::init(const AbstractShift &prevAShift,
                                             const Stretch &stretch,
                                             const shared_ptr<RCArc> &pArc,
                                             int indResource) {
  if (isPreprocessed_) return nullptr;
  return initExpander<HardForbiddenPatternExpander,
                      HardForbiddenPatternResource>(
      prevAShift, stretch, pArc, *this, indResource);
}

// TODO(JO): If at the beginning of the horizon, we should actually use an
//  initial consumption of the resource I guess
// TODO(JO): Also, we once again get an error if looking at backwards paths
double SoftForbiddenPatternResource::computeBaseCost(
    const Stretch &stretch,
    const PAbstractShift &pPrevAShift) const {
  // count the number of times the pattern appears in the stretch
  double baseCost(0);
  auto itD1 = stretch.pDays().begin();
  auto itS1 = stretch.pShifts().begin();

  // look from the first day of the stretch if the previous shift is the
  // first one of the pattern
  if (pattern_.pAShifts_[0]->includes(*pPrevAShift)
      && pattern_.pADays_[0]->includes(*(*itD1)->previous())) {
    auto itSpat = pattern_.pAShifts_.begin() + 1;
    auto itDpat = pattern_.pADays_.begin() + 1;
    auto itS2 = itS1;
    auto itD2 = itD1;
    int conso = 1;
    for (; itS2 < stretch.pShifts().end();
           ++itS2, ++itD2, ++itSpat, ++itDpat) {
      if ((*itSpat)->includes(**itS2) && (*itDpat)->includes(**itD2)) {
        conso += 1;
        if (conso == pattern_.size_) {
          // case where the whole forbidden pattern appears in a stretch
          // if the resource is soft, we must update arc base cost and
          // continue with the computation of endConsumption
          baseCost += cost_;
          break;
        }
      } else {
        break;
      }
    }
  }
  for (; itS1 != stretch.pShifts().end(); ++itS1, ++itD1) {
    // start checking from here only if it matches the first day-shift of the
    // pattern
    if (!pattern_.pAShifts_[0]->includes(**itS1)
        || !pattern_.pADays_[0]->includes(**itD1))
      continue;
    auto itSpat = pattern_.pAShifts_.begin() + 1;
    auto itDpat = pattern_.pADays_.begin() + 1;
    auto itS2 = itS1 + 1;
    auto itD2 = itD1 + 1;
    int conso = 1;
    for (; itS2 < stretch.pShifts().end();
           ++itS2, ++itD2, ++itSpat, ++itDpat) {
      if ((*itSpat)->includes(**itS2) && (*itDpat)->includes(**itD2)) {
        conso += 1;
        if (conso == pattern_.size_) {
          // case where the whole forbidden pattern appears in a stretch
          // if the resource is soft, we must update arc base cost and
          // continue with the computation of endConsumption
          baseCost += cost_;
          break;
        }
      } else {
        break;
      }
    }
  }
  return baseCost;
}

// TODO(AL): ensure later that domination could be improved.
//  Beware of the indices and conso = 0
DominationStatus SoftForbiddenPatternResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  int conso1 = pL1->getConsumption(id_);
  int conso2 = pL2->getConsumption(id_);
  // first, treat the  cases where the pattern cannot appear in the
  // propagation of L1 if it does not appear in that of L2
  if ((conso1 == 0) || (conso2 == conso1)) {
    return DOMINATED;
  }
  // compute the size of repeatStartPattern_[c2] \ repeatStartPattern_[c1]
  int nRepeat = 1;
  if (conso2 > conso1) {
    // when the consumption of pL2 is greater than that in pL1, we must check
    // if the piece of the pattern recorded in pL2 is not repeated at the end
    // of the piece recorded in pL1
    // this test must be done differently in forward and backward propagation,
    // which can be tested with the presence of an inArc in the label
    if (pL1->getInArc() != nullptr) {
      for (const auto &c : repeatStartPattern_[conso1]) {
        if (c == conso2) break;
        nRepeat++;
      }
    } else {
      for (const auto &c : repeatEndPattern_[conso1]) {
        if (c == conso2) break;
        nRepeat++;
      }
    }
  }

  // if no cost -> same as hard domination
  if (cost == nullptr) return NOT_DOMINATED;

  // in every other case, the pattern can appear in the propagation of L1,
  // but not in that of L2
  *cost += this->cost_ * nRepeat;
  return DOMINATED;
}

DominationStatus HardForbiddenPatternResource::dominates(
    RCLabel *pL1, RCLabel *pL2, double *cost) const {
  int conso1 = pL1->getConsumption(id_);
  int conso2 = pL2->getConsumption(id_);
  // first, treat the  cases where the pattern cannot appear in the
  // propagation of L1 if it does not appear in that of L2
  if ((conso1 == 0) || (conso2 == conso1))
    return DOMINATED;

  if (conso2 > conso1) {
    // when the consumption of pL2 is greater than that in pL1, we must check
    // if the piece of the pattern recorded in pL2 is not repeated at the end
    // of the piece recorded in pL1
    // this test must be done differently in forward and backward propagation,
    // which can be tested with the presence of an inArc in the label
    if (pL1->getInArc() != nullptr) {
      for (const auto &c : repeatStartPattern_[conso1])
        if (c == conso2) return DOMINATED;
    } else {
      for (const auto &c : repeatEndPattern_[conso1])
        if (c == conso2) return DOMINATED;
    }
  }

  // in every other case, the pattern can appear in the propagation of L1,
  // but not in that of L2
  return NOT_DOMINATED;
}

bool SoftForbiddenPatternResource::merge(const ResourceValues &vForward,
                                         const ResourceValues &vBack,
                                         ResourceValues *vMerged,
                                         const PRCLabel &pLMerged) {
  int forwardConso = vForward.consumption;
  int backwardConso = vBack.consumption;
  // if the pattern didn't start backward, nothing happen
  if (backwardConso == 0) {
    vMerged->consumption = forwardConso;
    return true;
  }
  // check if the forbidden pattern appears in the merged label, this is the
  // case only if the sum of backward and forwards consumptions exactly
  // matches the length of the pattern
  int missingConso = getUb() - backwardConso;
  if (missingConso == forwardConso) {
    vMerged->consumption = this->getUb();
    pLMerged->addBaseCost(this->cost_);
    return true;
  }
  for (const auto &conso : repeatStartPattern_[missingConso]) {
    if (forwardConso == conso) {
      vMerged->consumption = this->getUb();
      pLMerged->addBaseCost(this->cost_);
      return true;
    }
  }
  vMerged->consumption = 0;
  return true;
}

bool HardForbiddenPatternResource::merge(const ResourceValues &vForward,
                                         const ResourceValues &vBack,
                                         ResourceValues *vMerged,
                                         const PRCLabel &pLMerged) {
  int forwardConso = vForward.consumption;
  int backwardConso = vBack.consumption;
  // check if the forbidden pattern appears in the merged label, this is the
  // case only if the sum of backward and forwards consumptions exactly
  // matches the length of the pattern
  int missingConso = getUb() - backwardConso;
  for (const auto& conso : repeatStartPattern_[missingConso]) {
    if (forwardConso == conso) return false;
  }
  vMerged->consumption = 0;
  return true;
}

bool SoftForbiddenPatternExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  vChild->consumption = forwardConsumption_[vChild->consumption];

  // if the arc is a long stretch, the pattern may be completed and another
  // can be started after: so we set the consumption to the consumption at
  // the end of the stretch
  if (vChild->consumption == resource_.getUb()) {
    pLChild->addBaseCost(resource_.getCost());
    vChild->consumption = endConsumption_;
  }
  return true;
}

bool SoftForbiddenPatternExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  vChild->consumption = backwardConsumption_[vChild->consumption];

  // if the arc is a long stretch, the pattern may be completed and another
  // can be started after: so we set the consumption to the consumption at
  // the start of the stretch
  if (vChild->consumption == resource_.getUb()) {
    pLChild->addBaseCost(resource_.getCost());
    vChild->consumption = startConsumption_;
  }
  return true;
}

bool HardForbiddenPatternExpander::expand(const PRCLabel &pLChild,
                                          ResourceValues *vChild) {
  vChild->consumption = forwardConsumption_[vChild->consumption];

  return vChild->consumption < resource_.getUb();
}

bool HardForbiddenPatternExpander::expandBack(const PRCLabel &pLChild,
                                              ResourceValues *vChild) {
  vChild->consumption = backwardConsumption_[vChild->consumption];

  return vChild->consumption < resource_.getUb();
}

void SoftForbiddenPatternResource::preprocess(const PRCGraph &pRCGraph) {
  if (pattern_.size_ == 2) {
    // as a generic preprocessing, we only want to preprocess the forbidden
    // pattern including two shifts: this requires to modify the base costs of
    // the corresponding arcs
    for (const PRCArc &pA : pRCGraph->pArcs()) {
      double cost = 0;
      preprocess(pA, &cost);
      pA->addBaseCost(cost);
    }
    // mark the resource as preprocessed to ignore in the RCSPP
    isPreprocessed_ = true;
  }
}

bool SoftForbiddenPatternResource::preprocess(
    const PRCArc &pA, double *cost) {
  *cost = 0;
  if (pattern_.size_ > 2) return true;

  // if the resource is soft we update the base cost of the arc for each
  // occurrence of the forbidden sequence
  Stretch st = pA->stretch;

  // shifts and days of the pattern
  PAbstractShift firstShift = pattern_.pAShifts_[0];
  PAbstractDay firstDay = pattern_.pADays_[0];
  PAbstractShift secondShift = pattern_.pAShifts_[1];
  PAbstractDay secondDay = pattern_.pADays_[1];

  // day and shift at the origin of the arc
  PAbstractShift pSPrev = pA->origin->pAShift;
  PDay pDPrev = pA->origin->pDay;
  // true if the first element of the pattern is found on the previous
  // (day,shift)
  bool isStarted = false;
  if (firstShift->includes(*pSPrev) && firstDay->includes(*pDPrev))
    isStarted = true;

  auto itS = st.pShifts().begin();
  for (const auto& pD : st.pDays()) {
    if (isStarted && secondShift->includes(**itS)
        && secondDay->includes(*pD)) {
      // here, the forbidden pattern appears on the stretch :update the
      // base cost
      *cost += getCost();
    }
    // check if the shift is the first of the pattern before going to
    // next day
    if (firstShift->includes(**itS) && firstDay->includes(*pD))
      isStarted = true;
    else
      isStarted = false;
    ++itS;
  }
#ifdef NS_DEBUG
  double c = computeBaseCost(pA->stretch, pA->origin->pAShift);
  if (abs(c-*cost) > 1e-5) {
    std::cout << pA->stretch.toString() << std::endl;
    Tools::throwError("ForbiddenPattern: wrong cost.");
  }
#endif
  return true;
}

void HardForbiddenPatternResource::preprocess(const PRCGraph &pRCGraph) {
  if (pattern_.size_ == 2) {
    // as a generic preprocessing, we only want to preprocess the forbidden
    // pattern including two shifts: this requires to modify the base costs of
    // the corresponding arcs
    for (const PRCArc &pA : pRCGraph->pArcs()) {
      double cost = 0;
      if (!preprocess(pA, &cost))
        pRCGraph->forbidArc(pA, true);
    }
    // mark the resource as preprocessed to ignore in the RCSPP
    isPreprocessed_ = true;
  }
}

bool HardForbiddenPatternResource::preprocess(
    const PRCArc &pA, double *cost) {
  if (pattern_.size_ > 2) return true;

  // if the resource is soft we update the base cost of the arc for each
  // occurrence of the forbidden sequence
  Stretch st = pA->stretch;

  // shifts and days of the pattern
  PAbstractShift firstShift = pattern_.pAShifts_[0];
  PAbstractDay firstDay = pattern_.pADays_[0];
  PAbstractShift secondShift = pattern_.pAShifts_[1];
  PAbstractDay secondDay = pattern_.pADays_[1];

  // day and shift at the origin of the arc
  PAbstractShift pSPrev = pA->origin->pAShift;
  PDay pDPrev = pA->origin->pDay;
  // true if the first element of the pattern is found on the previous
  // (day,shift)
  bool isStarted = false;
  if (firstShift->includes(*pSPrev) && firstDay->includes(*pDPrev))
    isStarted = true;

  auto itS = st.pShifts().begin();
  for (const auto& pD : st.pDays()) {
    if (isStarted && secondShift->includes(**itS)
        && secondDay->includes(*pD)) {
      // here, the forbidden pattern appears on the stretch ->
      // mark arc as forbidden
      return false;
    }
    // check if the shift is the first of the pattern before going to
    // next day
    if (firstShift->includes(**itS) || firstDay->includes(*pD))
      isStarted = true;
    else
      isStarted = false;
  }
  return true;
}
