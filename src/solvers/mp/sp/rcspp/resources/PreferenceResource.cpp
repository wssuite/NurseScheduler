/*
 * Copyright (C) 2021 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "PreferenceResource.h"

void SoftPreferenceResource::preprocess(const PRCGraph &pRCGraph) {
  for (const PRCArc& pA : pRCGraph->pArcs()) {
    double cost = 0;
    preprocess(pA, &cost);
    pA->addBaseCost(cost);
  }
  isPreprocessed_ = true;
}

bool SoftPreferenceResource::preprocess(const PRCArc& pA, double *cost) {
  *cost = 0;
  auto itS = pA->stretch.pShifts().begin();
  for (const auto& pDay : pA->stretch.pDays()) {
    if (pADay_->includes(*pDay))
      *cost += wish_.cost(*itS);
    itS++;
  }
  return true;
}

PExpander SoftPreferenceResource::init(const AbstractShift &prevAShift,
                                       const Stretch &stretch,
                                       const shared_ptr<RCArc> &pArc,
                                       int indResource) {
  /*Tools::throwError("This resource should never be initialized since it is "
                    "not a resource of the RCSPP");*/
  return nullptr;
}
