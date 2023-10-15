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

void PreferenceResource::preprocess(const PRCGraph &pRCGraph) {
  for (const PRCArc& pA : pRCGraph->pArcs()) {
    double cost = 0;
    if (preprocess(pA, &cost)) {
      pA->addBaseCost(cost);
    } else {
      // forbid arc forever. Must be changed manually if needed
      pRCGraph->forbidArc(pA, true);
    }
  }
  isPreprocessed_ = true;
}

bool PreferenceResource::preprocess(const PRCArc& pA, double *cost) {
  *cost = 0;
  auto itS = pA->stretch.pShifts().begin();
  for (const auto& pDay : pA->stretch.pDays()) {
    if (!pEndShift_->equals(**itS) && pADay_->includes(*pDay)) {
      *cost += wish_.cost(*itS);
      if (wish_.forbid(*itS))
        return false;
    }
    itS++;
  }

  // check also the previous day to ensure to forbid
  // all the arcs leaving a forbidden day
  if (pADay_->includes(*pA->origin->pDay)) {
    for (const auto &pS : pA->origin->pAShift->pIncludedShifts())
      if (!wish_.forbid(pS))
        return true;
    // all included shifts are forbidden
    return false;
  }

  return true;
}
