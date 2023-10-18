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

#include "AlternativeShiftResource.h"


void AlternativeShiftResource::preprocess(const PRCGraph &pRCGraph) {
  for (const PRCArc &pA : pRCGraph->pArcs()) {
    double cost = 0;
    preprocess(pA, &cost);
    pA->addBaseCost(cost);
    if (cost >= LARGE_SCORE)
      pRCGraph->forbidArc(pA);
  }
  isPreprocessed_ = true;
}

bool AlternativeShiftResource::preprocess(const PRCArc& pA, double *cost) {
  *cost = 0;
  for (const auto &pS : pA->stretch.pShifts())
    if (isAlternativeShift(pS))
      *cost += cost_;
  return true;
}

PExpander AlternativeShiftResource::init(const AbstractShift &prevAShift,
                                       const Stretch &stretch,
                                       const shared_ptr<RCArc> &pArc,
                                       int indResource) {
  /*Tools::throwError("This resource should never be initialized since it is "
                    "not a resource of the RCSPP");*/
  return nullptr;
}
