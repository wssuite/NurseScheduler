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

#include "FreeDaysAfterShiftResource.h"

// The enumeration creates one arc from each node with shift pAShift_ to the
// rest nodes at days +2,...,+nbFreeDays_ if they do not exist.
void SoftFreeDaysAfterShiftResource::enumerate(const PRCGraph &pRCGraph,
                                               bool forceEnum)  {
  auto it = std::find_if(
      pRCGraph->pShifts().begin(), pRCGraph->pShifts().end(),
      [](const PShift &pS) {
        return pS->isRest();
      });
  if (it == pRCGraph->pShifts().end())
    Tools::throwError("There is no rest shift in the RCGraph "
                      "for the resource SoftFreeDaysAfterShiftResource.");
  PShift pRestShift = *it;
  totalNbDays_ = pRCGraph->nDays();
  for (const PRCNode &pN : pRCGraph->pNodes()) {
    // only consider the nodes that correspond to a shift included in
    // pAShift_
    if (!pAShift_->includes(*pN->pAShift)) continue;

    // no enumeration needed if at a sink node
    if (pN->type == SINK_NODE) continue;

    // create the missing arcs
    // the arc corresponding to exactly one rest should exist, we thus add
    // only the arcs corresponding to more than one rest
    Stretch stretch(pN->pDay->next(), pRestShift);
    for (int nbRestDays = 2 ; nbRestDays <= nbFreeDays_ ; ++nbRestDays) {
      if (stretch.lastDayId() + 1 >= totalNbDays_) break;
      stretch.pushBack(pRestShift);
      for (const auto& pTarget :
           pRCGraph->pNodesPerDayId(stretch.lastDayId()))
        if (!pTarget->isWorkNode())
          PRCArc pA = pRCGraph->addSingleArc(pN, pTarget, stretch);
    }
  }
}

void SoftFreeDaysAfterShiftResource::preprocess(const PRCGraph &pRCGraph) {
  // no need to enumerate here, we just have to update the base costs of the
  // arcs by violating those that stop the succession of identical pAShift_
  for (const PRCArc &pA : pRCGraph->pArcs()) {
    double cost = 0;
    preprocess(pA, &cost);
    pA->addBaseCost(cost);
  }
  isPreprocessed_ = true;
}

bool SoftFreeDaysAfterShiftResource::preprocess(
    const PRCArc& pA, double *cost) {
  *cost = 0;
  bool startCounting = false;  // true when pAShift_ is met
  int count = 0;  // nb of free days after a pAShift_
  if (pAShift_->includes(*pA->origin->pAShift)) startCounting = true;
  auto itS = pA->stretch.pShifts().begin();
  for (const auto& pS : pA->stretch.pShifts()) {
    if (startCounting) {
      if (pS->isRest()) {
        count++;
      } else if ((count >= 1) || !(pAShift_->includes(*pS))) {
        if (count < nbFreeDays_) {
          *cost += cost_;
        }
        startCounting = false;
        count = 0;
      }
    } else if (pAShift_->includes(*pS)) {
      startCounting = true;
    }
  }
  // pay cost if stretch ends before last day and not enough free days
  // has been counted
  if ((count >= 1) && (count < nbFreeDays_)
      && (pA->stretch.lastDayId() + 1 < totalNbDays_) )
    *cost += cost_;
  return true;
}

PExpander SoftFreeDaysAfterShiftResource::init(const AbstractShift &prevAShift,
                                               const Stretch &stretch,
                                               const shared_ptr<RCArc> &pArc,
                                               int indResource) {
//  Tools::throwError("This resource should never be initialized since it is "
//                    "not a resource of the RCSPP");
  return nullptr;
}

PExpander FreeDaysAfterShiftResource::init(const AbstractShift &prevAShift,
                                           const Stretch &stretch,
                                           const shared_ptr<RCArc> &pArc,
                                           int indResource) {
//  Tools::throwError("This resource should never be initialized since it is "
//                    "not a resource of the RCSPP");
  return nullptr;
}
