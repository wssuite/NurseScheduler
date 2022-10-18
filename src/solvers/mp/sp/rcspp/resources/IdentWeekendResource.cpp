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

#include "IdentWeekendResource.h"

PExpander SoftIdentWeekendResource::init(const AbstractShift &prevAShift,
                                         const Stretch &stretch,
                                         const shared_ptr<RCArc> &pArc,
                                         int indResource) {
/* Tools::throwError("This resource should never be initialized since it is "
                    "not a resource of the RCSPP");*/
  return nullptr;
}

void SoftIdentWeekendResource::preprocess(const PRCGraph &pRCGraph) {
  // no need to enumerate here, we just have to update the base costs of the
  // arcs by violating those that stop the succession of identical pAShift_
  for (const PRCArc &pA : pRCGraph->pArcs()) {
    double cost = 0;
    preprocess(pA, &cost);
    pA->addBaseCost(cost);
  }
  isPreprocessed_ = true;
}

bool SoftIdentWeekendResource::preprocess(const PRCArc& pA, double *cost) {
  *cost = 0;
  PAbstractShift pSPrev = pA->origin->pAShift;
  PDay pDPrev = pA->stretch.pDay(0)->previous();
  bool inWeekend = weekend_.isWeekend(pDPrev);
  // count the number of days since the beginning of the weekend
  double c = !inWeekend ? 0 : weekend_.positionInWeekend(*pDPrev) + 1;
  // the constraint is that the pShiftCmp_ must be true nbWeekendDays_,
  // if the nurse does only d days, she must pay
  // nbWeekendDays_ - d times the penalty.
  // Therefore, each time there is a change :
  // we pay the next weekend days (current one included) for the previous
  // shift if worked.
  // we pay the previous weekend days (current one excluded) for the current
  // shift if worked.
  auto itS = pA->stretch.pShifts().begin();
  for (const auto& pD : pA->stretch.pDays()) {
    if (inWeekend) {
      // if still in weekend
      if (weekend_.isWeekend(pD)) {
        // check if a change of type
        if (!pShiftCmp_->equals(*pSPrev, **itS)) {
          // if worked, we pay the next weekend days (current one included)
          // for the previous shift.
          if (pSPrev->isWork())
            *cost += (weekend_.nDays() - c) * cost_;
          // if worked, we pay the previous weekend days (current one excluded)
          // for the current shift.
          if ((*itS)->isWork())
            *cost += c * cost_;
        }
        c++;
      } else {
        // reset flag and counter
        inWeekend = false;
        c = 0;
      }
    } else if (weekend_.isWeekend(pD)) {
      inWeekend = true;
      c = 1;
    }
//      if (isSameType_ && !pSPrev->includes(**itS)) {
//        // the constraint is that the same shift type must be worked
//        // nbWeekendDays_, if the nurse does only d days, she must pay
//        // nbWeekendDays_ - d times the penalty
//        if (pSPrev->isWork()) {
//          int nbMissingDays = nbWeekendDays_ -
//              pD->positionInWeekend(firstWeekendDay_, lastWeekendDay_);
//          baseCost += nbMissingDays * cost_;
//        }
//        if ((*itS)->isWork()) {
//          int nbMissingDays =
//              pD->positionInWeekend(firstWeekendDay_, lastWeekendDay_);
//          baseCost += nbMissingDays * cost_;
//        }
//      }
//      // if the resource is not on shift with the same type, it requires
//      // only that the weekend be completely worked or rested. So we
//      // penalize only the transitions from work to rest and from rest to
//      // work
//      if (!isSameType_) {
//        if (pSPrev->isWork() && (*itS)->isRest()) {
//          int nbMissingDays = nbWeekendDays_ -
//           (pDPrev->positionInWeekend(firstWeekendDay_, lastWeekendDay_) + 1);
//          baseCost += nbMissingDays * cost_;
//        } else if (pSPrev->isRest() && (*itS)->isWork()) {
//          int nbMissingDays =
//              pD->positionInWeekend(firstWeekendDay_, lastWeekendDay_);
//          baseCost += nbMissingDays * cost_;
//        }
//      }
//    }
    pSPrev = *itS;
    itS++;
  }
  return true;
}

