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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;


/**
 * Resource corresponding to the soft min/max constraints on the total
 * number of worked weekends
 * Plain expanders are used, because precomputation is more complex than
 * useful for total weekend resources
 */
class ConsWeekend : public ConsShift {
 public:
  ConsWeekend(PAbstractShift pShift,
              BoundedResource* pR,
              DayOfWeek firstWeekendDay, DayOfWeek lastWeekendDay,
              bool endOnLastDay,
              bool cyclic) :
      ConsShift(pShift, endOnLastDay, cyclic),
      pR_(pR), weekend_(firstWeekendDay, lastWeekendDay) {}

  void computeConsumption(
      const Stretch &stretch, ResourceValues *vChild,
      const std::function<bool(ResourceValues *vChild)>& processReset) const;

  void computeConsumptionBack(
      const Stretch &stretch, ResourceValues *vChild,
      const std::function<bool(ResourceValues *vChild)>& processReset) const;

  // basic getters
  BoundedResource *pResource() const { return pR_; }

  // wrappers for weekend Tools methods
  int nWeekendsInStretch(const Stretch& stretch) const {
    std::pair<int, int> firstLastDays = pR_->getFirstLastDays(stretch);
    return weekend_.nWeekendsInInterval(firstLastDays.first,
                                        firstLastDays.second);
  }

  const Weekend& weekend() const { return weekend_; }

 protected:
  BoundedResource *pR_;
  Weekend weekend_;
};
/**
 * Resource corresponding to the soft min/max constraints on the number of
 * consecutive shifts on weekends: it may be used for any kind of abstract shift
 */
class SoftConsWeekendShiftResource :
    public SoftBoundedResource, public ConsWeekend {
 public:
  SoftConsWeekendShiftResource(
      int lb, int ub, double lbCost, double ubCost,
      const PAbstractShift &pShift, int totalNbDays,
      DayOfWeek firstWeekendDayId = SATURDAY,
      DayOfWeek lastWeekendDayId = SUNDAY,
      bool endOnLastDay = true,
      bool cyclic = false) :
      SoftBoundedResource("Soft Weekend Cons "+ pShift->name,
                          lb, ub, lbCost, ubCost),
      ConsWeekend(pShift, this, firstWeekendDayId, lastWeekendDayId,
                  endOnLastDay, cyclic) {
    totalNbDays_ = totalNbDays;
    costType_ = CONS_WEEKEND_SHIFTS_COST;
  }

  BaseResource* clone() const override {
    return new SoftConsWeekendShiftResource(
        lb_, ub_, lbCost_, ubCost_, pAShift_, totalNbDays_,
        weekend_.firstWeekendDay().getDayOfWeek(),
        weekend_.lastWeekendDay().getDayOfWeek(),
        lastDayEndsSequence_, cyclic_);
  }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  int getConsumption(const State &initialState) const override;

  // WARNING: the resource is always active as
  // the consumption is carried out from one weekend to the next
//  // the resource needs to be checked for dominance only on nodes
//  // corresponding to the one checked in this constraint
//  bool isActive(int dayId, const PAbstractShift &pAShift) const override {
//    return weekend_.isWeekend(dayId) && pAShift_->includes(*pAShift);
//  }

  bool isAnyWorkShiftResource() const override { return
        pAShift_->isAnyWork();
  }

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

class HardConsWeekendShiftResource :
    public HardBoundedResource, public ConsWeekend {
 public:
  HardConsWeekendShiftResource(
      int lb, int ub, const PAbstractShift &pShift,
      int totalNbDays,
      DayOfWeek firstWeekendDay = SATURDAY,
      DayOfWeek lastWeekendDay = SUNDAY,
      bool endOnLastDay = true,
      bool cyclic = false) :
      HardBoundedResource("Hard Weekend Cons "+pShift->name, lb, ub),
      ConsWeekend(pShift, this, firstWeekendDay, lastWeekendDay,
                  endOnLastDay, cyclic) {
    totalNbDays_ = totalNbDays;
  }

  BaseResource* clone() const override {
    return new HardConsWeekendShiftResource(
        lb_, ub_, pAShift_, totalNbDays_,
        weekend_.firstWeekendDay().getDayOfWeek(),
        weekend_.lastWeekendDay().getDayOfWeek(),
        lastDayEndsSequence_, cyclic_);
  }

  int getConsumption(const State &initialState) const override;

  bool isAnyWorkShiftResource() const override {
    return pAShift_->isAnyWork();
  }

  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost) const override;

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return true; };

 protected:
  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/*
 * Expanders for the  resources
 */

struct ConsWeekendShiftExpander : public Expander {
  ConsWeekendShiftExpander(int rId, CostType type, const Stretch& stretch,
                           bool cyclic,
                           int nWeekendsBefore, int nWeekendsAfter) :
      Expander(rId, type),
      nWeekendsBefore(nWeekendsBefore),
      nWeekendsAfter(nWeekendsAfter),
      stretch_(stretch),
      cyclic(cyclic) {}

 protected:
  // number of weekends after the end of the stretch,
  // if only a portion of a weekend is left (e.g., only sunday),
  // it still counts this weekend
  const int nWeekendsAfter, nWeekendsBefore;
  const Stretch stretch_;
  bool cyclic;  // true of solving the cyclic version
};


struct SoftConsWeekendShiftExpander : public ConsWeekendShiftExpander {
  SoftConsWeekendShiftExpander(int indResource,
                               const SoftConsWeekendShiftResource& resource,
                               const Stretch& stretch,
                               bool cyclic,
                               int nWeekendsBefore,
                               int nWeekendsAfter):
      ConsWeekendShiftExpander(indResource,
                               CONS_WEEKEND_SHIFTS_COST,
                               std::move(stretch),
                               cyclic,
                               nWeekendsBefore,
                               nWeekendsAfter),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const SoftConsWeekendShiftResource& resource_;
};

struct HardConsWeekendShiftExpander : public ConsWeekendShiftExpander {
  HardConsWeekendShiftExpander(int indResource,
                               const HardConsWeekendShiftResource& resource,
                               const Stretch& stretch,
                               bool cyclic,
                               int nWeekendsBefore,
                               int nWeekendsAfter):
      ConsWeekendShiftExpander(indResource,
                               NO_COST,
                               stretch,
                               cyclic,
                               nWeekendsBefore,
                               nWeekendsAfter),
      resource_(resource) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;

  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const HardConsWeekendShiftResource& resource_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_CONSWEEKENDSHIFTRESOURCE_H_
