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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FORBIDDENPATTERNRESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FORBIDDENPATTERNRESOURCE_H_

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data/Shift.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::vector;

/**
 * Resource corresponding to the soft min/max constraints on the number of
 * consecutive shifts: it may be used for any kind of abstract shift
 */
class ForbiddenPatternResource : public Resource {
 public:
  explicit ForbiddenPatternResource(
      std::string name, Pattern  pattern) :
      Resource(std::move(name)), pattern_(std::move(pattern)) {
    preprocessPattern();
  }

  int getConsumption(const State &initialState) const override { return 0;}

  // the resource needs to be checked for dominance only on nodes
  // corresponding to the one checked in this constraint
  bool isActive(int dayId, const PAbstractShift &pAShift) const override {
    return true;
  }

  int getUb() const {
    return pattern_.size_;
  }

  virtual double getCost() const {return 0;}

  const Pattern& getPattern() const { return pattern_;}

  const std::vector<int>& repeatStartPattern(int i) const {
    return repeatStartPattern_[i];
  }

  const std::vector<int>& repeatEndPattern(int i) const {
    return repeatEndPattern_[i];
  }

 protected:
  void preprocessPattern();

  const Pattern pattern_;
  /*
 * List the subpatterns that appear at the start or end of the pattern and
 * which are repeated several times in the pattern; this will be useful
 * to detect forbidden patterns when expanding or back-expanding a label in
 * the label-setting algorithm
 */
  // FOR FORWARD EXPANSION:
  // if vector repeatStartPattern_[i] includes index j (> i), this means that
  // the starting subpattern from index 0 to i is repeated from index j-i to j
  // example: consider pattern Day-Late-Day-Late:
  // starting subpattern Day appear twice and ends at indices 0 and 2, so
  // repeatStartPattern_[0] = [2]
  // starting subpattern Day-Late appear twice, ending at indices 1 and 3,
  // so repeatStartPattern_[1] = [3]
  // starting subpattern Day-Late-Day appears only once,  so
  // repeatStartPattern_[2] = []
  vector<vector<int>> repeatStartPattern_;
  // FOR BACKWARD EXPANSION:
  // if vector repeatEndPattern_[i] includes index j (> i), this means that
  // the ending subpattern from index size_-1-i to size_-1 appears also from
  // index size-1-j to size_-1-j+i
  // example: consider pattern Day-Late-Day-Late once again:
  // ending subpattern Late is appears twice, starting at indices 3-0 and
  // 3-2, so repeatEndPattern_[0] = [2]
  // ending subpattern Day-Late appears twice, starting at indices 3-1 and
  // 3-3, so repeatEndPattern_[1] = [3]
  // ending subpattern Late-Day-Late appears only once, so
  // repeatEndPattern_[2] = []
  vector<vector<int>> repeatEndPattern_;
};

class SoftForbiddenPatternResource : public ForbiddenPatternResource {
 public:
  explicit SoftForbiddenPatternResource(
      Pattern pattern,
      double cost,
      const std::string& _name = "Soft forbidden pattern cons") :
      ForbiddenPatternResource(_name, std::move(pattern)), cost_(cost) {
    costType_ = FORBIDDEN_PATTERN_COST;
  }

  BaseResource* clone() const override {
    return new SoftForbiddenPatternResource(
        pattern_, cost_, name);
  }

  bool isHard() const override {return false;}

  double getCost() const override {return cost_;}

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;

  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost) const override;

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  bool isInRosterMaster() const override { return false; };
  // WARNING: we will hit an issue when the pattern includes rest shifts
  // in the middle for the rotation decomposition ; this is not handled at all
  bool isInRotationMaster() const override {
    for (int k = 0; k < pattern_.size_; k++)
      if (pattern_.pAShifts_[k]->isRest() && k > 0 && k < pattern_.size_ - 1) {
        std::cerr << "ForbiddenPatternResource cannot handle a rest shift "
                     "in the middle of the pattern if using a "
                     "rotation decomposition." << std::endl;
        return true;
      }
    return false;
  };

  double computeBaseCost(const Stretch &stretch,
                         const PAbstractShift &pPrevAShift) const;

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  const double cost_;
};

class HardForbiddenPatternResource : public ForbiddenPatternResource {
 public:
  explicit HardForbiddenPatternResource(
      Pattern  pattern,
      const std::string& name = "Hard forbidden pattern cons") :
      ForbiddenPatternResource(name, std::move(pattern)) {
    costType_ = NO_COST;
  }

  BaseResource* clone() const override {
    return new HardForbiddenPatternResource(pattern_, name);
  }

  bool isHard() const override {return true;}

  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;

  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost) const override;

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;

  bool isInRosterMaster() const override { return false; };
  // TODO(AL): we will hit an issue when the pattern includes rest shifts in
  //  the rotation decomposition ; this is not handled at all
  bool isInRotationMaster() const override {
    for (const auto& pS : pattern_.pAShifts_)
      if (pS->isRest())
        Tools::throwError("ForbiddenPatternResource cannot handle a rest shift "
                          "in the pattern if using a rotation decomposition.");
    return false;
  };

 protected:
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;
};

/*
 * Expander for the  resources
 */

struct SoftForbiddenPatternExpander : public Expander {
  SoftForbiddenPatternExpander(int indResource,
                               const SoftForbiddenPatternResource& resource,
                               vector<int> forwardConsumption,
                               vector<int> backwardConsumption,
                               int endConsumption,
                               int startConsumption):
      Expander(indResource, FORBIDDEN_PATTERN_COST),
      resource_(resource),
      forwardConsumption_(std::move(forwardConsumption)),
      backwardConsumption_(std::move(backwardConsumption)),
      endConsumption_(endConsumption),
      startConsumption_(startConsumption) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const SoftForbiddenPatternResource& resource_;
  // consumption on the arc for each possible value of consumption in a label
  // in forward label-setting
  vector<int> forwardConsumption_;
  // consumption on the arc for each possible value of consumption in a label
  // in backward label-setting
  vector<int> backwardConsumption_;
  // below, the values are useful only to get the consumption of a label if
  // the pattern is completed while expanding along this arc
  // a. resource consumption at the end of the stretch:
  int endConsumption_;
  // b. resource consumtion at the start of the stretch:
  int startConsumption_;
};

struct HardForbiddenPatternExpander : public Expander {
  HardForbiddenPatternExpander(int indResource,
                               const HardForbiddenPatternResource& resource,
                               vector<int> forwardConsumption,
                               vector<int> backwardConsumption,
                               int end = 0,
                               int start = 0):
      Expander(indResource, NO_COST),
      resource_(resource),
      forwardConsumption_(std::move(forwardConsumption)),
      backwardConsumption_(std::move(backwardConsumption)) {}

  bool expand(const PRCLabel &pLChild, ResourceValues *vChild) override;
  bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) override;

 private:
  const HardForbiddenPatternResource& resource_;
  // consumption on the arc for each possible value of consumption in a label
  // in forward label-setting
  vector<int> forwardConsumption_;
  // consumption on the arc for each possible value of consumption in a label
  // in backward label-setting
  vector<int> backwardConsumption_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_FORBIDDENPATTERNRESOURCE_H_
