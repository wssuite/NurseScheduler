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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_PREFERENCERESOURCE_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_PREFERENCERESOURCE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data/Nurse.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

class PreferenceResource : public Resource {
 public:
  explicit PreferenceResource(
      const std::string &_name,
      PAbstractDay pADay,
      Wish wish) :
      Resource(_name.empty() ? "Pref. on day "+std::to_string(pADay->getId())+
               " for "+wish.toString() : _name),
      pADay_(std::move(pADay)),
      wish_(std::move(wish)) {}

  int getConsumption(const State &initialState) const override {return 0;}

  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override = 0;

  bool isInRosterMaster() const override { return false; };
  bool isInRotationMaster() const override { return false; };

 protected:
  const PAbstractDay pADay_;
  const Wish wish_;
};

class SoftPreferenceResource : public PreferenceResource {
 public:
  explicit SoftPreferenceResource(
      PAbstractDay pADay, Wish wish, const std::string &_name = "") :
      PreferenceResource(_name, std::move(pADay), std::move(wish)) {
    costType_ = PREFERENCE_COST;
  }

  BaseResource* clone() const override {
    return new SoftPreferenceResource(pADay_, wish_);
  }

  // getters
  bool isHard() const override {return false;}

  // initialize the expander on a given arc
  PExpander init(const AbstractShift &prevAShift,
                 const Stretch &stretch,
                 const shared_ptr<RCArc> &pArc,
                 int indResource) override;

  // add the cost of preference violation to all the arcs of the input graph
  void preprocess(const PRCGraph &pRCGraph) override;
  bool preprocess(const PRCArc& pA, double *cost) override;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RESOURCES_PREFERENCERESOURCE_H_
