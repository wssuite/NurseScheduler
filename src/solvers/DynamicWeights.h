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

#ifndef SRC_SOLVERS_DYNAMICWEIGHTS_H_
#define SRC_SOLVERS_DYNAMICWEIGHTS_H_

#include <map>
#include <memory>
#include <vector>
#include <string>
#include <utility>

#include "data/Scenario.h"
#include "Parameters.h"


class LiveNurse;
typedef std::shared_ptr<LiveNurse> PLiveNurse;

class DynamicWeights {
 public:
  DynamicWeights(PScenario pScenario, std::vector<PLiveNurse> pLiveNurses);

  // Compute the weights o the violation of the min/max number of working days
  // For now, the update depends only on the initial states and on the contract
  // of the nurses, on the number of days on the demand, on the number of weeks
  // already treated and on the number of weeks left
  // The required data on the nurses is mostly computed in preprocessTheNurses
  void boundsAndWeights(WeightStrategy strategy, int horizon);

  void computeWeightsTotalShiftsForStochastic(int horizon);

  void computeWeightsTotalShiftsForPrimalDual(
      WeightStrategy strategy, int horizon);

  int version() const {
    return version_;
  }

  const vector<double> &getMinTotalShiftsAvg() const {
    return minTotalShiftsAvg_;
  }
  const vector<double> &getMaxTotalShiftsAvg() const {
    return maxTotalShiftsAvg_;
  }
  const vector<double> &getWeightTotalShiftsAvg() const {
    return weightTotalShiftsAvg_;
  }
  const vector<double> &getMaxTotalWeekendsAvg() const {
    return maxTotalWeekendsAvg_;
  }
  const vector<double> &getWeightTotalWeekendsAvg() const {
    return weightTotalWeekendsAvg_;
  }
  const vector<double> &getMinTotalShiftsContractAvg() const {
    return minTotalShiftsContractAvg_;
  }
  const vector<double> &getMaxTotalShiftsContractAvg() const {
    return maxTotalShiftsContractAvg_;
  }
  const vector<double> &getMaxTotalWeekendsContractAvg() const {
    return maxTotalWeekendsContractAvg_;
  }
  const vector<double> &getWeightTotalShiftsContractAvg() const {
    return weightTotalShiftsContractAvg_;
  }
  const vector<double> &getWeightTotalWeekendsContractAvg() const {
    return weightTotalWeekendsContractAvg_;
  }

 protected:
  PScenario pScenario_;
  std::vector<PLiveNurse> pLiveNurses_;

  // clear all vectors
  void clear();

  // increment this flag everytime a set of modifications is performed
  int version_;

  // Interval inside of which there is no penalty for the total number of
  // working days (for each nurse)
  // This interval is computed from the max/min number of working days averaged
  // over the number of remaining weeks
  std::vector<double> minTotalShiftsAvg_;
  std::vector<double> maxTotalShiftsAvg_;

  // Penalties for values outside of [minTotalShiftsAvg_,maxTotalShiftsAvg_]
  std::vector<double> weightTotalShiftsAvg_;

  // Number of worked weekends below which there is no penalty for the
  // total number of working weekends
  // This interval is computed from the max number of working weekends averaged
  // over the number of remaining weeks
  std::vector<double> maxTotalWeekendsAvg_;

  // Penalties for the number of working weekends on the current period
  // (for each nurse)
  std::vector<double> weightTotalWeekendsAvg_;

  // Number of min, max and weekends allowed for all nurses under the same
  // contract
  std::vector<double> minTotalShiftsContractAvg_,
      maxTotalShiftsContractAvg_,
      maxTotalWeekendsContractAvg_;
  // Penalties for exceeding the average number of shifts allowed for all the
  // nurses under the same contract
  std::vector<double> weightTotalShiftsContractAvg_,
      weightTotalWeekendsContractAvg_;
};

#endif  // SRC_SOLVERS_DYNAMICWEIGHTS_H_
