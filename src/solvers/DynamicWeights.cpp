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

#include "DynamicWeights.h"

#include <algorithm>

#include "solvers/Solver.h"
#include "tools/Tools.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

using std::vector;

DynamicWeights::DynamicWeights(
    PScenario pScenario, std::vector<PLiveNurse> pLiveNurses):
    pScenario_(std::move(pScenario)), pLiveNurses_(std::move(pLiveNurses)),
    version_(0) {
  clear();  // fill all vectors with a 0 for each nurse
}

//-----------------------------------------------------------------------------
// Compute the weights of the violation of the min/max number of working days
// For now, the update depends only on the initial states and on the contract
// of the nurses, on the number of days on the demand, on the number of weeks
// already treated and on the number of weeks left
// The required data on the nurses is mostly computed in preprocessTheNurses
//-----------------------------------------------------------------------------
void DynamicWeights::boundsAndWeights(WeightStrategy strategy, int horizon) {
  // clear all vectors and increment the version depending on the strategy.
  // the version is incremented if any member is modified.
  clear();
  version_++;
  switch (strategy) {
    case MAX :
    case MEAN:
      computeWeightsTotalShiftsForPrimalDual(strategy, horizon);
      break;

    case RANDOMMEANMAX:
      if (Tools::randomInt(0, 1) == 0)
        computeWeightsTotalShiftsForPrimalDual(MEAN, horizon);
      else
        computeWeightsTotalShiftsForPrimalDual(MAX, horizon);
      break;

    case NO_STRAT:
      computeWeightsTotalShiftsForStochastic(horizon);
      break;

    default: Tools::throwError("Weight/bound strategy not defined");
      break;
  }
}

void DynamicWeights::computeWeightsTotalShiftsForStochastic(int horizon) {
  // The important value to infer the importance of respecting the strict
  // constraints on the total number of working days/weekends is the remaining
  // number of days/weekends after the demand currently treated
  int remainingWeekends = pScenario_->nWeeks() - pScenario_->thisWeek();
  int remainingDays = 7 * remainingWeekends - horizon;
  double factorRemainingDays = remainingDays / (7.0 * pScenario_->nWeeks());

  // portion of these remaining weekends in this demand
  int nbWeeks = horizon / 7;
  double factorRemainingWeekends =
      static_cast<double>(nbWeeks) / remainingWeekends;

  // Compute the non-penalized intervals and the associated penalties
  for (const auto &pN : pLiveNurses_) {
    // compute the interval that must be respected to have a chance of
    // not paying penalties in the future
    const auto &pTotalDuration = pN->totalShiftDurationResource_;
    pTotalDuration->setLb(pN->minWorkDaysNoPenaltyTotalDays_);
    pTotalDuration->setLbCost(pScenario_->weights().totalShifts);
    pTotalDuration->setUb(pN->maxWorkDaysNoPenaltyTotalDays_);
    pTotalDuration->setUbCost(pScenario_->weights().totalShifts);

    // penalize each worked weekend with primal-dual costs
    const auto &pTotalWeekend = pN->totalWeekendResource_;
    pTotalWeekend->setUb(0);
    int remainingWeekendsToWork =
        pN->maxTotalWeekends() - pN->pStateIni_->totalWeekendsWorked_;

    // no need to penalize twice with maxTotalWeekendsAvg_
    if (remainingWeekendsToWork <= 0) {
      pTotalWeekend->setUbCost(pScenario_->weights().totalWeekends);
      maxTotalWeekendsAvg_[pN->num_] = 0;
      weightTotalWeekendsAvg_[pN->num_] = 0;
    } else {
      pTotalWeekend->setUbCost(
          pScenario_->weights().totalWeekends / pN->maxTotalWeekends()
              * pN->pStateIni_->totalWeekendsWorked_);
      // compute the proportion of weekends that can be worked in this demand
      // without exceeding the max in the future.
      // round with a certain probability to the floor or the ceil
      double numberOfAuthorizedWeekend =
          remainingWeekendsToWork * factorRemainingWeekends;
      maxTotalWeekendsAvg_[pN->num_] =
          Tools::roundWithProbability(numberOfAuthorizedWeekend);
      weightTotalWeekendsAvg_[pN->num_] = pScenario_->weights().totalWeekends;
    }

    double factorMarginOnAvg =
        0.25 * factorRemainingDays * 7.0 / horizon;
    factorMarginOnAvg = (horizon % 14 == 0) ? .0 : factorMarginOnAvg;
    minTotalShiftsAvg_[pN->num_] =
        (1.0 - factorMarginOnAvg) * pN->minAvgWorkDaysNoPenaltyTotalDays_;
    maxTotalShiftsAvg_[pN->num_] =
        (1.0 + factorMarginOnAvg) * pN->maxAvgWorkDaysNoPenaltyTotalDays_;
    weightTotalShiftsAvg_[pN->num_] = pScenario_->weights().totalShifts;
  }
}

void DynamicWeights::computeWeightsTotalShiftsForPrimalDual(
    WeightStrategy strategy, int horizon) {
  vector<double> maxPrimalDualCostForContractDays(pScenario_->nContracts());
  vector<double> maxPrimalDualCostForContractWE(pScenario_->nContracts());
  vector<double> meanPrimalDualCostForContractDays(pScenario_->nContracts());
  vector<double> meanPrimalDualCostForContractWE(pScenario_->nContracts());
  vector<int> nbNursesPerContract(pScenario_->nContracts());

  // initialize the non-penalized intervals and the associated penalties
  // for each contract
  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    minTotalShiftsContractAvg_[p] = 0;
    maxTotalShiftsContractAvg_[p] = 0;
    weightTotalShiftsContractAvg_[p] = pScenario_->weights().totalShifts;
    maxTotalWeekendsContractAvg_[p] = 0;
    weightTotalWeekendsContractAvg_[p] = pScenario_->weights().totalWeekends;
  }

  // Compute the non-penalized intervals and the associated penalties
  // for each contract
  for (const auto &pN : pLiveNurses_) {
    // Deactivate minimum constraints
    const auto &pTotalDuration = pN->totalShiftDurationResource_;
    pTotalDuration->setLb(0);
    pTotalDuration->setLbCost(0);

    // Activate maximum constraints (work) with primal-dual cost
    // Remark: these costs are always nonnegative here because,
    // for the competition, nurses are always in shortage.
    // If this does not apply -> find another formulation
    pTotalDuration->setUb(0);
    // Primal-dual cost of max working days
    double w = pScenario_->weights().totalShifts / pN->maxTotalShifts()
        * pN->pStateIni_->totalTimeWorked_;
    // Must not be higher that WEIGHT
    if (w > pScenario_->weights().totalShifts)
      w = pScenario_->weights().totalShifts;
    pTotalDuration->setUbCost(w);

    // Activate maximum constraints (WE) with primal-dual cost
    const auto &pTotalWeekend = pN->totalWeekendResource_;
    pTotalWeekend->setUb(0);
    w = pScenario_->weights().totalWeekends / pN->maxTotalWeekends()
        * pN->pStateIni_->totalWeekendsWorked_;
    if (w > pScenario_->weights().totalWeekends)
      w = pScenario_->weights().totalWeekends;
    pTotalWeekend->setUbCost(w);

    // Update average days/weekends in contract
    int p = pN->pContract_->id_;
    minTotalShiftsContractAvg_[p] +=
        std::max(0, pN->minTotalShifts() - pN->pStateIni_->totalTimeWorked_);
    maxTotalShiftsContractAvg_[p] +=
        std::max(0, pN->maxTotalShifts() - pN->pStateIni_->totalTimeWorked_);
    maxTotalWeekendsContractAvg_[p] += std::max(
        0, pN->maxTotalWeekends() - pN->pStateIni_->totalWeekendsWorked_);

    // Check for maximum primal-dual costs
    maxPrimalDualCostForContractDays[p] = std::max(
        maxPrimalDualCostForContractDays[p], pTotalDuration->getUbCost());
    maxPrimalDualCostForContractWE[p] = std::max(
        maxPrimalDualCostForContractWE[p], pTotalWeekend->getUbCost());
    meanPrimalDualCostForContractDays[p] += pTotalDuration->getUbCost();
    meanPrimalDualCostForContractWE[p] += pTotalWeekend->getUbCost();
    nbNursesPerContract[p]++;
  }

  // round the min/max values of the interval associated to the contract
  int remainingWeeks = pScenario_->nWeeks() - pScenario_->thisWeek();
  int remainingDays = 7 * remainingWeeks;
  int nWeekendsInHorizon = (horizon + 1) / 7;
  double ratioDays = horizon * (1.0 / remainingDays);
  double ratioWeekends = nWeekendsInHorizon * (1.0 / remainingWeeks);

  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    minTotalShiftsContractAvg_[p] =
        Tools::roundWithProbability(minTotalShiftsContractAvg_[p] * ratioDays);
    maxTotalShiftsContractAvg_[p] =
        Tools::roundWithProbability(maxTotalShiftsContractAvg_[p] * ratioDays);
    maxTotalWeekendsContractAvg_[p] = Tools::roundWithProbability(
        maxTotalWeekendsContractAvg_[p] * ratioWeekends);

    meanPrimalDualCostForContractDays[p] /= nbNursesPerContract[p];
    meanPrimalDualCostForContractWE[p] /= nbNursesPerContract[p];

    // to penalize some contracts vs the others
    PContract pContract = pScenario_->pContract(p);
    int lengthStintMin =
        pContract->minConsDaysWork_ + pContract->maxConsDaysOff_;
    // compute the number of days worked, if the nurse works the minimum days
    // without penalties
    int minWorkDaysNoPenalty = pContract->minConsDaysWork_ *
        (7 * pScenario_->nWeeks()) / lengthStintMin;
    // compute the ratio between the min work days without penalties and the
    // maximum number of shifts
    double ratioMinMax =
        minWorkDaysNoPenalty * 1.0 / pContract->maxTotalShifts_;
    double weightContract = 1 + ratioMinMax / 10;
    weightTotalShiftsContractAvg_[p] *= weightContract;

    if (strategy == MAX) {
      weightTotalShiftsContractAvg_[p] -= maxPrimalDualCostForContractDays[p];
      weightTotalWeekendsContractAvg_[p] -= maxPrimalDualCostForContractWE[p];
    } else if (strategy == MEAN) {
      weightTotalShiftsContractAvg_[p] -= meanPrimalDualCostForContractDays[p];
      weightTotalWeekendsContractAvg_[p] -= meanPrimalDualCostForContractWE[p];
    }

    if (false) {
      std::cout << "# " << std::endl;
      std::cout << "##################################################"
                << std::endl;
      std::cout << "# " << (*(pScenario_->pContract(p)))
                << std::endl;
      std::cout << "# min/max : " << std::endl;
      std::cout << "#    | min total shifts: " << minTotalShiftsContractAvg_[p]
                << std::endl;
      std::cout << "#    | max total shifts: " << maxTotalShiftsContractAvg_[p]
                << std::endl;
      std::cout << "#    | max total we    : "
                << maxTotalWeekendsContractAvg_[p] << std::endl;
      std::cout << "# costs   : " << std::endl;
      std::cout << "#    | cost shift: " << weightTotalShiftsContractAvg_[p]
                << std::endl;
      std::cout << "#    | cost we   : " << weightTotalWeekendsContractAvg_[p]
                << std::endl;
      std::cout << "##################################################"
                << std::endl;
      std::cout << "# " << std::endl;
    }
  }
}

void DynamicWeights::clear() {
  // clear the vectors that are about to be filled
  minTotalShiftsAvg_ = vector<double>(pScenario_->nNurses());
  maxTotalShiftsAvg_ = vector<double>(pScenario_->nNurses());
  weightTotalShiftsAvg_ = vector<double>(pScenario_->nNurses());
  maxTotalWeekendsAvg_ = vector<double>(pScenario_->nNurses());
  weightTotalWeekendsAvg_ = vector<double>(pScenario_->nNurses());

  minTotalShiftsContractAvg_ = vector<double>(pScenario_->nContracts());
  maxTotalShiftsContractAvg_ = vector<double>(pScenario_->nContracts());
  weightTotalShiftsContractAvg_ = vector<double>(pScenario_->nContracts());
  maxTotalWeekendsContractAvg_ = vector<double>(pScenario_->nContracts());
  weightTotalWeekendsContractAvg_ = vector<double>(pScenario_->nContracts());
}
