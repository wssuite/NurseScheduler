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

#include "Demand.h"

#include <cmath>
#include <algorithm>
#include <set>
#include <sstream>
#include <utility>

#include "data/Day.h"
#include "data/Scenario.h"

//-----------------------------------------------------------------------------
//
// C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

// constructor and destructor
//
Demand::Demand(int nbDays,
               int firstDayId,
               int nbShifts,
               int nbSkills,
               std::string name,
               vector3D<int> demand,
               DemandCtr ctr,
               double cost) :
    name_(std::move(name)),
    nDays_(nbDays),
    firstDayId_(firstDayId),
    nShifts_(nbShifts),
    nSkills_(nbSkills),
    demand_(std::move(demand)),
    cost_(cost),
    ctr_(ctr),
    demandTotal_(0) {
  // run the preprocessing
  this->preprocess();
}

Demand::~Demand() = default;

// compute all the potentially helpful attributes of a demand
// this includes the total demand per skill, per shift,
//
void Demand::preprocess() {
  // initialize the preprocessed vectors
  demandTotal_ = 0;
  Tools::initVector(&demandPerDay_, nDays_, 0);
  Tools::initVector(&demandPerShift_, nShifts_, 0);
  Tools::initVector(&demandPerSkill_, nSkills_, 0);
  Tools::initVector(&demandHighestPerSkill_, nSkills_, 0);

  for (int day = 0; day < nDays_; day++) {
    for (int shift = 1; shift < nShifts_; shift++) {
      for (int skill = 0; skill < nSkills_; skill++) {
        // update the total demand
        demandTotal_ += demand_[day][shift][skill];
        // update the demand per day
        demandPerDay_[day] += demand_[day][shift][skill];
        // update the demand per shift
        demandPerShift_[shift] += demand_[day][shift][skill];
        // update the demand per skill
        demandPerSkill_[skill] += demand_[day][shift][skill];
        // update the demand per day
        demandHighestPerSkill_[skill] =
            std::max(demand_[day][shift][skill], demandHighestPerSkill_[skill]);
      }
    }
  }
  isPreprocessed_ = true;
}

// add another week demand at the end of the current one
// update all the parameters
void Demand::pushBack(const PDemand& pDemand) {
  // check if same scenario
  if ((nShifts_ != pDemand->nShifts_) || (nSkills_ != pDemand->nSkills_)) {
    std::string error = "Demands are not compatible";
    Tools::throwError(error.c_str());
  }

  /*
   * Build new demand
   */

  // number of days covered by the demand and index of the first day
  //
  nDays_ += pDemand->nDays_;

  // pushes back the second demand on the first
  for (const vector2D<int> &vector : pDemand->demand_)
    demand_.push_back(vector);

  // run the preprocessing
  this->preprocess();
}

// Returns a new demand that appends pDemand to the current one
PDemand Demand::append(const PDemand& pDemand) const {
  PDemand pNewDemand = std::make_shared<Demand>(*this);
  pNewDemand->pushBack(pDemand);
  return pNewDemand;
}

int Demand::findMaxOptimalGap() const {
  if (costs_.empty())
    return isSoftCost(cost_) ? cost_ : 0;

  std::set<int> allGaps;
  for (const auto &costs2 : costs_)
    for (const auto &costs1 : costs2)
      for (double c : costs1)
        if (c > 1e-3 && isSoftCost(c))
          allGaps.insert(static_cast<int>(c));
  return Tools::gcd(allGaps);
}

// modify the demand by randomly swapping the demand of nnSwaps days
//
void Demand::swapDays(int nbSwaps) {
  for (int i = 0; i < nbSwaps; i++) {
    int day1 = Tools::randomInt(0, nDays_ - 1);
    int day2 = Tools::randomInt(0, nDays_ - 1);

    // save the demand on day 1
    vector2D<int> demandTmp = demand_[day1];

    // make the modification in the demand
    demand_[day1] = demand_[day2];
    demand_[day2] = demandTmp;
  }
}

// modify the demand by randomly swapping the demand of nbSwaps shifts
// the swapped shifts necessarily correspond to the same skill
//
void Demand::swapShifts(int nbSwaps) {
  for (int i = 0; i < nbSwaps; i++) {
    int sk = Tools::randomInt(0, nSkills_ - 1);
    int day1 = Tools::randomInt(0, nDays_ - 1);
    int sh1 =
        Tools::randomInt(1, nShifts_ - 1);  // make sure shift 0 is not taken
    int day2 = Tools::randomInt(0, nDays_ - 1);
    int sh2 =
        Tools::randomInt(1, nShifts_ - 1);  // make sure shift 0 is not taken

    // save the demand on day1/shift1
    int demandTmp;
    demandTmp = demand_[day1][sh1][sk];

    // make the modification in the demand
    demand_[day1][sh1][sk] = demand_[day2][sh2][sk];
    demand_[day2][sh2][sk] = demandTmp;
  }
}

// perturb the demand by adding demand in a number of shifts randomly chosen
// in the interval [minPerturb,maxPerturb]
// if the generate number is negative, then shifts are removed
// the perturbed shifts are also randomly chosen
// for a given skill the demand on a shift cannot become greater than the
// largest demand observed on the week
//
void Demand::perturbShifts(int minPerturb, int maxPerturb) {
  // generate the number of perturbations
  int nbPerturb = Tools::randomInt(minPerturb, maxPerturb);
  int valPerturb = (nbPerturb >= 0) ? 1 : -1;
  nbPerturb = std::abs(nbPerturb);

  // presolve the demand to find the highest demand per skill
  if (!isPreprocessed_) this->preprocess();

  int coTrials = 0;  // used to avoid infinite loop
  for (int i = 0; i < std::abs(nbPerturb); i++) {
    // draw the particular shift whose demand should be perturbed
    // select only a demand that is below the highest demand of the week for
    // this shift if demand should be added (this constraint is added to avoid
    // non feasibility due to the perturbation)
    bool isAtUpperBound = true;
    int day, sh, sk;
    while (isAtUpperBound && coTrials < 10 * nbPerturb) {
      day = Tools::randomInt(0, nDays_ - 1);
      sh = Tools::randomInt(1, nShifts_ - 1);
      sk = Tools::randomInt(0, nSkills_ - 1);
      isAtUpperBound = (valPerturb >= 0)
          && (demand_[day][sh][sk] >= demandHighestPerSkill_[sk]);
      coTrials++;
    }
    if (coTrials >= 10 * nbPerturb) break;

    // perturb the demand
    demand_[day][sh][sk] += valPerturb;
  }
}

// copy the input demand and apply a perturbation to generate random demand
//
PDemand Demand::randomPerturbation() {
  PDemand pDemand = std::make_shared<Demand>(*this);

  // three different types of perturbations are made
  // the order does not seem to be important
  pDemand->swapDays(nDays_ / 2);
  pDemand->swapShifts(nDays_ * nSkills_);
  pDemand->perturbShifts(-nDays_, nDays_);

  // get the main characteristics of the new demand
  pDemand->preprocess();

  return pDemand;
}

// Keep the preferences relative to the days in [begin,end)
PDemand Demand::keep(int begin, int end) {
  PDemand pDemand = std::make_shared<Demand>(*this);

  pDemand->demand_.erase(
          pDemand->demand_.begin() + end, pDemand->demand_.end());

  if (begin > 0)
    pDemand->demand_.erase(pDemand->demand_.begin(),
                           pDemand->demand_.begin() + begin);

  // update the number of days and indicate that this particular demand has not
  // been preprocessed
  pDemand->nDays_ = end - begin;
  pDemand->isPreprocessed_ = false;

  return pDemand;
}

// shorten the demand by keeping only the nbDays first days
//
void Demand::keepFirstNDays(int nbDays) {
  demand_.erase(demand_.begin() + nbDays, demand_.end());

  // update the number of days and indicate that this particular demand has not
  // been preprocessed
  nDays_ = nbDays;
  isPreprocessed_ = false;
}

// display the demand, and include the preprocessed information if the input
// boolean is set to true
//
std::string Demand::toString(bool withPreprocessedInfo) {
  std::stringstream rep;

  rep << "# " << std::endl;
  rep << "# DEMAND" << std::endl;

  // describe the demand being written
  //
  rep << "# " << std::endl;
  rep << "# Name of the demand: " << name_ << std::endl;

  rep << "# The demand refers to " << nSkills_ << " skills for ";
  rep << nShifts_ - 1 << " shifts per day on " << nDays_ << " days"
      << std::endl;
  rep << "# The demand is ";
  if (isGE()) rep << "a greater or equal";
  else if (isLE()) rep << "a lower or equal";
  else  rep << "an equality";
  rep << " constraint." << std::endl;

  // write the number of staff required per shift for each skill
  //
  rep << "#\t\t\t\t\t";
  for (int dayId = 0; dayId < 7; dayId++) {
    rep << " " << Day::toDayOfWeekShortName(dayId) << "\t";
  }
  rep << "# " << std::endl;
  for (int sh = 0; sh < nShifts_; sh++) {
    for (int sk = 0; sk < nSkills_; sk++) {
      // string str = "#   " + Tools::intToShift_[sh] + " "
      // + Scenario::intToSkill_[sk] + " ";
      // rep << str;
      // if(str.length() < 16) rep << "\t";
      rep << "#\tShift " << sh << " Skill " << sk << " " << "\t";
      for (int day = 0; day < 7; day++) {
        rep << "\t(" << demand_[day][sh][sk];
      }
      rep << std::endl;
    }
    rep << "# " << std::endl;
  }

  // write the preprocessed indicators if presolve was run
  //
  if (withPreprocessedInfo) {
    if (!isPreprocessed_)
      Tools::throwError("Trying to write the preprocessed information "
                        "of a demand that was not preprocessed!");
    // enumerate the global indicators
    //
    rep << "# " << std::endl;
    rep << "# PREPROCESSED DATA ON THE DEMAND" << std::endl;

    rep << "# Total demand = " << demandTotal_ << std::endl;

    rep << "# " << std::endl;
    rep << "# Demand per day" << std::endl;
    for (int i = 0; i < nDays_; i++) {
      rep << "#\t\t" << Day::toDayOfWeekShortName(i + firstDayId_) <<
          " (" << i + firstDayId_ << ")" << ": " << demandPerDay_[i]
          << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Demand per shift" << std::endl;
    rep << "#\t\tShift 0 is rest;" << std::endl;
    for (int i = 1; i < nShifts_; i++) {
      rep << "#\t\tShift " << i << ": " << demandPerShift_[i] << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Demand per skill" << std::endl;
    for (int i = 0; i < nSkills_; i++) {
      rep << "#\t\tSkill " << i << ": " << demandPerSkill_[i] << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Highest demand per skill for one shift" << std::endl;
    for (int i = 0; i < nSkills_; i++) {
      rep << "#\t\tSkill " << i << ": " << demandHighestPerSkill_[i]
          << std::endl;
    }
  }

  return rep.str();
}

// compute the total duration needed for the associated demand
int Demand::getTotalDuration(const PScenario &pScenario) const {
  int duration = 0;
  for (const PShift &pS : pScenario->pShifts())
    duration += getTotalDemand(pS->id) * pS->duration;
  return duration;
}

int Demand::getTotalDemand() const {
  int n = 0;
  for (int sh=0; sh < nShifts_; ++sh)
    n += getTotalDemand(sh);
  return n;
}

int Demand::getTotalDemand(int shift) const {
  int n = 0;
  for (int sk=0; sk < nSkills_; ++sk)
    n += getTotalDemand(shift, sk);
  return n;
}

int Demand::getTotalDemand(int shift, int skill) const {
  int n = 0;
  for (const auto & demandPerDay : demand_)
    n += demandPerDay[shift][skill];
  return n;
}

double Demand::coverageCost(int coverage, int day, int sh, int sk) const {
  int diff = coverage - demand_[day][sh][sk];
  double c = cost(day, sh, sk);
  if (isEQ())
    return c * abs(diff);
  if (isGE())
    return c * std::max(0, -diff);
  return c * std::max(0, diff);
}
