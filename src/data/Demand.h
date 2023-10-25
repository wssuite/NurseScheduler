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

#ifndef SRC_DATA_DEMAND_H_
#define SRC_DATA_DEMAND_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tools/Tools.h"
#include "data/Shift.h"

class Scenario;

typedef std::shared_ptr<Scenario> PScenario;

class Demand;
typedef std::shared_ptr<Demand> PDemand;

enum DemandCtr {
    D_LE, D_EQ, D_GE
};

//-----------------------------------------------------------------------------
//
// C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

class Demand {
 public:
  // generic constructor and destructor
  Demand() : nDays_(0), firstDayId_(0), nShifts_(0), nSkills_(0),
             ctr_(D_GE), cost_(INFEAS_COST) {}
  Demand(int nbDays, int firstDayId, int nbShifts, int nbSkills,
         std::string name, vector3D<int> demand, DemandCtr ctr = D_GE,
         double cost = INFEAS_COST);
  ~Demand();

  // constant attributes of the demand
  //
 public:
  // name of the demand
  //
  std::string name_;

  // number of days covered by the demand and index of the first day
  //
  int nDays_;
  const int firstDayId_;

  // number of shifts per day and number of skills to cover
  const int nShifts_, nSkills_;

  // demand for each day, shift and skill
  vector3D<int> demand_;

  // default cost
  double cost_;

  // cost for each day, shift and skill
  // if empty, cost is used
  vector3D<double> costs_;

  // type
  DemandCtr ctr_ = D_GE;

 public:
  // total demand
  //
  int demandTotal_{};

  // preprocessed attributes aggregating the information of the demand
  //
  bool isPreprocessed_{};

  // total demand per skill
  //
  std::vector<int> demandPerSkill_;

  // total demand per shift
  //
  std::vector<int> demandPerShift_;

  // total demand per day
  //
  std::vector<int> demandPerDay_;

  // highest demands per skill over the considered period
  //
  std::vector<int> demandHighestPerSkill_;

 protected:
  // modify the demand by randomly swapping the demand of nnSwaps days
  //
  void swapDays(int nbSwaps);

  // modify the demand by randomly swapping the demand of nbSwaps shifts
  // the swapped shifts necessarily correspond to the same skill
  //
  void swapShifts(int nbSwaps);

  // perturb the demand by adding demand in a number of shifts randomly chosen
  // in the interval [minPerturb,maxPerturb]
  // the perturbed shifts are also randomly chosen
  // for a given skill the demand on a shift cannot become greater than the
  // largest demand observed on the week
  //
  void perturbShifts(int minPerturb, int maxPerturb);

 public:
  // compute all the potentially helpful attributes of a demand
  // this includes the total demand per skill, per shift,
  void preprocess();

  const vector3D<int> & demand() {
    return demand_;
  }

  // add another week demand at the end of the current one
  // update all the parameters
  void pushBack(const PDemand& pDemand);

  // Returns a new demand that appends pDemand to the current one
  PDemand append(const PDemand& pDemand) const;

  // return the maximum gap to be optimal for the demand
  int findMaxOptimalGap() const;

  // display the demand, and include the preprocessed information if the input
  // boolean is set to true
  //
  std::string toString(bool withPreprocessedInfo);

  // copy the input demand and apply a perturbation to generate random demand
  //
  PDemand randomPerturbation();

  // shorten the demand by keeping only the nbDays in [begin, end)
  // return a new demand
  PDemand keep(int begin, int end);

  // shorten the demand by keeping only the nbDays first days
  //
  void keepFirstNDays(int nbDays);

  // compute the total duration needed for the associated demand
  int getTotalDuration(const PScenario &pScenario) const;
  int getTotalDemand() const;
  int getTotalDemand(int shift) const;
  int getTotalDemand(int shift, int skill) const;

  bool isGE() const {
    return ctr_ == D_GE;
  }

  bool isLE() const {
    return ctr_ == D_LE;
  }

  bool isEQ() const {
    return ctr_ == D_EQ;
  }

  double cost(int day, int shift, int skill) const {
    if (costs_.empty()) return cost_;
    return costs_.at(day).at(shift).at(skill);
  }

  void setCost(double cost) {
    costs_.clear();
    cost_ = cost;
  }

  void setCosts(vector3D<double> costs) {
    costs_ = std::move(costs);
  }

  double coverageCost(int coverage, int day, int sh, int sk) const;
};

#endif  // SRC_DATA_DEMAND_H_
