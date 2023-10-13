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
#include <vector>

#include "tools/Tools.h"
#include "data/Shift.h"

class Scenario;

typedef std::shared_ptr<Scenario> PScenario;

class Demand;
typedef std::shared_ptr<Demand> PDemand;

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
             isOptDemand_(true) {}
  Demand(int nbDays, int firstDayId, int nbShifts, int nbSkills,
         std::string name, vector3D<int> minDemand);
  Demand(int nbDays, int firstDayId, int nbShifts, int nbSkills,
         std::string name, vector3D<int> minDemand, vector3D<int> optDemand);
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

  // minimum and optimal demand for each day, shift and skill
  //
  vector3D<int> minDemand_;
  bool isOptDemand_;
  vector3D<int> optDemand_;

 public:
  // total demand in the minimal and optimal demands
  //
  int minTotal_{}, optTotal_{};

  // preprocessed attributes aggregating the information of the demand
  //
  bool isPreprocessed_{};

  // total demand per skill in the minimal and optimal demands
  //
  std::vector<int> minPerSkill_, optPerSkill_;

  // total demand per shift in the minimal and optimal demands
  //
  std::vector<int> minPerShift_, optPerShift_;

  // total demand per day in the minimal and optimal demands
  //
  std::vector<int> minPerDay_, optPerDay_;

  // highest demands per skill over the considered period
  //
  std::vector<int> minHighestPerSkill_, optHighestPerSkill_;

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
  void preprocessMinDemand();
  void preprocessOptDemand();

  // add another week demand at the end of the current one
  // update all the parameters
  void pushBack(const PDemand& pDemand);

  // Returns a new demand that appends pDemand to the current one
  PDemand append(const PDemand& pDemand) const;

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
  int getTotalMinDuration(const PScenario &pScenario) const;
  int getTotalMinDemand() const;
  int getTotalMinDemand(int shift) const;
  int getTotalMinDemand(int shift, int skill) const;

  int getTotalDemand(
      int shift, int skill, const vector3D<int>& demand) const;
};

#endif  // SRC_DATA_DEMAND_H_
