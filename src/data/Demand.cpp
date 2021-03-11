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

#include <math.h>

#include <sstream>
#include <algorithm>

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
               int firstDay,
               int nbShifts,
               int nbSkills,
               std::string name,
               vector3D<int> minDemand,
               vector3D<int> optDemand) : name_(name),
                                          nDays_(nbDays),
                                          firstDay_(firstDay),
                                          nShifts_(nbShifts),
                                          nSkills_(nbSkills),
                                          minDemand_(minDemand),
                                          optDemand_(optDemand),
                                          minTotal_(0),
                                          optTotal_(0),
                                          isPreprocessed_(false) {
  // run the preprocessing
  this->preprocessDemand();
}

Demand::~Demand() {}

// compute all the potentially helpful attributes of a demand
// this includes the total demand per skill, per shift,
//
void Demand::preprocessDemand() {
  // initialize the preprocessed vectors
  Tools::initVector(&minPerDay_, nDays_, 0);
  Tools::initVector(&optPerDay_, nDays_, 0);
  Tools::initVector(&minPerShift_, nShifts_, 0);
  Tools::initVector(&optPerShift_, nShifts_, 0);
  Tools::initVector(&minPerSkill_, nSkills_, 0);
  Tools::initVector(&optPerSkill_, nSkills_, 0);
  Tools::initVector(&minHighestPerSkill_, nSkills_, 0);
  Tools::initVector(&optHighestPerSkill_, nSkills_, 0);

  for (int day = 0; day < nDays_; day++) {
    for (int shift = 1; shift < nShifts_; shift++) {
      for (int skill = 0; skill < nSkills_; skill++) {
        // update the total demand
        minTotal_ += minDemand_[day][shift][skill];
        optTotal_ += optDemand_[day][shift][skill];

        // update the demand per day
        minPerDay_[day] += minDemand_[day][shift][skill];
        optPerDay_[day] += optDemand_[day][shift][skill];

        // update the demand per shift
        minPerShift_[shift] += minDemand_[day][shift][skill];
        optPerShift_[shift] += optDemand_[day][shift][skill];

        // update the demand per skill
        minPerSkill_[skill] += minDemand_[day][shift][skill];
        optPerSkill_[skill] += optDemand_[day][shift][skill];

        // update the demand per day
        minHighestPerSkill_[skill] =
            std::max(minDemand_[day][shift][skill], minHighestPerSkill_[skill]);
        optHighestPerSkill_[skill] =
            std::max(optDemand_[day][shift][skill], optHighestPerSkill_[skill]);
      }
    }
  }
  isPreprocessed_ = true;
}

// add another week demand at the end of the current one
// update all the parameters
void Demand::pushBack(PDemand pDemand) {
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
  for (const vector2D<int> &vector : pDemand->minDemand_)
    minDemand_.push_back(vector);
  for (const vector2D<int> &vector : pDemand->optDemand_)
    optDemand_.push_back(vector);

  // run the preprocessing
  this->preprocessDemand();
}

// Returns a new demand that appends pDemand to the current one
PDemand Demand::append(PDemand pDemand) {
  PDemand bigDemand = std::make_shared<Demand>(*this);

  // check if same scenario
  if ((nShifts_ != pDemand->nShifts_) || (nSkills_ != pDemand->nSkills_)) {
    std::string error = "Demands are not compatible";
    Tools::throwError(error.c_str());
  }

  // number of days covered by the demand and index of the first day
  //
  bigDemand->nDays_ += pDemand->nDays_;

  // pushes back the second demand on the first
  for (const vector2D<int> &vector : pDemand->minDemand_) {
    bigDemand->minDemand_.push_back(vector);
  }

  for (const vector2D<int> &vector : pDemand->optDemand_) {
    bigDemand->optDemand_.push_back(vector);
  }

  // run the preprocessing
  bigDemand->preprocessDemand();

  return bigDemand;
}

// modify the demand by randomly swapping the demand of nnSwaps days
//
void Demand::swapDays(int nbSwaps) {
  for (int i = 0; i < nbSwaps; i++) {
    int day1 = Tools::randomInt(0, nDays_ - 1);
    int day2 = Tools::randomInt(0, nDays_ - 1);

    // save the demand on day 1
    vector2D<int> minDemandTmp = minDemand_[day1],
        optDemandTmp = optDemand_[day1];

    // make the modification in the demand
    minDemand_[day1] = minDemand_[day2];
    optDemand_[day1] = minDemand_[day2];

    minDemand_[day2] = minDemandTmp;
    optDemand_[day2] = optDemandTmp;
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
    int minDemandTmp, optDemandTmp;
    minDemandTmp = minDemand_[day1][sh1][sk];
    optDemandTmp = optDemand_[day1][sh1][sk];

    // make the modification in the demand
    minDemand_[day1][sh1][sk] = minDemand_[day2][sh2][sk];
    optDemand_[day1][sh1][sk] = optDemand_[day2][sh2][sk];

    minDemand_[day2][sh2][sk] = minDemandTmp;
    optDemand_[day2][sh2

    ][sk] = optDemandTmp;
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

  // preprocess the demand to find the highest demand per skill
  if (!isPreprocessed_) this->preprocessDemand();

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
      isAtUpperBound = (valPerturb >= 0) ? (minDemand_[day][sh][sk]
          >= minHighestPerSkill_[sk]) : false;
      coTrials++;
    }
    if (coTrials >= 10 * nbPerturb) break;

    // perturb the demand
    minDemand_[day][sh][sk] += valPerturb;
    optDemand_[day][sh][sk] += valPerturb;
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
  pDemand->preprocessDemand();

  return pDemand;
}

// Keep the preferences relative to the days in [begin,end)
PDemand Demand::keep(int begin, int end) {
  PDemand pDemand = std::make_shared<Demand>(*this);

  pDemand->minDemand_.erase(pDemand->minDemand_.begin() + end,
                            pDemand->minDemand_.end());
  pDemand->optDemand_.erase(pDemand->optDemand_.begin() + end,
                            pDemand->optDemand_.end());

  if (begin > 0) {
    pDemand->minDemand_.erase(pDemand->minDemand_.begin(),
                              pDemand->minDemand_.begin() + begin);
    pDemand->optDemand_.erase(pDemand->optDemand_.begin(),
                              pDemand->optDemand_.begin() + begin);
  }

  // update the number of days and indicate that this particular demand has not
  // been preprocessed
  pDemand->nDays_ = end - begin;
  pDemand->isPreprocessed_ = false;

  return pDemand;
}

// shorten the demand by keeping only the nbDays first days
//
void Demand::keepFirstNDays(int nbDays) {
  minDemand_.erase(minDemand_.begin() + nbDays, minDemand_.end());
  optDemand_.erase(optDemand_.begin() + nbDays, optDemand_.end());

  // update the number of days and indicate that this particular demand has not
  // been preprocessed
  nDays_ = nbDays;
  isPreprocessed_ = false;
}

// shorten the demand by removing the nbDays first days
//
void Demand::removeFirstNDays(int nbDays) {
  minDemand_.erase(minDemand_.begin(), minDemand_.begin() + nbDays);
  optDemand_.erase(optDemand_.begin(), optDemand_.begin() + nbDays);

  // update the number of days and indicate that this particular demand has not
  // been preprocessed
  nDays_ = nDays_ - nbDays;
  isPreprocessed_ = false;
}

// remove a list of skills from the demand
//
void Demand::removeSkills(std::vector<int> skills) {
  // make sure that the indices of skills are ordered
  std::stable_sort(skills.begin(), skills.end());

  // erase the skills for each day/shift
  for (int day = 0; day < nDays_; day++) {
    for (int shift = 0; shift < nShifts_; shift++) {
      for (int skill = skills.size() - 1; skill >= 0; skill--) {
        minDemand_[day][shift].erase(minDemand_[day][shift].begin() + skill);
        optDemand_[day][shift].erase(optDemand_[day][shift].begin() + skill);
      }
    }
  }
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
  rep << std::endl;

  // write the number of staff required per shift for each skill
  //
  rep << "#\t\t\t\t\t";
  for (int dayId = 0; dayId < 7; dayId++) {
    rep << " " << Tools::intToDay(dayId) << "\t";
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
        rep << "\t(" << minDemand_[day][sh][sk] << ","
            << optDemand_[day][sh][sk] << ")";
      }
      rep << std::endl;
    }
    rep << "# " << std::endl;
  }

  // write the preprocessed indicators if preprocess was run
  //
  if (withPreprocessedInfo) {
    if (!isPreprocessed_)
      Tools::throwError("Trying to write the preprocessed information "
                        "of a demand that was not preprocessed!");
    // enumerate the global indicators
    //
    rep << "# " << std::endl;
    rep << "# PREPROCESSED DATA ON THE DEMAND" << std::endl;

    rep << "# Total minimum demand = " << minTotal_ << std::endl;
    rep << "# Total optimal demand = " << optTotal_ << std::endl;

    rep << "# " << std::endl;
    rep << "# Demand per day" << std::endl;
    for (int i = 0; i < nDays_; i++) {
      rep << "#\t\t" << Tools::intToDay(i + firstDay_) << " (" << i + firstDay_
          << ")" << ": ";
      rep << "minimum = " << minPerDay_[i] << " ; optimal = " << optPerDay_[i];
      rep << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Demand per shift" << std::endl;
    rep << "#\t\tShift 0 is rest;" << std::endl;
    for (int i = 1; i < nShifts_; i++) {
      rep << "#\t\tShift " << i << ": ";
      rep << "minimum = " << minPerShift_[i] << " ; optimal = "
          << optPerShift_[i];
      rep << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Demand per skill" << std::endl;
    for (int i = 0; i < nSkills_; i++) {
      rep << "#\t\tSkill " << i << ": ";
      rep << "minimum = " << minPerSkill_[i] << " ; optimal = "
          << optPerSkill_[i];
      rep << std::endl;
    }

    rep << "# " << std::endl;
    rep << "# Highest demand per skill for one shift" << std::endl;
    for (int i = 0; i < nSkills_; i++) {
      rep << "#\t\tSkill " << i << ": ";
      rep << "minimum = " << minHighestPerSkill_[i] << " ; optimal = "
          << optHighestPerSkill_[i];
      rep << std::endl;
    }
  }

  return rep.str();
}
