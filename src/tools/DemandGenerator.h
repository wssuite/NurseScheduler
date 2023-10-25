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

#ifndef SRC_TOOLS_DEMANDGENERATOR_H_
#define SRC_TOOLS_DEMANDGENERATOR_H_

#include <map>
#include <string>
#include <vector>

#include "data/Scenario.h"
#include "data/Demand.h"


//-----------------------------------------------------------------------------
//
// C l a s s  D e m a n d G e n e r a t o r
//
// Generate and store random demand scenarios for stochastic solution of the
// problem
//
//-----------------------------------------------------------------------------

class DemandGenerator {
 public:
  // default constructor and destructor
  DemandGenerator(int nbDemands,
                  int nbDays,
                  vector2D<PDemand> demands,
                  PScenario pScenario) :
      nbDemandsToGenerate_(nbDemands),
      nbDaysInGeneratedDemands_(nbDays),
      demandHistory_(demands),
      pScenario_(pScenario),
      rdm_(Tools::getANewRandomGenerator()) {
  }
  ~DemandGenerator();

 public:
  // check the feasibility of a demand scenario
  bool checkDemandFeasibility(PDemand pDemand);

  // generate nbScenarios_ through perturbations of the demand history
  vector2D<PDemand> generatePerturbedDemands();

  // generate 1 demand through perturbations of the demand history
  vector<PDemand> generateSinglePerturbatedDemands(
          bool checkFeasibility = true);

 protected:
  // number of demand scenarios that should be generated
  int nbDemandsToGenerate_;

  // number of days that must be considered in each generated demand
  int nbDaysInGeneratedDemands_;

  // demand history from which the random scenarios should be generated
  vector2D<PDemand> demandHistory_;

  // nurse rostering scenario under study.
  // this attribute is necessary to check the feasibility
  // of the generated demands.
  PScenario pScenario_;

  // random generator
  std::minstd_rand rdm_;
};

#endif  // SRC_TOOLS_DEMANDGENERATOR_H_
