/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, Samuel Rosat, and
 * contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 *  full license detail.
 */

#ifndef SRC_SOLVERS_SOLVERINPUT_H_
#define SRC_SOLVERS_SOLVERINPUT_H_

#include <vector>

#include "tools/MyTools.h"
#include "data/Nurse.h"
#include "data/Roster.h"
#include "data/Scenario.h"


//-----------------------------------------------------------------------------
//
//  C l a s s   S o l v e r I n p u t
//
//  All the information to make an instance
//
//-----------------------------------------------------------------------------

class SolverInput {
 public:
  // Constructor and Destructor
  SolverInput();
  ~SolverInput();

// Attributes are public because they are const
 protected:
  // Recall the "const" attributes as pointers :
  // Nurses and Scenario informations
  //
  PScenario pScenario_;
  std::vector<Nurse> *pTheNurses_;

  // Pointers to the minimum and optimum demand for each day, shift and skill
  //
  vector3D<int> *pMinDemand_, pOptDemand_;

  // pointer to the preferences of the nurses nurse
  // (that vector must be of same length and in the same order as the nurses)
  PPreferences pPreferences_;

  // pointer to the state of each nurse at the beginning of the time horizon
  //
  std::vector<State> *pInitState_;
};

#endif  // SRC_SOLVERS_SOLVERINPUT_H_
