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

#ifndef SRC_SOLVERS_MP_SP_ROSTERSP_H_
#define SRC_SOLVERS_MP_SP_ROSTERSP_H_

#include <vector>

#include "solvers/mp/sp/SubProblem.h"

class RosterSP : public SubProblem {
 public:
  // Constructor that correctly sets the resource (time + bounds),
  // but NOT THE COST
  RosterSP(PScenario scenario,
           int nbDays,
           PConstContract contract,
           std::vector<State> *pInitState);

  virtual ~RosterSP();

 protected:
  //----------------------------------------------------------------
  //
  // Update of the costs / network for solve function
  //
  //----------------------------------------------------------------

  // FUNCTIONS -- SOLVE
  // return a function that will post process any path found by the RC graph
  RCSPPSolver *initRCSSPSolver() override;

  // Creates all nodes of the rcspp (including resource window)
  void createNodes() override;

  // override creation of arcs source -> principal
  void createArcsSourceToPrincipal() override;

  // override creation of arcs principal -> principal
  void createArcsPrincipalToPrincipal() override;

  // override creation of pricing arcs principal -> sink
  void createArcsPrincipalToSink() override;

  // override updateArcCosts to take into account start and end work costs
  // when going fom/to a rest shift
  void updateArcCosts() override;
};

#endif  // SRC_SOLVERS_MP_SP_ROSTERSP_H_
