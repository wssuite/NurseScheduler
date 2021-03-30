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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROSTERSP_H_
#define SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROSTERSP_H_

#include <vector>

#include "SubProblem.h"

namespace boostRCSPP {

class RosterSP : public SubProblem {
 public:
  // Constructor that correctly sets the resource (time + bounds),
  // but NOT THE COST
  RosterSP(PScenario scenario,
           int nDays,
           PLiveNurse pNurse,
           SubProblemParam param);

  virtual ~RosterSP();

  double startWorkCost(int a) const override;

  double shiftCost(int a) const override;

  double endWorkCost(int a) const override;

  void computeCost(MasterProblem *, RCSolution *rcSol) const override;

 protected:
  //----------------------------------------------------------------
  //
  // Update of the costs / network for solve function
  //
  //----------------------------------------------------------------

  // FUNCTIONS -- SOLVE
  // return a function that will post process any path found by the RC graph
  BoostRCSPPSolver *initRCSSPSolver() override;

  // Creates all nodes of the rcspp (including resource window)
  void createNodes() override;

  // override creation of arcs source -> principal
  void createArcsSourceToPrincipal() override;

  // override creation of arcs principal -> principal
  void createArcsPrincipalToPrincipal() override;

  // override creation of pricing arcs principal -> sink
  void createArcsPrincipalToSink() override;

  // override updateArcDualCosts to take into account start and end work costs
  // when going fom/to a rest shift
  void updateArcDualCosts() override;
};

}  // namespace boostRCSPP

#endif  // SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROSTERSP_H_
