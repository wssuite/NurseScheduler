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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROTATIONSP_H_
#define SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROTATIONSP_H_

#include <vector>

#include "SubProblem.h"

namespace boostRCSPP {

class RotationSP : public SubProblem {
 public:
  RotationSP(PScenario scenario,
             int nDays,
             PLiveNurse pNurse,
             SubProblemParam param);


  virtual ~RotationSP();

  double startWorkCost(int a) const override;

  double shiftCost(int a, const PAbstractShift &prevS) const override;

  double endWorkCost(int a) const override;

  void computeCost(MasterProblem *, RCSolution *rcSol) const override;

 protected:
  // Index: (shiftType, day) of origin
  vector2D<int> arcsPrincipalToSink;

  // Creates all nodes of the rcspp (including resource window)
  void createNodes() override;

  // override creation of arcs source -> principal
  void createArcsSourceToPrincipal() override;

  // override creation of arcs principal -> principal
  void createArcsPrincipalToPrincipal() override;

  // override creation of pricing arcs principal -> sinks
  void createArcsPrincipalToSink() override;

  // override updateArcDualCosts to take into account end work costs when going
  // to a daily sink
  void updateArcDualCosts() override;
};

}  // namespace boostRCSPP

#endif  // SRC_SOLVERS_MP_SP_RCSPP_BOOST_ROTATIONSP_H_
