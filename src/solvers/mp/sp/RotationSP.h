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

#ifndef SRC_SOLVERS_MP_SP_ROTATIONSP_H_
#define SRC_SOLVERS_MP_SP_ROTATIONSP_H_

#include <vector>

#include "solvers/mp/sp/SubProblem.h"

class RotationSP : public SubProblem {
 public:
  RotationSP() = default;

  RotationSP(PScenario scenario,
             int nbDays,
             PConstContract contract,
             std::vector<State> *pInitState);

  virtual ~RotationSP();

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

  // override updateArcCosts to take into account end work costs when going
  // to a daily sink
  void updateArcCosts() override;
};

#endif  // SRC_SOLVERS_MP_SP_ROTATIONSP_H_
