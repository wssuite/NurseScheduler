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

#ifndef SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_
#define SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_

#include <memory>
#include <vector>

#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/Solver.h"

#include "OsiSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_solution.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_vector.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_var.hpp"  // NOLINT (suppress cpplint error)

struct HeuristicMIP {
  ~HeuristicMIP();
  std::recursive_mutex mutex_;
  OsiSolverInterface *pInitialSolver_, *pSolver_ = nullptr;
  BCP_solution_generic * sol_ = nullptr;
  bool deleted = false;

  void clone() {
    // get a copy of the solver
    pSolver_ = pInitialSolver_->clone(true);
    pSolver_->messageHandler()->logLevel(
        pInitialSolver_->messageHandler()->logLevel());
  }
};

class BcpHeuristics {
 public:
  BcpHeuristics(BcpModeler *pModel, SolverType type, int verbosity = 0);

  ~BcpHeuristics();

  BCP_solution_generic * rounding(OsiSolverInterface *solver,
                                  const BCP_vec<BCP_var*> &vars);

  BCP_solution_generic * mip_solve(OsiSolverInterface *solver,
                                   const BCP_vec<BCP_var*> &vars,
                                   SolverType type = CBC);

  void stop();

 protected:
  BCP_solution_generic * buildSolution(
      OsiSolverInterface *solver, const std::vector<BcpColumn*> &columns) const;

  BcpModeler *pModel_;
  shared_ptr<HeuristicMIP> mip_;
};

#endif  // SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_
