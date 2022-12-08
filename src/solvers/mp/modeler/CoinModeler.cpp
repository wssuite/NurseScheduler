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

#include "solvers/mp/modeler/CoinModeler.h"

#ifdef USE_GUROBI
#include <mutex>  // NOLINT (suppress cpplint error)

OsiGrbSolverInterface *pGlobalGRBSolver = nullptr;
std::mutex m;
#endif

OsiSolverInterface * getNewSolver(SolverType type) {
  switch (type) {
    case CLP: return new OsiClpSolverInterface();
    case Gurobi: {
#ifdef USE_GUROBI
      // create a global environment, so only one global env is created
      // very useful when running on a server with a license server
      // -> the token is retrieve only once at the first call and kept while
      // the program is alive
      std::lock_guard<std::mutex> l(m);
      if (pGlobalGRBSolver == nullptr) {
        pGlobalGRBSolver = new OsiGrbSolverInterface();
        // delete the global environment at exit
        atexit([]() {
          delete pGlobalGRBSolver;
        });
      }
      return new OsiGrbSolverInterface();
#endif
    }
    case Cplex:
#ifdef USE_CPLEX
      return new OsiCpxSolverInterface();
#endif
    case CBC:
#ifdef USE_CBC
    {
      auto pSolver = new OsiCbcSolverInterface();
      pSolver->getModelPtr()->setLogLevel(0);
      return pSolver;
    }
#endif
    case FirstAvailable: {
      SolverType type2 = getFirstSolverTypeAvailable();
      return getNewSolver(type2);
    }
    default: break;
  }
  return nullptr;
}
