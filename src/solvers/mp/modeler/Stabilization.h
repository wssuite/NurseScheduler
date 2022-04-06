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

#ifndef SRC_SOLVERS_MP_MODELER_STABILIZATION_H_
#define SRC_SOLVERS_MP_MODELER_STABILIZATION_H_

#include <vector>

#include "solvers/Solver.h"

#include "OsiSolverInterface.hpp"  // NOLINT (suppress cpplint error)


class Modeler;
class MyVar;
class MyCons;

class Stabilization {
  //---------------------------------------------------------------------------
  //
  // STAB: Methods required to implement stabilization in the column generation
  //
  // Ref: Lübbecke, Marco E., and Jacques Desrosiers.
  // "Selected topics in column generation."
  // Operations research 53.6 (2005): 1007-1023.
  //
  //---------------------------------------------------------------------------
 public:
  explicit Stabilization(Modeler *pModel);

  // Add stabilization variables z for the box [b_, b+] with the penalties c
  // if getting outside of the box:
  // dual = obj += - c_+ z_+ - c_- z_-, s.t.: b_- - z_-<= Pi <= b_+ + z_+
  // primal = obj += -b_- y_- + b_+ y_+, s.t.: y_- <= c_-, y_+ <= c_+
  // if primal constraint is <= -> create just minus var
  // if primal constraint is >= -> create just plus var
  // WARNING: they are inactive at the beginning
  //
  void initAllStabVariables();

// STAB
// Update the stabilization variables based on the dual solution
// 1- When the dual lays inside the box:
//     - increase the penalty of the duals (the bound for the primal)
//     - decrease the radius of the duals (the cost for the primal).
// 2- When the dual lays outside the box:
//     - decrease the penalty of the duals (the bound for the primal)
//     - increase the radius of the duals (the cost for the primal).
// When a dual solution (of the original problem) of better quality
// is obtained, recenter the box.
// The issue here is that the  dual solution is not available as the lagrangian
// bound needs to be computed (and available) and all sub problems need to
// have been solved to optimality.
// Instead, the solution is recenter when asked (recenter=true).
// Currently, the box is recentered when no more columns are generated.
  void stabUpdate(OsiSolverInterface *solver, bool recenter = true);

  // STAB
  // initialize the stabilization variables and center them on the current duals
  void stabInitializeBoundAndCost(OsiSolverInterface *solver);

  // STAB
  // deactivate the stabilization variables
  void stabDeactivateBoundAndCost(OsiSolverInterface *solver);

  // STAB
  // Check the stopping criterion of the relaxation solution specific to the
  // the stabilization
  // The point is that current solution can be infeasible if  stabilization
  // variables are non zero
  bool stabCheckStoppingCriterion() const;

  // STAB
  // return the current cost of the stabilization variables
  double getStabCost() const;

  const SolverParam & param() const;

  double epsilon() const;

 private:
  Modeler *pModel_;
  std::vector<MyVar *> stabVariablesPlus_, stabVariablesMinus_;
  // constraints associated with each two stab variables
  std::vector<MyCons *> stabConstraints_;
  std::vector<double> stabBoxCenters_;

 protected:
  // STAB
  // Multiply the upper bound of the input variable by the input factor
  void multiplyUbInSolver(MyVar *pVar,
                          OsiSolverInterface *solver,
                          double factor);
  // Set the bound of the input variable to the input value
  void updateVarUbInSolver(MyVar *pVar,
                           OsiSolverInterface *solver,
                           double value);

  // STAB
  // Set the cost of the input variable to the input value
  void updateVarCostInSolver(MyVar *pVar,
                             OsiSolverInterface *solver,
                             double value);

 private:
  void addStabVariables(const char *name,
                        MyCons *cons,
                        bool LECons,
                        bool GECons);
};

#endif  // SRC_SOLVERS_MP_MODELER_STABILIZATION_H_
