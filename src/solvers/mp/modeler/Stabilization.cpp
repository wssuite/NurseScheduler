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

#include "Stabilization.h"

#include <string>

#include "solvers/mp/modeler/Modeler.h"


Stabilization::Stabilization(Modeler *pModel) : pModel_(pModel) {}

// STAB: Add stabilization variable
void Stabilization::initAllStabVariables() {
  for (MyCons *con : pModel_->getCoreCons()) {
    std::string name = "stab_" + std::string(con->name_);
    addStabVariables(name.c_str(), con,
                     con->getRhs() < pModel_->getInfinity(),
                     con->getLhs() > -pModel_->getInfinity());
  }
}

// Add stabilization variables z for the box [b_, b+] with the penalties c
// if getting outside of the box:
// dual = obj += - c_+ z_+ - c_- z_-, s.t.: b_- - z_-<= Pi <= b_+ + z_+
// primal = obj += - b_- y_- + b_+ y_+, s.t.: y_- <= c_-, y_+ <= c_+
// if primal constraint is <= -> dual <= 0 -> just need the LB of the box
//                            -> create just minus var
// if primal constraint is >= -> dual >= 0 -> just need the UB of the box
//                            -> create just plus var
// WARNING: they are inactive at the beginning
//
void Stabilization::addStabVariables(
    const char *name,
    MyCons *cons,
    bool LECons,
    bool GECons) {
  MyVar *var;
  char n[255];

  // The lower side of the box
  if (LECons) {
    snprintf(n, sizeof(n), "%s_minus", name);
    pModel_->createPositiveVar(&var, n, LARGE_SCORE, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, -1);
    stabVariablesMinus_.push_back(var);
  } else {
    stabVariablesMinus_.push_back(nullptr);
  }

  // The  upper side of the box
  if (GECons) {
    snprintf(n, sizeof(n), "%s_plus", name);
    pModel_->createPositiveVar(&var, n, LARGE_SCORE, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, 1);
    stabVariablesPlus_.push_back(var);
  } else {
    stabVariablesPlus_.push_back(nullptr);
  }

  stabConstraints_.push_back(cons);
  stabBoxCenters_.push_back(0);
}

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
void Stabilization::stabUpdate(OsiSolverInterface *solver, bool recenter) {
  // stabilization variables corresponding to the cover constraints
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    double center = stabBoxCenters_[i],
        radius = varPlus ? varPlus->getCost() - center :
                 center + varMinus->getCost();
    double dual = pModel_->getDual(stabConstraints_[i], true);
    double penaltyFactor = 1;

    // if dual within the box, decrease radius and increase cost
    if (dual > center + radius + epsilon() &&
        dual < center - radius - epsilon()) {
      // decrease radius
      if (param().isStabUpdateBoxRadius_)
        radius /= param().stabBoxRadiusFactor_;
      // increase penalty
      if (param().isStabUpdatePenalty_)
        penaltyFactor *= param().stabBoxRadiusFactor_;
    } else {
      // increase radius
      if (param().isStabUpdateBoxRadius_)
        radius *= param().stabPenaltyFactor_;
      // decrease penalty
      if (param().isStabUpdatePenalty_)
        penaltyFactor /= param().stabPenaltyFactor_;
    }

    if (radius > param().stabBoxRadiusMax_)
      radius = param().stabBoxRadiusMax_;

    // recenter box
    if (recenter) {
      center = dual;
      stabBoxCenters_[i] = dual;
    }

    // update box
    if (varPlus) updateVarCostInSolver(varPlus, solver, center + radius);
    if (varMinus) updateVarCostInSolver(varMinus, solver, -center + radius);
    // update penalty
    if (param().isStabUpdatePenalty_) {
      if (varPlus) multiplyUbInSolver(varPlus, solver, penaltyFactor);
      if (varMinus) multiplyUbInSolver(varMinus, solver, penaltyFactor);
    }
  }
}

// STAB
// activate the stabilization variables and center them on the current duals
void Stabilization::stabInitializeBoundAndCost(OsiSolverInterface *solver) {
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    double center = pModel_->getDual(stabConstraints_[i], true);
    stabBoxCenters_[i] = center;
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver,
                            center + param().stabBoxRadiusIni_);
      updateVarUbInSolver(varPlus, solver, param().stabPenaltyIni_);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver,
                            -center + param().stabBoxRadiusIni_);
      updateVarUbInSolver(varMinus, solver, param().stabPenaltyIni_);
    }
  }
}

// STAB
// deactivate the stabilization variables
void Stabilization::stabDeactivateBoundAndCost(OsiSolverInterface *solver) {
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    stabBoxCenters_[i] = 0;
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver, LARGE_SCORE);
      updateVarUbInSolver(varPlus, solver, 0);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver, LARGE_SCORE);
      updateVarUbInSolver(varMinus, solver, 0);
    }
  }
}

// STAB
// Stop when the stabilization variables are all null
bool Stabilization::stabCheckStoppingCriterion() const {
  if (!pModel_->getParameters().isStabilization_)
    return true;

  for (MyVar *var : stabVariablesPlus_)
    if (var && pModel_->getVarValue(var) > epsilon())
      return false;
  for (MyVar *var : stabVariablesMinus_)
    if (var && pModel_->getVarValue(var) > epsilon())
      return false;

  return true;
}

// STAB
// return the current cost of the stabilization variables
double Stabilization::getStabCost() const {
  if (!pModel_->getParameters().isStabilization_)
    return 0;
  return pModel_->getTotalCost(stabVariablesPlus_) +
      pModel_->getTotalCost(stabVariablesMinus_);
}

// STAB
// Multiply the upper bound of the input variable by the input factor
void Stabilization::multiplyUbInSolver(MyVar *pVar,
                                       OsiSolverInterface *solver,
                                       double factor) {
  int varind = pVar->getIndex();
  double ub = pVar->getUB();

  if (ub != solver->getColUpper()[varind]) {
    Tools::throwError("multiplyUbInSolver: the upper bound stored in the "
                      "variable is not the same as that in the solver!");
  }

  ub *= factor;
  if (ub > param().stabBoxBoundMax_) ub = param().stabBoxBoundMax_;

  solver->setColUpper(varind, ub);
  pVar->setUB(ub);
}

// STAB
// Set the bound of the input variable to the input value
void Stabilization::updateVarUbInSolver(MyVar *pVar,
                                        OsiSolverInterface *solver,
                                        double value) {
  int varind = pVar->getIndex();
  double ub = pVar->getUB();

  if (ub != solver->getColUpper()[varind]) {
    Tools::throwError("updateVarUbInSolver: the upper bound stored in the "
                      "variable is not the same as that in the solver!");
  }

  if (value > param().stabBoxBoundMax_) value = param().stabBoxBoundMax_;

  solver->setColUpper(varind, value);
  pVar->setUB(value);
}

// STAB
// Set the cost of the input variable to the input value
void Stabilization::updateVarCostInSolver(MyVar *pVar,
                                          OsiSolverInterface *solver,
                                          double value) {
  int varind = pVar->getIndex();
  double cost = pVar->getCost();

  if (cost != solver->getObjCoefficients()[varind]) {
    Tools::throwError("updateVarCostInSolver: the cost stored in the variable "
                      "is not the same as that in the solver!");
  }

  if (value > param().stabPenaltyMax_) value = param().stabPenaltyMax_;
  else if (-value > param().stabPenaltyMax_) value = -param().stabPenaltyMax_;

  solver->setObjCoeff(varind, value);
  pVar->setCost(value);
}

const SolverParam & Stabilization::param() const {
  return pModel_->getParameters();
}

double Stabilization::epsilon() const { return param().epsilon_; }
