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
void Stabilization::initStabVariables(const std::vector<MyCons *> & cons) {
  for (MyCons *con : cons) {
    if (con == nullptr) continue;
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
    pModel_->createIntVar(&var, n, param().stabPenaltyIni_, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, -1);
    stabVariablesMinus_.push_back(var);
  } else {
    stabVariablesMinus_.push_back(nullptr);
  }

  // The  upper side of the box
  if (GECons) {
    snprintf(n, sizeof(n), "%s_plus", name);
    pModel_->createIntVar(&var, n, param().stabPenaltyIni_, {}, 0, 0);
    pModel_->addCoefLinear(cons, var, 1);
    stabVariablesPlus_.push_back(var);
  } else {
    stabVariablesPlus_.push_back(nullptr);
  }

  stabConstraints_.push_back(cons);
  stabBoxCenters_.push_back(0);
}

// Update the stabilization variables based on the dual solution
// 1- When all dual lays inside the box:
//     - decrease the radius of the duals (the cost for the primal).
// 2- When the dual lays outside the box and column generation has ended:
//     - decrease the penalty of the duals (the bound for the primal)
// The issue here is that the dual solution is not available as the lagrangian
// bound needs to be computed (and available) and all sub problems need to
// have been solved to optimality.
// The solution is recenter at each iteration.
void Stabilization::stabUpdate(
    OsiSolverInterface *solver, bool columnGenerated) {
  // stabilization variables corresponding to the cover constraints
  std::vector<double> duals = pModel_->getDuals(stabConstraints_, true);
  double penaltyFactor = 1;
  // check if all duals within the box -> decrease box radius
  if (stabCheckStoppingCriterion()) {
    if (param().isStabUpdateBoxRadius_)
      stabRadius_ /= param().stabBoxRadiusFactor_;
  } else if (!columnGenerated) {
    // otherwise, decrease penalties if column generation has ended
    penaltyFactor /= param().stabPenaltyFactor_;
  }

  // update centers, radius and penalties
  stabBoxCenters_ = duals;
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    // update box
    double center = stabBoxCenters_[i];
    MyVar *varPlus = stabVariablesPlus_[i],
            *varMinus = stabVariablesMinus_[i];
    if (varPlus) updateVarCostInSolver(varPlus, solver, center + stabRadius_);
    if (varMinus)
      updateVarCostInSolver(varMinus, solver, -center + stabRadius_);
    // update penalty
    if (penaltyFactor < 1 - epsilon()) {
      if (varPlus) multiplyUbInSolver(varPlus, solver, penaltyFactor);
      if (varMinus) multiplyUbInSolver(varMinus, solver, penaltyFactor);
    }
  }
}

// compute the difference between the current dual variables and
// the previous ones
double Stabilization::computeStabAvgDifference(
    const std::vector<double> &duals) {
  double stabAvgDiff = 0;
  for (int i = 0; i < stabConstraints_.size(); ++i)
    stabAvgDiff += std::abs(duals[i] - stabBoxCenters_[i]);
  stabAvgDiff /= stabConstraints_.size();
//  std::cout << "Avg stab diff: " << stabAvgDiff << std::endl;
  return stabAvgDiff;
}


// STAB
// activate the stabilization variables and center them on the current duals
void Stabilization::stabInitializeBoundAndCost(OsiSolverInterface *solver) {
  // check if ub is lower than 1, so no integer solution can be found
  // with active stab variables
  if (param().stabPenaltyIni_ < -epsilon() ||
      param().stabPenaltyIni_ > 1 - epsilon())
    Tools::throwError("Stabilization penalty (%.2f) should be within [0, 1[ "
                      "such that no integer solution contains active "
                      "stabilization variables.", param().stabPenaltyIni_);
  // compute current duals average differences
  std::vector<double> duals = pModel_->getDuals(stabConstraints_, true);
  double avg = computeStabAvgDifference(duals);
  stabRadius_ = avg / param().stabBoxRadiusIniRatio_;
  stabBoxCenters_ = duals;
  for (int i = 0; i < stabConstraints_.size(); ++i) {
    MyVar *varPlus = stabVariablesPlus_[i],
        *varMinus = stabVariablesMinus_[i];
    double center = stabBoxCenters_[i];
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver, center + stabRadius_);
      updateVarUbInSolver(varPlus, solver, param().stabPenaltyIni_);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver, -center + stabRadius_);
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
    if (varPlus) {
      updateVarCostInSolver(varPlus, solver, INFEAS_COST);
      updateVarUbInSolver(varPlus, solver, 0);
    }
    if (varMinus) {
      updateVarCostInSolver(varMinus, solver, INFEAS_COST);
      updateVarUbInSolver(varMinus, solver, 0);
    }
  }
}

// STAB
// Stop when the stabilization variables are all null
bool Stabilization::stabCheckStoppingCriterion() const {
  if (!param().isStabilization_)
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
  if (!param().isStabilization_)
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
    Tools::throwException("multiplyUbInSolver: the upper bound stored in the "
                          "variable is not the same as that in the solver!");
  }

  ub *= factor;
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
    Tools::throwException("updateVarUbInSolver: the upper bound stored in the "
                          "variable is not the same as that in the solver!");
  }

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
    Tools::throwException(
        "updateVarCostInSolver: the cost stored in the variable "
        "is not the same as that in the solver!");
  }

  solver->setObjCoeff(varind, value);
  pVar->setCost(value);
}

const SolverParam & Stabilization::param() const {
  return pModel_->getParameters();
}

// epsilon needs to be like a zero here
double Stabilization::epsilon() const { return 1e-9; }
