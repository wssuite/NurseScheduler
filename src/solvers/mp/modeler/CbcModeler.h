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

#ifndef SRC_SOLVERS_MP_MODELER_CBCMODELER_H_
#define SRC_SOLVERS_MP_MODELER_CBCMODELER_H_

#include <map>
#include <vector>
#include <string>

#include "solvers/mp/modeler/CoinModeler.h"
#include "tools/Tools.h"

// Using CLP as the solver
#include "CbcModel.hpp"  // NOLINT (suppress cpplint error)
#include "OsiClpSolverInterface.hpp"  // NOLINT (suppress cpplint error)

class CbcModeler : public CoinModeler {
 public:
  // default constructor
  CbcModeler() : CoinModeler(), primalValues_(0), objVal_(0), model_(NULL) {}

  explicit CbcModeler(const char *name) :
      CoinModeler(), primalValues_(0), objVal_(0), model_(NULL) {}

  // useful constructor
  CbcModeler(const vector<CoinVar *> &coreVars,
             const vector<CoinVar *> &columnVars,
             const vector<CoinCons *> &cons);

  CbcModeler(const vector<CoinVar *> &coreVars,
             const vector<CoinVar *> &columnVars,
             const vector<CoinCons *> &cons,
             OsiSolverInterface *osiClp);

  ~CbcModeler() {
    if (model_ != NULL) delete model_;
  }

  // initialize the vectors columnVars, coreVars and cons
  void initializeVectors(const vector<CoinVar *> &coreVars,
                         const vector<CoinVar *> &columnVars,
                         const vector<CoinCons *> &cons);

  // solve the model
  int solve(bool relaxation = false);

  // Add a pricer
  int addObjPricer(MyPricer *pPricer) {
    Tools::throwError("There is no pricer "
                      "if Cbc is used to solve the problem!");
    return 0;
  }

  /*
   * Create variable:
   *    var is a pointer to the pointer of the variable
   *    var_name is the name of the variable
   *    lhs, rhs are the lower and upper bound of the variable
   *    vartype is the type of the variable:
   *    VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
   */
  virtual int createVar(MyVar **var,
                        const char *var_name,
                        int index,
                        double objCoeff,
                        double lb,
                        double ub,
                        VarType vartype,
                        const std::vector<double> &pattern,
                        double score);

  virtual int createColumnVar(MyVar **var,
                              const char *var_name,
                              int index,
                              double objCoeff,
                              const std::vector<double> &pattern,
                              double dualObj,
                              double lb,
                              double ub,
                              VarType vartype,
                              double score);

  /*
   * Create linear constraint:
   *    con is a pointer to the pointer of the constraint
   *    con_name is the name of the constraint
   *    lhs, rhs are the lower and upper bound of the constraint
   *    nonZeroVars is the number of non-zero coefficients to add to the constraint
   */
  virtual int createCoinConsLinear(CoinCons **con,
                                   const char *con_name,
                                   int index,
                                   double lhs,
                                   double rhs);

  /*
   * Create the Clp solver and assign the result to the Cbc model
   * The method can only be called after the creation of all the variables and
   * linear constraints
   *
  */
  int setModel();

//  /*
//   * Set a high priority on all the variable returned by the branching rule
//   */
//  void setBranchingRule();

  /*
   * Get the primal value
   */
  virtual double getVarValue(MyVar *var);

  /*
   * Get the dual variables
   */
  double getDual(MyObject *cons, bool transformed = false) {
    Tools::throwError("There is no dual solution if the problem "
                      "is solved with integrality constraints!");
    return 0;
  }

  double getObjective() { return objVal_; }

  /**************
   * Parameters *
   *************/
  int setVerbosity(int v) { verbosity_ = v; return 1; }

  /**************
   * Outputs *
   *************/

  virtual int printStats();

  virtual int printBestSol();

  virtual int writeProblem(string fileName);

  virtual int writeLP(string fileName);

/************************************
* Class own methods and parameters *
************************************/
  void setSolution();

 protected:
  // Cbc model
  CbcModel *model_;

  // OsiClpSolverInterface if available
  OsiSolverInterface *pOsiSolver_;

  // results
  double objVal_;
  double *primalValues_;
};

#endif  // SRC_SOLVERS_MP_MODELER_CBCMODELER_H_
