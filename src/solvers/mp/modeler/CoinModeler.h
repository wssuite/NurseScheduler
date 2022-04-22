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

#ifndef SRC_SOLVERS_MP_MODELER_COINMODELER_H_
#define SRC_SOLVERS_MP_MODELER_COINMODELER_H_

#include <map>
#include <vector>
#include <string>

#include "solvers/mp/modeler/Modeler.h"

/* Coin includes */
#include "CoinPackedMatrix.hpp"  // NOLINT (suppress cpplint error)
#include "OsiClpSolverInterface.hpp"  // NOLINT (suppress cpplint error)

#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_CBC
#include "OsiCbcSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

/*
 * My Constraints
 */
// Coin cons, just a virtual class
struct CoinCons : public MyCons {
 public:
  CoinCons(const char *name, int index, double lhs, double rhs) :
      MyCons(name, index, lhs, rhs) {}

  explicit CoinCons(const CoinCons &cons) :
      MyCons(cons) {}

  virtual ~CoinCons() {}
};

/*
 * My Variables
 */
// Coin var, just a virtual class
struct CoinVar : public MyVar {
  CoinVar(const char *name,
          int index,
          double cost,
          VarType type,
          double lb,
          double ub,
          const std::vector<double> &pattern = {},
          double dualCost = 99999,
          const std::vector<int> &indexRows = {},
          const std::vector<double> &coeffs = {}) :
      MyVar(name, index, cost, type, lb, ub, pattern),
      dualCost_(dualCost),
      indexRows_(indexRows),
      coeffs_(coeffs) {}

  explicit CoinVar(const CoinVar &var) :
      MyVar(var),
      dualCost_(var.dualCost_),
      indexRows_(var.indexRows_),
      coeffs_(var.coeffs_) {}

  virtual ~CoinVar() {}

  /*
  * Setters and Getters
  */

  void addRow(int index, double coeff) {
    indexRows_.push_back(index);
    coeffs_.push_back(coeff);
  }

  void toString(const std::vector<MyCons *> &cons) const {
    std::cout << name_ << ":";
    for (unsigned int i = 0; i < indexRows_.size(); ++i)
      std::cout << " " << cons[indexRows_[i]]->name_ << ":" << coeffs_[i];
    std::cout << std::endl;
  }

  int getNbRows() const { return indexRows_.size(); }

  const std::vector<int> &getIndexRows() const { return indexRows_; }

  int getIndexRow(int i) const { return indexRows_[i]; }

  const std::vector<double> &getCoeffRows() const { return coeffs_; }

  double getCoeffRow(int i) const { return coeffs_[i]; }

 protected:
  double dualCost_;  // dualCost of the variable
  // index of the rows of the matrix where the variable has non-zero coefficient
  std::vector<int> indexRows_;
  std::vector<double> coeffs_;  // value of these coefficients
};

static const std::map<VarType, char> typesToCoin = {
    {VARTYPE_CONTINUOUS, '\0'},
    {VARTYPE_BINARY, '\1'},
    {VARTYPE_INTEGER, '\2'}
};

class CoinModeler : public Modeler {
 public:
  CoinModeler() :
      Modeler() {}
  virtual ~CoinModeler() {}

  /*
   * Create linear constraint:
   *    con is a pointer to the pointer of the constraint
   *    con_name is the name of the constraint
   *    lhs, rhs are the lower and upper bound of the constraint
   *    nonZeroVars is the number of non-zero coefficients to add
   *    to the constraint
   *    vars is an array of pointers to the variables to add to the constraints
   *    (with non-zero coefficient)
   *    coeffs is the array of coefficient to add to the constraints
   */
  // has to be implement to create the good cons according to the coin modeler
  // chosen (BCP, CBC ...)
  // They do not store the coefficients of the variables
  // (the variables store their own coefficients)
  virtual int createCoinConsLinear(MyCons **con,
                                   const char *con_name,
                                   int index,
                                   double lhs,
                                   double rhs) = 0;

  virtual int createConsLinear(MyCons **con,
                               const char *con_name,
                               int index,
                               double lhs,
                               double rhs,
                               std::vector<MyVar *> vars = {},
                               std::vector<double> coeffs = {}) {
    createCoinConsLinear(con, con_name, index, lhs, rhs);
    for (unsigned int i = 0; i < vars.size(); ++i)
      addCoefLinear(*con, vars[i], coeffs[i]);
    return 1;
  }

  // the cut stores the coefficients and indices of the variables contained
  // within the cut
  // Indeed, the cuts can be removed and it's much easier this way
  virtual int createCoinCutLinear(MyCons **con,
                                  const char *con_name,
                                  int index,
                                  double lhs,
                                  double rhs,
                                  const std::vector<int> &indexVars,
                                  const std::vector<double> &coeffs) = 0;

  int createCutLinear(MyCons **con,
                      const char *con_name,
                      int index,
                      double lhs,
                      double rhs,
                      std::vector<MyVar *> vars,
                      std::vector<double> coeffs) {
    std::vector<int> indexes;
    for (MyVar *var : vars)
      indexes.push_back(var->getIndex());
    return createCoinCutLinear(con, con_name, index, lhs, rhs, indexes, coeffs);
  }

 public:
  /*
  * Add variables to constraints
  */
  int addCoefLinear(MyCons *cons,
                    MyVar *var,
                    double coeff,
                    bool transformed = false) {
    CoinVar *var2 = dynamic_cast<CoinVar *>(var);
    var2->addRow(cons->getIndex(), coeff);
    return 1;
  }

  /*
  * Build the CoinedPackMatrix corresponding to the problem
  * The indices need to be a sequence from 0 to N without any hole.
  */

  /* build the problem */
  CoinPackedMatrix * buildCoinMatrix() {
    // define nb rows and col
    const int corenum = coreVars_.size();
    // define a matrix as a vector of tuples (row_index, col_index, coeff)
    std::vector<int> row_indices, col_indices;
    std::vector<double> coeffs;
    int index_max = 0;
    for (int i = 0; i < corenum; ++i) {
      CoinVar *var = dynamic_cast<CoinVar *>(coreVars_[i]);

      // build the tuples of the matrix
      int col_index = var->getIndex();
      if (col_index > index_max) index_max = col_index;
      for (int j = 0; j < var->getNbRows(); ++j) {
        row_indices.push_back(var->getIndexRow(j));
        col_indices.push_back(col_index);
        coeffs.push_back(var->getCoeffRow(j));
      }
    }

    if (index_max != (corenum - 1)) {
      std::cerr << "There are " << corenum
                << " core variables, but the maximum index is " << index_max
                << std::endl;
      Tools::throwError(
          "The indices of the variables are not forming a sequence 0..N");
    }

    // initialize the matrix
    auto matrix = new CoinPackedMatrix(false,
                                       &(row_indices[0]),
                                       &(col_indices[0]),
                                       &(coeffs[0]),
                                       row_indices.size());

    return matrix;
  }

  /*
  * Get the primal value
  */
  virtual double getVarValue(MyVar *var) const { return 0; }

  /*
  * Get the dual variables
  */
  virtual double getDual(MyCons *cons,
                         bool transformed = false) const { return 0; }

  /*
  * Get the reduced cost
  */
  virtual double getReducedCost(MyVar *var) const { return 0; }

  /**************
  * Parameters *
  *************/
  virtual int setVerbosity(int v) { return 0; }

  /**************
  * Outputs *
  *************/

  virtual int nbSolutions() const { return 0; }

  virtual double getObjective() const { return pTree_->getBestUB(); }
  virtual double getObjective(int index) const { return LARGE_SCORE; }

  virtual int writeProblem(std::string fileName) const { return 0; }

  virtual int writeLP(std::string fileName) const { return 0; }

  virtual void toString(MyObject *obj) const {
    CoinVar *var = dynamic_cast<CoinVar *>(obj);
    if (var) var->toString(coreCons_);
    else
      Modeler::toString(obj);
  }
};

static OsiSolverInterface * getNewSolver(SolverType type) {
  switch (type) {
    case CLP: return new OsiClpSolverInterface();
    case Gurobi:
#ifdef USE_GUROBI
      return new OsiGrbSolverInterface();
#endif
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
    default: break;
  }
  return nullptr;
}

#endif  // SRC_SOLVERS_MP_MODELER_COINMODELER_H_

