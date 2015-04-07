/*
 * CbcModeler.h
 *
 *  Created on: April 7, 2015
 *      Author: jeremy omer
 */

#ifndef SRC_CBCMODELER_H_
#define SRC_CBCMODELER_H_

#include "CoinModeler.h"
#include "CbcModel.hpp"
#include "MyTools.h"

// Using CLP as the solver
#include "OsiClpSolverInterface.hpp"

class CbcModeler: public CoinModeler {
public:
 CbcModeler(const char* name):
    CoinModeler(), primalValues_(0), objVal_(0)
{ }
 ~CbcModeler() { }

 //solve the model
 int solve(bool relaxation = false);

 //Add a pricer
 int addObjPricer(MyPricer* pPricer){
    Tools::throwError("There is no pricer if Cbc is used to solve the problem!");
 }

 /*
  * Create variable:
  *    var is a pointer to the pointer of the variable
  *    var_name is the name of the variable
  *    lhs, rhs are the lower and upper bound of the variable
  *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
  */
 int createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub);

 int createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub) { return 1;}

 /*
  * Create linear constraint:
  *    con is a pointer to the pointer of the constraint
  *    con_name is the name of the constraint
  *    lhs, rhs are the lower and upper bound of the constraint
  *    nonZeroVars is the number of non-zero coefficients to add to the constraint
  */

 int createCoinConsLinear(CoinCons** con, const char* con_name, int index, double lhs, double rhs);

 /*
  * Create the Clp solver and assign the result the Cbc model
  * The method can only be called after the creation of all the variables and
  * linear constraints
  *
 */
 int setModel();

 /*
  * Get the primal value
  */

 double getVarValue(MyObject* var);

 /*
  * Get the dual variables
  */

 double getDual(MyObject* cons, bool transformed = false) {
   Tools::throwError("There is no dual solution if the problem is solved with integrality constraints!");
 }

 /**************
  * Parameters *
  *************/
 int setVerbosity(int v) {model_.setLogLevel(v);}

 /**************
  * Outputs *
  *************/

 int printStats();

 int printBestSol();

 int writeProblem(string fileName);

 int writeLP(string fileName);

/************************************
* Class own methods and parameters *
************************************/

  void setSolution();


protected:
  // Cbc model
  //
  CbcModel model_;

  //results
  //
  double objVal_;
  double *primalValues_;
};

#endif
