/*
 * ScipModeler.h:
 *    Tools to create a SCIP problem
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#ifndef SRC_SCIPMODELER_H_
#define SRC_SCIPMODELER_H_

/* standard library includes */
#include <cfloat>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "objscip/objpricer.h"

/* namespace usage */
using namespace std;
using namespace scip;

class ScipModeler {
public:
   ScipModeler(const char* name){
      initializeSCIP(name);
   };

   ~ScipModeler(){
      deleteSCIP();
   };

   int initializeSCIP(const char* name){
      /* initialize SCIP environment */
      SCIP_CALL( SCIPcreate(&scip_) );

      /* include default plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(scip_) );

      /* create empty problem */
      SCIP_CALL( SCIPcreateProb(scip_, name, 0, 0, 0, 0, 0, 0, 0) );
   }

   int deleteSCIP(){
      SCIP_CALL( SCIPfree(&scip_) );
      BMScheckEmptyMemory();
   }

   //solve the model
   int solve(bool relaxation = false){
      if(relaxation)
         SCIP_CALL( SCIPsolve(scip_) );
      else
         SCIP_CALL( SCIPsolve(scip_) );
   }

   //Add a pricer
   int addObjPricer(ObjPricer* pricer, const char* name){
      /* include the pricer */
      SCIP_CALL( SCIPincludeObjPricer(scip_, pricer, true) );
      /* activate the pricer */
      SCIP_CALL( SCIPactivatePricer(scip_, SCIPfindPricer(scip_, name)) );
   }

   /*
    * Add variables
    */

   int createVar(SCIP_VAR* var, const char* var_name, double lhs, double rhs, double objCoeff, SCIP_VARTYPE vartype){
      if(rhs==DBL_MAX)
         SCIP_CALL( SCIPcreateVar(scip_, &var, var_name, lhs, SCIPinfinity(scip_), objCoeff, vartype,
            true, false, 0, 0, 0, 0, 0) );
      else
         SCIP_CALL( SCIPcreateVar(scip_, &var, var_name, lhs, rhs, objCoeff, vartype,
            true, false, 0, 0, 0, 0, 0) );
   }

   void createPositiveVar(SCIP_VAR* var, const char* var_name, double objCoeff, double rhs = DBL_MAX){
         createVar(var, var_name, 0.0, rhs, objCoeff, SCIP_VARTYPE_CONTINUOUS);
   }

   void createIntVar(SCIP_VAR* var, const char* var_name, double objCoeff, double rhs = DBL_MAX){
      createVar(var, var_name, 0, rhs, objCoeff, SCIP_VARTYPE_INTEGER);
   }

   void createBoolVar(SCIP_VAR* var, const char* var_name, double objCoeff, double rhs = DBL_MAX){
      createVar(var, var_name, 0.0, SCIPinfinity(scip_), objCoeff, SCIP_VARTYPE_BINARY);
   }

   /*
    * Add linear constraints:
    *    scip is a pointer to the problem
    *    con is a pointer to the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
    *    coeffs is the array of coefficient to add to the constraints
    */

   int createConsLinear(SCIP_CONS* con, const char* con_name, double lhs, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      SCIP_CALL( SCIPcreateConsLinear(scip_, &con, con_name, nonZeroVars, vars, coeffs, lhs, rhs,
         true, false, true, true, true, false, true, false, false, false) );
   }

   //Add a lower or equal constraint
   void createLEConsLinear(SCIP_CONS* con, const char* con_name, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, nonZeroVars, vars, coeffs);
   }

   //Add a greater or equal constraint
   void createGEConsLinear(SCIP_CONS* con, const char* con_name, double lhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createConsLinear(con, con_name, lhs, SCIPinfinity(scip_), nonZeroVars, vars, coeffs);
   }

   //Add an equality constraint
   void createEQConsLinear(SCIP_CONS* con, const char* con_name, double eq,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createConsLinear(con, con_name, eq, eq, nonZeroVars, vars, coeffs);
   }

   //Add final linear constraints
   int createFinalConsLinear(SCIP_CONS* con, const char* con_name, double lhs, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      SCIP_CALL( SCIPcreateConsLinear(scip_, &con, con_name, nonZeroVars, vars, coeffs, lhs, rhs,
         true, false, true, true, true, false, false, false, false, false) );
   }

   void createFinalLEConsLinear(SCIP_CONS* con, const char* con_name, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createFinalConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, nonZeroVars, vars, coeffs);
   }

   void createFinalGEConsLinear(SCIP_CONS* con, const char* con_name, double lhs,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createFinalConsLinear(con, con_name, lhs, SCIPinfinity(scip_), nonZeroVars, vars, coeffs);
   }

   void createFinalEQConsLinear(SCIP_CONS* con, const char* con_name, double eq,
      int nonZeroVars = 0, SCIP_VAR** vars = NULL, double* coeffs = NULL){
      createFinalConsLinear(con, con_name, eq, eq, nonZeroVars, vars, coeffs);
   }

   /**************
    * Statistics *
    *************/
   int printStats(){
      SCIP_CALL( SCIPprintStatistics(scip_, NULL) );
      SCIP_CALL( SCIPprintBestSol(scip_, NULL, FALSE) );
   }

private:
   //SCIP pointer
   SCIP* scip_;
};


#endif /* SRC_SCIPMODELER_H_ */
