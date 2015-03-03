/*
 * ScipModeler.h:
 *    Tools to create a SCIP problem.
 *    In general, when you want to create a new scip object,
 *    you have to give a pointer to the pointer of the object.
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
   }

   ~ScipModeler(){ }

   int initializeSCIP(const char* name){
      /* initialize SCIP environment */
      SCIP_CALL( SCIPcreate(&scip_) );

      /* include default plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(scip_) );

      /* create empty problem */
      SCIP_CALL( SCIPcreateProb(scip_, name, 0, 0, 0, 0, 0, 0, 0) );
   }

   //delete all the model built by scip
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
   int addObjPricer(ObjPricer* pricer){
      /* include the pricer */
      SCIP_CALL( SCIPincludeObjPricer(scip_, pricer, true) );
      /* activate the pricer */
      SCIP_CALL( SCIPactivatePricer(scip_, SCIPfindPricer(scip_, pricer->scip_name_)) );
   }

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
    */

   int createVar(SCIP_VAR** var, const char* var_name, double objCoeff, double lhs, double rhs, SCIP_VARTYPE vartype, double score){
      if(rhs==DBL_MAX)
         SCIP_CALL( SCIPcreateVar(scip_, var, var_name, lhs, SCIPinfinity(scip_), objCoeff, vartype,
            true, false, 0, 0, 0, 0, 0) );
      else
         SCIP_CALL( SCIPcreateVar(scip_, var, var_name, lhs, rhs, objCoeff, vartype,
            true, false, 0, 0, 0, 0, 0) );

      if(score > 0)
         SCIP_CALL( SCIPaddPricedVar(scip_, *var, score) );
      else SCIP_CALL( SCIPaddVar(scip_, *var) );
   }

   void createPositiveVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0, double rhs = DBL_MAX){
      createVar(var, var_name, objCoeff, 0.0, rhs, SCIP_VARTYPE_CONTINUOUS, score);
   }

   void createIntVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0, double rhs = DBL_MAX){
      createVar(var, var_name, objCoeff, 0, rhs, SCIP_VARTYPE_INTEGER, score);
   }

   void createBinaryVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0){
      createVar(var, var_name, objCoeff, 0.0, 1.0, SCIP_VARTYPE_BINARY, score);
   }

   /*
    * Create linear constraint:
    *    con is a pointer to the pointer of the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
    *    coeffs is the array of coefficient to add to the constraints
    */

   int createConsLinear(SCIP_CONS** con, const char* con_name, double lhs, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      SCIP_CALL( SCIPcreateConsLinear(scip_, con, con_name, nonZeroVars, vars, coeffs, lhs, rhs,
         true, false, true, true, true, false, true, false, false, false) );
      SCIP_CALL( SCIPaddCons(scip_, *con) );
   }

   //Add a lower or equal constraint
   void createLEConsLinear(SCIP_CONS** con, const char* con_name, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, nonZeroVars, vars, coeffs);
   }

   //Add a greater or equal constraint
   void createGEConsLinear(SCIP_CONS** con, const char* con_name, double lhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createConsLinear(con, con_name, lhs, SCIPinfinity(scip_), nonZeroVars, vars, coeffs);
   }

   //Add an equality constraint
   void createEQConsLinear(SCIP_CONS** con, const char* con_name, double eq,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createConsLinear(con, con_name, eq, eq, nonZeroVars, vars, coeffs);
   }

   //Add final linear constraints
   int createFinalConsLinear(SCIP_CONS** con, const char* con_name, double lhs, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      SCIP_CALL( SCIPcreateConsLinear(scip_, con, con_name, nonZeroVars, vars, coeffs, lhs, rhs,
         true, false, true, true, true, false, false, false, false, false) );
      SCIP_CALL( SCIPaddCons(scip_, *con) );
   }

   void createFinalLEConsLinear(SCIP_CONS** con, const char* con_name, double rhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createFinalConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, nonZeroVars, vars, coeffs);
   }

   void createFinalGEConsLinear(SCIP_CONS** con, const char* con_name, double lhs,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createFinalConsLinear(con, con_name, lhs, SCIPinfinity(scip_), nonZeroVars, vars, coeffs);
   }

   void createFinalEQConsLinear(SCIP_CONS** con, const char* con_name, double eq,
      int nonZeroVars = 0, SCIP_VAR** vars = {}, double* coeffs = {}){
      createFinalConsLinear(con, con_name, eq, eq, nonZeroVars, vars, coeffs);
   }

   /*
    * Add variables to constraints
    */

   int addCoefLinear(SCIP_CONS* cons, SCIP_VAR* var, double coeff){
      SCIP_CALL( SCIPaddCoefLinear(scip_, cons, var, coeff) );
   }

   /*
    * Add new Column to the SCIP problem
    */

   void createColumn(SCIP_VAR** var, const char* var_name, double objCoeff, SCIP_VARTYPE vartype,
      int nonZeroCons = 0, SCIP_CONS** cons = {}, double* coeffs = {}, bool transformed = false, double score = 0){
      switch(vartype){
      case SCIP_VARTYPE_BINARY:
         createBinaryVar(var, var_name, objCoeff, score);
         break;
      case SCIP_VARTYPE_INTEGER:
         createIntVar(var, var_name, objCoeff, score);
         break;
      default:
         createPositiveVar(var, var_name, objCoeff, score);
         break;
      }

      for(int i=0; i<nonZeroCons; i++){
         if(transformed)
            getTransformedCons(cons[i], &(cons[i]));
         addCoefLinear(cons[i], *var, coeffs[i]);
      }
   }

   void createPositiveColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      int nonZeroCons = 0, SCIP_CONS** cons = {}, double* coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_CONTINUOUS, nonZeroCons, cons, coeffs, transformed, score);
   }

   void createBinaryColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      int nonZeroCons = 0, SCIP_CONS** cons = {}, double* coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_BINARY, nonZeroCons, cons, coeffs, transformed, score);
   }

   void createIntColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      int nonZeroCons = 0, SCIP_CONS** cons = {}, double* coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_INTEGER, nonZeroCons, cons, coeffs, transformed, score);
   }

   /*
    * Get the transformed variables and constraints
    *
    *  Because SCIP transforms the original problem in preprocessing, we need to get the references to
    *  the variables and constraints in the transformed problem from the references in the original
    *  problem.
    */

   int getTransformedVar(SCIP_VAR* var, SCIP_VAR** var2){
      SCIP_CALL( SCIPgetTransformedVar(scip_, var, var2) );
   }

   int getTransformedCons(SCIP_CONS* cons, SCIP_CONS** cons2){
      SCIP_CALL( SCIPgetTransformedCons(scip_, cons, cons2) );
   }

   /*
    * get the primal values
    */

   SCIP_SOL* getBestSol(){
      return SCIPgetBestSol(scip_);
   }

   double getVarValue(SCIP_SOL* sol, SCIP_VAR* var){
      SCIPgetSolVal(scip_, sol, var);
   }

   vector<double> getVarValues(SCIP_SOL* sol, vector<SCIP_VAR*> vars){
      vector<double> values(vars.size());
      for(int i=0; i<vars.size(); ++i)
         values[i] = getVarValue(sol, vars[i]);
      return values;
   }

   /*
    * Get the dual variables
    */

   double getDual(SCIP_CONS* cons, bool transformed = false){
      if(transformed)
         getTransformedCons(cons, &cons);
      SCIPgetDualsolLinear(scip_, cons);
   }

   vector<double> getDuals(vector<SCIP_CONS*> cons, bool transformed = false){
      vector<double> dualValues(cons.size());
      for(int i=0; i<cons.size(); ++i)
         dualValues[i] = getDual(cons[i], transformed);
      return dualValues;
   }

   /**************
    * Parameters *
    *************/
   int setVerbosity(int v){
      SCIP_CALL( SCIPsetIntParam(scip_, "display/verblevel", v) );
      /* SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); */
   }

   /**************
    * Outputs *
    *************/
   int printStats(){
      SCIP_CALL( SCIPprintStatistics(scip_, NULL) );
   }

   int printBestSol(){
      SCIP_CALL( SCIPprintBestSol(scip_, NULL, FALSE) );
   }

   int writeProblem(string fileName){
      SCIP_CALL( SCIPwriteOrigProblem(scip_, fileName.c_str(), "lp", FALSE) );
   }

   int writeLP(string fileName){
      SCIP_CALL( SCIPwriteLP(scip_, fileName.c_str()) );
   }

   SCIP* getScip(){
      return scip_;
   }

private:
   //SCIP pointer
   SCIP* scip_;
};


#endif /* SRC_SCIPMODELER_H_ */
