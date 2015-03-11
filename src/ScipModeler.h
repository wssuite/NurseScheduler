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
#include <typeinfo>

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

   ~ScipModeler(){
      deleteSCIP();
   }

   //solve the model
   inline int solve(bool relaxation = false){
      if(relaxation)
         SCIP_CALL( SCIPsolve(scip_) );
      else
         SCIP_CALL( SCIPsolve(scip_) );
   }

   //Add a pricer
   inline int addObjPricer(ObjPricer* pricer){
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

   inline int createVar(SCIP_VAR** var, const char* var_name, double objCoeff, double lhs, double rhs, SCIP_VARTYPE vartype, double score){
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

   inline void createPositiveVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0, double rhs = DBL_MAX){
      createVar(var, var_name, objCoeff, 0.0, rhs, SCIP_VARTYPE_CONTINUOUS, score);
   }

   inline void createIntVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0, double rhs = DBL_MAX){
      createVar(var, var_name, objCoeff, 0, rhs, SCIP_VARTYPE_INTEGER, score);
   }

   inline void createBinaryVar(SCIP_VAR** var, const char* var_name, double objCoeff, double score = 0){
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

   inline int createConsLinear(SCIP_CONS** con, const char* con_name, double lhs, double rhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      SCIP_CALL( SCIPcreateConsLinear(scip_, con, con_name, vars.size(), &(vars)[0], &(coeffs)[0], lhs, rhs,
         true, false, true, true, true, false, true, false, false, false) );
      SCIP_CALL( SCIPaddCons(scip_, *con) );
   }

   //Add a lower or equal constraint
   inline void createLEConsLinear(SCIP_CONS** con, const char* con_name, double rhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, vars, coeffs);
   }

   //Add a greater or equal constraint
   inline void createGEConsLinear(SCIP_CONS** con, const char* con_name, double lhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(con, con_name, lhs, SCIPinfinity(scip_), vars, coeffs);
   }

   //Add an equality constraint
   inline void createEQConsLinear(SCIP_CONS** con, const char* con_name, double eq,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(con, con_name, eq, eq, vars, coeffs);
   }

   //Add final linear constraints
   inline int createFinalConsLinear(SCIP_CONS** con, const char* con_name, double lhs, double rhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      SCIP_CALL( SCIPcreateConsLinear(scip_, con, con_name, vars.size(), &(vars)[0], &(coeffs)[0], lhs, rhs,
         true, false, true, true, true, false, false, false, false, false) );
      SCIP_CALL( SCIPaddCons(scip_, *con) );
   }

   inline void createFinalLEConsLinear(SCIP_CONS** con, const char* con_name, double rhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(con, con_name, -SCIPinfinity(scip_), rhs, vars, coeffs);
   }

   inline void createFinalGEConsLinear(SCIP_CONS** con, const char* con_name, double lhs,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(con, con_name, lhs, SCIPinfinity(scip_), vars, coeffs);
   }

   inline void createFinalEQConsLinear(SCIP_CONS** con, const char* con_name, double eq,
      vector<SCIP_VAR*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(con, con_name, eq, eq, vars, coeffs);
   }

   /*
    * Add variables to constraints
    */

   inline int addCoefLinear(SCIP_CONS* cons, SCIP_VAR* var, double coeff){
      SCIP_CALL( SCIPaddCoefLinear(scip_, cons, var, coeff) );
   }

   /*
    * Add new Column to the SCIP problem
    */

   inline void createColumn(SCIP_VAR** var, const char* var_name, double objCoeff, SCIP_VARTYPE vartype,
      vector<SCIP_CONS*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
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

      for(int i=0; i<cons.size(); i++){
         SCIP_CONS* con = cons[i];
         if(transformed)
            getTransformedCons(con, &con);
         addCoefLinear(con, *var, coeffs[i]);
      }
   }

   inline void createPositiveColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      vector<SCIP_CONS*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_CONTINUOUS, cons, coeffs, transformed, score);
   }

   inline void createBinaryColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      vector<SCIP_CONS*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_BINARY, cons, coeffs, transformed, score);
   }

   inline void createIntColumn(SCIP_VAR** var, const char* var_name, double objCoeff,
      vector<SCIP_CONS*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, SCIP_VARTYPE_INTEGER, cons, coeffs, transformed, score);
   }

   /*
    * Get the transformed variables and constraints
    *
    *  Because SCIP transforms the original problem in preprocessing, we need to get the references to
    *  the variables and constraints in the transformed problem from the references in the original
    *  problem.
    */

   inline int getTransformedVar(SCIP_VAR* var, SCIP_VAR** var2){
      SCIP_CALL( SCIPgetTransformedVar(scip_, var, var2) );
   }

   inline int getTransformedCons(SCIP_CONS* cons, SCIP_CONS** cons2){
      SCIP_CALL( SCIPgetTransformedCons(scip_, cons, cons2) );
   }

   /*
    * get the primal values
    */

   inline SCIP_SOL* getBestSol(){
      return SCIPgetBestSol(scip_);
   }

   inline double getVarValue(SCIP_SOL* sol, SCIP_VAR* var){
      SCIPgetSolVal(scip_, sol, var);
   }

   inline vector<double> getVarValues(SCIP_SOL* sol, vector<SCIP_VAR*> vars){
      vector<double> values(vars.size());
      for(int i=0; i<vars.size(); ++i)
         values[i] = getVarValue(sol, vars[i]);
      return values;
   }

   /*
    * Get the dual variables
    */

   inline double getDual(SCIP_CONS* cons, bool transformed = false){
      if(transformed)
         getTransformedCons(cons, &cons);
      SCIPgetDualsolLinear(scip_, cons);
   }

   inline vector<double> getDuals(vector<SCIP_CONS*> cons, bool transformed = false){
      vector<double> dualValues(cons.size());
      for(int i=0; i<cons.size(); ++i)
         dualValues[i] = getDual(cons[i], transformed);
      return dualValues;
   }

   /**************
    * Parameters *
    *************/
   inline int setVerbosity(int v){
      SCIP_CALL( SCIPsetIntParam(scip_, "display/verblevel", v) );
      /* SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); */
   }

   /**************
    * Outputs *
    *************/

   //compute the total cost of SCIP_VAR* in the solution sol*
   inline double getTotalCost(SCIP_SOL* sol, SCIP_VAR* var){
      double value = getVarValue(sol, var);
      return value *  var->branchfactor;
   }

   //compute the total cost of a vector of SCIP_VAR* in the solution sol*
   template<typename T> inline double getTotalCost(SCIP_SOL* sol, map<SCIP_VAR*, T> map){
      double value = 0 ;
      for(pair<SCIP_VAR*, T> var: map)
         value += getTotalCost(sol, var.first);
      return value;
   }

   //compute the total cost of a multiple vectors of SCIP_VAR* in the solution sol*
   template<typename V> inline double getTotalCost(SCIP_SOL* sol, vector<V> vector){
      double value = 0 ;
      for(V vect: vector)
         value += getTotalCost(sol, vect);
      return value;
   }

   inline int printStats(){
      SCIP_CALL( SCIPprintStatistics(scip_, NULL) );
   }

   inline int printBestSol(){
      SCIP_CALL( SCIPprintBestSol(scip_, NULL, FALSE) );
   }

   inline int writeProblem(string fileName){
      SCIP_CALL( SCIPwriteOrigProblem(scip_, fileName.c_str(), "lp", FALSE) );
   }

   inline int writeLP(string fileName){
      SCIP_CALL( SCIPwriteLP(scip_, fileName.c_str()) );
   }

   /**************
    * Getters *
    *************/

   inline SCIP* getScip(){
      return scip_;
   }

private:

   inline int initializeSCIP(const char* name){
      /* initialize SCIP environment */
      SCIP_CALL( SCIPcreate(&scip_) );

      /* include default plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(scip_) );

      /* create empty problem */
      SCIP_CALL( SCIPcreateProb(scip_, name, 0, 0, 0, 0, 0, 0, 0) );
   }

   //delete all the model built by scip
   inline int deleteSCIP(){
      SCIP_CALL( SCIPfree(&scip_) );
      BMScheckEmptyMemory();
   }

   //SCIP pointer
   SCIP* scip_;
};


#endif /* SRC_SCIPMODELER_H_ */
