/*
 * ScipModeler.h:
 *    Tools to create a SCIP problem.
 *    In general, when you want to create a new scip object,
 *    you have to give a pointer to the pointer of the object.
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#ifndef SRC_MODELER_H_
#define SRC_MODELER_H_

/* standard library includes */
#include <cfloat>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>

#include "MyTools.h"

/* namespace usage */
using namespace std;

/*
 * My Variables
 */
struct MyVar {
   MyVar(){ }
   ~MyVar(){ }
};

/*
 * My Constraints
 */
struct MyCons {
   MyCons(){ }
   ~MyCons(){ }
};

/*
 * My Pricer
 */
struct MyPricer {
   MyPricer(){ }
   ~MyPricer(){ }
};

/*
 * My Branching Rule
 */
struct MyRule {
   MyRule(){ }
   ~MyRule(){ }
};

/*
 * Var types
 */
enum VarType {VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY};

class Modeler {
public:

   Modeler(){ }

   virtual ~Modeler(){ }

   //solve the model
   virtual int solve(bool relaxation = false);

   //Add a pricer
   virtual int addObjPricer(MyPricer pPricer);

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
    */

   virtual void createVar(MyVar* var, const char* var_name, double objCoeff, double lb, double ub, VarType vartype, double score);

   inline void createPositiveVar(MyVar* var, const char* var_name, double objCoeff, double score = 0, double ub = DBL_MAX){
      return createVar(var, var_name, objCoeff, 0.0, ub, VARTYPE_CONTINUOUS, score);
   }

   inline void createIntVar(MyVar* var, const char* var_name, double objCoeff, double score = 0, double ub = DBL_MAX){
      return createVar(var, var_name, objCoeff, 0, ub, VARTYPE_INTEGER, score);
   }

   inline void createBinaryVar(MyVar* var, const char* var_name, double objCoeff, double score = 0){
      return createVar(var, var_name, objCoeff, 0.0, 1.0, VARTYPE_BINARY, score);
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

   virtual int createConsLinear(MyCons* cons, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {});

   //Add a lower or equal constraint
   inline void createLEConsLinear(MyCons* cons, const char* con_name, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, DBL_MIN, rhs, vars, coeffs);
   }

   //Add a greater or equal constraint
   inline void createGEConsLinear(MyCons* cons, const char* con_name, double lhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, lhs, DBL_MAX, vars, coeffs);
   }

   //Add an equality constraint
   inline void createEQConsLinear(MyCons* cons, const char* con_name, double eq,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, eq, eq, vars, coeffs);
   }

   //Add final linear constraints
   virtual int createFinalConsLinear(MyCons* cons, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {});

   inline void createFinalLEConsLinear(MyCons* cons, const char* con_name, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, DBL_MIN, rhs, vars, coeffs);
   }

   inline void createFinalGEConsLinear(MyCons* cons, const char* con_name, double lhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, lhs, DBL_MAX, vars, coeffs);
   }

   inline void createFinalEQConsLinear(MyCons* cons, const char* con_name, double eq,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, eq, eq, vars, coeffs);
   }

   /*
    * Add variables to constraints
    */

   virtual int addCoefLinear(MyCons cons, MyVar var, double coeff);

   /*
    * Add new Column to the problem
    */

   inline void createColumn(MyVar* var, const char* var_name, double objCoeff, VarType vartype,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      switch(vartype){
      case VARTYPE_BINARY:
         createBinaryVar(var, var_name, objCoeff, score);
         break;
      case VARTYPE_INTEGER:
         createIntVar(var, var_name, objCoeff, score);
         break;
      default:
         createPositiveVar(var, var_name, objCoeff, score);
         break;
      }

      for(int i=0; i<cons.size(); i++){
         MyCons con = *(cons[i]);
         addCoefLinear(con, *var, coeffs[i]);
      }
   }

   inline void createPositiveColumn(MyVar* var, const char* var_name, double objCoeff,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, VARTYPE_CONTINUOUS, cons, coeffs, transformed, score);
   }

   inline void createBinaryColumn(MyVar* var, const char* var_name, double objCoeff,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, VARTYPE_BINARY, cons, coeffs, transformed, score);
   }

   inline void createIntColumn(MyVar* var, const char* var_name, double objCoeff,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, VARTYPE_INTEGER, cons, coeffs, transformed, score);
   }

   /*
    * get the primal values
    */

   virtual double getVarValue(MyVar var);

   inline vector<double> getVarValues(vector<MyVar*> vars){
      vector<double> values(vars.size());
      for(int i=0; i<vars.size(); ++i)
         values[i] = getVarValue(*(vars[i]));
      return values;
   }

   /*
    * Get the dual variables
    */

   virtual  double getDual(MyCons cons, bool transformed = false);

   inline vector<double> getDuals(vector<MyCons*> cons, bool transformed = false){
      vector<double> dualValues(cons.size());
      for(int i=0; i<cons.size(); ++i)
         dualValues[i] = getDual(*(cons[i]), transformed);
      return dualValues;
   }

   /**************
    * Parameters *
    *************/
   virtual int setVerbosity(int v);

   /**************
    * Outputs *
    *************/

   //compute the total cost of MyVar in the solution
   virtual double getTotalCost(MyVar* var);

   //compute the total cost of a vector of MyVar in the solution
   template<typename T>  inline double getTotalCost(map<MyVar*, T> map0){
      double value = 0 ;
      for(pair<MyVar*, T> var: map0)
         value += getTotalCost(var.first);
      return value;
   }

   //compute the total cost of a multiple vectors of MyVar in the solution
   template<typename V> inline double getTotalCost(vector<V> vector){
      double value = 0 ;
      for(V vect: vector)
         value += getTotalCost(vect);
      return value;
   }

   virtual int printStats();

   virtual int printBestSol();

   virtual int writeProblem(string fileName);

   virtual int writeLP(string fileName);

   /**************
    * Getters *
    *************/

   template<typename M> M getModel(){
      string error = "This template has not been implemented.";
      Tools::throwError(error.c_str());
   }
};


#endif /* SRC_MODELER_H_ */
