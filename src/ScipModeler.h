/*
 * ScipModeler.h
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#ifndef SRC_SCIPMODELER_H_
#define SRC_SCIPMODELER_H_

#include "Modeler.h"

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "objscip/objpricer.h"

using namespace scip;

/*
 * My Variables
 */
class ScipVar: public MyObject {
public:
   ScipVar(SCIP_VAR* var):MyObject(){ var_=var; }
   ~ScipVar(){ }
   void get(SCIP_VAR** var){ *var = var_; }
   void set(SCIP_VAR* var){ var_ = var; }
protected:
   SCIP_VAR* var_;
};

/*
 * My Constraints
 */
struct ScipCons: public MyObject  {
public:
   ScipCons(SCIP_CONS* cons):MyObject(){ cons_=cons; }
   ~ScipCons(){ }
   void get(SCIP_CONS** cons){ *cons = cons_; }
   void set(SCIP_CONS* cons){ cons_ = cons; }
protected:
   SCIP_CONS* cons_;
};

/*
 * My Pricer
 */
struct ScipPricer: public MyObject  {
public:
   ScipPricer(ObjPricer* pricer):MyObject(){ pricer_=pricer; }
   ~ScipPricer(){ }
   void get(ObjPricer** pricer){ *pricer = pricer_; }
   void set(ObjPricer* pricer){ pricer_ = pricer; }
protected:
   ObjPricer* pricer_;
};

/*
 * My Branching Rule
 */
//struct ScipRule: public MyRule  {
//   ScipRule(){ }
//   ~ScipRule(){ }
//   SCIP_VAR* rule_;
//};

class ScipModeler: public Modeler
{
public:
   ScipModeler(const char* name);
   ~ScipModeler();

   //solve the model
   int solve(bool relaxation = false);

   //Add a pricer
   int addObjPricer(MyObject* pPricer);

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
    */
   int createVar(MyObject** var, const char* var_name, double objCoeff, double lb, double ub, VarType vartype, double score);

   /*
    * Create linear constraint:
    *    con is a pointer to the pointer of the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
    *    coeffs is the array of coefficient to add to the constraints
    */

   int createConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
      vector<MyObject*> vars = {}, vector<double> coeffs = {});

   //Add final linear constraints
   int createFinalConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
      vector<MyObject*> vars = {}, vector<double> coeffs = {});

   /*
    * Add variables to constraints
    */

   int addCoefLinear(MyObject* cons, MyObject* var, double coeff, bool transformed=false);

   /*
    * Get the transformed variables and constraints
    *
    *  Because SCIP transforms the original problem in preprocessing, we need to get the references to
    *  the variables and constraints in the transformed problem from the references in the original
    *  problem.
    */

   int getTransformedVar(SCIP_VAR* var, SCIP_VAR** var2);

   int getTransformedCons(SCIP_CONS* cons, SCIP_CONS** cons2);

   double getVarValue(MyObject* var);

   /*
    * Get the dual variables
    */

   double getDual(MyObject* cons, bool transformed = false);

   /**************
    * Parameters *
    *************/
   int setVerbosity(int v);

   /**************
    * Outputs *
    *************/

   //compute the total cost of SCIP_VAR* in the solution sol*
   double getTotalCost(MyObject* var);

   int printStats();

   int printBestSol();

   int writeProblem(string fileName);

   int writeLP(string fileName);

   /**************
    * Getters *
    *************/

   ScipModeler* getModel();

   Scip* getScip();


private:
   //initialize the model
   int initializeSCIP(const char* name);
   //delete all the model built by scip
   int deleteSCIP();
   //get the primal values
   SCIP_SOL* getBestSol();

   //SCIP pointer
   SCIP* scip_;
};

#endif /* SRC_SCIPMODELER_H_ */
