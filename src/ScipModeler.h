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
struct ScipVar: public MyVar {
   ScipVar(SCIP_VAR* var){ var_=var; }
   ~ScipVar(){ }
   SCIP_VAR* var_;
};

/*
 * My Constraints
 */
struct ScipCons: public MyCons  {
   ScipCons(SCIP_CONS* cons){ cons_=cons; }
   ~ScipCons(){ }
   SCIP_CONS* cons_;
};

/*
 * My Pricer
 */
struct ScipPricer: public MyPricer  {
   ScipPricer(ObjPricer* pricer){ pricer_=pricer; }
   ~ScipPricer(){ }
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
   int addObjPricer(MyPricer pPricer);

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
    */
   int createVar(MyVar* var, const char* var_name, double objCoeff, double lb, double ub, SCIP_VARTYPE vartype, double score);

   /*
    * Create linear constraint:
    *    con is a pointer to the pointer of the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
    *    coeffs is the array of coefficient to add to the constraints
    */

   int createConsLinear(MyCons* con, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {});

   //Add final linear constraints
   int createFinalConsLinear(MyCons* con, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {});

   /*
    * Add variables to constraints
    */

   int addCoefLinear(MyCons cons, MyVar var, double coeff);

   /*
    * Get the transformed variables and constraints
    *
    *  Because SCIP transforms the original problem in preprocessing, we need to get the references to
    *  the variables and constraints in the transformed problem from the references in the original
    *  problem.
    */

   int getTransformedVar(SCIP_VAR* var, SCIP_VAR** var2);

   int getTransformedCons(SCIP_CONS* cons, SCIP_CONS** cons2);

   double getVarValue(MyVar var);

   /*
    * Get the dual variables
    */

   double getDual(MyCons cons, bool transformed = false);

   /**************
    * Parameters *
    *************/
   virtual int setVerbosity(int v);

   /**************
    * Outputs *
    *************/

   //compute the total cost of SCIP_VAR* in the solution sol*
   double getTotalCost(MyVar* var);

   virtual int printStats();

   virtual int printBestSol();

   virtual int writeProblem(string fileName);

   virtual int writeLP(string fileName);

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
