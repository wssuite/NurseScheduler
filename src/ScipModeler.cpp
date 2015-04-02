/*
 * MyScipModeler.cpp
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#include "ScipModeler.h"

   /** reduced cost pricing method of variable pricer for feasible LPs */
   /** Pricing of additional variables if LP is feasible.
    *
    *  - get the values of the dual variables you need
    *  - construct the reduced-cost arc lengths from these values
    *  - find the shortest admissible tour with respect to these lengths
    *  - if this tour has negative reduced cost, add it to the LP
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
    *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
    */
   SCIP_DECL_PRICERREDCOST(ScipPricer::scip_redcost)
   {
      /* call pricing routine and set result pointer, see above*/
      *result = SCIP_SUCCESS;
      if(!pPricer_->pricing())
         *result = SCIP_DIDNOTRUN;

      return SCIP_OKAY;
   }

ScipModeler::ScipModeler(const char* name): Modeler()
{
   initializeSCIP(name);
}

ScipModeler::~ScipModeler()
{
   deleteSCIP();
}

//solve the model
int ScipModeler::solve(bool relaxation){
   if(relaxation)
      SCIP_CALL( SCIPsolve(scip_) );
   else
      SCIP_CALL( SCIPsolve(scip_) );
}

//Add a pricer
int ScipModeler::addObjPricer(MyPricer* pPricer){
   ScipPricer* pricer2 = new ScipPricer(pPricer, scip_);
   /* include the pricer and delete it (because of the true) */
   SCIP_CALL( SCIPincludeObjPricer(scip_, pricer2, true) );
   /* activate the pricer */
   SCIP_CALL( SCIPactivatePricer(scip_, SCIPfindPricer(scip_, pricer2->scip_name_)) );
}

/*
 * Create variable:
 *    var is a pointer to the pointer of the variable
 *    var_name is the name of the variable
 *    lhs, rhs are the lower and upper bound of the variable
 *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
 */
int ScipModeler::createVar(MyObject** var, const char* var_name, double objCoeff, double lb, double ub, VarType vartype, double score){
   if(lb==DBL_MIN)
      lb = -SCIPinfinity(scip_);
   if(ub==DBL_MAX)
      ub = SCIPinfinity(scip_);

   SCIP_VARTYPE type;
   switch(vartype){
   case VARTYPE_BINARY:
      type = SCIP_VARTYPE_BINARY;
      break;
   case VARTYPE_INTEGER:
      type = SCIP_VARTYPE_INTEGER;
      break;
   default:
      type = SCIP_VARTYPE_CONTINUOUS;
      break;
   }

   SCIP_VAR* var2;
   SCIP_CALL( SCIPcreateVar(scip_, &(var2), var_name, lb, ub, objCoeff, type,
      true, false, 0, 0, 0, 0, 0) );

   if(score > 0)
      SCIP_CALL( SCIPaddPricedVar(scip_, var2, score) );
   else SCIP_CALL( SCIPaddVar(scip_, var2) );

   *var = new ScipVar(var2);
   objects_.push_back(*var);
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

int ScipModeler::createConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
   vector<MyObject*> vars, vector<double> coeffs){
   vector<SCIP_VAR*> vars2;
   for(MyObject* var: vars){
      ScipVar* var2 = (ScipVar*) var;
      vars2.push_back(var2->var_);
   }

   SCIP_CONS* con2;
   SCIP_CALL( SCIPcreateConsLinear(scip_, &( con2 ), con_name, vars2.size(), &(vars2)[0], &(coeffs)[0], lhs, rhs,
      true, false, true, true, true, false, true, false, false, false) );
   SCIP_CALL( SCIPaddCons(scip_, con2) );

   *con = new ScipCons(con2);
   objects_.push_back(*con);
}

//Add final linear constraints
int ScipModeler::createFinalConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
   vector<MyObject*> vars, vector<double> coeffs){
   vector<SCIP_VAR*> vars2;
   for(MyObject* var: vars){
      ScipVar* var2 = (ScipVar*) var;
      vars2.push_back(var2->var_);
   }

   SCIP_CONS* con2;
   SCIP_CALL( SCIPcreateConsLinear(scip_, &(con2), con_name, vars2.size(), &(vars2)[0], &(coeffs)[0], lhs, rhs,
      true, false, true, true, true, false, false, false, false, false) );
   SCIP_CALL( SCIPaddCons(scip_, con2) );

   *con = new ScipCons(con2);
   objects_.push_back(*con);
}

/*
 * Add variables to constraints
 */

int ScipModeler::addCoefLinear(MyObject* cons, MyObject* var, double coeff, bool transformed){
   ScipCons* cons2 = (ScipCons*) cons;
   SCIP_CONS* cons3 = cons2->cons_;
   ScipVar* var2 = (ScipVar*) var;
   SCIP_VAR* var3 = var2->var_;

   if(transformed)
      getTransformedCons(cons3, &(cons3));
   SCIP_CALL( SCIPaddCoefLinear(scip_, cons3, var3, coeff) );
}

/*
 * Get the transformed variables and constraints
 *
 *  Because SCIP transforms the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */

int ScipModeler::getTransformedVar(SCIP_VAR* var, SCIP_VAR** var2){
   SCIP_CALL( SCIPgetTransformedVar(scip_, var, var2 ) );
}

int ScipModeler::getTransformedCons(SCIP_CONS* cons, SCIP_CONS** cons2){
   SCIP_CALL( SCIPgetTransformedCons(scip_, cons, cons2 ) );
}

/*
 * get the primal values
 */

SCIP_SOL* ScipModeler::getBestSol(){
   return SCIPgetBestSol(scip_);
}

double ScipModeler::getVarValue(MyObject* var){
   ScipVar* var2 = (ScipVar*) var;
   SCIP_VAR* var3 = var2->var_;
   SCIP_SOL* sol = getBestSol();

   SCIPgetSolVal(scip_, sol, var3);
}

/*
 * Get the dual variables
 */

double ScipModeler::getDual(MyObject* cons, bool transformed){
   ScipCons* cons2 = (ScipCons*) cons;
   SCIP_CONS* cons3 = cons2->cons_;

   if(transformed)
      getTransformedCons(cons3, &(cons3));
   SCIPgetDualsolLinear(scip_, cons3);
}

/**************
 * Parameters *
 *************/
int ScipModeler::setVerbosity(int v){
   SCIP_CALL( SCIPsetIntParam(scip_, "display/verblevel", v) );
   /* SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); */
}

/**************
 * Outputs *
 *************/

//compute the total cost of SCIP_VAR* in the solution sol*
double ScipModeler::getTotalCost(MyObject* var){
   ScipVar* var2 = (ScipVar*) var;
   SCIP_VAR* var3 = var2->var_;

   double value = getVarValue(var);
   return value *  var3->branchfactor;
}

int ScipModeler::printStats(){
   SCIP_CALL( SCIPprintStatistics(scip_, NULL) );
}

int ScipModeler::printBestSol(){
   SCIP_CALL( SCIPprintBestSol(scip_, NULL, FALSE) );
}

int ScipModeler::writeProblem(string fileName){
   SCIP_CALL( SCIPwriteOrigProblem(scip_, fileName.c_str(), "lp", FALSE) );
}

int ScipModeler::writeLP(string fileName){
   SCIP_CALL( SCIPwriteLP(scip_, fileName.c_str()) );
}

/*
 * Private methods
 */


int ScipModeler::initializeSCIP(const char* name){
   /* initialize SCIP environment */
   SCIP_CALL( SCIPcreate(&scip_) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip_) );

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip_, name, 0, 0, 0, 0, 0, 0, 0) );
}

//delete all the model built by scip
int ScipModeler::deleteSCIP(){
   SCIP_CALL( SCIPfree(&scip_) );
   BMScheckEmptyMemory();
}
