/*
 * MyScipModeler.cpp
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#include "ScipModeler.h"

ScipModeler::ScipModeler(const char* name)
{
   initializeSCIP(name);
}

ScipModeler::~ScipModeler()
{
   deleteSCIP();
}

//solve the model
int ScipModeler::solve(bool relaxation = false){
   if(relaxation)
      SCIP_CALL( SCIPsolve(scip_) );
   else
      SCIP_CALL( SCIPsolve(scip_) );
}

//Add a pricer
int ScipModeler::addObjPricer(MyPricer pricer){
   ObjPricer* pricer2 = ((ScipPricer)pricer).pricer_;
   /* include the pricer */
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
int ScipModeler::createVar(MyVar* var, const char* var_name, double objCoeff, double lb, double ub, SCIP_VARTYPE vartype, double score){
   if(lb==DBL_MIN)
      lb = -SCIPinfinity(scip_);
   if(ub==DBL_MAX)
      ub = SCIPinfinity(scip_);

   SCIP_VAR* var2;
   SCIP_CALL( SCIPcreateVar(scip_, &(var2), var_name, lb, ub, objCoeff, vartype,
      true, false, 0, 0, 0, 0, 0) );

   if(score > 0)
      SCIP_CALL( SCIPaddPricedVar(scip_, var2, score) );
   else SCIP_CALL( SCIPaddVar(scip_, var2) );

   var = new ScipVar(var2);
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

int ScipModeler::createConsLinear(MyCons* con, const char* con_name, double lhs, double rhs,
   vector<MyVar*> vars = {}, vector<double> coeffs = {}){
   vector<SCIP_VAR*> vars2;
   for(ScipVar* var: vars)
      vars2.push_back(var->var_);
   SCIP_CONS* con2 = ((ScipCons*)con)->cons_;
   SCIP_CALL( SCIPcreateConsLinear(scip_, &( con2 ), con_name, vars2.size(), &(vars2)[0], &(coeffs)[0], lhs, rhs,
      true, false, true, true, true, false, true, false, false, false) );
   SCIP_CALL( SCIPaddCons(scip_, con2) );

   con = new ScipCons(con2);
}

//Add final linear constraints
int ScipModeler::createFinalConsLinear(MyCons* con, const char* con_name, double lhs, double rhs,
   vector<MyVar*> vars = {}, vector<double> coeffs = {}){
   vector<SCIP_VAR*> vars2;
   for(ScipVar* var: vars)
      vars2.push_back(var->var_);
   SCIP_CONS* con2 = ((ScipCons*)con)->cons_;
   SCIP_CALL( SCIPcreateConsLinear(scip_, &(con2), con_name, vars2.size(), &(vars2)[0], &(coeffs)[0], lhs, rhs,
      true, false, true, true, true, false, false, false, false, false) );
   SCIP_CALL( SCIPaddCons(scip_, con2) );

   con = new ScipCons(con2);
}

/*
 * Add variables to constraints
 */

int ScipModeler::addCoefLinear(MyCons cons, MyVar var, double coeff){
   getTransformedCons(((ScipCons)cons).cons_, &(((ScipCons)cons).cons_));
   SCIP_CALL( SCIPaddCoefLinear(scip_, ((ScipCons)cons).cons_, ((ScipVar)var).var_, coeff) );
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

double ScipModeler::getVarValue(MyVar var){
   SCIPgetSolVal(scip_, getBestSol(), ((ScipVar)var).var_);
}

/*
 * Get the dual variables
 */

double ScipModeler::getDual(MyCons cons, bool transformed = false){
   if(transformed)
      getTransformedCons(((ScipCons)cons).cons_, &(((ScipCons)cons).cons_));
   SCIPgetDualsolLinear(scip_, ((ScipCons)cons).cons_);
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
double ScipModeler::getTotalCost(MyVar* var){
   double value = getVarValue(*var);
   return value *  ((ScipVar*)var)->var_->branchfactor;
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

/**************
 * Getters *
 *************/

ScipModeler* ScipModeler::getModel(){
   return this;
}

Scip* ScipModeler::getScip(){
   return scip_;
}


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
