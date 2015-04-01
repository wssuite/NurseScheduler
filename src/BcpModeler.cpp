/*
 * MyBcpModeler.cpp
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#include "BcpModeler.h"
#include "OsiClpSolverInterface.hpp"
#include "CoinTime.hpp"
#include "BCP_lp.hpp"

/*
 * BCP_lp_user methods
 */

//Convert a set of variables into corresponding columns for the current LP relaxation.
void BcpLpModel::vars_to_cols(const BCP_vec<BCP_cut*>& cuts, // on what to expand
   BCP_vec<BCP_var*>& vars,       // what to expand
   BCP_vec<BCP_col*>& cols,       // the expanded cols
   // few things that the user can use for lifting vars if allowed
   const BCP_lp_result& lpres,
   BCP_object_origin origin, bool allow_multiple)
{
   TransformVarsToColumns(vars, cols);
}

//vars = are just the giver vars
//cols is the vector where the new columns will be stored
void BcpLpModel::TransformVarsToColumns(BCP_vec<BCP_var*>& vars, BCP_vec<BCP_col*>& cols){
   const int varnum = vars.size();
   if (varnum == 0)
      return;
   cols.reserve(varnum);

   for (int i = 0; i < varnum; ++i) {
      CoinVar* var = dynamic_cast<CoinVar*>(vars[i]);
      if(!var)
         Tools::throwError("Bad variable casting.");
      const int size = var->indexRows_.size();
      int* indices = &(var->indexRows_[0]);
      double* coeffs = &(var->coeffs_[0]);
      cols.unchecked_push_back(
         new BCP_col(size, indices, coeffs, var->cost_, var->lb_, var->ub_));
   }
}

//Generate variables within the LP process.
void BcpLpModel::generate_vars_in_lp(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars, const BCP_vec<BCP_cut*>& cuts, const bool before_fathom,
   BCP_vec<BCP_var*>& new_vars, BCP_vec<BCP_col*>& new_cols)
{
   pModel_->setLPSol(lpres);
   pModel_->pricing(0);

   //check if new columns add been added since the last time
   //if there are some, add all of them in new_vars
   int size = pModel_->getNbColumns();
   if ( (size != nbCurrentColumnVarsBeforePricing_) || ! before_fathom) {
      new_vars.reserve(size-nbCurrentColumnVarsBeforePricing_); //reserve the memory for the new columns
      for(int i=nbCurrentColumnVarsBeforePricing_; i<size; ++i)
         new_vars.push_back((BcpColumn*)pModel_->getColumns()[i]);
//      TransformVarsToColumns(new_vars, new_cols);
      nbCurrentColumnVarsBeforePricing_ = size;
      return;
   }

   // must be before fathoming. we need vars with red cost below the
   // negative of (lpobj-ub)/ks_num otherwise we can really fathom.
   //    const double rc_bound =
   //   (lpres.dualTolerance() + (lpres.objval() - upper_bound()))/kss.ks_num;
   //    generate_vars(lpres, vars, rc_bound, new_vars);
}

/*
 * BcpBranchingTree
 */

// setting the base
//Create the core of the problem by filling out the last three arguments.
void BcpBranchingTree::initialize_core(BCP_vec<BCP_var_core*>& vars,
   BCP_vec<BCP_cut_core*>& cuts, BCP_lp_relax*& matrix){
   //define nb rows and col
   const int rownum = pModel_->getCons().size();
   const int colnum = pModel_->getCoreVars().size();

   // bounds and objective
   double lb[colnum], ub[colnum], obj[colnum], rhs[rownum], lhs[rownum];
   //copy of the core variables
   vars.reserve(colnum);
   for(int i=0; i<colnum; ++i){
      BcpCoreVar* var = dynamic_cast<BcpCoreVar*>(pModel_->getCoreVars()[i]);
      if(!var)
         Tools::throwError("Bad variable casting.");
      vars.push_back(var);
      lb[i] = var->lb_;
      ub[i] = var->ub_;
      obj[i] = var->cost_;
   }

   //copy of the core cuts
   cuts.reserve(rownum);
   for(int i=0; i<rownum; ++i){
      BcpCoreCons* cut = dynamic_cast<BcpCoreCons*>(pModel_->getCons()[i]);
      if(!cut)
         Tools::throwError("Bad constraint casting.");
      cuts.push_back(cut);
      rhs[i] = cut->rhs_;
      lhs[i] = cut->lhs_;
   }

   matrix = new BCP_lp_relax;
   matrix->copyOf(pModel_->buildCoinMatrix(true), obj, lb, ub, rhs, lhs);
}

// create the root node
//Create the set of extra variables and cuts that should be added
//to the formulation in the root node.
void BcpBranchingTree::create_root(BCP_vec<BCP_var*>& added_vars,
   BCP_vec<BCP_cut*>& added_cuts,
   BCP_user_data*& user_data){

   added_vars.reserve(nbInitialColumnVars_);
   for(int i=0; i<nbInitialColumnVars_; ++i)
      added_vars.unchecked_push_back((BcpColumn*)pModel_->getColumns()[i]);
}

/*
 * BcpModeler
 */

//solve the model
int BcpModeler::solve(bool relaxatione){
   BcpInitialize bcp(this);
   char* argv[0];
   return bcp_main(0, argv, &bcp);
}

/*
 * Create core variable:
 *    var is a pointer to the pointer of the variable
 *    var_name is the name of the variable
 *    lb, ub are the lower and upper bound of the variable
 *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
 */
int BcpModeler::createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub){
   *var = new BcpCoreVar(var_name, index, objCoeff, vartype, lb, ub);
   return 1;
}

int BcpModeler::createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub){
   *var = new BcpColumn(var_name, index, objCoeff, vartype, lb, ub);
   return 1;
}


/*
 * Create linear constraint:
 *    con is a pointer to the pointer of the constraint
 *    con_name is the name of the constraint
 *    lhs, rhs are the lower and upper bound of the constraint
 *    nonZeroVars is the number of non-zero coefficients to add to the constraint
 */

int BcpModeler::createCoinConsLinear(CoinCons** con, const char* con_name, int index, double lhs, double rhs){
   *con = new BcpCoreCons(con_name, index, lhs, rhs);
   return 1;
}

/*
 * get the primal values
 */

double BcpModeler::getVarValue(MyObject* var){
   CoinVar* var2 = (CoinVar*) var;
   double value = 0;

   BcpCoreVar* core_var = dynamic_cast<BcpCoreVar*>(var2);
   if (core_var) {
   }

   BcpColumn* column_var = dynamic_cast<BcpColumn*>(var2);
   if (column_var) {
   }

   return value;
}

/*
 * Get the dual variables
 */

double BcpModeler::getDual(MyObject* cons, bool transformed){
   BcpCoreCons* cons2 = (BcpCoreCons*) cons;
   return 0;
}

/**************
 * Parameters *
 *************/
int BcpModeler::setVerbosity(int v){
}

/**************
 * Outputs *
 *************/

int BcpModeler::printStats(){
}

int BcpModeler::printBestSol(){
}

int BcpModeler::writeProblem(string fileName){
}

int BcpModeler::writeLP(string fileName){
//   OsiClpSolverInterface solver = ;
//   solver.writeLp(fileName.c_str(), "lp", 1e-5, 10, 5);
}
