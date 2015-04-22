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

//Modify parameters of the LP solver before optimization.
//This method provides an opportunity for the user to change parameters of the LP solver before optimization in the LP solver starts.
//The second argument indicates whether the optimization is a "regular" optimization or it will take place in strong branching.
//Default: empty method.
void BcpLpModel::modify_lp_parameters ( OsiSolverInterface* lp, const int changeType, bool in_strong_branching){
   for(pair<BCP_lp_par::chr_params, bool> entry: pModel_->getLpParameters())
      set_param(entry.first, entry.second);
}

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

      //Copy the vectors var->getIndexRows() and var->getCoeffRows() in arrays
      const int size = var->getNbRows();

      //create a new array which will be deleted by ~BCP_col()
      int* indexRows = new int[size];
      vector<int> index = var->getIndexRows();
      std::copy(index.begin(), index.end(), indexRows);

      //create a new array which will be deleted by ~BCP_col()
      double* coeffRows= new double[size];
      vector<double> coeffs = var->getCoeffRows();
      std::copy(coeffs.begin(), coeffs.end(), coeffRows);

      cols.unchecked_push_back(
         new BCP_col(size, indexRows, coeffRows, var->getCost(), var->getLB(), var->getUB()));
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
   if ( size != nbCurrentColumnVarsBeforePricing_ ) { //|| ! before_fathom
      new_vars.reserve(size-nbCurrentColumnVarsBeforePricing_); //reserve the memory for the new columns
      for(int i=nbCurrentColumnVarsBeforePricing_; i<size; ++i){
         BcpColumn* var = dynamic_cast<BcpColumn*>(pModel_->getColumns()[i]);
         if(!var)
            Tools::throwError("Bad variable casting.");
         //create a new BcpColumn which will be deleted by BCP
         new_vars.unchecked_push_back(new BcpColumn(*var));
      }
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
 * BCP_DoNotBranch_Fathomed: The node should be fathomed without even trying to branch.
 * BCP_DoNotBranch: BCP should continue to work on this node.
 * BCP_DoBranch: branch on one of the candidates cands
 *
 */
BCP_branching_decision BcpLpModel::select_branching_candidates(const BCP_lp_result& lpres, //the result of the most recent LP optimization.
   const BCP_vec<BCP_var*> &  vars, //the variables in the current formulation.
   const BCP_vec< BCP_cut*> &  cuts, //the cuts in the current formulation.
   const BCP_lp_var_pool& local_var_pool, //the local pool that holds variables with negative reduced cost.
   //In case of continuing with the node the best so many variables will be added to the formulation (those with the most negative reduced cost).
   const BCP_lp_cut_pool& local_cut_pool, //the local pool that holds violated cuts.
   //In case of continuing with the node the best so many cuts will be added to the formulation (the most violated ones).
   BCP_vec<BCP_lp_branching_object*>&  cands, //the generated branching candidates.
   bool force_branch) //indicate whether to force branching regardless of the size of the local cut/var pools{
{
   pModel_->setLPSol(lpres);

   //if some variables have been generated, do not branch
   if(local_var_pool.size() > 0)
      return BCP_DoNotBranch;

   //manage integrality on the skill allocation variables
   vector<CoinVar*> integerCoreVariables;
   for(CoinVar* coreVar: pModel_->getCoreVars()){
      if(coreVar->is_integer() && !pModel_->isInteger(coreVar)){
//         integerCoreVariables.push_back(coreVar);
         appendCoreIntegerVars(coreVar, cands);
         return BCP_DoBranch;
      }
   }
//   if(integerCoreVariables.size() > 0){
//      for(CoinVar* coreVar: integerCoreVariables)
//         appendCoreIntegerVars(coreVar, cands);
//      return BCP_DoBranch;
//   }

   //fixing candidates
   vector<MyObject*> fixingCandidates;
   pModel_->logical_fixing(fixingCandidates);

   //fix if some candidates
   if(fixingCandidates.size()>0){
      appendNewFixingVars(fixingCandidates, cands);
      return BCP_DoBranch;
   }

   //branching candidates
   vector<MyObject*> branchingCandidates;
   pModel_->branching_candidates(branchingCandidates);

   //branch if some candidates
   if(branchingCandidates.size() > 0){
      appendNewBranchingVars(branchingCandidates, cands);
      return BCP_DoBranch;
   }

   //otherwise fathomed
   return BCP_DoNotBranch_Fathomed;
}

void BcpLpModel::appendCoreIntegerVars(CoinVar* coreVar, BCP_vec<BCP_lp_branching_object*>&  cands){
   BCP_vec<int> vpos; //positions of the variables
   BCP_vec<double> vbd; // old bound and then new one for each variable

   vpos.push_back(coreVar->getIndex());
   double value = pModel_->getVarValue(coreVar);
   vbd.push_back(coreVar->getLB()); // old lower bound
   vbd.push_back(floor(value)); // new upper bound
   vbd.push_back(ceil(value)); // new lower bound
   vbd.push_back(coreVar->getUB()); // new lower bound


   cands.push_back(new  BCP_lp_branching_object(2, //just one children where
      //all the columns with positions in vpos are fixed to 1
      0, 0, /* vars/cuts_to_add */
      &vpos, 0, &vbd, 0, /* forced parts: position and bounds (old bound and then new one) */
      0, 0, 0, 0 /* implied parts */));
}

void BcpLpModel::appendNewFixingVars(vector<MyObject*> columns, BCP_vec<BCP_lp_branching_object*>&  cands){
   BCP_vec<int> vpos; //positions of the variables
   BCP_vec<double> vbd; // old bound and then new one for each variable

   for (MyObject* var: columns) {
      BcpColumn* col = dynamic_cast<BcpColumn*>(var);
      if(col){
         vpos.push_back(col->getIndex());
         //fix the old bound (0) and the new bound (1)
         vbd.push_back(0); // old lower bound
         vbd.push_back(1); // new lower bound
      }
      else
         Tools::throwError("The variable fixed is not a column.");
   }

   cands.push_back(new  BCP_lp_branching_object(1, //just one children where
      //all the columns with positions in vpos are fixed to 1
      0, 0, /* vars/cuts_to_add */
      &vpos, 0, &vbd, 0, /* forced parts: position and bounds (old bound and then new one) */
      0, 0, 0, 0 /* implied parts */));
}

void BcpLpModel::appendNewBranchingVars(vector<MyObject*> columns, BCP_vec<BCP_lp_branching_object*>&  cands){
   const int nbChildren = columns.size();
   BCP_vec<int> vpos; //positions of the variables
   BCP_vec<double> vbd; // old bound and then new one for each variable and for each children
   //this vector is filled is this order:
   //for the first child: old and new bounds for all the variables in vpos
   //then for the second child: old and new bounds for all the variables in vpos
   //....

   int child = 0;
   for (MyObject* var: columns) {
      BcpColumn* col = dynamic_cast<BcpColumn*>(var);
      if(col)
         vpos.push_back(col->getIndex());
      else
         Tools::throwError("The variable fixed is not a column.");
      //set the new bound of col to 1 and all the others columns keep 0
      for(int i=0; i<nbChildren; ++i){
         //fix the old bound (0)
         vbd.push_back(0);
         //if i==child, the new bound==1, otherwise 0
         vbd.push_back( (i==child) ? 1: 0 );
      }
      ++child;
   }

   cands.push_back(new  BCP_lp_branching_object(nbChildren, 0, 0, /* vars/cuts_to_add */
      &vpos, 0, &vbd, 0, /* forced parts */
      0, 0, 0, 0 /* implied parts */));
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
      //create a new BcpCoreVar which will be deleted by BCP
      vars.push_back(new BcpCoreVar(*var));
      lb[i] = var->getLB();
      ub[i] = var->getUB();
      obj[i] = var->getCost();
   }

   //copy of the core cuts
   cuts.reserve(rownum);
   for(int i=0; i<rownum; ++i){
      BcpCoreCons* cut = dynamic_cast<BcpCoreCons*>(pModel_->getCons()[i]);
      if(!cut)
         Tools::throwError("Bad constraint casting.");
      //create a new BcpCoreCons which will be deleted by BCP
      cuts.push_back(new BcpCoreCons(*cut));
      lhs[i] = cut->getLhs();
      rhs[i] = cut->getRhs();
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
   for(int i=0; i<nbInitialColumnVars_; ++i){
      BcpColumn* var = dynamic_cast<BcpColumn*>(pModel_->getColumns()[i]);
      if(!var)
         Tools::throwError("Bad variable casting.");
      //create a new BcpColumn which will be deleted by BCP
      added_vars.unchecked_push_back(new BcpColumn(*var));
   }

}

void BcpBranchingTree::display_feasible_solution(const BCP_solution* sol){
   // store the solution
   pModel_->setBestUb(sol->objective_value());

   BCP_solution_generic* sol2 = (BCP_solution_generic*) sol;
   vector<double> primal(pModel_->getNbVars());
   for(int i=0; i<sol2->_vars.size(); ++i){
      int index = sol2->_vars[i]->bcpind();
      primal[index] = sol2->_values[i];
   }
   pModel_->setPrimal(primal);
}

// various initializations before a new phase (e.g., pricing strategy)
void BcpBranchingTree::init_new_phase(int phase, BCP_column_generation& colgen, CoinSearchTreeBase*& candidates) {
   colgen = BCP_GenerateColumns;
   for(pair<BCP_tm_par::chr_params, bool> entry: pModel_->getTmParameters())
      set_param(entry.first, entry.second);
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
   objects_.push_back(*var);
   return 1;
}

int BcpModeler::createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, double dualObj, VarType vartype, double lb, double ub){
   *var = new BcpColumn(var_name, index, objCoeff, dualObj, vartype, lb, ub);
   objects_.push_back(*var);
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
   objects_.push_back(*con);
   return 1;
}

/*
 * get the primal values
 */

double BcpModeler::getVarValue(MyObject* var){
   CoinVar* var2 = (CoinVar*) var;
   if(primalValues_.size() ==0 )
      Tools::throwError("Primal solution has been initialized.");
   return primalValues_[var2->getIndex()];
}

/*
 * Get the dual variables
 */

double BcpModeler::getDual(MyObject* cons, bool transformed){
   CoinCons* cons2 = (CoinCons*) cons;
   if(dualValues_.size() == 0)
      Tools::throwError("Dual solution has been initialized.");
   return dualValues_[cons2->getIndex()];
}

/**************
 * Parameters *
 *************/
int BcpModeler::setVerbosity(int v){

   if(v>=1){
      tm_parameters[BCP_tm_par::VerbosityShutUp] = 1;
      tm_parameters[BCP_tm_par::TmVerb_First] = 1;
      tm_parameters[BCP_tm_par::TmVerb_BestFeasibleSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_FinalStatistics] = 1;
      tm_parameters[BCP_tm_par::TmVerb_Last] = 1;

      lp_parameters[BCP_lp_par::LpVerb_IterationCount] = 1; // Print the "Starting iteration x" line. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_GeneratedVarCount] = 1; // Print the number of variables generated during this iteration. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_FinalRelaxedSolution] = 1; // Turn on the user hook "display_lp_solution" for the last LP relaxation solved at a search tree node. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_Last] = 1; // Just a marker for the last LpVerb
   }

   if(v>=2){
      tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolutionValue] = 1;
      tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_NewPhaseStart] = 1;

      lp_parameters[BCP_lp_par::LpVerb_ReportVarGenTimeout] = 1; // Print information if receiving variables is timed out. (BCP_lp_generate_vars)
      lp_parameters[BCP_lp_par::LpVerb_ReportLocalVarPoolSize] = 1; // Similar as above for variables. (BCP_lp_generate_vars)
      lp_parameters[BCP_lp_par::LpVerb_AddedVarCount] = 1; // Print the number of variables added from the local variable pool in the current iteration. (BCP_lp_main_loop)
   }

   if(v>=3){
      tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolutionValue] = 1;
      tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_PrunedNodeInfo] = 1;

      lp_parameters[BCP_lp_par::LpVerb_ChildrenInfo] = 1; // After a branching object is selected print what happens to the presolved children (e.g., fathomed). (BCP_print_brobj_stat)
      lp_parameters[BCP_lp_par::LpVerb_ColumnGenerationInfo] = 1; // Print the number of variables generated before resolving the Lp ir fathoming a node. (BCP_lp_fathom)
      lp_parameters[BCP_lp_par::LpVerb_FathomInfo] = 1; // Print information related to fathoming. (BCP_lp_main_loop, BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
      lp_parameters[BCP_lp_par::LpVerb_RelaxedSolution] = 1; // Turn on the user hook "display_lp_solution". (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_LpSolutionValue] = 1; // Print the size of the problem matrix and the LP solution value after resolving the LP. (BCP_lp_main_loop)
   }

   if(v>=4){
      tm_parameters[BCP_tm_par::TmVerb_TimeOfImprovingSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_TrimmedNum] = 1;
      tm_parameters[BCP_tm_par::TmVerb_ReportDefault] = 1;

      lp_parameters[BCP_lp_par::LpVerb_MatrixCompression] = 1; // Print the number of columns and rows that were deleted during matrix compression. (BCP_lp_delete_cols_and_rows)
      lp_parameters[BCP_lp_par::LpVerb_PresolvePositions] = 1; // Print detailed information about all the branching candidates during strong branching. LpVerb_PresolveResult must be set for this parameter to have an effect. (BCP_lp_perform_strong_branching)
      lp_parameters[BCP_lp_par::LpVerb_PresolveResult] = 1; // Print information on the presolved branching candidates during strong branching. (BCP_lp_perform_strong_branching)
      lp_parameters[BCP_lp_par::LpVerb_ProcessedNodeIndex] = 1; // Print the "Processing NODE x on LEVEL y" line. (BCP_lp-main_loop)
      lp_parameters[BCP_lp_par::LpVerb_VarTightening] = 1; // Print the number of variables whose bounds have been changed by reduced cost fixing or logical fixing. (BCP_lp_fix_vars)
      lp_parameters[BCP_lp_par::LpVerb_RowEffectivenessCount] = 1; // Print the number of ineffective rows in the current problem. The definition of what rows are considered ineffective is determined by the paramter IneffectiveConstraints. (BCP_lp_adjust_row_effectiveness)
      lp_parameters[BCP_lp_par::LpVerb_StrongBranchPositions] = 1; // Print detailed information on the branching candidate selected by strong branching. LpVerb_StrongBranchResult must be set fo this parameter to have an effect. (BCP_print_brobj_stat)
      lp_parameters[BCP_lp_par::LpVerb_StrongBranchResult] = 1; // Print information on the branching candidate selected by strong branching. (BCP_print_brobj_stat)
   }

   return 1;
}

/**************
 * Outputs *
 *************/

int BcpModeler::printStats(){
}

int BcpModeler::writeProblem(string fileName){
}

int BcpModeler::writeLP(string fileName){
   //   OsiClpSolverInterface solver = ;
   //   solver.writeLp(fileName.c_str(), "lp", 1e-5, 10, 5);
}
