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
#include "BCP_lp_node.hpp"
#include "CbcModeler.h"
#include "RotationPricer.h"

/*
 * BCP_lp_user methods
 */

BcpLpModel::BcpLpModel(BcpModeler* pModel):
pModel_(pModel),nbCurrentColumnVarsBeforePricing_(pModel->getNbColumns()),
lpIteration_(0), last_node(-1), heuristicHasBeenRun_(false)
{ }

//Initialize the lp parameters and the OsiSolver
OsiSolverInterface* BcpLpModel::initialize_solver_interface(){
   for(pair<BCP_lp_par::chr_params, bool> entry: pModel_->getLpParameters())
      set_param(entry.first, entry.second);
   OsiClpSolverInterface * clp = new OsiClpSolverInterface();
   int verbosity = max(0, pModel_->getVerbosity()-1);
   clp->messageHandler()->setLogLevel(verbosity);
   return clp;
}

//Try to generate a heuristic solution (or return one generated during cut/variable generation.
//Return a pointer to the generated solution or return a NULL pointer.
BCP_solution* BcpLpModel::generate_heuristic_solution(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts){

   BCP_solution_generic* sol = NULL;

   //if heuristic has already been run in these node or
   //it has not been long enough since the last run or
   //the objective of the sub-problem is too negative
   if(heuristicHasBeenRun_ || current_index()%10 != 0 || pModel_->getLastMinDualCost() < -1)
      return sol;

   heuristicHasBeenRun_ = true;

   //copy the solver
   OsiSolverInterface* solver = getLpProblemPointer()->lp_solver;

   //store the basis
   const CoinWarmStart* ws = solver->getWarmStart();

//   // prepare for heuristic branching
//   solver->markHotStart();

   //define different size
   const int size = vars.size(), coreSize = pModel_->getCoreVars().size(), nbColumnsActive = size - coreSize;

   //store lower bounds
   map<int, double> indexColLbChanched;

   //while the solution is feasible
   solver->resolve();
   while( solver->isProvenOptimal() ){

      //find the best not integer columns
      vector<pair<int,double>> candidates;
      for(int i=coreSize; i<size; ++i){
         double value = solver->getColSolution()[i];
         if(value < EPSILON || value > 1 - EPSILON)
            continue;
         candidates.push_back(pair<int,double>(i, 1-value));
      }

      stable_sort(candidates.begin(), candidates.end(), compareCol);

      //if we have found a column
      if(candidates.size() > 0){
         double valueLeft = .99;
         for(pair<int,double>& p: candidates){
            if(p.second > valueLeft)
               break;
            if(p.second > .2)
               valueLeft -= p.second;
            indexColLbChanched.insert( pair<int, double>(p.first, solver->getColLower()[p.first]) );
            solver->setColLower(p.first, 1);
         }
         solver->resolve();
      }
      //else the solution is integer, create a BCP_solution_generic to return
      else{
         sol = new BCP_solution_generic();
         for(int i=0; i<size; ++i)
            if(solver->getColSolution()[i] > EPSILON){
               //create new var that will be deleted by the solution sol
               if(i<coreSize){
                  BcpCoreVar* var0 = dynamic_cast<BcpCoreVar*>(pModel_->getCoreVars()[i]);
                  sol->add_entry(new BcpCoreVar(*var0), solver->getColSolution()[i]);
               }
               else{
                  BcpColumn* var0 = dynamic_cast<BcpColumn*>(vars[i]);
                  sol->add_entry(new BcpColumn(*var0), solver->getColSolution()[i]);
               }
            }

         break;
      }
   }

   //restore bounds
   for(pair<int, double> p: indexColLbChanched)
      solver->setColLower(p.first, p.second);

//   // indicate to the lp solver that the heuristic branching is done
//   solver->unmarkHotStart();
   solver->setWarmStart(ws);

   delete ws;

   return sol;
}

bool BcpLpModel::compareCol(const pair<int,double>& p1, const pair<int,double>& p2){
   return (p1.second < p2.second);
}

//Modify parameters of the LP solver before optimization.
//This method provides an opportunity for the user to change parameters of the LP solver before optimization in the LP solver starts.
//The second argument indicates whether the optimization is a "regular" optimization or it will take place in strong branching.
//Default: empty method.
void BcpLpModel::modify_lp_parameters ( OsiSolverInterface* lp, const int changeType, bool in_strong_branching){
   if(current_index() != last_node){
      last_node = current_index();
      printSummaryLine();
   }
}

//print in cout a line summary of the current solver state
void BcpLpModel::printSummaryLine(const BCP_vec<BCP_var*>& vars){

   if(pModel_->getVerbosity() > 0){

//      double lower_bound = (getLpProblemPointer()->node->true_lower_bound < DBL_MIN) ? pModel_->myMax :
//         getLpProblemPointer()->node->true_lower_bound;

      if( vars.size() == 0 ){
         printf("BCP: %13s %5s | %10s %10s %10s | %8s %10s %12s %10s | %10s %5s %5s \n",
            "Node", "Lvl", "BestUB", "RootLB", "BestLB","#It",  "Obj", "#Frac", "#Active", "ObjSP", "#SP", "#Col");
         printf("BCP: %5d / %5d %5d | %10.0f %10.2f %10.2f | %8s %10s %12s %10s | %10s %5s %5s \n",
            current_index(), pModel_->getTreeSize(), current_level(),
            pModel_->getBestUB(), pModel_->getRootLB(), pModel_->getBestLB(),
            "-", "-", "-", "-", "-", "-", "-");
      }

      else{
         /* compute number of fractional columns */
         int frac = 0;
         int non_zero = 0;
         int nbCol = 0;
         for(int i=0; i<nbCurrentColumnVarsBeforePricing_; ++i){
            double value = pModel_->getVarValue(pModel_->getColumns()[i]);
            if(value < EPSILON)
               continue;
            non_zero ++;
            if(value < 1 - EPSILON)
               frac++;
         }

         int nbColGenerated = pModel_->getNbColumns() - nbCurrentColumnVarsBeforePricing_;

         printf("BCP: %5d / %5d %5d | %10.0f %10.2f %10.2f | %8d %10.2f %5d / %4d %10d | %10.2f %5d %5d  \n",
            current_index(), pModel_->getTreeSize(), current_level(),
            pModel_->getBestUB(), pModel_->getRootLB(), pModel_->getBestLB(),
            lpIteration_, pModel_->getLastObj(), frac, non_zero, vars.size() - pModel_->getCoreVars().size(),
            pModel_->getLastMinDualCost(), pModel_->getLastNbSubProblemsSolved(), nbColGenerated);
      }
   }
}

//stop this node or BCP
bool BcpLpModel::doStop(){
   pModel_->doStop();

   //fathom if the true lower bound greater than current upper bound
   if(pModel_->getBestUB() - getLpProblemPointer()->node->true_lower_bound <
         pModel_->getAbsoluteGap() - EPSILON)
      return true;

   return false;
}

//This method provides an opportunity for the user to tighten the bounds of variables.
//The method is invoked after reduced cost fixing. The results are returned in the last two parameters.
//Parameters:
//lpres    the result of the most recent LP optimization,
//vars  the variables in the current formulation,
//status   the stati of the variables as known to the system,
//var_bound_changes_since_logical_fixing    the number of variables whose bounds have changed (by reduced cost fixing) since the most recent invocation of this method that has actually forced changes returned something in the last two arguments,
//changed_pos    the positions of the variables whose bounds should be changed
//new_bd   the new bounds (lb/ub pairs) of these variables.
///////////////////////////////////////////////////////////////////////////////////////
//WARNING: JUST FIX BOUNDS THAT ARE COMPULSORY IN AN INTEGER POINT OF VIEW
//AS THE FINAL LP BOUND OF THE NODE IS GIVEN TO ITS CHILDREN,
//FIXING WRONG BOUNDS CAN LEAD TO A FALSE LOWER BOUND, THUS FATHOMING WRONG NODES.
///////////////////////////////////////////////////////////////////////////////////////
void BcpLpModel::logical_fixing (const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts,
   const BCP_vec<BCP_obj_status>& var_status,
   const BCP_vec<BCP_obj_status>& cut_status,
   const int var_bound_changes_since_logical_fixing,
   BCP_vec<int>& changed_pos,
   BCP_vec<double>& new_bd){

}

// Restoring feasibility.
//This method is invoked before fathoming a search tree node that has been found infeasible and
//the variable pricing did not generate any new variables.
void BcpLpModel::restore_feasibility(const BCP_lp_result& lpres,
   const std::vector<double*> dual_rays,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts,
   BCP_vec<BCP_var*>& vars_to_add,
   BCP_vec<BCP_col*>& cols_to_add){
   //dive is finished
   //dive_ = false;
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
      vector<int>& index = var->getIndexRows();
      copy(index.begin(), index.end(), indexRows);

      //create a new array which will be deleted by ~BCP_col()
      double* coeffRows = new double[size];
      vector<double>& coeffs = var->getCoeffRows();
      copy(coeffs.begin(), coeffs.end(), coeffRows);

      cols.unchecked_push_back(
         new BCP_col(size, indexRows, coeffRows, var->getCost(), var->getLB(), var->getUB()) );
   }
}

//Generate variables within the LP process.
void BcpLpModel::generate_vars_in_lp(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars, const BCP_vec<BCP_cut*>& cuts, const bool before_fathom,
   BCP_vec<BCP_var*>& new_vars, BCP_vec<BCP_col*>& new_cols)
{
   if(doStop())
      return;

   ++lpIteration_;
   pModel_->setLPSol(lpres, vars);
   pModel_->pricing(0, before_fathom);

   /* Print a line summary of the solver state */
   printSummaryLine(vars);

   //check if new columns add been added since the last time
   //if there are some, add all of them in new_vars
   int size = pModel_->getNbColumns();
   if ( size != nbCurrentColumnVarsBeforePricing_ ) { //|| ! before_fathom
      new_vars.reserve(size-nbCurrentColumnVarsBeforePricing_); //reserve the memory for the new columns
      for(int i=nbCurrentColumnVarsBeforePricing_; i<size; ++i){
         BcpColumn* var = dynamic_cast<BcpColumn*>(pModel_->getColumns()[i]);
         //create a new BcpColumn which will be deleted by BCP
         new_vars.unchecked_push_back(new BcpColumn(*var));
      }
      nbCurrentColumnVarsBeforePricing_ = size;
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
	//update node
	if(local_var_pool.size() == 0)
		pModel_->updateNodeLB(lpres.objval());

   //stop this process for BCP or the node
   if(doStop()){
      return BCP_DoNotBranch_Fathomed;

   }

   //if some variables have been generated, do not branch
   if(local_var_pool.size() > 0)
      return BCP_DoNotBranch;

   //update true_lower_bound, as we reach the end of the column generation
   getLpProblemPointer()->node->true_lower_bound = lpres.objval();
   heuristicHasBeenRun_ = false;

   //fathom if greater than current upper bound
   if(pModel_->getBestUB() - lpres.objval() < pModel_->getAbsoluteGap() - EPSILON)
      return BCP_DoNotBranch_Fathomed;

   //branching candidates: numberOfNursesByPosition_, rest on a day, ...
   vector<MyObject*> branchingCandidates;
   pModel_->branching_candidates(branchingCandidates);

   //fixing candidates: branch on columns greater than BRANCHLB
   vector<MyObject*> fixingCandidates;
   pModel_->logical_fixing(fixingCandidates);

   if(branchingCandidates.size() == 0)
      return BCP_DoNotBranch_Fathomed;

//   CoinVar* integerCoreVar = dynamic_cast<CoinVar*>(branchingCandidates[0]);
//   appendNewBranchingVarsOnNumberOfNurses(integerCoreVar, fixingCandidates, vars, cands);

   appendNewBranchingVarsOnRest(cuts.size(), branchingCandidates, fixingCandidates, vars, cands);

   return BCP_DoBranch;
}

//Decide what to do with the children of the selected branching object.
//Fill out the _child_action field in best. This will specify for every child what to do with it.
//Possible values for each individual child are BCP_PruneChild, BCP_ReturnChild and BCP_KeepChild.
//There can be at most child with this last action specified.
//It means that in case of diving this child will be processed by this LP process as the next search tree node.
//Default: Every action is BCP_ReturnChild.
//However, if BCP dives then one child will be mark with BCP_KeepChild. The decision which child to keep is based on the ChildPreference parameter in BCP_lp_par.
//Also, if a child has a presolved lower bound that is higher than the current upper bound then that child is mark as BCP_FathomChild.
void BcpLpModel::set_actions_for_children(BCP_presolved_lp_brobj* best){
   best->action()[0] = BCP_KeepChild;
}

void BcpLpModel::appendNewBranchingVarsOnNumberOfNurses(CoinVar* integerCoreVar, vector<MyObject*>& columns,
   const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands){
//   const int nbChildren = 2+columns.size();
   const int nbChildren = (columns.size() > 0 ) ? 3 : 2;

   BCP_vec<int> vpos; //positions of the variables
   BCP_vec<double> vbd; // old bound and then new one for each variable and for each children
   //this vector is filled is this order:
   //for the first child: old and new bounds for all the variables in vpos
   //then for the second child: old and new bounds for all the variables in vpos
   //....

   // set the positions vector and the first node for the columns */
   vector<int> currentIndex = buildBranchingColumns(integerCoreVar, columns, vars, vpos, vbd);

   /*
    * BOUNDS VECTOR
    */

   //keep same bound for the columns
   for(int j=0; j<columns.size(); ++j){
      vbd.push_back(vars[currentIndex[j]]->lb());
      vbd.push_back(vars[currentIndex[j]]->ub());
   }
   //branch on floor(value) of integerCoreVar
   double value = pModel_->getVarValue(integerCoreVar);
   vbd.push_back(integerCoreVar->getLB()); // old lower bound
   vbd.push_back(floor(value)); // new upper bound
   /* update tree */
   pModel_->pushBackNewNode(integerCoreVar, integerCoreVar->getLB(),floor(value));

   //keep same bound for the columns
   for(int j=0; j<columns.size(); ++j){
      vbd.push_back(vars[currentIndex[j]]->lb());
      vbd.push_back(vars[currentIndex[j]]->ub());
   }
   //branch on ceil(value) of integerCoreVar
   vbd.push_back(ceil(value)); // new lower bound
   vbd.push_back(integerCoreVar->getUB()); // new lower bound
   /* update tree */
   pModel_->pushBackNewNode(integerCoreVar, ceil(value),integerCoreVar->getUB());

   /* create the new candidate */
   cands.push_back(new  BCP_lp_branching_object(nbChildren, 0, 0, /* vars/cuts_to_add */
      &vpos, 0, &vbd, 0, /* forced parts */
      0, 0, 0, 0 /* implied parts */));
}


   //Branch on the core integer var: the rest arcs on a day for a nurse
   //just 2 children
   //Try also to fix to 1 some columns
   void BcpLpModel::appendNewBranchingVarsOnRest(int nbCuts, vector<MyObject*>& coreVars, vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands){

      const int nbChildren = (columns.size() > 0 ) ? 3 : 2;

      BCP_vec<BCP_cut*> new_cuts; //add a branching cut for the set of arcs
      BCP_vec<int> vpos; //positions of the variables
      BCP_vec<double> vbd; // old bound and then new one for each variable and for each children
      BCP_vec<int> cpos; //positions of the cuts
      BCP_vec<double> cbd; //bounds of the cuts

      // set the positions vector and the first node for the columns */
      vector<int> currentIndex = buildBranchingColumns(0, columns, vars, vpos, vbd);

      //keep same bounds for the columns for two resting branching nodes
      for(int j=0; j<columns.size(); ++j){
         vbd.push_back(vars[currentIndex[j]]->lb());
         vbd.push_back(vars[currentIndex[j]]->ub());
      }
      for(int j=0; j<columns.size(); ++j){
         vbd.push_back(vars[currentIndex[j]]->lb());
         vbd.push_back(vars[currentIndex[j]]->ub());
      }

      /* create the day off */
      pair<LiveNurse*, int> dayOff = pModel_->getLastBranchingRest();

      //creating the branching cut
      char name[50];
      sprintf(name, "RestBranchingCons_N%d_%d", dayOff.first->id_, dayOff.second);
      vector<int> indexes;
      vector<double> coeffs;
      for(MyObject* var: coreVars){
         CoinVar* var2 = dynamic_cast<CoinVar*>(var);
         indexes.push_back(var2->getIndex());
         coeffs.push_back(1);
      }
      //create a new BcpBranchCons which will be deleted by BCP
      BcpBranchCons* cons = new BcpBranchCons(name, pModel_->getCons().size()+pModel_->getBranchingCons().size(),
         0, 1, indexes, coeffs);
      new_cuts.push_back(cons);
      //create our own new BcpBranchCons
      pModel_->pushBackBranchingCons(new BcpBranchCons(*cons));

      /* bounds of the cut */
      cpos.push_back(nbCuts);
      if(columns.size() > 0){
         cbd.push_back(0); cbd.push_back(1);
      }
      //push the node rest on day dayOff.second
      cbd.push_back(1); cbd.push_back(1);
      /* update tree */
      pModel_->pushBackNewNode(dayOff.first, dayOff.second, true, coreVars);
      //push the node work on day dayOff.second
      cbd.push_back(0); cbd.push_back(0);
      /* update tree */
      pModel_->pushBackNewNode(dayOff.first, dayOff.second, false, coreVars);

      //creates the nodes

      cands.push_back(new  BCP_lp_branching_object(nbChildren, 0, &new_cuts, /* vars/cuts_to_add */
         &vpos, &cpos, &vbd, &cbd, /* forced parts */
         0, 0, 0, 0 /* implied parts */));
   }

   //Try also to fix to 1 some columns
   void BcpLpModel::appendNewBranchingVarsOnColumns(vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands){

      const int nbChildren = 1;
      BCP_vec<int> vpos; //positions of the variables
      BCP_vec<double> vbd; // old bound and then new one for each variable and for each children

      buildBranchingColumns(0, columns, vars, vpos, vbd);

      cands.push_back(new  BCP_lp_branching_object(nbChildren, 0, 0, /* vars/cuts_to_add */
         &vpos, 0, &vbd, 0, /* forced parts */
         0, 0, 0, 0 /* implied parts */));
   }

   //build the vector of the branching candidates for the columns
   //return the indexes of the columns in the current formulation
   vector<int> BcpLpModel::buildBranchingColumns(CoinVar* var, vector<MyObject*>& columns, const BCP_vec<BCP_var*>&  vars,
      BCP_vec<int>& vpos, BCP_vec<double>& vbd){
      //current index in the BCP formulation
      vector<int> currentIndex(columns.size());
      for(int j=0; j<columns.size(); ++j){
         bool found = false;
         BcpColumn* col = dynamic_cast<BcpColumn*>(columns[j]);
         for(int i=pModel_->getCoreVars().size(); i<vars.size(); ++i){
            BcpColumn* col2 = dynamic_cast<BcpColumn*>(vars[i]);
            if(col->getIndex() == col2->getIndex()){
               currentIndex[j] = i;
               found = true;
               break;
            }
         }
         if(!found)
            Tools::throwError("The column has not been found.");
      }

      //positions vector
      for(int j=0; j<columns.size(); ++j)
         vpos.push_back(currentIndex[j]);
      if(var) vpos.push_back(var->getIndex());

      //branch on all the columns
      if(columns.size() > 0){
         //set all the columns bounds to 1
         for(int l=0; l<columns.size(); ++l){
            vbd.push_back(1);
            vbd.push_back(vars[currentIndex[l]]->ub());
         }
         //don't change bounds for coreVars
         if(var) {
            vbd.push_back(var->getLB());
            vbd.push_back(var->getUB());
         }
      }

      /* update tree */
      if(columns.size() > 0)
         pModel_->pushBackNewNode(columns);

      return currentIndex;
   }

//Here, we generate a cut to branch on a set of variables
void BcpLpModel::cuts_to_rows(const BCP_vec<BCP_var*>& vars, //the variables currently in the relaxation (IN)
BCP_vec<BCP_cut*>& cuts, //the cuts to be converted (IN/OUT)
BCP_vec<BCP_row*>& rows, //the rows into which the cuts are converted (OUT)
const BCP_lp_result&   lpres, //solution to the current LP relaxation (IN)
BCP_object_origin origin, //where the cuts come from (IN)
bool  allow_multiple) //whether multiple expansion, i.e., lifting, is allowed (IN)
{
   rows.reserve(cuts.size());
   for(BCP_cut* cut: cuts){
      BcpBranchCons* branchingCut = dynamic_cast<BcpBranchCons*>(cut);
      if(!cut)
         Tools::throwError("Should be a branching cut.");

      //create new arrays which will be deleted by ~BCP_row()
      const int size = branchingCut->getIndexCols().size();
      int* elementIndices = new int[size];
      copy(branchingCut->getIndexCols().begin(), branchingCut->getIndexCols().end(), elementIndices);
      double* elementValues = new double[size];
      copy(branchingCut->getCoeffCols().begin(), branchingCut->getCoeffCols().end(), elementValues);
      rows.unchecked_push_back( new BCP_row(size, elementIndices, elementValues, branchingCut->getLhs(), branchingCut->getRhs()) );
   }
}

/*
 * BcpBranchingTree
 */

BcpBranchingTree::BcpBranchingTree(BcpModeler* pModel):
      pModel_(pModel) , nbInitialColumnVars_(pModel->getNbColumns()), minGap_(0.0)
{}

// setting the base
//Create the core of the problem by filling out the last three arguments.
void BcpBranchingTree::initialize_core(BCP_vec<BCP_var_core*>& vars,
   BCP_vec<BCP_cut_core*>& cuts, BCP_lp_relax*& matrix){
   // initialize tm parameters
   set_param(BCP_tm_par::TmVerb_SingleLineInfoFrequency, pModel_->getFrequency());
   //always dive
   set_param(BCP_tm_par::UnconditionalDiveProbability, 1);
   if (pModel_->getMaxSolvingtime() < DBL_MAX)
      set_param(BCP_tm_par::MaxRunTime, pModel_->getMaxSolvingtime());
   for(pair<BCP_tm_par::chr_params, bool> entry: pModel_->getTmParameters())
      set_param(entry.first, entry.second);

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
   pModel_->setBestUB(sol->objective_value());

   //store the solution
   pModel_->addBcpSol(sol);
}

// various initializations before a new phase (e.g., branching strategy)
void BcpBranchingTree::init_new_phase(int phase, BCP_column_generation& colgen, CoinSearchTreeBase*& candidates) {
   colgen = BCP_GenerateColumns;
//   CoinSearchTreeCompareHighestIncrease::pModel = pModel_;
   switch(pModel_->getSearchStrategy()){
         case BestFirstSearch:
            set_param(BCP_tm_par::TreeSearchStrategy, 0);
            candidates = new MyCoinSearchTree<CoinSearchTreeCompareBest>(pModel_);
            break;
         case BreadthFirstSearch:
            set_param(BCP_tm_par::TreeSearchStrategy, 1);
            candidates = new MyCoinSearchTree<CoinSearchTreeCompareBreadth>(pModel_);
            break;
         case DepthFirstSearch:
            set_param(BCP_tm_par::TreeSearchStrategy, 2);
            candidates = new MyCoinSearchTree<CoinSearchTreeCompareDepth>(pModel_);
            break;
         default:
            candidates = new MyCoinSearchTree<CoinSearchTreeCompareBest>(pModel_);
            break;
         }
}

/*
 * BcpModeler
 */
BcpModeler::BcpModeler(const char* name):
   CoinModeler(), currentNode_(0), tree_size_(1), nb_nodes_last_incumbent_(0), diveDepth_(0), diveLenght_(myMax),
   primalValues_(0), dualValues_(0), reducedCosts_(0), lhsValues_(0),
   best_lb_in_root(myMax), best_lb(DBL_MAX), lastNbSubProblemsSolved_(0), lastMinDualCost_(0)
{
   //create the root
   pushBackNewNode();
}

//solve the model
int BcpModeler::solve(bool relaxatione){
   BcpInitialize bcp(this);
   char* argv[0];
   try{
      return bcp_main(0, argv, &bcp);
   }catch(BcpStop& e) { }

   return 0;
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
 * Set the solution
 */

void BcpModeler::setLPSol(const BCP_lp_result& lpres, const BCP_vec<BCP_var*>&  vars){
   obj_history_.push_back(lpres.objval());

   //copy the new arrays in the vectors for the core vars
   const int nbCoreVar = coreVars_.size();
   const int nbColVar = columnVars_.size();
   const int nbCons = cons_.size();

   //clear all
   primalValues_.clear();
   dualValues_.clear();
   reducedCosts_.clear();
   lhsValues_.clear();

   //assign value for core variables
   primalValues_.assign(lpres.x(), lpres.x()+nbCoreVar);
   dualValues_.assign(lpres.pi(), lpres.pi()+nbCons);
   reducedCosts_.assign(lpres.dj(), lpres.dj()+nbCoreVar);
   lhsValues_.assign(lpres.lhs(), lpres.lhs()+nbCons);

   //reserve some space for the columns and fill it with 0
   double zeroArray[nbColVar];
   CoinFillN(zeroArray, nbColVar, 0.0);
   primalValues_.insert(primalValues_.end(), zeroArray, zeroArray+nbColVar);
   reducedCosts_.insert(reducedCosts_.end(), zeroArray, zeroArray+nbColVar);
   //loop through the variables and link the good columns together
   for(int i=nbCoreVar; i<vars.size(); ++i){
      BCP_var* var0 = vars[i];
      BcpColumn* var = dynamic_cast<BcpColumn*>(var0);
      primalValues_[var->getIndex()] = lpres.x()[i];
      reducedCosts_[var->getIndex()] = lpres.dj()[i];
   }
}

void BcpModeler::addBcpSol(const BCP_solution* sol){
   //create a solution which is not going to delete the vars at the end
   BCP_solution_generic mySol(false);
   BCP_solution_generic* sol2 = (BCP_solution_generic*) sol;

   int coreSize = coreVars_.size();
   for(int i=0; i<sol2->_vars.size(); ++i){
      BcpColumn* col = dynamic_cast<BcpColumn*>(sol2->_vars[i]);
      if(col){
         BcpColumn* myCol = (BcpColumn*) columnVars_[col->getIndex()-coreSize];
         mySol.add_entry(myCol, sol2->_values[i]);
      }
      else{
         BcpCoreVar* myVar = (BcpCoreVar*) coreVars_[sol2->_vars[i]->bcpind()];
         mySol.add_entry(myVar, sol2->_values[i]);
      }
   }

   bcpSolutions_.push_back(mySol);
}

bool BcpModeler::loadBestSol(){
   int i=0, index = -1;
   double bestObj = DBL_MAX;
   for(BCP_solution_generic& sol: bcpSolutions_){
      if(sol.objective_value() < bestObj){
         bestObj = sol.objective_value();
         index = i;
      }
      ++i;
   }
   if(index == -1)
      return false;

   loadBcpSol(index);
   return true;
}

void BcpModeler::loadBcpSol(int index){
   BCP_solution_generic& sol = bcpSolutions_[index];
   const int size1 = sol._vars.size(), size2 = getNbVars();
   vector<double> primal(size2);
   for(int i=0; i<sol._vars.size(); ++i){
      int index = sol._vars[i]->bcpind();
      primal[index] = sol._values[i];
   }
   setPrimal(primal);
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
   verbosity_ = v;

   if(v>=1){
      tm_parameters[BCP_tm_par::VerbosityShutUp] = 1;
      tm_parameters[BCP_tm_par::TmVerb_First] = 1;
      tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolutionValue] = 1;
      tm_parameters[BCP_tm_par::TmVerb_NewPhaseStart] = 1;
      tm_parameters[BCP_tm_par::TmVerb_Last] = 1;

      lp_parameters[BCP_lp_par::LpVerb_Last] = 1; // Just a marker for the last LpVerb
   }

   if(v>=2){
      TmVerb_SingleLineInfoFrequency = 1;

      lp_parameters[BCP_lp_par::LpVerb_IterationCount] = 1; // Print the "Starting iteration x" line. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_GeneratedVarCount] = 1; // Print the number of variables generated during this iteration. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_ReportVarGenTimeout] = 1; // Print information if receiving variables is timed out. (BCP_lp_generate_vars)
      lp_parameters[BCP_lp_par::LpVerb_ReportLocalVarPoolSize] = 1; // Similar as above for variables. (BCP_lp_generate_vars)
      lp_parameters[BCP_lp_par::LpVerb_AddedVarCount] = 1; // Print the number of variables added from the local variable pool in the current iteration. (BCP_lp_main_loop)
   }

   if(v>=3){
      tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolutionValue] = 1;
      tm_parameters[BCP_tm_par::TmVerb_PrunedNodeInfo] = 1;

      lp_parameters[BCP_lp_par::LpVerb_ChildrenInfo] = 1; // After a branching object is selected print what happens to the presolved children (e.g., fathomed). (BCP_print_brobj_stat)
      lp_parameters[BCP_lp_par::LpVerb_ColumnGenerationInfo] = 1; // Print the number of variables generated before resolving the Lp ir fathoming a node. (BCP_lp_fathom)
      lp_parameters[BCP_lp_par::LpVerb_FathomInfo] = 1; // Print information related to fathoming. (BCP_lp_main_loop, BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
      lp_parameters[BCP_lp_par::LpVerb_RelaxedSolution] = 1; // Turn on the user hook "display_lp_solution". (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::LpVerb_LpSolutionValue] = 1; // Print the size of the problem matrix and the LP solution value after resolving the LP. (BCP_lp_main_loop)
   }

   if(v>=4){
      tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_TimeOfImprovingSolution] = 1;
      tm_parameters[BCP_tm_par::TmVerb_TrimmedNum] = 1;
      tm_parameters[BCP_tm_par::ReportWhenDefaultIsExecuted] = 1;

      lp_parameters[BCP_lp_par::LpVerb_FinalRelaxedSolution] = 1; // Turn on the user hook "display_lp_solution" for the last LP relaxation solved at a search tree node. (BCP_lp_main_loop)
      lp_parameters[BCP_lp_par::ReportWhenDefaultIsExecuted] = 1; // Print out a message when the default version of an overridable method is executed. Default: 1.
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

//Check if BCP must stop
bool BcpModeler::doStop(){
   //continue if doesn't have a lb
   if(getBestLB() == myMax)
      return false;

   //check relative gap
   if(best_ub - getBestLB() < minRelativeGap_ * getBestLB() - EPSILON){
      char error[100];
      sprintf(error, "Stopped: relative gap < %.2f.", minRelativeGap_);
      throw BcpStop(error);
   }
   if(best_ub - getBestLB() < relativeGap_ * getBestLB() - EPSILON){
      //if the relative gap is small enough and if same incumbent since the last dive, stop
      if(nb_nodes_last_incumbent_ > diveLenght_){
         char error[100];
         sprintf(error, "Stopped: relative gap < %.2f and more than %d nodes without new incumbent.",
            relativeGap_, diveLenght_);
         throw BcpStop(error);
      }
   }
   else
      //if the relative gap is too big, wait 2 dives before stopping
      if(nb_nodes_last_incumbent_ > 2*diveLenght_){
         char error[100];
         sprintf(error, "Stopped: relative gap > %.2f and more than %d nodes without new incumbent.",
            relativeGap_, diveLenght_);
         throw BcpStop(error);
      }

   return false;
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
