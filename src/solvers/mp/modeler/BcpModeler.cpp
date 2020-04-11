/*
 * MyBcpModeler.cpp
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#include "solvers/mp/modeler/BcpModeler.h"
#include "OsiClpSolverInterface.hpp"
#include "CoinTime.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_node.hpp"
#include "solvers/mp/RotationPricer.h"
#include "solvers/mp/TreeManager.h"
#include "solvers/MasterProblem.h"
#include <string>

#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"
#endif

#ifdef USE_CBC
#include "CbcModeler.h"
#endif

/*
 * BCP_lp_user methods
 */

BcpLpModel::BcpLpModel(BcpModeler* pModel):
pModel_(pModel), lpIteration_(0), last_node(-1), heuristicHasBeenRun_(false),nbNodesSinceLastHeuristic_(0), nbGeneratedColumns_(0)
{
   // Initialization of nb_dives_to_wait_before_branching_on_columns_
   for(int i=4; i<1000000; i*=2)
      nb_dives_to_wait_before_branching_on_columns_.push_back(i);
}

//Initialize the lp parameters and the OsiSolver
OsiSolverInterface* BcpLpModel::initialize_solver_interface(){
	for(pair<BCP_lp_par::chr_params, bool> entry: pModel_->getLpParameters())
	set_param(entry.first, entry.second);
	set_param(pModel_->strong_branching.first, pModel_->strong_branching.second);

	OsiSolverInterface* solver = nullptr;
	switch(pModel_->getLPSolverType()){
		case CLP:
		solver = new OsiClpSolverInterface();
		break;
		case Gurobi:
#ifdef USE_GUROBI
		solver = new OsiGrbSolverInterface();
#else
    Tools::throwError("BCP has not been built with Gurobi.");
#endif
		break;
		case Cplex:
#ifdef USE_CPLEX
		solver = new OsiCpxSolverInterface();
#else
    Tools::throwError("BCP has not been built with Cplex.");
#endif
		break;
		default:
		Tools::throwError("The LP solver requested is not implemented.");
	}

	int verbosity = max(0, pModel_->getVerbosity()-1);
	solver->messageHandler()->setLogLevel(verbosity);
	return solver;
}

//Try to generate a heuristic solution (or return one generated during cut/variable generation.
//Return a pointer to the generated solution or return a NULL pointer.
BCP_solution* BcpLpModel::generate_heuristic_solution(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts){

   BCP_solution_generic* sol = NULL;

   //if no integer solution is needed, don't run the heuristic
   if( pModel_->getParameters().performHeuristicAfterXNode_==-1 || pModel_->getParameters().stopAfterXSolution_ == 0)
      return sol;

   //if heuristic has already been run in these node or
   //it has not been long enough since the last run or
   //the objective of the sub-problem is too negative
	if (heuristicHasBeenRun_ || nbNodesSinceLastHeuristic_<pModel_->getParameters().performHeuristicAfterXNode_)
		return sol;

	// use a criterion on the fraction of the solution that is already integer or
	// on the level in the tree (the latter criterion should depend on the
	// instance though)
	if (pModel_->getParameters().heuristicMinIntegerPercent_ > 0) {
		//count the fraction of current solution that is integer
		//
		MasterProblem* pMaster = pModel_->getMaster();
		double fractionInteger = pMaster->computeFractionOfIntegerInCurrentSolution();

		if (pModel_->getParameters().printFractionOfInteger_){
			std::cout << "FRACTION OF INTEGER ROSTERS= " << fractionInteger << std::endl;
		}

		if (100.0*fractionInteger < pModel_->getParameters().heuristicMinIntegerPercent_) {
			return sol;
		}
		std::cout << "RUN THE HEURISTIC" << std::endl;
	}

   heuristicHasBeenRun_ = true;
	nbNodesSinceLastHeuristic_=0;

   //copy the solver
   OsiSolverInterface* solver = getLpProblemPointer()->lp_solver;

   //store the basis
   const CoinWarmStart* ws = solver->getWarmStart();

   //   // prepare for heuristic branching
   //   solver->markHotStart();

   //define different size
   const int size = vars.size(), coreSize = pModel_->getCoreVars().size();

   //store lower bounds
   map<int, double> indexColLbChanged;

	// REMARK
	// another possible, but more costly heuristic is simply to solve the MIP
	// involving all the columns currently in the problem, for this simply run
	// solver->branchAndBound()

   //while the solution is feasible
   solver->resolve();
   while( solver->isProvenOptimal() ){

      //find the best not integer columns
      vector<pair<int,double>> candidates;
      for(int i=coreSize; i<size; ++i){
         double value = solver->getColSolution()[i];
         if(value < EPSILON || solver->getColLower()[i]==1) // value > 1 - EPSILON)
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
            indexColLbChanged.insert( pair<int, double>(p.first, solver->getColLower()[p.first]) );
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
	for(pair<int, double> p: indexColLbChanged)
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

		//print a line as it is the first iteration of this node
		if (pModel_->getParameters().printBcpSummary_|| pModel_->getVerbosity() > 0) {
			cout << "======================================================================================================================================" << endl;
			printSummaryLine();
		}

		if (pModel_->getParameters().printRelaxationLp_) {
			lp->writeLp("outfiles/test");
		}

		// modify dual tolerance // DBG
		// double dualTol = std::min(0.1,-pModel_->getParameters().sp_max_reduced_cost_bound_+EPSILON);
		// lp->setDblParam( OsiDualTolerance,dualTol);
		if (pModel_->getLPSolverType() == Cplex) {
#ifdef USE_CPLEX
      // do not let cplex use more than one thread, otherwise it will, and it
      // it will also keep resetting the parameter to the default value for some
      // unexplained reason
			if (OsiCpxSolverInterface* pSolverTmp=dynamic_cast<OsiCpxSolverInterface*>(lp)) {
			     CPXsetintparam(pSolverTmp->getEnvironmentPtr(),CPX_PARAM_THREADS,1);
			}
      // modifiy the value of cplex random seed otherwise the results will never
      // be reproduced
      if (OsiCpxSolverInterface* pSolverTmp=dynamic_cast<OsiCpxSolverInterface*>(lp)) {
        CPXsetintparam(pSolverTmp->getEnvironmentPtr(),CPXPARAM_RandomSeed, 1);
        std::cout << "Cplex random seed set to " << 1 << std::endl;
      }
      // use barrier optimization instead of dual simplex for reoptimization
      // after branching,because the solution parameter is just too degenerate
      if (changeType%2 == 1) {
        if (OsiCpxSolverInterface* pSolverTmp=dynamic_cast<OsiCpxSolverInterface*>(lp)) {
          CPXLPptr cpxlp = pSolverTmp->getLpPtr( OsiCpxSolverInterface::FREECACHED_RESULTS );
          CPXbaropt(pSolverTmp->getEnvironmentPtr(),cpxlp);
        }
      }
#endif
		}

		// STAB
		// set the cost and upper bounds of every core variables back to their
		// values at the time of solution
		if (pModel_->getParameters().isStabilization_) {
			for (CoinVar* pVar:pModel_->getCoreVars()) {
				int varind = pVar->getIndex();
				double cost = lp->getObjCoefficients()[varind];
				double ub = lp->getColUpper()[varind];
				pVar->setCost(cost);
				pVar->setUB(ub);
			}
		}

		// DGN
		pModel_->setNbDegenerateIt(0);
	}
	else {
		// DBG
		// double dualTol = std::min(0.1,-pModel_->getParameters().sp_max_reduced_cost_bound_+EPSILON);
		// lp->setDblParam( OsiDualTolerance,dualTol);
	}
}

//print in cout a line summary of the current solver state
void BcpLpModel::printSummaryLine(const BCP_vec<BCP_var*>& vars){

   FILE * pFile;
   pFile = pModel_->logfile().empty() ? stdout : fopen (pModel_->logfile().c_str(),"a");

   if(pModel_->getVerbosity() > 0){

		// DBG
      //      double lower_bound = (getLpProblemPointer()->node->true_lower_bound < DBL_MIN) ? pModel_->LARGE_SCORE :
      //         getLpProblemPointer()->node->true_lower_bound;

      if( vars.size() == 0 ){
         fprintf(pFile,"BCP: %13s %5s | %10s %10s %10s | %8s %10s %12s %10s | %10s %5s %5s \n",
            "Node", "Lvl", "BestUB", "RootLB", "BestLB","#It",  "Obj", "#Frac", "#Active", "ObjSP", "#SP", "#Col");
         fprintf(pFile,"BCP: %5d / %5d %5d | %10.0f %10.2f %10.2f | %8s %10s %12s %10s | %10s %5s %5s \n",
            current_index(), pModel_->getTreeSize(), current_level(),
            pModel_->getObjective(), pModel_->getRootLB(), pModel_->getBestLB(),
            "-", "-", "-", "-", "-", "-", "-");
      }

      else{
         /* compute number of fractional columns */
         int frac = 0;
         int non_zero = 0;
         for(MyVar* var: pModel_->getActiveColumns()){
            double value = pModel_->getVarValue(var);
            if(value < EPSILON)
               continue;
            non_zero ++;
            if(value < 1 - EPSILON)
               frac++;
         }

         fprintf(pFile,"BCP: %5d / %5d %5d | %10.0f %10.2f %10.2f | %8d %10.2f %5d / %4d %10ld | %10.2f %5d %5d  \n",
            current_index(), pModel_->getTreeSize(), current_level(),
            pModel_->getObjective(), pModel_->getRootLB(), pModel_->getBestLB(),
            lpIteration_, pModel_->getLastObj(), frac, non_zero, vars.size() - pModel_->getCoreVars().size(),
            pModel_->getLastMinDualCost(), pModel_->getLastNbSubProblemsSolved(), nbGeneratedColumns_);
      }
   }

   if (!pModel_->logfile().empty()) fclose(pFile);
}

//stop this node or BCP
bool BcpLpModel::doStop(){
   pModel_->doStop();

   //fathom if the true lower bound greater than current upper bound
   if(pModel_->getObjective() - getLpProblemPointer()->node->true_lower_bound <
      pModel_->getParameters().absoluteGap_ - EPSILON)
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
   // DBG
   // pModel_->checkActiveColumns(vars);
}

////////////////////////////////////////////////////////////////////////////////
// WARNING: THIS METHOD NEEDS TO BE REIMPLEMENTED AS EMPTY, BECAUSE IT CAUSES
// CYCLING IN THE COLUMN GENERATION WHEN THE LOWER BOUND IS TOO CLOSE FROM THE
// UPPER BOUND
////////////////////////////////////////////////////////////////////////////////

void BcpLpModel::reduced_cost_fixing(const double* dj, const double* x, const double gap,
		BCP_vec<BCP_var*>& vars, int& newly_changed) {
		BCP_lp_user::reduced_cost_fixing(dj,x,gap,vars,newly_changed);
		// DBG
		// std::cout << "NO REDUCED COST FIXING" << std::endl;
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
   // DBG dive_ = false;
   if (pModel_->getParameters().printBranchStats_) {
      std::cout << "DIVING STOPPED WITH INFEASIBLE SOLUTION" << std::endl;
      pModel_->printStats();
   }
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


// 	static int  cpt = 0;
// 	char  nom[1024];
// 	sprintf(nom, "/tmp/kkk%03d", cpt++);
// pModel_->pBcp_->getBcpLpModel()->getLpProblemPointer()->lp_solver->writeLp(nom);
//  cout << "LP WRITE: " << nom << endl;

}

//Generate variables within the LP process.
void BcpLpModel::generate_vars_in_lp(const BCP_lp_result& lpres,
	const BCP_vec<BCP_var*>& vars, const BCP_vec<BCP_cut*>& cuts, const bool before_fathom,
	BCP_vec<BCP_var*>& new_vars, BCP_vec<BCP_col*>& new_cols)
{
	if(doStop())
	return;


	// STAB
	// Detect when the coluln generation is stalling
	// If stabilization is used this will determine when the stabilization costs
	// are updated
	// Otherwise, an option can be set on to stop column generation after a given
	// number of degenerate iterations
	//
	bool isStall = (lpres.objval() >= pModel_->getLastObj()-EPSILON) && (lpres.objval() <= pModel_->getLastObj()+EPSILON);
	if (isStall) {
		pModel_->incrementNbDegenerateIt();
		if (current_index() > 0
		&& lpres.objval() <= pModel_->getObjective()-pModel_->getParameters().absoluteGap_+EPSILON
		&& pModel_->getNbDegenerateIt() == pModel_->getParameters().stopAfterXDegenerateIt_ ) {
			// DBG
			std::cout << "BRANCH BECAUSE COLUMN GENERATION IS STALLING" << std::endl;
			return;
		}
	}
	else {
		pModel_->setNbDegenerateIt(0);
	}

	// DBG: not sure to understand, what if the last objective value of the
	// relaxation is exactly equal to current LB?-> I changed +EPSION to -EPSILON
	// Stop the algorithm if the last objective value of the relaxation is
	// smaller than current LB
	//
	if(current_index() > 0 && lpres.objval() < pModel_->getCurrentLB() - EPSILON) {
		return;
	}

	// update the total number of LP solutions (from the beginning)
	++lpIteration_;

	// call the rotation pricer to find columns that should be added to the LP
	//
	pModel_->setLPSol(lpres, vars, lpIteration_);
	double maxReducedCost = pModel_->getParameters().sp_max_reduced_cost_bound_; // max reduced cost of a rotation that would be added to MP (a tolerance is substracted in the SP)
	vector<MyVar*> generatedColumns = pModel_->pricing(maxReducedCost, before_fathom);
	nbGeneratedColumns_ = generatedColumns.size();

	// Print a line summary of the solver state
	pModel_->setCurrentTreeLevel(current_level());
	if (pModel_->getParameters().printBcpSummary_) {
		printSummaryLine(vars);
	}

	// STAB: compute the Lagrangian bound
	// It can also be used in general to fathom nodes when the the Lagrangian
	// bound is larger than the best UB
	//
	double lagLb = -LARGE_SCORE;
	bool isImproveQuality = false;
	MasterProblem* pMaster = pModel_->getMaster();
	if ( (current_index() > 0 && pModel_->getParameters().isLagrangianFathom_)
		|| pModel_->getParameters().isStabilization_) {
		lagLb = pModel_->getBestLB();

		if (pModel_->getParameters().isStabilization_) {
			lagLb=pMaster->computeLagrangianBound(lpres.objval(),pModel_->getLastMinDualCost());
		}
		else if (pModel_->getParameters().isLagrangianFathom_) {
			//&& pModel_->getLastNbSubProblemsSolved() >= pMaster->getNbNurses()) {
			lagLb=pMaster->computeLagrangianBound(lpres.objval(),pModel_->getLastMinDualCost());
		}

		isImproveQuality = pModel_->updateNodeLagLB(lagLb);

		// LAGLB: fathom if Lagrangian bound greater than current upper bound
		if(pModel_->getParameters().isLagrangianFathom_
		&& pModel_->getObjective() - pModel_->getNodeLastLagLB() < pModel_->getParameters().absoluteGap_ - EPSILON){
			nbGeneratedColumns_ = 0;
			for(MyVar* var: generatedColumns){
				BcpColumn* col = dynamic_cast<BcpColumn*>(var);
				delete col;
			}
			generatedColumns.clear();
			std::cout << "Forcibly fathom, because Lagrangian bound is exceeded." << std::endl;
			return;
		}
	}
	//check if new columns add been added since the last time
	//if there are some, add all of them in new_vars
	//
	new_vars.reserve(nbGeneratedColumns_); //reserve the memory for the new columns
	for(MyVar* var: generatedColumns){
		BcpColumn* col = dynamic_cast<BcpColumn*>(var);
		//create a new BcpColumn which will be deleted by BCP
		new_vars.unchecked_push_back(col);
		col->addActiveIteration(lpIteration_); //initialize the counter of active iteration for this new variable
	}

	// STAB: this is where we have an opportunity to update the costs and bounds
	// of the stabilization variables
	if (pModel_->getParameters().isStabilization_) {
		this->stabUpdateBoundAndCost(isStall,isImproveQuality);
	}


// 	static int  cpt = 0;
// 	char  nom[1024];
// 	sprintf(nom, "/tmp/fff%03d", cpt++);
// pModel_->pBcp_->getBcpLpModel()->getLpProblemPointer()->lp_solver->writeLp(nom);
//  cout << "LP WRITE: " << nom << endl;

	//	//debug
	//	pModel_->checkActiveColumns(new_vars);

	// must be before fathoming. we need vars with red cost below the
	// negative of (lpobj-ub)/ks_num otherwise we can really fathom.
	//    const double rc_bound =
	//   (lpres.dualTolerance() + (lpres.objval() - upper_bound()))/kss.ks_num;
	//    generate_vars(lpres, vars, rc_bound, new_vars);
}

// STAB
// Update all the bounds and costs of the stabilization variables
// The costs are updated only when solution process is stalling or no column
// is generated
// The bounds are updated whenever the quality of the lagrangian bound varies
//
bool BcpLpModel::stabUpdateBoundAndCost(bool isStall,bool isImproveQuality) {
	SolverParam param = pModel_->getParameters();
	MasterProblem* pMaster = pModel_->getMaster();
	OsiSolverInterface* solver = getLpProblemPointer()->lp_solver;

	if (param.isStabUpdateBounds_) {
		if (isImproveQuality) {
			pMaster->stabUpdateBound(solver, 1.0/param.stabBoundFactor_);
		}
		else {
			pMaster->stabUpdateBound(solver, param.stabBoundFactor_);
		}
	}
	if (param.isStabUpdateCost_) {
		if (!nbGeneratedColumns_||isStall) {
			pMaster->stabUpdateCost(solver, param.stabCostMargin_);
		}
	}
	return false;
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
	//if some variables have been generated, do not branch
	if(local_var_pool.size() > 0 ) {
		return BCP_DoNotBranch;
	}

	// STAB: Do not branch if some stabilization variables are positive
	if (!pModel_->getMaster()->stabCheckStoppingCriterion()) {
		return BCP_DoNotBranch;
	}

	//update node
	pModel_->updateNodeLB(lpres.objval());

	//update true_lower_bound, as we reach the end of the column generation
	getLpProblemPointer()->node->true_lower_bound = lpres.objval();
	heuristicHasBeenRun_ = true;
	nbNodesSinceLastHeuristic_++;
	pModel_->incrementNbNodes();

	pModel_->setLPSol(lpres, vars, lpIteration_);

	//if root and a variable with the obj LARGE_SCORE is positive -> INFEASIBLE
	// otherwise, record the root solution for future use
	//
	if(current_index() == 0) {
		for(MyVar* col: pModel_->getActiveColumns()){
			if(((CoinVar*)col)->getCost() >= 1000.0 && pModel_->getVarValue(col) > EPSILON)
			throw InfeasibleStop("Feasibility columns are still present in the solution");
		}
		pModel_->recordLpSol();
		if (pModel_->gettimeFirstRoot() < EPSILON) {
			pModel_->settimeFirstRoot(CoinWallclockTime()-start_time());
		}
	}

	//stop this process for BCP or the node
	if(doStop()){
		return BCP_DoNotBranch_Fathomed;
	}

	//fathom if greater than current upper bound
	if(pModel_->getObjective() - lpres.objval() < pModel_->getParameters().absoluteGap_ - EPSILON){
		return BCP_DoNotBranch_Fathomed;
	}

	//LAGLB: fathom if Lagrangian bound greater than current upper bound
	if(pModel_->getParameters().isLagrangianFathom_
	&& pModel_->getObjective() - pModel_->getNodeLastLagLB() < pModel_->getParameters().absoluteGap_ - EPSILON){
		// DBG
		std::cout << "Fathom node with Largangian bound" << std::endl;
		return BCP_DoNotBranch_Fathomed;
	}

	//print the current solution
	if (pModel_->getParameters().printRelaxationSol_) {
		pModel_->getParameters().saveFunction_->printCurrentSol();
	}

	// if not currently diving with the same rule as in the heuristic,
	// try and run the heuristic
	//
	if (pModel_->getParameters().performHeuristicAfterXNode_ > -1 &&
		(!pModel_->is_columns_node() || !pModel_->getParameters().branchColumnUntilValue_) ) {
		heuristicHasBeenRun_ = false;
		generate_heuristic_solution(lpres, vars, cuts);
		heuristicHasBeenRun_ = true;
	}

	//continue the dive
	MyBranchingCandidate candidate;
	if(pModel_->is_columns_node()){
		pModel_->column_candidates(candidate);
		buildCandidate(candidate, vars, cuts, cands);
		pModel_->updateDive();
		return BCP_DoBranch;
	}

	//do we branch on columns ?
	if(pModel_->getNbDives() < pModel_->getParameters().nbDiveIfBranchOnColumns_) {
		pModel_->column_candidates(candidate);
	}

	//branching candidates: numberOfNursesByPosition_, rest on a day, ...
	bool generate = pModel_->branching_candidates(candidate);

	//after a given number of nodes since last dive, prepare to go fro a new dive
	if(pModel_->getNbDives() >= nb_dives_to_wait_before_branching_on_columns_.front()){
		nb_dives_to_wait_before_branching_on_columns_.pop_front();
		pModel_->resetNbNodesSinceLastDive();
	}

	if(!generate) {
		return BCP_DoNotBranch_Fathomed;
	}

	buildCandidate(candidate, vars, cuts, cands);

	return BCP_DoBranch;
}

void BcpLpModel::buildCandidate(const MyBranchingCandidate& candidate, const BCP_vec<BCP_var*>&  vars, const BCP_vec< BCP_cut*> &  cuts, BCP_vec<BCP_lp_branching_object*>&  cands){
   BCP_vec<BCP_var*> new_vars; //add a branching cut for the set of arcs
   BCP_vec<BCP_cut*> new_cuts; //add a branching cut for the set of arcs
   BCP_vec<int> vpos; //positions of the variables
   BCP_vec<double> vbd; // new bounds for each variable and for each children
   BCP_vec<int> cpos; //positions of the cuts
   BCP_vec<double> cbd; //bounds of the cuts

   /*
    * Branching variables
    */

   BCP_vec<double> generalVarBounds;
   int nbNewVar = 0;
   for(MyVar* var: candidate.getBranchingVars()){
      //search the var in the vars
      BcpColumn* col = dynamic_cast<BcpColumn*>(var);
      int index = var->getIndex();
      if(col){
         index = pModel_->getIndexCol(index);
         if(index == -1){ // new var
            index = vars.size()+nbNewVar;
            ++nbNewVar;
         }
      }
      vpos.push_back(index);
   }

   for(MyVar* newVar: candidate.getNewBranchingVars())
      new_vars.push_back(dynamic_cast<BCP_var*>(newVar));

   //bounds
   for(const MyBranchingNode& node: candidate.getChildren()){
      //		cout << "Node " << ++l << " -- modified bounds:" << endl;
      vector<double>::const_iterator lbIt = node.getLb().begin();
      vector<MyVar*>::const_iterator varIt = candidate.getBranchingVars().begin();
      int i = 0;
      for(vector<double>::const_iterator ubIt = node.getUb().begin(); ubIt != node.getUb().end(); ++ubIt){
         vbd.push_back(*lbIt);
         vbd.push_back(*ubIt);
         			//debug
         			// Rotation rot((*varIt)->getPattern());
         			// if(*lbIt != (*varIt)->getLB())
         			// 	cout << "Var " << vpos[i] << "(" << (*varIt)->getIndex() << ")" << ": LB=" << *lbIt << " - " << rot.toString() << endl;
         			// if(*ubIt != (*varIt)->getUB())
         			// 	cout << "Var " << vpos[i] << "(" << (*varIt)->getIndex() << ")" << ": UB=" << *ubIt << " - " << rot.toString() << endl;
         ++lbIt;
         ++varIt;
         ++i;
      }
   }

   /*
    * Branching cuts
    */

   int nbNewCut = 0;
   for(MyCons* cons: candidate.getBranchingCons()){
      //search the var in the vars
      BcpBranchCons* cut = dynamic_cast<BcpBranchCons*>(cons);
      int index = cons->getIndex();
      if(cut){
         index = cuts.size()+nbNewCut;
         ++nbNewCut;
      }
      cpos.push_back(index);
   }

   for(MyCons* newCut: candidate.getNewBranchingCons())
      new_cuts.push_back(dynamic_cast<BCP_cut*>(newCut));

   //bounds
   for(const MyBranchingNode& node: candidate.getChildren()){
      vector<double>::const_iterator lhsIt = node.getLhs().begin();
      for(vector<double>::const_iterator rhsIt = node.getRhs().begin(); rhsIt != node.getRhs().end(); ++rhsIt){
         cbd.push_back(*lhsIt);
         cbd.push_back(*rhsIt);
         ++lhsIt;
      }
   }

   cands.push_back(new  BCP_lp_branching_object(candidate.getChildren().size(), &new_vars, &new_cuts, /* vars/cuts_to_add */
      &vpos, &cpos, &vbd, &cbd, /* forced parts */
      0, 0, 0, 0 /* implied parts */));
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
   //	if(pModel_->continueDiving()) best->action()[0] = BCP_KeepChild;
   //	else pModel_->addCurrentNodeToStack();
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

void BcpLpModel::select_vars_to_delete(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts,
   const bool before_fathom,
   BCP_vec<int>& deletable)
{
	// BCP_lp_user::select_vars_to_delete(lpres, vars, cuts, before_fathom, deletable);


	if (before_fathom && getLpProblemPointer()->param(BCP_lp_par::NoCompressionAtFathom))
	return;
	const int varnum = vars.size();
	deletable.reserve(varnum);
	for (int i = getLpProblemPointer()->core->varnum(); i < varnum; ++i) {
		BCP_var *var = vars[i];
		if (var->is_to_be_removed()
		||	(! var->is_non_removable() && var->lb() == 0 && var->ub() == 0))
		{
			deletable.unchecked_push_back(i);
		}
	}

	if(pModel_->getParameters().printBranchStats_ && before_fathom){
		std::cout << "ABOUT TO FATHOM CURRENT NODE" << std::endl;
		pModel_->printStats();
	}
	// DBG
			for (unsigned int i = pModel_->getCoreVars().size(); i < vars.size(); ++i) {
				BcpColumn* var = dynamic_cast<BcpColumn*>(vars[i]);

				if (var->lb() == 0 && var->ub() == 0) {
					deletable.unchecked_push_back(i);
					continue;
				}

				int inactive_iteration = lpIteration_ - var->getLastActive();
				double activity_rate = var->getActiveCount() * 1.0 / (lpIteration_ - var->getIterationCreation());
				if(inactive_iteration > min_inactive_iteration && activity_rate < max_activity_rate)
					deletable.unchecked_push_back(i);
			}
}

/*
 * BcpBranchingTree
 */

BcpBranchingTree::BcpBranchingTree(BcpModeler* pModel):
	pModel_(pModel) , nbInitialColumnVars_(pModel->getActiveColumns().size())
{}

// setting the base
//Create the core of the problem by filling out the last three arguments.
void BcpBranchingTree::initialize_core(BCP_vec<BCP_var_core*>& vars,
   BCP_vec<BCP_cut_core*>& cuts, BCP_lp_relax*& matrix){
   // initialize tm parameters
   set_param(BCP_tm_par::TmVerb_SingleLineInfoFrequency, pModel_->getFrequency());
   //always dive
   set_param(BCP_tm_par::UnconditionalDiveProbability, 1);
   set_param(BCP_tm_par::MaxRunTime, LARGE_TIME);//pModel_->getParameters().maxSolvingTimeSeconds_);
   for(pair<BCP_tm_par::chr_params, bool> entry: pModel_->getTmParameters())
      set_param(entry.first, entry.second);

   //define nb rows and col
   const int rownum = pModel_->getCons().size();
   const int colnum = pModel_->getCoreVars().size();

   // bounds and objective
   double* lb, *ub, *obj, *rhs, *lhs;
   lb = (double*) malloc(colnum*sizeof(double));
   ub = (double*) malloc(colnum*sizeof(double));
   obj = (double*) malloc(colnum*sizeof(double));
   rhs = (double*) malloc(rownum*sizeof(double));
   lhs = (double*) malloc(rownum*sizeof(double));

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
   matrix->copyOf(pModel_->buildCoinMatrix(), obj, lb, ub, lhs, rhs);

   if (lb) free(lb);
   if (ub) free(ub);
   if (obj) free(obj);
   if (lhs) free(lhs);
   if (rhs) free(rhs);
}

// create the root node
//Create the set of extra variables and cuts that should be added
//to the formulation in the root node.
void BcpBranchingTree::create_root(BCP_vec<BCP_var*>& added_vars,
   BCP_vec<BCP_cut*>& added_cuts,
   BCP_user_data*& user_data){

   added_vars.reserve(pModel_->getActiveColumns().size());
   for(MyVar* col: pModel_->getActiveColumns()){
		 BcpColumn* var = dynamic_cast<BcpColumn*>(col);
		 if(!var)
			 Tools::throwError("Bad variable casting.");
		 // add BcpColumn which will be deleted by BCP - the variables in activeColumns need to be owned by no other model
		 added_vars.unchecked_push_back(var);
   }
	pModel_->clearActiveColumns();
}

void BcpBranchingTree::display_feasible_solution(const BCP_solution* sol){
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


//-----------------------------------------------------------------------------
//
//  C l a s s   B c p M o d e l e r
//
// Specific class of modeler designed to use BCP.
// In particular, may methods defined in the BCP library need to be implemented
// here
//
//-----------------------------------------------------------------------------


BcpModeler::BcpModeler(MasterProblem* pMaster, const char* name, LPSolverType type):
CoinModeler(), pMaster_(pMaster), pBcp_(0), primalValues_(0), dualValues_(0), reducedCosts_(0), lhsValues_(0),
lastNbSubProblemsSolved_(0), lastMinDualCost_(0), nbNodes_(0), LPSolverType_(type) {
  pBcp_ = new BcpInitialize(this);
}

//destroy all the column in the solutions
BcpModeler::~BcpModeler() {
   for(BCP_solution_generic& sol: bcpSolutions_){
      int size = sol._vars.size();
      for(int i=coreVars_.size(); i<size; ++i) delete sol._vars[i];
   }
   for(MyCons* cons: branchingCons_)
      if(cons){
         delete cons;
         cons = 0;
      }
   branchingCons_.clear();
   for(MyVar* var: columnsInSolutions_)
      if(var){
         delete var;
         var = 0;
      }
  delete  pBcp_;
}

//solve the model
int BcpModeler::solve(bool relaxation){
   //create the root
   pTree_->pushBackNewNode();
   
   //solve
   char** argv=NULL;

   int value=-1;
	try{
		value = bcp_main(0, argv, pBcp_);
		getMaster()->setStatus(OPTIMAL);
	}catch(OptimalStop& e) {
		getMaster()->setStatus(OPTIMAL);
	}catch(FeasibleStop& e) {
		getMaster()->setStatus(FEASIBLE);
	}catch(InfeasibleStop& e) {
		getMaster()->setStatus(INFEASIBLE);
	}catch(TimeoutStop& e) {
		getMaster()->setStatus(TIME_LIMIT);
	}

	// retrieve the statistics of the solution
	timeStats_.add(pBcp_->getBcpLpModel()->getLpProblemPointer()->stat);
	nbLpIterations_ = pBcp_->getBcpLpModel()->getNbLpIterations();

	static int  cpt = 0;
	char  nom[1024];
	sprintf(nom, "ccc%03d", cpt++);
pBcp_->getBcpLpModel()->getLpProblemPointer()->lp_solver->writeLp(nom);

   // clear tree
   pTree_->clear();
   treeMapping_.clear();
   return value;
}

//reinitialize all parameters and clear vectors
void BcpModeler::reset(bool rollingHorizon) {
   pTree_->reset();
   lastNbSubProblemsSolved_=0;
   lastMinDualCost_=0;
   solHasChanged_ = false;

   obj_history_.clear();
   primalValues_.clear();
   dualValues_.clear();
   reducedCosts_.clear();
   lhsValues_.clear();

	// delete the best solutions that were if solving with rolling horizon
   if (rollingHorizon) {
		unsigned int index = getBestSolIndex();
		for(unsigned int ind = 0; ind < bcpSolutions_.size(); ind++){
			if (ind == index) continue;
			BCP_solution_generic sol = bcpSolutions_[ind];
			int size = sol._vars.size();
			for(int i=coreVars_.size(); i<size; ++i) delete sol._vars[i];
		}
      bcpSolutions_.clear();
   }

   //create the root
   pTree_->pushBackNewNode();
}



/*
 * Create core variable:
 *    var is a pointer to the pointer of the variable
 *    var_name is the name of the variable
 *    lb, ub are the lower and upper bound of the variable
 *    vartype is the type of the variable: SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY
 */
int BcpModeler::createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub, const vector<double>& pattern){
   *var = new BcpCoreVar(var_name, index, objCoeff, vartype, lb, ub, pattern);
   objects_.push_back(*var);
   return 1;
}

int BcpModeler::createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, const vector<double>& pattern, double dualObj, VarType vartype, double lb, double ub){
   *var = new BcpColumn(var_name, index, objCoeff, pattern, dualObj, vartype, lb, ub);
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
   return 0;
}

void BcpModeler::createCutLinear(MyCons** cons, const char* con_name, double lhs, double rhs,
   vector<MyVar*> vars, vector<double> coeffs){
   vector<int> indexes;
   for(MyVar* var: vars)
      indexes.push_back(var->getIndex());
   BcpBranchCons* cons2 = new BcpBranchCons(con_name, cons_.size()+branchingCons_.size(),	lhs, rhs, indexes, coeffs);
   branchingCons_.push_back(cons2);
   *cons = new BcpBranchCons(*cons2);
}

/*
 * Set the solution: this is where we set the active column variables
 */

 void BcpModeler::setLPSol(const BCP_lp_result& lpres, const BCP_vec<BCP_var*>&  vars, int lpIteration){
	 obj_history_.push_back(lpres.objval());
	 solHasChanged_ = false;

	 //copy the new arrays in the vectors for the core vars
	 const int nbCoreVar = coreVars_.size();
	 const int nbVar = vars.size();
	 const int nbCons = cons_.size();

	 //clear all
	 primalValues_.clear();
	 dualValues_.clear();
	 reducedCosts_.clear();
	 lhsValues_.clear();

	 //assign value for core variables
	 primalValues_.assign(lpres.x(), lpres.x()+nbVar);
	 dualValues_.assign(lpres.pi(), lpres.pi()+nbCons);
	 reducedCosts_.assign(lpres.dj(), lpres.dj()+nbVar);
	 lhsValues_.assign(lpres.lhs(), lpres.lhs()+nbCons);

	clearActiveColumns();
	 for(int i=nbCoreVar; i<nbVar; ++i){
		 BcpColumn* var = dynamic_cast<BcpColumn*>(vars[i]);
		 if(primalValues_[i] > EPSILON)
		 var->addActiveIteration(lpIteration); //update the different counters
		 addActiveColumn(var, i);
	 }
	 //	//debug
	 //	checkActiveColumns(vars);
 }

void BcpModeler::checkActiveColumns(const BCP_vec<BCP_var*>&  vars){
   ShiftNode* shiftNode = dynamic_cast<ShiftNode*>(pTree_->getCurrentNode());
   if(shiftNode == 0) return;

   for(unsigned int i=coreVars_.size(); i<vars.size(); ++i){
      BcpColumn* var = dynamic_cast<BcpColumn*>(vars[i]);
      Rotation rot(var->getPattern());
      if(var->getUB() == 0 || var->ub() == 0 || shiftNode->pNurse_->id_ != rot.nurseId_) continue;
      for(pair<int,int> p: rot.shifts_)
         if(shiftNode->day_ == p.first &&
            find(shiftNode->forbiddenShifts_.begin(), shiftNode->forbiddenShifts_.end(), p.second) != shiftNode->forbiddenShifts_.end()){
            cout << "problem: active column " << var->bcpind() << " with forbidden shift " << p.second << endl;
            cout << rot.toString() << endl;
            getchar();
         }
   }
}

void BcpModeler::addBcpSol(const BCP_solution* sol){
   //if no integer solution is needed, don't store the solutions
   if(parameters_.stopAfterXSolution_ == 0)
      return;

   //create a solution which is not going to delete the vars at the end (argument=false)
   BCP_solution_generic mySol(false);
   BCP_solution_generic* sol2 = (BCP_solution_generic*) sol;

	bool isArtificialSol=false;
   for(unsigned int i=0; i<sol2->_vars.size(); ++i){
      BcpColumn* col = dynamic_cast<BcpColumn*>(sol2->_vars[i]);
      if(col){
			// if ( (col->is_integer() && sol2->_values[i] < 1-EPSILON) || (sol2->_values[i] < EPSILON) )  {
			// 	std::cout << "Var type continuous = " << col->getVarType() << std::endl;
			// 	std::cout << "Column " << i << " has value " << sol2->_values[i] << std::endl;
			// 	continue;
			// }
         col = new BcpColumn(*col);
         mySol.add_entry(col, sol2->_values[i]);
         columnsInSolutions_.push_back(col);
      }
      else {
			mySol.add_entry((BcpCoreVar*)coreVars_[sol2->_vars[i]->bcpind()], sol2->_values[i]);
			if (((BcpCoreVar*)coreVars_[sol2->_vars[i]->bcpind()])->getCost()  > 1000.0) {
				if (sol2->_values[i] > EPSILON) isArtificialSol = true;
			}
		}
   }

	if (!isArtificialSol) {
   	bcpSolutions_.push_back(mySol);
	}

   // if the solution improves the upper bound, record the new upper bound and load the integer solution
   if(pTree_->getBestUB() > sol->objective_value() + EPSILON){
      pTree_->setBestUB(sol->objective_value());

      // print the solution in a text file
      if(parameters_.printEverySolution_){
         loadBcpSol(bcpSolutions_.size()-1);
         solHasChanged_ = true;
         parameters_.saveFunction_->save(parameters_.weekIndices_, parameters_.outfile_);
      }
   }

	if(pTree_->getBestLB()<1.0e5 && pTree_->getBestUB() - pTree_->getBestLB() < parameters_.absoluteGap_ - EPSILON){
		char error[100];
		sprintf(error, "Stopped: absolute gap < %.2f.", parameters_.absoluteGap_);
		throw OptimalStop(error);
	}
}

// Get the index of the best solution in the vector of solutions of BCP
//
int BcpModeler::getBestSolIndex() {
   int i=0, index = -1;
   double bestObj = LARGE_SCORE;
   for(BCP_solution_generic& sol: bcpSolutions_){
      if(sol.objective_value() < bestObj){
         bestObj = sol.objective_value();
         index = i;
      }
      ++i;
   }
   return index;
}

bool BcpModeler::loadBestSol(){
   int index = getBestSolIndex();

   if(index == -1)
      return false;

	this->loadInputSol(bcpSolutions_[index]);

   return true;
}

void BcpModeler::loadBcpSol(int index){
	BCP_solution_generic& sol = bcpSolutions_[index];

	this->loadInputSol(sol);
}

// Clear the active column, set the active columns with those in the input
// solution and set the primal values accordingly
//
void BcpModeler::loadInputSol(BCP_solution_generic& sol){
	const int size = sol._vars.size(), core_size = coreVars_.size();
	vector<double> primal(core_size);

	int ind = core_size;
	clearActiveColumns();

	for(int i=0; i<size; ++i){
		 // type of sol._vars[i] is either BcpColumn or BcpCoreVar, dynamic_cast will get the proper one
		 BcpColumn* col = dynamic_cast<BcpColumn*>(sol._vars[i]);
		 if(col){
			  addActiveColumn(col, ind++);
				primal.push_back(sol._values[i]);
		 }
		 else{
				BcpCoreVar* var = dynamic_cast<BcpCoreVar*>(sol._vars[i]);
				primal[var->getIndex()] = sol._values[i];
		 }
	}

	setPrimal(primal);
}

// DBG
// Set the value of the active columns with those in the best solution
void BcpModeler::setActiveColumnsValuesWithBestSol() {
   int index = getBestSolIndex();

   BCP_solution_generic& sol = bcpSolutions_[index];
   for(unsigned int i = 0; i < primalValues_.size(); i++){
      primalValues_[i] = 0.0;
   }

   for (unsigned int i=0; i < sol._vars.size(); i++) {
      // type of sol._vars[i] is either BcpColumn or BcpCoreVar, dynamic_cast will get the proper one
      BcpColumn* col = dynamic_cast<BcpColumn*>(sol._vars[i]);
      if (col) {
         setVarValue(col,sol._values[i]);
      }
      else{
         BcpCoreVar* var = dynamic_cast<BcpCoreVar*>(sol._vars[i]);
         setVarValue(var,sol._values[i]);
      }
   }
}

// Set every rotation to one : this is useful only when the active columns
// are only the rotations included in a provided initial solution
//
void BcpModeler::setEveryRotationToOne() {
	for (MyVar* pVar:activeColumnVars_) {
		pVar->setLB(1.0);
		pVar->setUB(1.0);
	}
}

/*
 * fix/unfix all the rotations variables starting from the input vector of days
 */
void BcpModeler::fixRotationsStartingFromDays(vector<bool> isFixDay) {

   // get the best solution currently in BCP
   int index = getBestSolIndex();
   BCP_solution_generic& sol = bcpSolutions_[index];

   const int size = sol._vars.size();

   for(int i=0; i<size; ++i){
      // type of sol._vars[i] is either BcpColumn or BcpCoreVar, dynamic_cast will get the proper one
      BcpColumn* col = dynamic_cast<BcpColumn*>(sol._vars[i]);
      if(col){
         double value = sol._values[i];
         if (isFixDay[col->getFirstDay()]) {
            if(value >= EPSILON && value <= 1 - EPSILON)
               Tools::throwError("BcpModeler::fixRotationsStartingFromDays: the value to be fixed is not integer!");
            col->setLB(1.0);
            col->setUB(1.0);
         }
      }
   }
}
void BcpModeler::unfixRotationsStartingFromDays(vector<bool> isUnfixDay) {
   for(MyVar* var: activeColumnVars_){
      if(isUnfixDay[var->getFirstDay()]) {
         var->setLB(0.0);
         switch(var->getVarType()){
         case VARTYPE_BINARY:
            var->setUB(1.0);
            break;
         default:
            var->setUB(infinity_);
            break;
         }
      }
   }
}

// fix/unfix all the rotations variables of the input nurses
void BcpModeler::fixRotationsOfNurses(vector<bool> isFixNurse) {

   // get the best solution currently in BCP

	if (!bcpSolutions_.empty()) {
	   int index = getBestSolIndex();
	   BCP_solution_generic& sol = bcpSolutions_[index];

	   const int size = sol._vars.size();

	   for(int i=0; i<size; ++i){
	      // type of sol._vars[i] is either BcpColumn or BcpCoreVar, dynamic_cast will get the proper one
	      BcpColumn* col = dynamic_cast<BcpColumn*>(sol._vars[i]);
	      if(col){
	         double value = sol._values[i];
	         if (isFixNurse[col->getNurseId()]) {
	            if(value >= EPSILON && value <= 1 - EPSILON)
	               Tools::throwError("BcpModeler::fixRotationsStartingFromDays: the value to be fixed is not integer!");
	            col->setLB(value);
	            col->setUB(value);
	         }
	      }
	   }
	}
	// DBG: WARNING THIS IMPLEMENTATION IS RISKY, IT ASSUMES THAT THERE IS NO
	// SOLUTION AVAILABLE ONLY WHEN THE SOLUTION HAS JUST BEEN LOADED
	else {
		for(MyVar* var: activeColumnVars_){
	      if(isFixNurse[var->getNurseId()]) {
	         var->setLB(1.0);
				var->setUB(1.0);
	      }
	   }
	}

}
void BcpModeler::unfixRotationsOfNurses(vector<bool> isUnfixNurse) {
   for(MyVar* var: activeColumnVars_){
      if(isUnfixNurse[var->getNurseId()]) {
         var->setLB(0.0);
         switch(var->getVarType()){
         case VARTYPE_BINARY:
            var->setUB(1.0);
            break;
         default:
            var->setUB(infinity_);
            break;
         }
      }
   }
}

/*
 * relax/unrelax all the rotations variables starting from the input vector of days
 */
void BcpModeler::relaxRotationsStartingFromDays(vector<bool> isRelaxDay) {
   for(MyVar* var: activeColumnVars_){
      if(isRelaxDay[var->getFirstDay()]) {
         var->setVarType(VARTYPE_CONTINUOUS);
      }
   }
}
void BcpModeler::unrelaxRotationsStartingFromDays(vector<bool> isUnrelaxDay) {
   for(MyVar* var: activeColumnVars_){
      if(isUnrelaxDay[var->getFirstDay()]) {
         var->setVarType(VARTYPE_INTEGER);
      }
   }
}


/*
 * get/set the primal values
 */

double BcpModeler::getVarValue(MyVar* var){
   if(primalValues_.size() ==0 )
      Tools::throwError("Primal solution has not been initialized.");

   unsigned int index = ((CoinVar*) var)->getIndex();
   //if a column, fetch index
   if(index>=coreVars_.size()){
      //if the column is not active, return 0
      if(columnsToIndex_.count(index)==0) return 0;
      index = columnsToIndex_[index];
   }
   return primalValues_[index];
}

void BcpModeler::setVarValue(MyVar* var, double value){
   if(primalValues_.size() ==0 )
      Tools::throwError("Primal solution has not been initialized.");

   unsigned int index = ((CoinVar*) var)->getIndex();
   //if a column, fetch index
   if(index>=coreVars_.size()){
      //if the column is not active, return 0
      if(columnsToIndex_.count(index)==0) return;
      index = columnsToIndex_[index];
   }
   primalValues_[index] = value;
}

/*
 * Get the dual variables
 */

double BcpModeler::getDual(MyCons* cons, bool transformed){
   CoinCons* cons2 = (CoinCons*) cons;
   if(dualValues_.size() == 0)
      Tools::throwError("Dual solution has been initialized.");
   return dualValues_[cons2->getIndex()];
}

/*
 * Get the reduced cost
 */

double BcpModeler::getReducedCost(MyVar* var){
   if(reducedCosts_.size() == 0)
      Tools::throwError("Reduced cost solution has been initialized.");

   unsigned int index = ((CoinVar*) var)->getIndex();
   //if a column, fetch index
   if(index>=coreVars_.size()){
      //if the column is not active, return 0
      if(columnsToIndex_.count(index)==0) return 0;
      index = columnsToIndex_[index];
   }
   return reducedCosts_[index];
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
      tm_parameters[BCP_tm_par::TmVerb_NewPhaseStart] = 0;
      tm_parameters[BCP_tm_par::TmVerb_Last] = 1;
		tm_parameters[BCP_tm_par::DebugLpProcesses] =1;

      lp_parameters[BCP_lp_par::LpVerb_Last] = 1; // Just a marker for the last LpVerb
		lp_parameters[BCP_lp_par::LpVerb_FathomInfo] = 1; // Print information related to fathoming. (BCP_lp_main_loop, BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
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
      lp_parameters[BCP_lp_par::LpVerb_LpSolutionValue] = 1; // Print the size of the problem matrix and the LP solution value after resolving the LP. (BCP_lp_main_loop)
   }

   if(v>=4){
      lp_parameters[BCP_lp_par::LpVerb_RelaxedSolution] = 1; // Turn on the user hook "display_lp_solution". (BCP_lp_main_loop)
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
   if(pTree_->getBestLB() >= LARGE_SCORE)
      return false;

   //check relative gap
	if(pTree_->getBestUB() - pTree_->getBestLB() < parameters_.absoluteGap_ - EPSILON){
      char error[100];
      sprintf(error, "Stopped: absolute gap < %.2f.", parameters_.absoluteGap_);
      throw OptimalStop(error);
   }
	else if (getMaster()->getTimerTotal()->dSinceStart() > getParameters().maxSolvingTimeSeconds_ ) {
		char error[100];
		std::cout << "Total time spent solving the problem " << getMaster()->getTimerTotal()->dSinceStart() << std::endl;
      sprintf(error, "Stopped: Time has run out.");
      throw TimeoutStop(error);
	}
	else if(parameters_.solveToOptimality_) {
	      return false;
	}

	//check the number of solutions
	if(nbSolutions() >= parameters_.stopAfterXSolution_){
		char error[100];
		sprintf(error, "Stopped: %d solutions have been founded", nbSolutions());
		throw FeasibleStop(error);
	}
	else if(pTree_->getBestUB() - pTree_->getBestLB() < parameters_.minRelativeGap_ * pTree_->getBestLB() - EPSILON){
		char error[100];
		sprintf(error, "Stopped: relative gap < %.2f.", parameters_.minRelativeGap_);
		throw FeasibleStop(error);
	}
   else if(pTree_->getBestUB() - pTree_->getBestLB() < parameters_.relativeGap_ * pTree_->getBestLB() - EPSILON){
      //if the relative gap is small enough and if same incumbent since the last dive, stop
      if(pTree_->getNbNodesSinceLastIncumbent() > parameters_.nbDiveIfMinGap_*pTree_->getDiveLength()){
         char error[100];
         sprintf(error, "Stopped: relative gap < %.2f and more than %d nodes without new incumbent.",
            parameters_.relativeGap_, parameters_.nbDiveIfMinGap_*pTree_->getDiveLength());
         throw FeasibleStop(error);
      }
   }
   else if (nbSolutions() > 0) {
	//if(pTree_->getBestUB() - pTree_->getBestLB() < 10.0 * pTree_->getBestLB()) {
      //if the relative gap is too big, wait 2 dives before stopping
      if(pTree_->getNbNodesSinceLastIncumbent() > parameters_.nbDiveIfRelGap_*pTree_->getDiveLength()){
         char error[100];
         sprintf(error, "Stopped: relative gap > %.2f and more than %d nodes without new incumbent.",
            parameters_.relativeGap_, parameters_.nbDiveIfRelGap_*pTree_->getDiveLength());
         throw FeasibleStop(error);
      }
	}

   return false;
}


/**************
 * Outputs *
 *************/

int BcpModeler::writeProblem(string fileName){
  // pBcp_->getBcpLpModel()->getLpProblemPointer()->lp_solver->writeLp(fileName.c_str());
  return 0;
}

int BcpModeler::writeLP(string fileName){
   //   OsiClpSolverInterface solver = ;
   //   solver.writeLp(fileName.c_str(), "lp", 1e-5, 10, 5);
   return 0;
}
