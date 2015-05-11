/*
 * ScipModeler.h
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#ifndef SRC_SCIPMODELER_H_
#define SRC_SCIPMODELER_H_

#include "CoinModeler.h"

/* BCP includes */
#include "BCP_enum.hpp"
#include "BCP_vector.hpp"
#include "BCP_var.hpp"
#include "BCP_cut.hpp"
#include "BCP_buffer.hpp"
#include "BCP_tm_user.hpp"
#include "BCP_lp_user.hpp"
#include "BCP_USER.hpp"
#include "BCP_solution.hpp"
#include "OsiClpSolverInterface.hpp"

/*
 * My Variables
 */

//these variables are not generated, they are always in the LP problem
struct BcpCoreVar: public CoinVar, public BCP_var_core {
   static BCP_var_t getBcpVarType(VarType vartype){
      BCP_var_t type;
      switch(vartype){
      case VARTYPE_BINARY:
         type = BCP_BinaryVar;
         break;
      case VARTYPE_INTEGER:
         type = BCP_IntegerVar;
         break;
      default:
         type = BCP_ContinuousVar;
         break;
      }
      return type;
   }

   BcpCoreVar(const char* name, int index, double cost, VarType type, double lb, double ub):
      CoinVar(name, index, cost, type, lb, ub),
      BCP_var_core(getBcpVarType(type), cost, lb, ub)
   {
      set_bcpind(index);
   }

   BcpCoreVar(const BcpCoreVar& var) :
      CoinVar(var), BCP_var_core(getBcpVarType(type_), cost_, lb_, ub_)
   {
      set_bcpind(var.index_);
   }

   ~BcpCoreVar(){ }
};

//these variables are generated during the process, they are not always in the LP problem
struct BcpColumn: public CoinVar, public BCP_var_algo{
   BcpColumn(const char* name, int index, double cost, double dualCost, VarType type, double lb, double ub):
      CoinVar(name, index, cost, type, lb, ub, dualCost),
      BCP_var_algo(BcpCoreVar::getBcpVarType(type), cost, lb, ub)
   {
      set_bcpind(index);
   }

   BcpColumn(const BcpColumn& var) :
      CoinVar(var), BCP_var_algo(BcpCoreVar::getBcpVarType(type_), cost_, lb_, ub_)
   {
      set_bcpind(var.index_);
   }

   ~BcpColumn(){ }
};

/*
 * My Constraints
 */
//these constraints are not generated, they are always in the LP problem
struct BcpCoreCons: public CoinCons, public BCP_cut_core{
   BcpCoreCons(const char* name, int index, double lhs, double rhs):
      CoinCons(name, index, lhs, rhs),
      BCP_cut_core(lhs, rhs)
   { }

   BcpCoreCons(const BcpCoreCons& cons) :
      CoinCons(cons), BCP_cut_core(lhs_, rhs_)
   { }

   ~BcpCoreCons(){ }
};

//these constraints are not generated, they are always in the LP problem
struct BcpBranchCons: public CoinCons, public BCP_cut_algo{
   BcpBranchCons(const char* name, int index, double lhs, double rhs,
      vector<int>& indexCols, vector<double>& coeffs):
      CoinCons(name, index, lhs, rhs),
      BCP_cut_algo(lhs, rhs),
      indexCols_(indexCols), coeffs_(coeffs)
   { }

   BcpBranchCons(const BcpBranchCons& cons) :
      CoinCons(cons), BCP_cut_algo(lhs_, rhs_), indexCols_(cons.indexCols_), coeffs_(cons.coeffs_)
   { }

   ~BcpBranchCons(){ }

   vector<int>& getIndexCols() { return indexCols_; }

   vector<double>& getCoeffCols() { return coeffs_; }

protected:
vector<int> indexCols_; //index of the cols of the matrix where the col has non-zero coefficient
vector<double> coeffs_; //value of these coefficients
};

/*
 * My Pricer
 */
//struct BcpPricer: public MyObject  {
//public:
//   BcpPricer(ObjPricer* pricer):MyObject(){ pricer_=pricer; }
//   ~BcpPricer(){ }
//   ObjPricer* pricer_;
//};

/*
 * My Branching Rule
 */
//struct ScipRule: public MyRule  {
//   ScipRule(){ }
//   ~ScipRule(){ }
//   SCIP_VAR* rule_;
//};

struct BcpNode{

   BcpNode(): index_(0), bestLB_(DBL_MAX), pParent_(0), pNurse_(0), day_(0), rest_(false), pNumberOfNurses_(0), lb_(-DBL_MAX), ub_(DBL_MAX) {}
   BcpNode(int index, BcpNode* pParent, vector<MyObject*>& columns):
      index_(index), bestLB_(pParent->bestLB_), pParent_(pParent),
      columns_(columns), pNurse_(0), day_(0), rest_(false),
      pNumberOfNurses_(0), lb_(-DBL_MAX), ub_(DBL_MAX) {}
   BcpNode(int index, BcpNode* pParent, Nurse* pNurse, int day, bool rest):
      index_(index), bestLB_(pParent->bestLB_), pParent_(pParent),
      pNurse_(pNurse), day_(day), rest_(rest),
      pNumberOfNurses_(0), lb_(-DBL_MAX), ub_(DBL_MAX) {}
   BcpNode(int index, BcpNode* pParent, CoinVar* var, double lb, double ub):
      index_(index), bestLB_(pParent->bestLB_), pParent_(pParent),
      pNurse_(0), day_(0), rest_(false),
      pNumberOfNurses_(var), lb_(lb), ub_(ub) {}
   virtual ~BcpNode() {}

   const int index_;
   double bestLB_;

   //parent
   BcpNode* pParent_;

   //children
   vector<BcpNode*> children_;

   //vector of the columns on which we have branched. Can be empty
   const vector<MyObject*> columns_;

   //nurse and day on which we have branched for rest or work. pNurse_ can be 0
   const Nurse* pNurse_;
   const int day_;
   const bool rest_;

   //number of nurse on which we have branched. pNumberOfNurses_ can be 0
   CoinVar* pNumberOfNurses_;
   double lb_, ub_;
};

class BcpModeler: public CoinModeler {
public:
   BcpModeler(const char* name);
   ~BcpModeler() {
      for(BcpNode* node: tree_) delete node;
      for(BcpBranchCons* cons: branchingCons_) delete cons;
   }

   //solve the model
   int solve(bool relaxation = false);

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
    */
   int createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub);

   int createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, double dualObj, VarType vartype, double lb, double ub);

   /*
    * Create linear constraint:
    *    con is a pointer to the pointer of the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    */

   int createCoinConsLinear(CoinCons** con, const char* con_name, int index, double lhs, double rhs);

   /*
    * Get the primal value
    */

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

   double getObjective(){ return best_ub; }

   double getRelaxedObjective() { return best_lb_in_root; }

   int printStats();

   int writeProblem(string fileName);

   int writeLP(string fileName);

   /*
    * Class own methods and parameters
    */

   void setLPSol(const BCP_lp_result& lpres, const BCP_vec<BCP_var*>&  vars);

   void addBcpSol(const BCP_solution* sol);

   void loadBestSol();

   void loadBcpSol(int index);

   inline void setPrimal(vector<double> primal){ primalValues_ = primal; }

   inline void setBestLb(double bestLBRoot){ best_lb_in_root = bestLBRoot; }

   inline double getBestLb(){ return best_lb_in_root; }

   inline int getFrequency() { return TmVerb_SingleLineInfoFrequency; }

   inline void setLastNbSubProblemsSolved(int lastNbSubProblemsSolved){ lastNbSubProblemsSolved_ = lastNbSubProblemsSolved; }

   inline int getLastNbSubProblemsSolved(){ return lastNbSubProblemsSolved_; }

   inline void setLastMinDualCost(double lastMinDualCost){ lastMinDualCost_ = lastMinDualCost; }

   inline double getLastMinDualCost(){ return lastMinDualCost_; }

   inline double getLastObj(){ return obj_history_[obj_history_.size()-1]; }

   inline void updateNodeLB(double lb){ currentNode_->bestLB_ = lb; }

   inline void pushBackNewNode(){
         tree_.push_back(new BcpNode);
      }

   inline void pushBackNewNode(CoinVar* var, double lb, double ub){
      BcpNode* node = new BcpNode(tree_.size(), currentNode_, var, lb, ub);
      tree_.push_back(node);
      currentNode_->children_.push_back(node);
   }

   inline void pushBackNewNode(Nurse* pNurse, int day, bool rest){
      BcpNode* node = new BcpNode(tree_.size(), currentNode_, pNurse, day, rest);
      tree_.push_back(node);
      currentNode_->children_.push_back(node);
   }

   inline void pushBackNewNode(vector<MyObject*>& columns){
      BcpNode* node = new BcpNode(tree_.size(), currentNode_, columns);
      tree_.push_back(node);
      currentNode_->children_.push_back(node);
   }

   inline  void addForbidenShifts(LiveNurse* pNurse, set<pair<int,int> >& forbidenShifts) {
      BcpNode* node = currentNode_;
      while(node->pParent_){
         if(node->pNurse_ == pNurse && node->rest_)
            for(int i=1; i<pNurse->pScenario_->nbShifts_; ++i)
               forbidenShifts.insert(pair<int,int>(node->day_, i));
         node = node->pParent_;
      }
   }

   inline void setCurrentNode(int index) { currentNode_ = tree_[index]; }

   inline void pushBackBranchingCons(BcpBranchCons* cons){ branchingCons_.push_back(cons); }

   inline vector<BcpBranchCons*>& getBranchingCons(){ return branchingCons_; }

   inline  map<BCP_tm_par::chr_params, bool>& getTmParameters(){ return tm_parameters; }

   inline  map<BCP_lp_par::chr_params, bool>& getLpParameters(){ return lp_parameters; }

protected:
   //branching tree
   vector<BcpNode*> tree_;
   //current node
   BcpNode* currentNode_;
   //best lb in root
   double best_lb_in_root;
   //results
   vector<double> obj_history_;
   vector<double> primalValues_, dualValues_, reducedCosts_, lhsValues_;
   //bcp branching cons
   vector<BcpBranchCons*> branchingCons_;
   //bcp solution
   vector<BCP_solution_generic> bcpSolutions_;

   /* stats */
   //number of sub problems solved on the last iteration of column generation
   int lastNbSubProblemsSolved_;
   //min dual cost for a rotation on the last iteration of column generation
   double lastMinDualCost_;

   /* Parameters */
   //At every this many search tree node provide a single line info on the progress of the search tree.
   //If <= 0 then never.
   //Default: 0.
   int TmVerb_SingleLineInfoFrequency = 0;

   /* Tree Manager verbosity parameters */
   map<BCP_tm_par::chr_params, bool> tm_parameters = {
      { BCP_tm_par::VerbosityShutUp, 0},
      { BCP_tm_par::TmVerb_First, 0},
      { BCP_tm_par::TmVerb_AllFeasibleSolutionValue, 0},
      { BCP_tm_par::TmVerb_AllFeasibleSolution, 0},
      { BCP_tm_par::TmVerb_BetterFeasibleSolutionValue, 0},
      { BCP_tm_par::TmVerb_BetterFeasibleSolution, 1}, //need this method to store every better feasible solution
      { BCP_tm_par::TmVerb_BestFeasibleSolution, 1},
      { BCP_tm_par::TmVerb_NewPhaseStart, 0},
      { BCP_tm_par::TmVerb_PrunedNodeInfo, 0},
      { BCP_tm_par::TmVerb_TimeOfImprovingSolution, 0},
      { BCP_tm_par::TmVerb_TrimmedNum, 0},
      { BCP_tm_par::TmVerb_FinalStatistics, 1}, //need this method to store the best feasible solution
      { BCP_tm_par::ReportWhenDefaultIsExecuted, 0},
      { BCP_tm_par::TmVerb_Last, 0}
   };

   /* LP verbosity parameters */
   map<BCP_lp_par::chr_params, bool> lp_parameters = {
      { BCP_lp_par::ReportWhenDefaultIsExecuted, 0},// Print out a message when the default version of an overridable method is executed.
      { BCP_lp_par::LpVerb_AddedCutCount, 0},// Print the number of cuts added from the local cut pool in the current iteration. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_CutsToCutPoolCount, 0},// Print the number of cuts sent from the LP to the cut pool. (BCP_lp_send_cuts_to_cp)
      { BCP_lp_par::LpVerb_ReportLocalCutPoolSize, 0},// Print the current number of cuts in the cut pool. This number is printed several times: before and after generating columns at the current iteration, after removing non-essential cuts, etc. (BCP_lp_generate_cuts)
      { BCP_lp_par::LpVerb_ReportCutGenTimeout, 0},// Print information if receiving cuts is timed out. (BCP_lp_generate_cuts)
      { BCP_lp_par::LpVerb_GeneratedCutCount, 0},// Print the number of cuts generated during this iteration (since the LP was resolved last time). (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_AddedVarCount, 0},// Print the number of variables added from the local variable pool in the curent iteration. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_ChildrenInfo, 0},// After a branching object is selected print what happens to the presolved children (e.g., fathomed). (BCP_print_brobj_stat)
      { BCP_lp_par::LpVerb_ColumnGenerationInfo, 0},// Print the number of variables generated before resolving the Lp ir fathoming a node. (BCP_lp_fathom)
      { BCP_lp_par::LpVerb_FathomInfo, 0},// Print information related to fathoming. (BCP_lp_main_loop, BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
      { BCP_lp_par::LpVerb_IterationCount, 0},// Print the "Starting iteration x" line. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_RelaxedSolution, 0},// Turn on the user hook "display_lp_solution". (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_FinalRelaxedSolution, 0},// Turn on the user hook "display_lp_solution" for the last LP relaxation solved at a search tree node. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_LpSolutionValue, 0},// Print the size of the problem matrix and the LP solution value after resolving the LP. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_MatrixCompression, 0},// Print the number of columns and rows that were deleted during matrix compression. (BCP_lp_delete_cols_and_rows)
      { BCP_lp_par::LpVerb_PresolvePositions, 0},// Print detailed information about all the branching candidates during strong branching. LpVerb_PresolveResult must be set for this parameter to have an effect. (BCP_lp_perform_strong_branching)
      { BCP_lp_par::LpVerb_PresolveResult, 0},// Print information on the presolved branching candidates during strong branching. (BCP_lp_perform_strong_branching)
      { BCP_lp_par::LpVerb_ProcessedNodeIndex, 0},// Print the "Processing NODE x on LEVEL y" line. (BCP_lp-main_loop)
      { BCP_lp_par::LpVerb_ReportVarGenTimeout, 0},// Print information if receiving variables is timed out. (BCP_lp_generate_vars)
      { BCP_lp_par::LpVerb_ReportLocalVarPoolSize, 0},// Similar as above for variables. (BCP_lp_generate_vars)
      { BCP_lp_par::LpVerb_VarTightening, 0},// Print the number of variables whose bounds have been changed by reduced cost fixing or logical fixing. (BCP_lp_fix_vars)
      { BCP_lp_par::LpVerb_RowEffectivenessCount, 0},// Print the number of ineffective rows in the current problem. The definition of what rows are considered ineffective is determined by the paramter IneffectiveConstraints. (BCP_lp_adjust_row_effectiveness)
      { BCP_lp_par::LpVerb_StrongBranchPositions, 0},// Print detailed information on the branching candidate selected by strong branching. LpVerb_StrongBranchResult must be set fo this parameter to have an effect. (BCP_print_brobj_stat)
      { BCP_lp_par::LpVerb_StrongBranchResult, 0},// Print information on the branching candidate selected by strong branching. (BCP_print_brobj_stat)
      { BCP_lp_par::LpVerb_GeneratedVarCount, 0},// Print the number of variables generated during this iteration. (BCP_lp_main_loop)
      { BCP_lp_par::LpVerb_Last, 0}// Just a marker for the last LpVerb
   };
};

/*
 * The class inherits the BCP_lp_user class from which the user can derive a problem specific class
 * to be used in the LP process.
 * In that derived class the user can store data to be used in the methods she overrides.
 * Also that is the object the user must return in the USER_initialize::lp_init() method.
 *
 * There are two kind of methods in the class.
 * The non-virtual methods are helper functions for the built-in defaults, but the user can use them as well.
 * The virtual methods execute steps in the BCP algorithm where the user might want to override the default behavior.
 */
class BcpLpModel: public BCP_lp_user {
public:
   BcpLpModel(BcpModeler* pModel);
   ~BcpLpModel() { }

   /*
    * BCP_lp_user methods
    */

   void unpack_module_data(BCP_buffer& buf) {    buf.unpack(pModel_); }

   OsiSolverInterface* initialize_solver_interface();

   //Initializing a new search tree node.
   //This method serves as hook for the user to do some preprocessing on a search tree node before the node is processed.
   //Also, logical fixing results can be returned in the last four parameters.
   //This might be very useful if the branching implies significant tightening.
//   void initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
//      const BCP_vec<BCP_cut*>& cuts,
//      const BCP_vec<BCP_obj_status>& var_status,integerCoreVariables
//      const BCP_vec<BCP_obj_status>& cut_status,
//      BCP_vec<int>& var_changed_pos,
//      BCP_vec<double>& var_new_bd,
//      BCP_vec<int>& cut_changed_pos,
//      BCP_vec<double>& cut_new_bd)
//   { }

   //Try to generate a heuristic solution (or return one generated during cut/variable generation.
   //Return a pointer to the generated solution or return a NULL pointer.
   BCP_solution* generate_heuristic_solution(const BCP_lp_result& lpres,
   const BCP_vec<BCP_var*>& vars,
   const BCP_vec<BCP_cut*>& cuts);

   static bool compareCol(const pair<int,double>& p1, const pair<int,double>& p2);

   //Modify parameters of the LP solver before optimization.
   //This method provides an opportunity for the user to change parameters of the LP solver before optimization in the LP solver starts.
   //The second argument indicates whether the optimization is a "regular" optimization or it will take place in strong branching.
   //Default: empty method.
   void modify_lp_parameters ( OsiSolverInterface* lp, const int changeType, bool in_strong_branching);

   //print in cout a line summary of the current solver state
   void printSummaryLine(const BCP_vec<BCP_var*>& vars = {});

   //This method provides an opportunity for the user to tighten the bounds of variables.
   //The method is invoked after reduced cost fixing. The results are returned in the last two parameters.
   //Parameters:
   //lpres    the result of the most recent LP optimization,
   //vars  the variables in the current formulation,
   //status   the stati of the variables as known to the system,
   //var_bound_changes_since_logical_fixing    the number of variables whose bounds have changed (by reduced cost fixing) since the most recent invocation of this method that has actually forced changes returned something in the last two arguments,
   //changed_pos    the positions of the variables whose bounds should be changed
   //new_bd   the new bounds (lb/ub pairs) of these variables.
   void logical_fixing (const BCP_lp_result& lpres,
      const BCP_vec<BCP_var*>& vars,
      const BCP_vec<BCP_cut*>& cuts,
      const BCP_vec<BCP_obj_status>& var_status,
      const BCP_vec<BCP_obj_status>& cut_status,
      const int var_bound_changes_since_logical_fixing,
      BCP_vec<int>& changed_pos,
      BCP_vec<double>& new_bd);

   // Restoring feasibility.
   //This method is invoked before fathoming a search tree node that has been found infeasible and
   //the variable pricing did not generate any new variables.
   void restore_feasibility(const BCP_lp_result& lpres,
      const std::vector<double*> dual_rays,
      const BCP_vec<BCP_var*>& vars,
      const BCP_vec<BCP_cut*>& cuts,
      BCP_vec<BCP_var*>& vars_to_add,
      BCP_vec<BCP_col*>& cols_to_add);

   //Convert a set of variables into corresponding columns for the current LP relaxation.
   void vars_to_cols(const BCP_vec<BCP_cut*>& cuts, // on what to expand
      BCP_vec<BCP_var*>& vars,       // what to expand
      BCP_vec<BCP_col*>& cols,       // the expanded cols
      // things that the user can use for lifting vars if allowed
      const BCP_lp_result& lpres,
      BCP_object_origin origin, bool allow_multiple);


   //Convert (and possibly lift) a set of cuts into corresponding rows for the current LP relaxation.
   //Converting means computing for each cut the coefficients corresponding to each variable and
   //creating BCP_row objects that can be added to the formulation.
   //This method has different purposes depending on the value of the last argument.
   //If multiple expansion is not allowed then the user must generate a unique row for each cut.
   //This unique row must always be the same for any given cut.
   //This kind of operation is needed so that an LP relaxation can be exactly recreated.
   //On the other hand if multiple expansion is allowed then the user has (almost) free reign over what
   //she returns. She can delete some of the cuts or append new ones (e.g., lifted ones) to the end.
   //The result of the LP relaxation and the origin of the cuts are there to help her to make a decision about what to do.
   //For example, she might want to lift cuts coming from the Cut Generator, but not those coming from the Cut Pool.
   //The only requirement is that when this method returns the number of cuts and rows must be the same and
   //the i-th row must be the unique row corresponding to the i-th cut.
   //Here, we generate a cut to branch on a set of variables
   void cuts_to_rows(const BCP_vec<BCP_var*>& vars, //the variables currently in the relaxation (IN)
   BCP_vec<BCP_cut*>& cuts, //the cuts to be converted (IN/OUT)
   BCP_vec<BCP_row*>& rows, //the rows into which the cuts are converted (OUT)
   const BCP_lp_result&   lpres, //solution to the current LP relaxation (IN)
   BCP_object_origin origin, //where the cuts come from (IN)
   bool  allow_multiple); //whether multiple expansion, i.e., lifting, is allowed (IN)


   //Generate variables within the LP process.
   void generate_vars_in_lp(const BCP_lp_result& lpres,
      const BCP_vec<BCP_var*>& vars, const BCP_vec<BCP_cut*>& cuts, const bool before_fathom,
      BCP_vec<BCP_var*>& new_vars, BCP_vec<BCP_col*>& new_cols);

   /*
    * BCP_DoNotBranch_Fathomed: The node should be fathomed without even trying to branch.
    * BCP_DoNotBranch: BCP should continue to work on this node.
    * BCP_DoBranch: branch on one of the candidates cands
    *
    */
   BCP_branching_decision select_branching_candidates(const BCP_lp_result& lpres, //the result of the most recent LP optimization.
      const BCP_vec<BCP_var*> &  vars, //the variables in the current formulation.
      const BCP_vec< BCP_cut*> &  cuts, //the cuts in the current formulation.
      const BCP_lp_var_pool& local_var_pool, //the local pool that holds variables with negative reduced cost.
      //In case of continuing with the node the best so many variables will be added to the formulation (those with the most negative reduced cost).
      const BCP_lp_cut_pool& local_cut_pool, //the local pool that holds violated cuts.
      //In case of continuing with the node the best so many cuts will be added to the formulation (the most violated ones).
      BCP_vec<BCP_lp_branching_object*>&  cands, //the generated branching candidates.
      bool force_branch = false); //indicate whether to force branching regardless of the size of the local cut/var pools

   //Decide what to do with the children of the selected branching object.
   //Fill out the _child_action field in best. This will specify for every child what to do with it.
   //Possible values for each individual child are BCP_PruneChild, BCP_ReturnChild and BCP_KeepChild.
   //There can be at most child with this last action specified.
   //It means that in case of diving this child will be processed by this LP process as the next search tree node.
   //Default: Every action is BCP_ReturnChild.
   //However, if BCP dives then one child will be mark with BCP_KeepChild. The decision which child to keep is based on the ChildPreference parameter in BCP_lp_par.
   //Also, if a child has a presolved lower bound that is higher than the current upper bound then that child is mark as BCP_FathomChild.
   void set_actions_for_children(BCP_presolved_lp_brobj* best);

   //allow to count the nodes
   void select_cuts_to_delete (const BCP_lp_result &lpres, const BCP_vec< BCP_var * > &vars, const BCP_vec< BCP_cut * > &cuts, const bool before_fathom, BCP_vec< int > &deletable);

protected:
   BcpModeler* pModel_;
   //number of column in the master problem before the pricing
   int nbCurrentColumnVarsBeforePricing_;
   //count the iteration
   int lpIteration_;
   //count the nodes
   int last_node, nb_nodes;
   //if heuristic has been run. To be sure to run the heuristic no more than one time per node
   bool heuristicHasBeenRun_;

   //vars = are just the giver vars
   //cols is the vector where the new columns will be stored
   void TransformVarsToColumns(BCP_vec<BCP_var*>& vars, BCP_vec<BCP_col*>& cols);

   //Branch on the core integer var: the number of nurses var
   //just 2 children
   //Try also to fix to 1 some columns
   void appendNewBranchingVarsOnNumberOfNurses(CoinVar* integerCoreVar, vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands);

   //Branch on the core integer var: the rest arcs on a day for a nurse
   //just 2 children
   //Try also to fix to 1 some columns
   void appendNewBranchingVarsOnRest(int nbCuts, vector<MyObject*>& coreVars, vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands);

   //Try also to fix to 1 some columns
   void appendNewBranchingVarsOnColumns(vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<BCP_lp_branching_object*>&  cands);

   //build the vector of the branching candidates for the columns
   //return the indexes of the columns in the current formulation
   vector<int> buildBranchingColumns(CoinVar* var, vector<MyObject*>& columns,
      const BCP_vec<BCP_var*>&  vars, BCP_vec<int>& vpos, BCP_vec<double>& vbd);
};

/*
 * The class inherits the BCP_tm_user class from which the user can derive a problem specific class
 * to be used in the TM process.
 * In that derived class the user can store data to be used in the methods she overrides.
 * Also that is the object the user must return in the USER_initialize::tm_init() method.
 *
 * There are two kind of methods in the class.
 * The non-virtual methods are helper functions for the built-in defaults, but the user can use them as well.
 * The virtual methods execute steps in the BCP algorithm where the user might want to override the default behavior.
 */
class BcpBranchingTree: public BCP_tm_user{
public:
   BcpBranchingTree(BcpModeler* pModel);
   ~BcpBranchingTree() { }

   // pack the modeler
   void pack_module_data(BCP_buffer& buf, BCP_process_t ptype){
      switch (ptype) {
      case BCP_ProcessType_LP:
         buf.pack(pModel_); // Pack a pointer; does not work for parallel machines
         break;
      default:
         abort();
      }
   }

   // unpack an MIP feasible solution
//   BCP_solution* unpack_feasible_solution(BCP_buffer& buf);

   //display a feasible solution
   void display_feasible_solution(const BCP_solution* sol);

   // setting the base
   //Create the core of the problem by filling out the last three arguments.
   void initialize_core(BCP_vec<BCP_var_core*>& vars,
      BCP_vec<BCP_cut_core*>& cuts,
      BCP_lp_relax*& matrix);

   // create the root node
   //Create the set of extra variables and cuts that should be added
   //to the formulation in the root node.
   void create_root(BCP_vec<BCP_var*>& added_vars,
      BCP_vec<BCP_cut*>& added_cuts,
      BCP_user_data*& user_data);

   // feasible solution displaying
   //void display_feasible_solution(const BCP_solution* soln);

   // various initializations before a new phase (e.g., pricing strategy)
   void init_new_phase(int phase, BCP_column_generation& colgen, CoinSearchTreeBase*& candidates);

   // set search strategy
   //Values: 0 (BCP_BestFirstSearch), 1 (BCP_BreadthFirstSearch), 2 (BCP_DepthFirstSearch).
   void set_search_strategy(){
      switch(pModel_->getSearchStrategy()){
      case BestFirstSearch:
         set_param(BCP_tm_par::TreeSearchStrategy, 0);
         break;
      case BreadthFirstSearch:
         set_param(BCP_tm_par::TreeSearchStrategy, 1);
         break;
      case DepthFirstSearch:
         set_param(BCP_tm_par::TreeSearchStrategy, 2);
         break;
      }
   }

protected:
   BcpModeler* pModel_;
   int nbInitialColumnVars_;
   double minGap_;
};

/*
 * Define the default behaviour for the pack/unpack methods
 */
class BcpPacker : public BCP_user_pack {
public:
   BcpPacker(BcpModeler* pModel): pModel_(pModel){}
   ~BcpPacker() {}

   void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) { buf.pack(((BcpColumn*)var)->getIndex()); }

   BCP_var_algo* unpack_var_algo(BCP_buffer& buf) {
      int index = 0;
      buf.unpack(index);
      int i = index - pModel_->getCoreVars().size();
      BcpColumn* var = (BcpColumn*) pModel_->getColumns()[i];
      if(index != var->getIndex())
         Tools::throwError("Bad column unpacked or packed.");
      return new BcpColumn(*var);
   }

   void pack_cut_algo(const BCP_cut_algo* cons, BCP_buffer& buf) { buf.pack(((BcpBranchCons*)cons)->getIndex()); }

   BCP_cut_algo* unpack_cut_algo(BCP_buffer& buf) {
         int index = 0;
         buf.unpack(index);
         int i = index - pModel_->getCons().size();
         BcpBranchCons* cons = pModel_->getBranchingCons()[i];
         if(index != cons->getIndex())
            Tools::throwError("Bad cons unpacked or packed.");
         return new BcpBranchCons(*cons);
      }

protected:
   BcpModeler* pModel_;
};

/*
 * This class inherits the initializer class the user has to provide.
 *
 * The user will have to return an instance of the initializer class when the BCP_user_init() function is invoked.
 * The member methods of that instance will be invoked to create the various objects (well, pointers to them)
 * that are used/controlled by the user during the course of a run.
 */
class BcpInitialize : public USER_initialize {
public:
   BcpInitialize(BcpModeler* pModel): pModel_(pModel), pTree(0), pLPModel(0), pPacker(0) { }
   ~BcpInitialize() { }

   BCP_tm_user* tm_init(BCP_tm_prob& p, const int argnum, const char * const * arglist) {
      int size = pModel_->getNbColumns();
      pTree = new BcpBranchingTree(pModel_);
      return pTree;
   }

   BCP_lp_user* lp_init(BCP_lp_prob& p) {
      pLPModel = new BcpLpModel(pModel_);
      return pLPModel;
   }

   BCP_user_pack* packer_init(BCP_user_class* p)
   {
      pPacker = new BcpPacker(pModel_);
      return pPacker;
   }

protected:
   BcpModeler* pModel_;
   BcpLpModel* pLPModel;
   BcpBranchingTree* pTree;
   BcpPacker* pPacker;
};

#endif /* SRC_SCIPMODELER_H_ */
