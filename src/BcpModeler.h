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
#include "BCP_parameters.hpp"
#include "BCP_enum.hpp"
#include "BCP_vector.hpp"
#include "BCP_var.hpp"
#include "BCP_cut.hpp"
#include "BCP_buffer.hpp"
#include "BCP_parameters.hpp"
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
   { }

   BcpCoreVar(const BcpCoreVar& var) :
      CoinVar(var), BCP_var_core(getBcpVarType(type_), cost_, lb_, ub_)
   { }

   ~BcpCoreVar(){ }
};

//these variables are generated during the process, they are not always in the LP problem
struct BcpColumn: public CoinVar, public BCP_var_algo{
   BcpColumn(const char* name, int index, double cost, VarType type, double lb, double ub):
      CoinVar(name, index, cost, type, lb, ub),
      BCP_var_algo(BcpCoreVar::getBcpVarType(type), cost, lb, ub)
   { }

   BcpColumn(const BcpColumn& var) :
      CoinVar(var), BCP_var_algo(BcpCoreVar::getBcpVarType(type_), cost_, lb_, ub_)
   { }

   //BCP takes care of unpacking the BCP_var_algo
   BcpColumn(BCP_buffer& buf) :
      CoinVar(), BCP_var_algo(BCP_ContinuousVar, 0, 0, 0)
   {
      buf.unpack(varName_).unpack(index_).unpack(type_).unpack(cost_).unpack(lb_).unpack(ub_)
                                                        .unpack(indexRows_).unpack(coeffs_);
   }

   ~BcpColumn(){ }

   //pack method of BCP for BcpColumn
   //BCP takes care of packing the BCP_var_algo
   void pack(BCP_buffer& buf) const {
      buf.pack(varName_).pack(index_).pack(type_).pack(cost_).pack(lb_).pack(ub_)
                                                        .pack(indexRows_).pack(coeffs_);
   }
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

class BcpModeler: public CoinModeler {
public:
   BcpModeler(const char* name):
      CoinModeler(), primalValues_(0), dualValues_(0), reducedCosts_(0), lhsValues_(0), best_lb_in_root(-DBL_MAX)
      { }
   ~BcpModeler() { }

   //solve the model
   int solve(bool relaxation = false);

   //Add a pricer
   int addObjPricer(MyPricer* pPricer){
      pPricer_ = pPricer;
      return 1;
   }

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
    */
   int createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub);

   int createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub);

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

   int printStats();

   int printBestSol();

   int writeProblem(string fileName);

   int writeLP(string fileName);

   /*
    * Class own methods and parameters
    */

   void setLPSol(const BCP_lp_result& lpres){
      obj_history_.push_back(lpres.objval());
      primalValues_ = const_cast<double*>(lpres.x());
      dualValues_ = const_cast<double*>(lpres.pi());
      reducedCosts_ = const_cast<double*>(lpres.dj());
      lhsValues_ = const_cast<double*>(lpres.lhs());
   }


protected:
   //best lb in root
   double best_lb_in_root;

   //results
   vector<double> obj_history_;
   double *primalValues_, *dualValues_, *reducedCosts_, *lhsValues_;
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
   BcpLpModel(BcpModeler* pModel):
      pModel_(pModel), nbCurrentColumnVarsBeforePricing_(pModel->getNbColumns())
{ }
   ~BcpLpModel() { }

   /*
    * BCP_lp_user methods
    */

   void unpack_module_data(BCP_buffer& buf) {    buf.unpack(pModel_); }

   void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf){ ((BcpColumn*)var)->pack(buf); }

   BCP_var_algo* unpack_var_algo(BCP_buffer& buf) { return new BcpColumn(buf); }

   OsiSolverInterface* initialize_solver_interface(){ return new OsiClpSolverInterface(); }

   //Initializing a new search tree node.
   //This method serves as hook for the user to do some preprocessing on a search tree node before the node is processed.
   //Also, logical fixing results can be returned in the last four parameters.
   //This might be very useful if the branching implies significant tightening.
   void initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
      const BCP_vec<BCP_cut*>& cuts,
      const BCP_vec<BCP_obj_status>& var_status,
      const BCP_vec<BCP_obj_status>& cut_status,
      BCP_vec<int>& var_changed_pos,
      BCP_vec<double>& var_new_bd,
      BCP_vec<int>& cut_changed_pos,
      BCP_vec<double>& cut_new_bd)
   { }

   //Convert a set of variables into corresponding columns for the current LP relaxation.
   void vars_to_cols(const BCP_vec<BCP_cut*>& cuts, // on what to expand
      BCP_vec<BCP_var*>& vars,       // what to expand
      BCP_vec<BCP_col*>& cols,       // the expanded cols
      // things that the user can use for lifting vars if allowed
      const BCP_lp_result& lpres,
      BCP_object_origin origin, bool allow_multiple);

   //Generate variables within the LP process.
   void generate_vars_in_lp(const BCP_lp_result& lpres,
      const BCP_vec<BCP_var*>& vars, const BCP_vec<BCP_cut*>& cuts, const bool before_fathom,
      BCP_vec<BCP_var*>& new_vars, BCP_vec<BCP_col*>& new_cols);

protected:
   BcpModeler* pModel_;
   int nbCurrentColumnVarsBeforePricing_;

   //vars = are just the giver vars
   //cols is the vector where the new columns will be stored
   void TransformVarsToColumns(BCP_vec<BCP_var*>& vars, BCP_vec<BCP_col*>& cols);
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
   BcpBranchingTree(BcpModeler* pModel):
      pModel_(pModel) ,clp_(0), nbInitialColumnVars_(pModel->getNbColumns())
{ }
   ~BcpBranchingTree() {
      if (clp_)
         delete clp_;
   }

   // pack the modeler
   void pack_module_data(BCP_buffer& buf, BCP_process_t ptype){ buf.pack(pModel_); }

   // unpack an MIP feasible solution
   //BCP_solution* unpack_feasible_solution(BCP_buffer& buf);

   void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) { ((BcpColumn*)var)->pack(buf); }

   BCP_var_algo* unpack_var_algo(BCP_buffer& buf) { return new BcpColumn(buf); }

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
   void init_new_phase(int phase, BCP_column_generation& colgen) { colgen = BCP_GenerateColumns; }

protected:
   BcpModeler* pModel_;
   OsiClpSolverInterface* clp_;
   int nbInitialColumnVars_;
};

/*
 * Define the default behaviour for the pack/unpack methods
 */
class BcpPacker : public BCP_user_pack {
public:
   BcpPacker() {}
   ~BcpPacker() {}

   void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) { ((BcpColumn*)var)->pack(buf); }
   BCP_var_algo* unpack_var_algo(BCP_buffer& buf) { return new BcpColumn(buf); }
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
   ~BcpInitialize() {
      if(pLPModel)
         delete pLPModel;
      if(pTree)
         delete pTree;
      if(pPacker)
         delete pPacker;
   }

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
      pPacker = new BcpPacker();
      return pPacker;
   }

protected:
   BcpModeler* pModel_;
   BcpLpModel* pLPModel;
   BcpBranchingTree* pTree;
   BcpPacker* pPacker;
};

#endif /* SRC_SCIPMODELER_H_ */
