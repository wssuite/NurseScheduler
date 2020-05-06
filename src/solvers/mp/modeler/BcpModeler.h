/*
 * BcpModeler.h
 *
 *  Created on: Mar 17, 2015
 *      Author: legraina
 */

#ifndef SRC_BCPMODELER_H_
#define SRC_BCPMODELER_H_

#include "solvers/mp/modeler/CoinModeler.h"
#include "solvers/mp/MasterProblem.h"

/* BCP includes */
#include "BCP_enum.hpp"
#include "BCP_vector.hpp"
#include "BCP_var.hpp"
#include "BCP_cut.hpp"
#include "BCP_buffer.hpp"
#include "BCP_tm_user.hpp"
#include "BCP_lp_user.hpp"
#include "BCP_lp.hpp"
#include "BCP_USER.hpp"
#include "BCP_solution.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinSearchTree.hpp"


//-----------------------------------------------------------------------------
//
//  C l a s s   B c p L p S o l
//
// Set of attributes of a solution of a relaxation
//
//-----------------------------------------------------------------------------

struct BcpLpSol {

	// constructor/destructor
	BcpLpSol() {}

	BcpLpSol(const std::vector<MyVar*>& activeColumnVars,
					 const std::map<int,int>& columnsToIndex,
					 const std::vector<double>& primalValues, const std::vector<double>& dualValues,
					 const std::vector<double>& reducedCosts, const std::vector<double>& lhsValues):
	activeColumnVars_(activeColumnVars),columnsToIndex_(columnsToIndex),
	primalValues_(primalValues),dualValues_(dualValues),reducedCosts_(reducedCosts),lhsValues_(lhsValues){}

	~BcpLpSol(){}

	// active columns in the solution
	std::vector<MyVar*> activeColumnVars_;
	std::map<int,int> columnsToIndex_;

	// primal and dual values of the solution
	std::vector<double> primalValues_;
	std::vector<double> dualValues_;
	std::vector<double> reducedCosts_;
	std::vector<double> lhsValues_;

	// set all the attributes of the solution
	void setLpSol(std::vector<MyVar*> activeColumnVars, std::map<int,int> columnsToIndex,
		std::vector<double> primalValues, std::vector<double> dualValues,
		std::vector<double> reducedCosts, std::vector<double> lhsValues) {
		activeColumnVars_ = activeColumnVars;
		columnsToIndex_ = columnsToIndex;
		primalValues_ = primalValues;
		dualValues_ = dualValues;
		reducedCosts_ = reducedCosts;
		lhsValues_ = lhsValues;
	}

};


//-----------------------------------------------------------------------------
//
//  C l a s s   B c p C o r e V a r
//
// These variables are not generated, they are always in the LP problem
//
//-----------------------------------------------------------------------------


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

	static VarType getMyVarType(BCP_var_t vartype){
		VarType type;
		switch(vartype){
		case BCP_BinaryVar:
			type = VARTYPE_BINARY;
			break;
		case BCP_IntegerVar:
			type = VARTYPE_INTEGER;
			break;
		default:
			type = VARTYPE_CONTINUOUS;
			break;
		}
		return type;
	}

	BcpCoreVar(const char* name, int index, double cost, VarType type, double lb, double ub,
	    const std::vector<double>& pattern = DEFAULT_PATTERN):
		CoinVar(name, index, cost, type, lb, ub, pattern),
		BCP_var_core(getBcpVarType(type), cost, lb, ub)
	{
		set_bcpind(index_);
	}

	BcpCoreVar(const BcpCoreVar& var) :
		CoinVar(var), BCP_var_core(getBcpVarType(type_), cost_, lb_, ub_)
	{
		set_bcpind(var.index_);
	}
};

//-----------------------------------------------------------------------------
//
//  C l a s s   B c p C o l u m n
//
// These variables are generated during the process, they are not always in the
// LP problem
//
//-----------------------------------------------------------------------------


struct BcpColumn: public CoinVar, public BCP_var_algo{
	BcpColumn(const char* name, int index, double cost, const std::vector<double>& pattern, double dualCost, VarType type, double lb, double ub,
			const std::vector<int>& indexRows = Tools::EMPTY_INT_VECTOR, const std::vector<double>& coeffs = Tools::EMPTY_DOUBLE_VECTOR):
				CoinVar(name, index, cost, type, lb, ub, pattern, dualCost, indexRows, coeffs),
				BCP_var_algo(BcpCoreVar::getBcpVarType(type), cost, lb, ub)
	{
		set_bcpind(index_);
	}

	BcpColumn(const BcpColumn& var) :
		CoinVar(var), BCP_var_algo(BcpCoreVar::getBcpVarType(type_), cost_, lb_, ub_)
	{
		set_bcpind(index_);
	}

	void pack(BCP_buffer& buf){
		buf.pack(name_);
      buf.pack(index_);
		buf.pack(_var_type);
		buf.pack(cost_);
		buf.pack(lb_);
		buf.pack(ub_);
		buf.pack(dualCost_);
		buf.pack(indexRows_);
		buf.pack(coeffs_);
		buf.pack(pattern_);
		buf.pack(iteration_creation_);
		buf.pack(active_count_);
		buf.pack(last_active_);
	}

	static BcpColumn* unpack(BCP_buffer& buf) {
	  char * name;
    buf.unpack(name);
    int index;
    buf.unpack(index);
    VarType type;
    buf.unpack(type);
    double cost, lb, ub, dualCost;
    buf.unpack(cost);
    buf.unpack(lb);
    buf.unpack(ub);
    buf.unpack(dualCost);
    std::vector<int> indexRows;
    buf.unpack(indexRows);
    std::vector<double> coeffs, pattern;
    buf.unpack(coeffs);
    buf.unpack(pattern);
    BcpColumn* col = new BcpColumn(name, index, cost, pattern, dualCost, type, lb, ub, indexRows, coeffs);
    buf.unpack(col->iteration_creation_);
    buf.unpack(col->active_count_);
    buf.unpack(col->last_active_);
    return col;
  }

	//use BCP bounds
	double getLB() const {
		return std::max(lb_,lb());
	}

	double getUB() const {
		return std::min(ub_,ub());
	}

	// use original bounds to set BCP_var_algo bounds
  void resetBounds() {
	  set_lb(lb_);
	  set_ub(ub_);
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

//these constraints are not generated, they are always in the LP problem
struct BcpBranchCons: public CoinCons, public BCP_cut_algo{
	BcpBranchCons(const char* name, int index, double lhs, double rhs,
                const std::vector<int>& indexCols, const std::vector<double>& coeffs):
				CoinCons(name, index, lhs, rhs),
				BCP_cut_algo(lhs, rhs),
				indexCols_(indexCols), coeffs_(coeffs)
	{ }

	BcpBranchCons(const BcpBranchCons& cons) :
		CoinCons(cons), BCP_cut_algo(lhs_, rhs_), indexCols_(cons.indexCols_), coeffs_(cons.coeffs_)
	{ }

    void pack(BCP_buffer& buf){
      buf.pack(name_);
      buf.pack(index_);
      buf.pack(lhs_);
      buf.pack(rhs_);
      buf.pack(indexCols_);
      buf.pack(coeffs_);
    }

    static BcpBranchCons* unpack(BCP_buffer& buf) {
      char *name;
      buf.unpack(name);
      int index;
      buf.unpack(index);
      double lhs, rhs;
      buf.unpack(lhs);
      buf.unpack(rhs);
      std::vector<int> indexCols;
      buf.unpack(indexCols);
      std::vector<double> coeffs;
      buf.unpack(coeffs);

      return new BcpBranchCons(name, index, lhs, rhs, indexCols, coeffs);
    }

    std::vector<int>& getIndexCols() { return indexCols_; }

    std::vector<double>& getCoeffCols() { return coeffs_; }

protected:
    std::vector<int> indexCols_; //index of the cols of the matrix where the col has non-zero coefficient
    std::vector<double> coeffs_; //value of these coefficients
};

// LP solver types
//
enum LPSolverType {CLP, Gurobi, Cplex};
static std::map<std::string,LPSolverType> LPSolverTypesByName =
	{{"CLP",CLP},{"Gurobi",Gurobi},{"Cplex",Cplex}};



//-----------------------------------------------------------------------------
//
//  C l a s s   B c p M o d e l e r
//
// Specific class of modeler designed to use BCP.
// In particular, may methods defined in the BCP library need to be implemented
// here
//
//-----------------------------------------------------------------------------
class BcpInitialize; // class forward declaration

class BcpModeler: public CoinModeler {
public:
	BcpModeler(MasterProblem* pMaster, const char* name, LPSolverType type = CLP);
	~BcpModeler();

	//solve the model
	int solve(bool relaxation = false);

	//Reset and clear solving parameters
	void reset();

	void addActiveColumn(MyVar* var, int index=-1){
		BcpColumn* col = dynamic_cast<BcpColumn*>(var);
		if(!col)
			std::cout << "error";
		Modeler::addActiveColumn(var, index);
		if(index>=0) columnsToIndex_[col->getIndex()] = index;
	}

	// Delete all the objects owned by this modeler and then call the parent function
  void clear();

	void clearActiveColumns() {
    Modeler::clearActiveColumns();
    columnsToIndex_.clear();
	}

protected:
	/*
	 * Create variable:
	 *    var is a pointer to the pointer of the variable
	 *    var_name is the name of the variable
	 *    lhs, rhs are the lower and upper bound of the variable
	 *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
	 */
  int createVar(MyVar** var, const char* var_name, int index, double objCoeff,
                double lb, double ub, VarType vartype, const std::vector<double>& pattern, double score) override;

  int createColumnVar(MyVar** var, const char* var_name, int index, double objCoeff,
        const std::vector<double>& pattern, double dualObj, double lb, double ub, VarType vartype, double score) override ;

	/*
	 * Create linear constraint:
	 *    con is a pointer to the pointer of the constraint
	 *    con_name is the name of the constraint
	 *    lhs, rhs are the lower and upper bound of the constraint
	 *    nonZeroVars is the number of non-zero coefficients to add to the constraint
	 */

	int createCoinConsLinear(MyCons** con, const char* con_name, int index, double lhs, double rhs);

  int createCoinCutLinear(MyCons** con, const char* con_name, int index, double lhs, double rhs,
                                    const std::vector<int>& indexVars, const std::vector<double>& coeffs);

	MyVar* copyColumn(MyVar* var) const;

  //  Delete all the objects owned by this modeler
  void deleteSolutions();

public:
	/*
	 * Get/set the primal value
	 */

	double getVarValue(MyVar* var) const;
	void setVarValue(MyVar* var, double value);

	/*
	 * Get the dual variables
	 */

	double getDual(MyCons* cons, bool transformed = false) const;

	/*
	 * Get the reduced cost
	 */

	double getReducedCost(MyVar* var) const;

	/**************
	 * Parameters *
	 *************/
	int setVerbosity(int v);
	inline void setParameters(const SolverParam& parameters, PrintSolution* func=nullptr){
		Modeler::setParameters(parameters, func);
		if (pPricer_) {
			pPricer_->initPricerParameters(parameters);
		}
	}

	/**************
	 * Outputs *
	 *************/

	int writeProblem(std::string fileName) const;

	int writeLP(std::string fileName) const;

	/*
	 * Class own methods and parameters
	 */

	void setLPSol(const BCP_lp_result& lpres, const BCP_vec<BCP_var*>&  vars, int lpIteration);

	int getIndexCol(int colIndex) {
		auto it = columnsToIndex_.find(colIndex);
		if(it == columnsToIndex_.end()) return -1;
		return it->second;
	}

	void addBcpSol(const BCP_solution* sol);

	// Get the index of the best solution in the vector of solutions of BCP
	int getBestSolIndex();

	// Clear the active column, set the active columns with those in the best solution,
	// and set the primal values accordingly
	bool loadBestSol();

	// Set the value of the active columns with those in the best solution
	void setActiveColumnsValuesWithBestSol();

	// Clear the active column, set the active columns with those in the solution with input index
	// and set the primal values accordingly
	void loadBcpSol(int index);

	// Clear the active column, set the active columns with those in the input
	// solution and set the primal values accordingly
	void loadInputSol(const BCP_solution_generic& sol);

	inline void setPrimal(const std::vector<double>& primal){ primalValues_ = primal; }

	inline int getFrequency() { return TmVerb_SingleLineInfoFrequency; }

	inline void setLastNbSubProblemsSolved(int lastNbSubProblemsSolved){ lastNbSubProblemsSolved_ = lastNbSubProblemsSolved; }

	inline int getLastNbSubProblemsSolved(){ return lastNbSubProblemsSolved_; }

	inline double getLastMinDualCost(){ return lastMinDualCost_; }
	inline void setLastMinDualCost(double lastMinDualCost){ lastMinDualCost_ = lastMinDualCost; }

	inline double getLastObj(){ return obj_history_.empty()? infinity_:obj_history_.back(); }

	// Set every rotation to one : this is useful only when the active columns
	// are only the rotations included in a provided initial solution
	//
	void setEveryRotationToOne();

	// fix/unfix all the rotations variables starting from the input vector of days
	void fixRotationsStartingFromDays(const std::vector<bool>& isFixDay);
	void unfixRotationsStartingFromDays(const std::vector<bool>& isUnfixDay);

	// fix/unfix all the rotations variables of the input nurses
	void fixRotationsOfNurses(const std::vector<bool>& isFixNurse);
	void unfixRotationsOfNurses(const std::vector<bool>& isUnfixNurse);

	// relax/unrelax the integrality of all the rotations variables starting from the input vector of days
	void relaxRotationsStartingFromDays(const std::vector<bool>& isRelaxDay);
	void unrelaxRotationsStartingFromDays(const std::vector<bool>& isUnRelaxDay);

	/*
	 * Manage the storage of our own tree
	 */

	inline void setCurrentNode(const CoinTreeSiblings* s) {
		/* the current node of this siblings is already taken as processed */
		int nodeIndex = s->size() - s->toProcess() - 1;
		pTree_->setCurrentNode(treeMapping_[s][nodeIndex]);

		/* if no more child in this siblings */
		if(s->toProcess() == 0){
			treeMapping_.erase(s);
			pTree_->eraseCurrentSibblings();
		}
	}

	inline MyNode* getCurrentNode() {return pTree_->getCurrentNode();}


	inline void addToMapping(const CoinTreeSiblings* s) {
		const int nbLeaves = s->size();
    treeMapping_[s] = pTree_->addToMapping(nbLeaves, s->currentNode()->getDepth());
	}

	inline MyNode* getNode(const CoinTreeSiblings* s) {
    int nodeIndex = s->size() - s->toProcess();
    return treeMapping_[s][nodeIndex];
  }

	/*
	 * Parameters getters
	 */

	inline LPSolverType getLPSolverType() { return LPSolverType_; };

	inline  std::map<BCP_tm_par::chr_params, bool>& getTmParameters(){ return tm_parameters; }

	inline  std::map<BCP_lp_par::chr_params, bool>& getLpParameters(){ return lp_parameters; }

	inline bool is_solution_changed() { return solHasChanged_; }

	inline int nbSolutions() const { return bcpSolutions_.size(); }

	inline double getObjective() const {return Modeler::getObjective(); }

	inline double getObjective(int index) const { return bcpSolutions_[index].objective_value(); }

	inline MasterProblem* getMaster() const {return pMaster_;}

	//check if Bcp stops
	bool doStop() const;

	//check the active rotations
	void checkActiveColumns(const BCP_vec<BCP_var*>& vars) const;

	std::pair<BCP_lp_par::int_params, int> strong_branching = std::pair<BCP_lp_par::int_params, int>(BCP_lp_par::MaxPresolveIter, -1); //disable strong branching

	// Get/set the value of the current level in the branch and bound tree
	//
	int getCurrentTreeLevel() const {return currentTreeLevel_;}
	void setCurrentTreeLevel(int level) {currentTreeLevel_ = level;}

	// STAB
	// Get/set the number of consecutive column generation degenerate iterations
	//
	int getNbDegenerateIt() const {return nbDegenerateIt_;}
	void setNbDegenerateIt(int nbIt) {nbDegenerateIt_=nbIt;}
	void incrementNbDegenerateIt() {nbDegenerateIt_++;}

	// LNS
	// Record the current solution of the relaxation
	//
	BcpLpSol recordLpSol() const {
		BcpLpSol currentSol(activeColumnVars_, columnsToIndex_, primalValues_, dualValues_, reducedCosts_, lhsValues_);
		return currentSol;
	}

	// Get the solution of the root node relaxation
	//
	const BcpLpSol& getRootSolution() const {return rootSolution_;}

	// Get/set statistics
	//
	BCP_lp_statistics getTimeStats() const {return timeStats_;}
	double gettimeFirstRoot() const {return timeFirstRoot_;}
	void settimeFirstRoot(double t) {timeFirstRoot_ = t;}
	int getNbLpIterations() const {return nbLpIterations_;}
	void incrementNbNodes() {nbNodes_++;}
	int getNbNodes() const {return nbNodes_-1;}

// protected:

	// solver that called the model
	MasterProblem* pMaster_;
  BcpInitialize* pBcp_;

protected:

  //mapping between the CoinTreeSiblings* and my BcpNode*
	//a sibblings contains a list of all its leaves CoinTreeNode
  std::map<const CoinTreeSiblings*, std::vector<MyNode*>> treeMapping_;
	//results
  std::vector<double> obj_history_;
    std::vector<double> primalValues_, dualValues_, reducedCosts_, lhsValues_;
    std::map<int,int> columnsToIndex_;
	bool solHasChanged_ = false; //reload solution ?
	//bcp solution
  std::vector<BCP_solution_generic> bcpSolutions_;
	// bcp solution of the root node
	BcpLpSol rootSolution_;


	/* stats */
	//number of sub problems solved on the last iteration of column generation
	int lastNbSubProblemsSolved_;
	//min dual cost for a pattern on the last iteration of column generation
	double lastMinDualCost_;
	//all column generation times as computed by BCP
	BCP_lp_statistics timeStats_;
	//other important timer
	double timeFirstRoot_=0.0;
	// number of iterations of column generation
	int nbLpIterations_=0;
	// number of branch and bound nodes explored in the tree
	int nbNodes_;

	/* Parameters */
	LPSolverType LPSolverType_;

	// STAB: number of consecutive degenerate column generation iterations
	int nbDegenerateIt_=0;

	//At every this many search tree node provide a single line info on the progress of the search tree.
	//If <= 0 then never.
	//Default: 0.
	int TmVerb_SingleLineInfoFrequency = 0;

	// Current level in the branch and bound tree
	//
	int currentTreeLevel_ = 0;


	/* Tree Manager verbosity parameters */
  std::map<BCP_tm_par::chr_params, bool> tm_parameters = {
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
			{ BCP_tm_par::TmVerb_Last, 0},
			{ BCP_tm_par::DebugLpProcesses,0}
	};

	/* LP verbosity parameters */
  std::map<BCP_lp_par::chr_params, bool> lp_parameters = {
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
			{ BCP_lp_par::LpVerb_Last, 0},// Just a marker for the last LpVerb
			{ BCP_lp_par::DoReducedCostFixingAtZero, 0}, // Deactivate reduced cost fixing
			{ BCP_lp_par::DoReducedCostFixingAtAnything, 0} // Deactivate reduced cost fixing
	};
};

//-----------------------------------------------------------------------------
//
//  C l a s s   B c p L p M o d e l
//
// The class inherits the BCP_lp_user class from which the user can derive a
// problem specific class to be used in the LP process.
// In that derived class the user can store data to be used in the methods she
// overrides.
//  Also that is the object the user must return in the
// USER_initialize::lp_init() method.
//
// There are two kind of methods in the class.
// The non-virtual methods are helper functions for the built-in defaults, but
// the user can use them as well.
// The virtual methods execute steps in the BCP algorithm where the user might
// want to override the default behavior.
//
//-----------------------------------------------------------------------------


class BcpLpModel: public BCP_lp_user {
public:
	BcpLpModel(BcpModeler* pModel);
	~BcpLpModel() {
	}

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

	static bool compareCol(const std::pair<int,double>& p1, const std::pair<int,double>& p2);

	//Modify parameters of the LP solver before optimization.
	//This method provides an opportunity for the user to change parameters of the LP solver before optimization in the LP solver starts.
	//The second argument indicates whether the optimization is a "regular" optimization or it will take place in strong branching.
	//Default: empty method.
	void modify_lp_parameters ( OsiSolverInterface* lp, const int changeType, bool in_strong_branching);

	//print in cout a line summary of the current solver state
	void printSummaryLine(const BCP_vec<BCP_var*>& vars = {});

	//stop this node or BCP
	bool doStop();

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

	// This method needs to be reimplemented as empty, because it causes cycling in
	// the column generation when the lower bound is too close from the upper bound
	//
	void reduced_cost_fixing(const double* dj, const double* x, const double gap,
			BCP_vec<BCP_var*>& vars, int& newly_changed);

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

	void select_vars_to_delete(const BCP_lp_result& lpres,
			const BCP_vec<BCP_var*>& vars,
			const BCP_vec<BCP_cut*>& cuts,
			const bool before_fathom,
			BCP_vec<int>& deletable);

    int writeLP(std::string fileName) {
      getLpProblemPointer()->lp_solver->writeLp(fileName.c_str());
      return 0;
    }

	// STAB
	// Update the bounds and/or costs of the stabilization variables
	bool stabUpdateBoundAndCost(bool isStall, bool isImproveQuality) ;

	//Parameters to be deleted
	int min_inactive_iteration = 50;
	double max_activity_rate = .1;

	// getters/setters
  BCP_lp_statistics getTimeStats() {
    return getLpProblemPointer()->stat;
  }

	int getNbLpIterations() const {return lpIteration_;}

protected:
	BcpModeler* pModel_;
	//count the iteration
	int currentNodelpIteration_, lpIteration_;
	//count the nodes
	int last_node;
	// stored if the current node corresponds to a backtracking
	bool backtracked_;
	//if heuristic has been run. To be sure to run the heuristic no more than one time per node
	bool heuristicHasBeenRun_;
	int nbNodesSinceLastHeuristic_;
	//number of generated columns
	int nbGeneratedColumns_;
	// Number of dives to wait before branching on columns again
  std::list<int> nb_dives_to_wait_before_branching_on_columns_;

	//vars = are just the giver vars
	//cols is the vector where the new columns will be stored
	void TransformVarsToColumns(BCP_vec<BCP_var*>& vars, BCP_vec<BCP_col*>& cols);

	//build the candidate from my candidate
	void buildCandidate(const MyBranchingCandidate& candidate, const BCP_vec<BCP_var*>&  vars, const BCP_vec< BCP_cut*> &  cuts, BCP_vec<BCP_lp_branching_object*>&  cands);

	// rerun the code use to test the integer feasibility of a solution and find why a solution is not feasible
  void find_infeasibility(const BCP_lp_result& lpres, //the result of the most recent LP optimization.
                          const BCP_vec<BCP_var*> &  vars);
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

	// various initializations before a new phase (e.g., pricing strategy)
	void init_new_phase(int phase, BCP_column_generation& colgen, CoinSearchTreeBase*& candidates);

	//override current method
	//WARNINGS: if not, BCP loads automatically a CoinSearchTree<CoinSearchTreeCompareDepth> in new solution()
	//   void
	//   BCP_tm_user::change_candidate_heap(CoinSearchTreeManager& candidates,
	//                  const bool new_solution)
	//   {
	//       if (new_solution) {
	//      candidates.newSolution(p->ub());
	//       } else {
	//      candidates.reevaluateSearchStrategy();
	//       }
	//   }
	void change_candidate_heap(CoinSearchTreeManager& candidates, const bool new_solution) {}

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
		default:
			break;
		}
	}

protected:
	BcpModeler* pModel_;
	int nbInitialColumnVars_;
};

template<class Comp> class MyCoinSearchTree: public CoinSearchTreeBase {
public:
	MyCoinSearchTree(BcpModeler* pModel): CoinSearchTreeBase(), pModel_(pModel), comp_() {}
	MyCoinSearchTree(const MyCoinSearchTree& t) :
		CoinSearchTreeBase(), comp_(), pModel_(t.pModel_) {
		candidateList_ = t.getCandidates();
		numInserted_ = t.numInserted();
		size_ = t.size();
	}
	virtual ~MyCoinSearchTree() {}
	const char* compName() const { return Comp::name(); }

protected:
	/*
	 * Tree: allow to update our own tree
	 */

	void realpop() {
		/* update the current node of the modeler */
		pModel_->setCurrentNode(this->candidateList_[0]);
		/* the siblings is now empty -> choose the next one */
		//copy the best candidate at the first place
		candidateList_[0] = candidateList_.back();
		//and remove the last item
		candidateList_.pop_back();
	}

	void fixTop() {
		/* update the current node of the modeler */
		pModel_->setCurrentNode(this->candidateList_[0]);
	}

	void realpush(CoinTreeSiblings* s) {
		//add the current node to the BcpModeler
		if(s->toProcess() > 0){
			pModel_->addToMapping(s);
			this->candidateList_.push_back(s);
		}

		/* update the quality  and then the candidateList */
		//update quality
		for(CoinTreeSiblings* s1: this->candidateList_){
			double q = pModel_->getNode(s1)->getQuality();
			s1->currentNode()->setQuality(q);
		}
		//reorder the list with the new quality
		//we dont use a heap
		std::stable_sort(candidateList_.begin(), candidateList_.end(), comp_);
	}

private:
	BcpModeler* pModel_;
	Comp comp_;
};

/*
 * Define the default behaviour for the pack/unpack methods
 */
class BcpPacker : public BCP_user_pack {
public:
	BcpPacker(BcpModeler* pModel): pModel_(pModel){}
	~BcpPacker() {}

	void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) {
	  ((BcpColumn*)var)->pack(buf);
	}

	BCP_var_algo* unpack_var_algo(BCP_buffer& buf) {
    return BcpColumn::unpack(buf);
  }

  void pack_cut_algo(const BCP_cut_algo* cons, BCP_buffer& buf) {
    ((BcpBranchCons*)cons)->pack(buf);
  }

  BCP_cut_algo* unpack_cut_algo(BCP_buffer& buf) {
    return BcpBranchCons::unpack(buf);
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
	BcpInitialize(BcpModeler* pModel): pModel_(pModel), pLpModel_(0), pTree_(0), pPacker_(0) { }
	~BcpInitialize() {
		if (pLpModel_->getLpProblemPointer()) delete pLpModel_->getLpProblemPointer();
		// if (pLpModel_->getOsiBabSolver()) delete pLpModel_->getOsiBabSolver();
		delete pTree_;
	}

	BCP_tm_user* tm_init(BCP_tm_prob& p, const int argnum, const char * const * arglist) {
		if(!pTree_) pTree_ = new BcpBranchingTree(pModel_);
		return pTree_;
	}

	BCP_lp_user* lp_init(BCP_lp_prob& p) {
		if(!pLpModel_) pLpModel_ = new BcpLpModel(pModel_);
		return pLpModel_;
	}

	BCP_user_pack* packer_init(BCP_user_class* p)
	{
		if(!pPacker_) pPacker_ = new BcpPacker(pModel_);
		return pPacker_;
	}

	BCP_lp_statistics getTimeStats() {
    return pLpModel_->getTimeStats();
	}

	int getNbLpIterations() const {
	  return pLpModel_->getNbLpIterations();
	}

  int writeLP(std::string fileName) const{
    if(pLpModel_) return pLpModel_->writeLP(fileName);
    std::cout << "WARNING: BCP cannot write the model as the LP solver has not been initialized." << std::endl;
    return 1;
  }

// protected:
	BcpModeler* pModel_;
	BcpLpModel* pLpModel_;
	BcpBranchingTree* pTree_;
	BcpPacker* pPacker_;
};

#endif /* SRC_BCPMODELER_H_ */
