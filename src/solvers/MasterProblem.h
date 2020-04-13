/*
* MasterProblem.h
*
*  Created on: 3 fevr. 2015
*      Author: samuel
*/

#ifndef MASTERPROBLEM_H_
#define MASTERPROBLEM_H_

/* Inheritance */
#include "solvers/Solver.h"

/* Tools include */
#include "tools/MyTools.h"

/* My includes */
#include "data/Nurse.h"
#include "solvers/mp/modeler/Modeler.h"
#include "OsiSolverInterface.hpp"

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------
enum CostType {TOTAL_COST, CONS_SHIFTS_COST, CONS_WORKED_DAYS_COST, COMPLETE_WEEKEND_COST, PREFERENCE_COST, INIT_REST_COST};

struct DualCosts{
public:

	DualCosts(std::vector< std::vector<double> > & workCosts, std::vector<double> & startWorkCosts, std::vector<double> & endWorkCosts, double workedWeekendCost, bool doNotCopy = true):
	workCosts_(workCosts), startWorkCosts_(startWorkCosts), endWorkCosts_(endWorkCosts), workedWeekendCost_(workedWeekendCost) {

		if(!doNotCopy){
			workCosts_.clear();
			workCosts_ = Tools::appendVectors(workCosts_, workCosts);
			startWorkCosts_.clear();
			startWorkCosts_ = Tools::appendVectors(startWorkCosts_, startWorkCosts);
			endWorkCosts_.clear();
			endWorkCosts_ = Tools::appendVectors(endWorkCosts_, endWorkCosts);
		}
	}

	// GETTERS
	//
	inline double dayShiftWorkCost(int day, int shift){return (workCosts_[day][shift-1]);}
	inline double startWorkCost(int day){return (startWorkCosts_[day]);}
	inline double endWorkCost(int day){return (endWorkCosts_[day]);}
	inline double workedWeekendCost(){return workedWeekendCost_;}


protected:

	// Indexed by : (day, shift) !! 0 = shift 1 !!
	vector2D<double> & workCosts_;

	// Indexed by : day
	std::vector<double> & startWorkCosts_;

	// Indexed by : day
	std::vector<double> & endWorkCosts_;

	// Reduced cost of the weekends
	double workedWeekendCost_;

};


struct Rotation {

	// Specific constructors and destructors
	//
	Rotation(std::map<int,int> shifts, int nurseId = -1, double cost = DBL_MAX, double dualCost = DBL_MAX) :
	shifts_(shifts), id_(s_count),nurseId_(nurseId), cost_(cost),
	consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
	dualCost_(dualCost), length_(shifts.size()), timeDuration_(shifts.size())
	{
		++s_count;
		firstDay_ = 999;
		for(std::map<int,int>::iterator itS = shifts.begin(); itS != shifts.end(); ++itS)
		if(itS->first < firstDay_) firstDay_ = itS->first;
	};

	Rotation(int firstDay, std::vector<int> shiftSuccession, int nurseId = -1, double cost = DBL_MAX, double dualCost = DBL_MAX) :
					id_(s_count),nurseId_(nurseId), cost_(cost),
	consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
	dualCost_(dualCost), firstDay_(firstDay), length_(shiftSuccession.size()), timeDuration_(shiftSuccession.size())
	{
		++s_count;
		for(int k=0; k<length_; k++) shifts_.insert(std::pair<int,int>( (firstDay+k) , shiftSuccession[k] ));
	}

	Rotation(std::vector<double> compactPattern) :
					id_(s_count),nurseId_((int)compactPattern[0]), cost_(DBL_MAX),
	consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
	dualCost_(DBL_MAX), firstDay_((int)compactPattern[1]), length_(compactPattern.size()-2), timeDuration_(compactPattern.size()-2)
	{
		++s_count;
		for(int k=0; k<length_; k++) shifts_.insert(std::pair<int,int>( (firstDay_+k) , (int)compactPattern[k+2] ));
	}

	Rotation(Rotation& rotation, int nurseId) :
	  shifts_(rotation.shifts_), id_(rotation.id_), nurseId_(nurseId), cost_(rotation.cost_),
	consShiftsCost_(rotation.consShiftsCost_), consDaysWorkedCost_(rotation.consDaysWorkedCost_),
	completeWeekendCost_(rotation.completeWeekendCost_), preferenceCost_(rotation.preferenceCost_), initRestCost_(rotation.initRestCost_),
					dualCost_(rotation.dualCost_), firstDay_(rotation.firstDay_), length_(rotation.length_),
					timeDuration_(rotation.timeDuration_)
	{
		if(rotation.nurseId_ != nurseId_){
			cost_ = DBL_MAX;
			dualCost_ = DBL_MAX;
		}
	}

	~Rotation(){}

	//count rotations
	//
	static unsigned int s_count;

	// Shifts to be performed
	//
	std::map<int,int> shifts_;

	//Id of the rotation
	//
	long id_;

	//the nurse
	//
	int nurseId_;

	// Cost
	//
	double cost_;
	double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_, initRestCost_ ;

	// Dual cost as found in the subproblem
	//
	double dualCost_;

	// First worked day
	//
	int firstDay_;

	// Duration
	//
	int length_;

	// Level of the branch and bound tree where the rotation has been generated
	//
	int treeLevel_=0;

        // Time duration (in hours)

        int timeDuration_;

	//compact the rotation in a vector
	const std::vector<double> getCompactPattern(){
		std::vector<double> compact;
		compact.push_back(nurseId_);
		compact.push_back(firstDay_);
		for(std::pair<int,int> p: shifts_) compact.push_back(p.second);
		return compact;
	}

	//Compute the cost of a rotation
	//
	void computeCost(Scenario* pScenario, Preferences* pPreferences, const std::vector<LiveNurse*>& liveNurses, int horizon);

	//Compute the dual cost of a rotation
	//
	void checkDualCost(DualCosts& costs);


        // calcule le nombre d'heures d'une rotation
        void computeTimeDuration(Scenario* pScenario){
	  timeDuration_=0;
	  for(std::pair<int,int> p: shifts_) {
	    timeDuration_ += pScenario->hoursToWork_[p.second];
	  }
	}

	std::string toString(int nbDays = -1, std::vector<int> shiftIDToShiftTypeID={}) const {
		if(nbDays == -1) nbDays = firstDay_+length_;
		std::stringstream rep;
		rep << "#   | ROTATION: N=" << nurseId_ << "  cost=" << cost_ << "  dualCost=" << dualCost_ << "  firstDay=" << firstDay_ << "  length=" << length_ << std::endl;
		rep << "#               |";
		std::vector<int> allTasks (nbDays);
		for(auto itTask = shifts_.begin(); itTask != shifts_.end(); ++itTask)
		    allTasks[itTask->first] = itTask->second;
		for(unsigned int i=0; i<allTasks.size(); i++){
			if(allTasks[i] < 1) rep << "\t|";
			else {
			    int t = allTasks[i];
			    if(t < shiftIDToShiftTypeID.size())
			        rep << shiftIDToShiftTypeID[t] << ":";
			    rep << allTasks[i]<< "|";
			}
		}
		rep << std::endl;
		return rep.str();
	}

	bool operator!=(const Rotation& rot2){
		if(nurseId_ != rot2.nurseId_) return true;
		if(firstDay_ != rot2.firstDay_) return true;
		if(length_ != rot2.length_) return true;
		for(const std::pair<int,int>& p: rot2.shifts_) {
			if(shifts_[p.first] != p.second) return true;
		}
		return false;
	}
	bool operator==(const Rotation& rot2){
		if(nurseId_ != rot2.nurseId_) return false;
		if(firstDay_ != rot2.firstDay_) return false;
		if(length_ != rot2.length_) return false;
		for(const std::pair<int,int>& p: rot2.shifts_) {
			if(shifts_[p.first] != p.second) return false;
		}
		return true;
	}

	//Compare rotations on index
	//
	static bool compareId(const Rotation& rot1, const Rotation& rot2);

	//Compare rotations on cost
	//
	static bool compareCost(const Rotation& rot1, const Rotation& rot2);

	//Compare rotations on dual cost
	//
	static bool compareDualCost(const Rotation& rot1, const Rotation& rot2);

	// Returns true if both rotations are disjoint PLUS ONE DAY INBETWEEN !!
	bool isDisjointWith(Rotation& rot2){
		return ((this->firstDay_+this->length_ < rot2.firstDay_ - 1)
		|| (rot2.firstDay_+rot2.length_ < this->firstDay_ - 1)
	);
}

// Returns true if both rotations are disjoint for the shift
bool isShiftDisjointWith(Rotation& rot2){
	if( (this->firstDay_+this->length_ < rot2.firstDay_) || (rot2.firstDay_+rot2.length_ < this->firstDay_) )
	return true;

	int commomFirstDay = std::max(firstDay_, rot2.firstDay_),
	    commomLastDay=std::min(firstDay_+length_, rot2.firstDay_+rot2.length_)-1;
	for(int k=commomFirstDay; k<=commomLastDay; ++k)
	if(shifts_[k] == rot2.shifts_[k]) return false;

	return true;
}
};


//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

enum MySolverType { S_SCIP, S_CLP, S_Gurobi, S_Cplex, S_CBC };
static std::map<std::string,MySolverType> MySolverTypesByName =
{{"CLP",S_CLP},{"Gurobi",S_Gurobi},{"Cplex",S_Cplex},{"CBC",S_CBC},{"SCIP",S_SCIP}};

class MasterProblem : public Solver, public PrintSolution{
	//allows RotationPricer to access all private arguments and methods of MasterProblem
	friend class RotationPricer;
	friend class DiveBranchingRule;
public:
	// Specific constructor and destructor
	MasterProblem(Scenario* pScenario, Demand* pDemand,
		Preferences* pPreferences, std::vector<State>* pInitState, MySolverType solver);
	~MasterProblem();

	//solve the rostering problem
	double solve(std::vector<Roster> solution = {});

	//solve the rostering problem or just the relaxation(root node)
	double solve(std::vector<Roster> solution, bool rebuild);

	// Solve with parameters
	double solve(SolverParam parameters, std::vector<Roster> solution = {});

	//Resolve the problem with another demand and keep the same preferences
	//
	double resolve(Demand* pDemand, SolverParam parameters, std::vector<Roster> solution = {});

	//get the pointer to the model
	Modeler* getModel(){
		return pModel_;
	}

	//override PrintSolution virtual method
	void save(std::vector<int>& weekIndices, std::string outdir);
	void printCurrentSol();

	//get a reference to the restsPerDay_ for a Nurse
	inline vector2D<MyVar*>& getRestsPerDay(Nurse* pNurse){
		return restsPerDay_[pNurse->id_];
	}

	// build the, possibly fractional, roster corresponding to the solution
	// currently stored in the model
	vector3D<double> getFractionalRoster() ;

	//------------------------------------------------
	// Solution with rolling horizon process
	//------------------------------------------------

	// relax/unrelax the integrality constraints of the variables corresponding to input days
	void relaxDays(std::vector<bool> isRelax);
	void unrelaxDays(std::vector<bool> isUnrelax);

	// fix/unfix all the variables corresponding to the input vector of days
	void fixDays(std::vector<bool> isFixDay);
	void unfixDays(std::vector<bool> isUnfixDay);

	// fix/unfix all the variables corresponding to the input vector of nurses
	void fixNurses(std::vector<bool> isFixNurse);
	void unfixNurses(std::vector<bool> isUnfixNurse);

	// Solve the problem with a method that allows for a warm start
	double rollingSolve(SolverParam parameters, int firstDay);

	// Special solve function for LNS
	// It is a priori the same as a regular, but it might be modified if needed
	double LNSSolve(SolverParam parameters);

	//---------------------------------------------------------------------------
	//
	// Methods required to implement stabilization in the column generation
	//
	//---------------------------------------------------------------------------

	// STAB
	// Multiply the upper bound of the input variable by the input factor
	void multiplyUbInSolver(MyVar* pVar, OsiSolverInterface* solver, double factor);
	// Set the bound of the input variable to the input value
	void updateVarUbInSolver(MyVar* pVar, OsiSolverInterface* solver, double value);

	// STAB
	// Set the cost of the input variable to the input value
	void updateVarCostInSolver(MyVar* pVar, OsiSolverInterface* solver, double value);

	// STAB
	// Update all the upper bounds of the stabilization variables by multiplying
	// them by an input factor
	void stabUpdateBound(OsiSolverInterface* solver, double factor);

	// STAB
	// Update all the costs of the stabilization variables to the values
	// corresponding dual variables with a small margin in input
	void stabUpdateCost(OsiSolverInterface* solver, double margin);

	// STAB
	// Check the stopping criterion of the relaxation solution specific to the
	// the stabilization
	// The point is that current solution can be infeasible if  stabilization
	// variables are non zero
	bool stabCheckStoppingCriterion();

	// STAB: compute the lagrangian bound
	//
	double computeLagrangianBound(double objVal,double sumRedCost);

	// STAB: reset the costs and bounds of the stabilization variables
	//
	void stabResetBoundAndCost(OsiSolverInterface* solver, SolverParam param);

	/*
	* Solving parameterdoubles
	*/
	const char* PB_NAME = "GenCol";
	int solvingTime;

public:
	// getter/setters
	//
  std::vector<MyVar*> getMinWorkedDaysVars() {return minWorkedDaysVars_;}
    std::vector<MyVar*> getMaxWorkedDaysVars() {return maxWorkedDaysVars_;}
    std::vector<MyVar*> getMaxWorkedWeekendVars() {return maxWorkedWeekendVars_;}
	vector3D<MyVar*> getOptDemandVars() {return optDemandVars_;}

private:
	Modeler* pModel_;
	vector2D<int> positionsPerSkill_;//link positions to skills
	vector2D<int> skillsPerPosition_;//link skills to positions
	MyPricer* pPricer_;//prices the rotations
	MyTree* pTree_;//store the tree information
	MyBranchingRule* pRule_; //choose the variables on which we should branch
	MySolverType solverType_; //which solver is used

    std::vector<MyVar*> initialStateVars_; //stores all the initial rotations finishing on the first day
	vector3D<MyVar*> restsPerDay_; //stores all the arcs that are resting on a day for each nurse

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // IMPORTANT:  DANS LA VERSION MERINIO, LES 'NURSES' SONT REMPLACEES PAR DES EMPLOYES ET LES 'WORKEDDAYS' PAR DES HEURES TRAVAILLEES
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	* Variables
	*/
	vector2D<MyVar*> columnVars_; //binary variables for the columns

	vector2D<MyVar*> restingVars_; //binary variables for the resting arcs in the rotation network
	vector3D<MyVar*> longRestingVars_; //binary variables for the resting arcs in the rotation network

    std::vector<MyVar*> minWorkedDaysVars_; //count the number of missing worked days per nurse
    std::vector<MyVar*> maxWorkedDaysVars_; //count the number of exceeding worked days per nurse
    std::vector<MyVar*> maxWorkedWeekendVars_; //count the number of exceeding worked weekends per nurse

    std::vector<MyVar*> minWorkedDaysAvgVars_; //count the number of missing worked days from average per nurse
    std::vector<MyVar*> maxWorkedDaysAvgVars_; // count the number of exceeding worked days from average per nurse
    std::vector<MyVar*> maxWorkedWeekendAvgVars_; //count the number of exceeding worked weekends from average per nurse

    std::vector<MyVar*> minWorkedDaysContractAvgVars_; //count the number of missing worked days from average per contract
    std::vector<MyVar*> maxWorkedDaysContractAvgVars_; // count the number of exceeding worked days from average per contract
    std::vector<MyVar*> maxWorkedWeekendContractAvgVars_; //count the number of exceeding worked weekends from average per contract

	vector3D<MyVar*> optDemandVars_; //count the number of missing nurse to reach the optimal
	vector3D<MyVar*> numberOfNursesByPositionVars_; // count the number of nurses by position on each day, shift
	vector4D<MyVar*> skillsAllocVars_; //makes the allocation of the skills

	/*
	* Constraints
	*/
	//transmission of the flow on the resting nodes
	//initialization of the flow constraint at the first position of each restFlowCons_[i] (i=nurse)
	vector2D<MyCons*> restFlowCons_;
	//transmission of the flow on the working nodes
	//end of the flow constraint at the last position of each workFlowCons_[i] (i=nurse)
	vector2D<MyCons*> workFlowCons_;

	std::vector<MyCons*> minWorkedDaysCons_; //count the number of missing worked days per nurse
    std::vector<MyCons*> maxWorkedDaysCons_; //count the number of exceeding worked days per nurse
    std::vector<MyCons*> maxWorkedWeekendCons_; //count the number of exceeding worked weekends per nurse
	MyCons* sumMaxWorkedWeekendCons_;	// count the total number of weekends that will be penalized


    std::vector<MyCons*> minWorkedDaysAvgCons_; //count the number of missing worked days from average per nurse
    std::vector<MyCons*> maxWorkedDaysAvgCons_; // count the number of exceeding worked days from average per nurse
    std::vector<MyCons*> maxWorkedWeekendAvgCons_; //count the number of exceeding worked weekends from average per nurse

    std::vector<MyCons*> minWorkedDaysContractAvgCons_; //count the number of missing worked days from average per contract
    std::vector<MyCons*> maxWorkedDaysContractAvgCons_; // count the number of exceeding worked days from average per contract
    std::vector<MyCons*> maxWorkedWeekendContractAvgCons_; //count the number of exceeding worked weekends from average per contract

	vector3D<MyCons*> minDemandCons_; //ensure a minimal coverage per day, per shift, per skill
	vector3D<MyCons*> optDemandCons_; //count the number of missing nurse to reach the optimal
	vector3D<MyCons*> numberOfNursesByPositionCons_; //ensure there are enough nurses for numberOfNursesByPositionVars_
	vector3D<MyCons*> feasibleSkillsAllocCons_; // ensures that each nurse works with the good skill

	// STAB
	// Stabilization variables for each constraint
	// Two variables are needed for equality constraints and one for inequalities
	// The constraints on average values are not stabilized yet
	// The position and allocation constraints do not require stabilization
	vector2D<MyVar*> stabRestFlowPlus_;
	vector2D<MyVar*> stabRestFlowMinus_;
	vector2D<MyVar*> stabWorkFlowPlus_;
	vector2D<MyVar*> stabWorkFlowMinus_;

    std::vector<MyVar*> stabMinWorkedDaysPlus_;
    std::vector<MyVar*> stabMaxWorkedDaysMinus_;
    std::vector<MyVar*> stabMaxWorkedWeekendMinus_;

	vector3D<MyVar*> stabMinDemandPlus_; //ensure a minimal coverage per day, per shift, per skill
	vector3D<MyVar*> stabOptDemandPlus_; //count the number of missing nurse to reach the optimal

	// vectors of booleans indicating whether some above constraints are present
	// in the model
  std::vector<bool> isMinWorkedDaysAvgCons_;
    std::vector<bool> isMaxWorkedDaysAvgCons_;
    std::vector<bool> isMaxWorkedWeekendAvgCons_;

    std::vector<bool> isMinWorkedDaysContractAvgCons_, isMaxWorkedDaysContractAvgCons_, isMaxWorkedWeekendContractAvgCons_;

	/*
	* Methods
	*/

	// Initialize the solver at construction
	void initializeSolver(MySolverType solverType);

	// Main method to build the rostering problem for a given input
	void build(SolverParam parameters);

	// Initialization of the master problem with/without solution
	void initialize(SolverParam parameters, std::vector<Roster> solution={});

	// Provide an initial solution to the solver
	void initializeSolution(std::vector<Roster> solution);

	//solve method to catch execption
	void solveWithCatch();

	//solve a solution in the output
	void storeSolution();

	//Create a new rotation variable
	//add the correct constraints and coefficients for the nurse i working on a rotation
	//if s=-1, the nurse works on all shifts
	MyVar* addRotation(Rotation& rotation, const char* baseName, bool coreVar = false);

	//compute and add the last rotation finishing on the day just before the first one
	Rotation computeInitStateRotation(LiveNurse* pNurse);

	//update the demand with a new one of the same size
	//change the rhs of the constraints minDemandCons_ and optDemandCons_
	void updateDemand(Demand* pDemand);

	//get the cost of all chosen rotations in solution sol for a certain CostType
	double getRotationCosts(CostType costType = TOTAL_COST, bool initStateRotation = false);

	//get the cost of all vars in the solution for a certain CostType
	double getRotationCosts(CostType costType, const std::vector<MyVar*>& vars);

	/* Build each set of constraints - Add also the coefficient of a column for each set */
	void buildRotationCons(SolverParam parameters);
	int addRotationConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int k, bool firstDay, bool lastDay);
	void buildMinMaxCons(SolverParam parameters);
	int addMinMaxConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int nbDays, int nbWeekends);
	void buildSkillsCoverageCons(SolverParam parameters);
	int addSkillsCoverageConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int k, int s=-1);

	/* Display functions */
  std::string costsConstrainstsToString();
    std::string allocationToString(bool printInteger = true);
    std::string coverageToString(bool printInteger = true);
    std::string workedWeekendsToString(bool printInteger = true);
};

#endif /* MASTERPROBLEM_H_ */
