/*
* MasterProblem.h
*
*  Created on: 3 fevr. 2015
*      Author: samuel
*/

#ifndef MASTERPROBLEM_H_
#define MASTERPROBLEM_H_

/* Inheritance */
#include <jmorecfg.h>
#include <solvers/mp/rcspp/RCGraph.h>
#include "solvers/Solver.h"

/* Tools include */
#include "tools/MyTools.h"

/* My includes */
#include "data/Nurse.h"
#include "solvers/mp/modeler/Modeler.h"
#include "OsiSolverInterface.hpp"

enum CostType {TOTAL_COST, CONS_SHIFTS_COST, CONS_WORKED_DAYS_COST,
    COMPLETE_WEEKEND_COST, PREFERENCE_COST, REST_COST};

struct DualCosts{
public:

  DualCosts(const vector2D<double> & workedShiftsCosts,
            const std::vector<double> & startWorkCosts,
            const std::vector<double> & endWorkCosts,
            double workedWeekendCost):
              workedShiftsCosts_(workedShiftsCosts), startWorkCosts_(startWorkCosts),
              endWorkCosts_(endWorkCosts), workedWeekendCost_(workedWeekendCost) {}

  // GETTERS
  //
  inline int nDays() const { return startWorkCosts_.size(); }
  inline double workedDayShiftCost(int day, int shift){return (workedShiftsCosts_[day][shift-1]);}
  inline double startWorkCost(int day){return (startWorkCosts_[day]);}
  inline double endWorkCost(int day){return (endWorkCosts_[day]);}
  inline double workedWeekendCost(){return workedWeekendCost_;}


protected:

  // Indexed by : (day, shift) !! 0 = shift 1 !!
  vector2D<double> workedShiftsCosts_;

  // Indexed by : day
  std::vector<double> startWorkCosts_;

  // Indexed by : day
  std::vector<double> endWorkCosts_;

  // Reduced cost of the weekends
  double workedWeekendCost_;

};

struct Pattern;
typedef std::shared_ptr<Pattern> PPattern;

struct Pattern {
    Pattern(int nurseId, int firstDay, int length):
      nurseId_(nurseId), firstDay_(firstDay), length_(length) {}
      Pattern(const std::vector<double>& pattern):
      nurseId_((int)pattern[0]), firstDay_((int)pattern[1]), length_((int)pattern[2]) {}
    virtual ~Pattern() {}

    virtual bool equals(PPattern pat) const {
      if (nurseId_ != pat->nurseId_) return false;
      if (firstDay_ != pat->firstDay_) return false;
      if (length_ != pat->length_) return false;
      for (int k=firstDay_; k<firstDay_+length_; ++k) {
        if (getShift(k) != pat->getShift(k)) return false;
      }
      return true;
    }

    // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest) !!
    virtual bool isDisjointWith(PPattern pat, bool needRest=true) const {
      return ((firstDay_+length_ < pat->firstDay_ - needRest)
              || (pat->firstDay_+pat->length_ < firstDay_ - needRest)
      );
    };

    // Returns true if both columns are disjoint PLUS ONE DAY INBETWEEN (needRest) !!
    virtual bool isShiftDisjointWith(PPattern pat, bool needRest=true) const {
      if(isDisjointWith(pat, needRest))
        return true;

      int commomFirstDay = std::max(firstDay_, pat->firstDay_),
          commomLastDay=std::min(firstDay_+length_, pat->firstDay_+pat->length_)-1;
      for(int k=commomFirstDay; k<=commomLastDay; ++k)
        if(getShift(k) == pat->getShift(k)) return false;

      return true;
    };

    virtual std::string toString(int nbDays = -1, std::vector<int> shiftIDToShiftTypeID={}) const {
      return "";
    }

    virtual int getShift(int day) const = 0;

    // need to be able to write the pattern as a vector and to create a new one from it
    virtual std::vector<double> getCompactPattern() const {
      std::vector<double> compact;
      compact.push_back(nurseId_);
      compact.push_back(firstDay_);
      compact.push_back(length_);
      return compact;
    }

    int nurseId_;
    int firstDay_;
    int length_;

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
    double solve(const SolverParam& parameters, std::vector<Roster> solution = {});

    //Resolve the problem with another demand and keep the same preferences
    //
    double resolve(Demand* pDemand, const SolverParam& parameters, std::vector<Roster> solution = {});

    // needs to be specialized: add a colum  to the master from a solution of the subproblem
    virtual MyVar* addColumn(int nurseId, const RCSolution& solution) = 0;

    // retrieve the object represented ny the  vector pattern
    virtual PPattern getPattern(const std::vector<double>& pattern) const = 0;

    // throw an error if pattern is already present as an active column
    void checkIfPatternAlreadyPresent(const std::vector<double>& pattern) const;

    virtual const std::vector<MyVar*>& getRestVarsPerDay(LiveNurse* pNurse, int day) const = 0;

    //get the pointer to the model
    Modeler* getModel(){
      return pModel_;
    }

    //override PrintSolution virtual method
    void save(std::vector<int>& weekIndices, std::string outdir);
    void printCurrentSol();

    // build the, possibly fractional, roster corresponding to the solution
    // currently stored in the model
    vector3D<double> getFractionalRoster() ;

    // build a DualCosts structure
    DualCosts buildDualCosts(LiveNurse* pNurse) const;

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
    double rollingSolve(const SolverParam& parameters, int firstDay);

    // Special solve function for LNS
    // It is a priori the same as a regular, but it might be modified if needed
    double LNSSolve(const SolverParam& parameters);

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
    virtual void stabUpdateBound(OsiSolverInterface* solver, double factor);

    // STAB
    // Update all the costs of the stabilization variables to the values
    // corresponding dual variables with a small margin in input
    virtual void stabUpdateCost(OsiSolverInterface* solver, double margin);

    // STAB
    // Check the stopping criterion of the relaxation solution specific to the
    // the stabilization
    // The point is that current solution can be infeasible if  stabilization
    // variables are non zero
    virtual bool stabCheckStoppingCriterion() const;

    // STAB: compute the lagrangian bound
    //
    virtual double computeLagrangianBound(double objVal,double sumRedCost) const;

    // STAB: reset the costs and bounds of the stabilization variables
    //
    virtual void stabResetBoundAndCost(OsiSolverInterface* solver,
        const SolverParam& parameters);

    /*
    * Solving parameterdoubles
    */
    const char* PB_NAME = "GenCol";
    int solvingTime;

    // getter/setters
    //
    std::vector<MyVar*> getMinWorkedDaysVars() {return minWorkedDaysVars_;}
    std::vector<MyVar*> getMaxWorkedDaysVars() {return maxWorkedDaysVars_;}
    std::vector<MyVar*> getMaxWorkedWeekendVars() {return maxWorkedWeekendVars_;}
    vector3D<MyVar*> getOptDemandVars() {return optDemandVars_;}

  protected:
    Modeler* pModel_;
    vector2D<int> positionsPerSkill_;//link positions to skills
    vector2D<int> skillsPerPosition_;//link skills to positions
    MyPricer* pPricer_;//prices the rotations
    MyTree* pTree_;//store the tree information
    MyBranchingRule* pRule_; //choose the variables on which we should branch
    MySolverType solverType_; //which solver is used

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // IMPORTANT:  WORKED DAYS CAN ALSO BE USED WORKED HOURS
    //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    * Variables
    */
    vector2D<MyVar*> columnVars_; //binary variables for the columns

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
    std::vector<MyVar*> stabMinWorkedDaysPlus_;
    std::vector<MyVar*> stabMaxWorkedDaysMinus_;
    std::vector<MyVar*> stabMaxWorkedWeekendMinus_;

    vector3D<MyVar*> stabMinDemandPlus_; //ensure a minimal coverage per day, per shift, per skill
    vector3D<MyVar*> stabOptDemandPlus_; //count the number of missing nurse to reach the optimal

    // vectors of booleans indicating whether some above constraints are present
    // in the model
    std::vector<bool> isMinWorkedDaysAvgCons_,
                      isMaxWorkedDaysAvgCons_,
                      isMaxWorkedWeekendAvgCons_,
                      isMinWorkedDaysContractAvgCons_,
                      isMaxWorkedDaysContractAvgCons_,
                      isMaxWorkedWeekendContractAvgCons_;

    /*
    * Methods
    */

    // Initialize the solver at construction
    void initializeSolver(MySolverType solverType);

    // Main method to build the rostering problem for a given input
    virtual void build(const SolverParam& parameters);

    // Initialization of the master problem with/without solution
    void initialize(const SolverParam& parameters, std::vector<Roster> solution={});

    // Provide an initial solution to the solver
    virtual void initializeSolution(const std::vector<Roster>& solution) = 0;

    //solve method to catch execption
    void solveWithCatch();

    //solve a solution in the output
    void storeSolution();

    // return the costs of all active columns associated to the type
    virtual double getColumnsCost(CostType costType, bool justHistoricalCosts) const = 0;

    //update the demand with a new one of the same size
    //change the rhs of the constraints minDemandCons_ and optDemandCons_
    void updateDemand(Demand* pDemand);

    /* Build each set of constraints - Add also the coefficient of a column for each set */
    void buildMinMaxCons(const SolverParam& parameters);
    int addMinMaxConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int nbDays, int nbWeekends);
    void buildSkillsCoverageCons(const SolverParam& parameters);
    int addSkillsCoverageConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int k, int s=-1);

    /* retrieve the dual values */
    virtual vector2D<double> getShiftsDualValues(LiveNurse*  pNurse) const;
    virtual std::vector<double> getStartWorkDualValues(LiveNurse* pNurse) const;
    virtual std::vector<double> getEndWorkDualValues(LiveNurse* pNurse) const;
    virtual double getWorkedWeekendDualValue(LiveNurse* pNurse) const;

    /* Display functions */
    std::string costsConstrainstsToString();
    std::string allocationToString(bool printInteger = true);
    std::string coverageToString(bool printInteger = true);
    std::string workedWeekendsToString(bool printInteger = true);
};

#endif /* MASTERPROBLEM_H_ */
