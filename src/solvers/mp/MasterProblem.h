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
    COMPLETE_WEEKEND_COST, PREFERENCE_COST, REST_COST,
    MIN_DAYS_COST, MAX_DAYS_COST, MAX_WEEKEND_COST};

struct DualCosts{
public:

  DualCosts(const vector2D<double> & workedShiftsCosts,
            const std::vector<double> & startWorkCosts,
            const std::vector<double> & endWorkCosts,
            double workedWeekendCost, double constant=0):
              workedShiftsCosts_(workedShiftsCosts), startWorkCosts_(startWorkCosts),
              endWorkCosts_(endWorkCosts), workedWeekendCost_(workedWeekendCost), constant_(constant) {}

  // GETTERS
  //
  inline int nDays() const { return startWorkCosts_.size(); }
  inline double workedDayShiftCost(int day, int shift) const{
    if(!shift) return 0;
    return (workedShiftsCosts_[day][shift-1]);
  }
  inline double startWorkCost(int day) const{return (startWorkCosts_[day]);}
  inline double endWorkCost(int day) const{return (endWorkCosts_[day]);}
  inline double workedWeekendCost() const{return workedWeekendCost_;}
  inline double constant() const { return constant_; }


protected:

  // Indexed by : (day, shift) !! 0 = shift 1 !!
  vector2D<double> workedShiftsCosts_;

  // Indexed by : day
  std::vector<double> startWorkCosts_;

  // Indexed by : day
  std::vector<double> endWorkCosts_;

  // Reduced cost of the weekends
  double workedWeekendCost_;

  // constant part of the reduced cost
  double constant_;

};

struct Pattern;
typedef std::shared_ptr<Pattern> PPattern;

struct Pattern {
    Pattern(int firstDay, int length, int nurseId = -1, double cost = DBL_MAX, double dualCost = DBL_MAX):
      nurseId_(nurseId), firstDay_(firstDay), length_(length), id_(s_count++),
      cost_(cost), dualCost_(dualCost) {}
    Pattern(const Pattern& pat, int nurseId):
        nurseId_(nurseId), firstDay_(pat.firstDay_), length_(pat.length_), id_(pat.id_),
        cost_(pat.cost_), dualCost_(pat.dualCost_) {
      if(pat.nurseId_ != nurseId_){
        cost_ = DBL_MAX;
        dualCost_ = DBL_MAX;
      }
    }
    Pattern(const std::vector<double>& pattern):
      nurseId_((int)pattern[0]), firstDay_((int)pattern[1]), length_((int)pattern[2]), id_(s_count++),
      cost_(DBL_MAX), dualCost_(DBL_MAX){}
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

    //Id of the rotation
    //
    const long id_;

    // Cost
    //
    double cost_;

    // Dual cost as found in the subproblem
    //
    double dualCost_;

    //Compare rotations on index
    //
    static bool compareId(Pattern* pat1, Pattern* pat2);

    //Compare rotations on cost
    //
    static bool compareCost(Pattern* pat1, Pattern* pat2);

    //Compare rotations on dual cost
    //
    static bool compareDualCost(Pattern* pat1, Pattern* pat2);

  private:
    //count rotations
    //
    static unsigned long s_count;

};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------

class MasterProblem : public Solver, public PrintSolution{
  public:
    // Specific constructor and destructor
    MasterProblem(PScenario pScenario, PDemand pDemand,
      PPreferences pPreferences, std::vector<State>* pInitState, MySolverType solver);
    ~MasterProblem();

    //solve the rostering problem
    double solve(std::vector<Roster> solution = {});

    //solve the rostering problem or just the relaxation(root node)
    double solve(std::vector<Roster> solution, bool rebuild);

    // Solve with parameters
    double solve(const SolverParam& parameters, std::vector<Roster> solution = {});

    //Resolve the problem with another demand and keep the same preferences
    //
    double resolve(PDemand pDemand, const SolverParam& parameters, std::vector<Roster> solution = {});

    // needs to be specialized: add a colum  to the master from a solution of the subproblem
    virtual MyVar* addColumn(int nurseId, const RCSolution& solution) = 0;

    // retrieve the object represented ny the  vector pattern
    virtual PPattern getPattern(const std::vector<double>& pattern) const = 0;

    // throw an error if pattern is already present as an active column
    void checkIfPatternAlreadyPresent(const std::vector<double>& pattern) const;

    virtual const std::vector<MyVar*>& getRestVarsPerDay(PLiveNurse pNurse, int day) const = 0;

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
    DualCosts buildDualCosts(PLiveNurse pNurse) const;

    // return the value V used to choose the number of columns on which to branch.
    // Choose as many columns as possible such that: sum (1 - value(column)) < V
    virtual double getBranchColumnValueMax() const {
      return std::max (1.0, pScenario_->nbWeeks_ / 2.0);
    }

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
    vector3D<MyVar*> optDemandVars_; //count the number of missing nurse to reach the optimal
    vector3D<MyVar*> numberOfNursesByPositionVars_; // count the number of nurses by position on each day, shift
    vector4D<MyVar*> skillsAllocVars_; //makes the allocation of the skills

    /*
    * Constraints
    */

    vector3D<MyCons*> minDemandCons_; //ensure a minimal coverage per day, per shift, per skill
    vector3D<MyCons*> optDemandCons_; //count the number of missing nurse to reach the optimal
    vector3D<MyCons*> numberOfNursesByPositionCons_; //ensure there are enough nurses for numberOfNursesByPositionVars_
    vector3D<MyCons*> feasibleSkillsAllocCons_; // ensures that each nurse works with the good skill

    // STAB
    // Stabilization variables for each constraint
    // Two variables are needed for equality constraints and one for inequalities
    // The constraints on average values are not stabilized yet
    // The position and allocation constraints do not require stabilization

    vector3D<MyVar*> stabMinDemandPlus_; //ensure a minimal coverage per day, per shift, per skill
    vector3D<MyVar*> stabOptDemandPlus_; //count the number of missing nurse to reach the optimal

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

    virtual double getMinDaysCost() const = 0;
    virtual double getMaxDaysCost() const = 0;
    virtual double getMaxWeekendCost() const = 0;

    //update the demand with a new one of the same size
    //change the rhs of the constraints minDemandCons_ and optDemandCons_
    void updateDemand(PDemand pDemand);

    /* Build each set of constraints - Add also the coefficient of a column for each set */
    void buildSkillsCoverageCons(const SolverParam& parameters);
    int addSkillsCoverageConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, const Pattern& pat) const;

    /* retrieve the dual values */
    virtual vector2D<double> getShiftsDualValues(PLiveNurse pNurse) const;
    virtual std::vector<double> getStartWorkDualValues(PLiveNurse pNurse) const;
    virtual std::vector<double> getEndWorkDualValues(PLiveNurse pNurse) const;
    virtual double getWorkedWeekendDualValue(PLiveNurse pNurse) const;
    virtual double getConstantDualvalue(PLiveNurse pNurse) const;

    /* Display functions */
    std::string costsConstrainstsToString();
    std::string allocationToString(bool printInteger = true);
    std::string coverageToString(bool printInteger = true);
};

#endif /* MASTERPROBLEM_H_ */
