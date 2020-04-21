//
// Created by antoine legrain on 2020-04-15.
//

#ifndef NURSESCHEDULER_ROTATIONMP_H
#define NURSESCHEDULER_ROTATIONMP_H


#include "MasterProblem.h"
#include "solvers/mp/rcspp/RCGraph.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

struct Rotation: Pattern {

    // Specific constructors and destructors
    //
    Rotation(std::map<int,int> shifts, int nurseId = -1, double cost = DBL_MAX, double dualCost = DBL_MAX) :
        Pattern(nurseId, -1, shifts.size()),
        shifts_(shifts), id_(s_count), cost_(cost),
        consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
        dualCost_(dualCost), timeDuration_(shifts.size())
    {
      ++s_count;
      firstDay_ = INT_MAX;
      for(auto itS = shifts.begin(); itS != shifts.end(); ++itS)
        if(itS->first < firstDay_) firstDay_ = itS->first;
    };

    Rotation(int firstDay, std::vector<int> shiftSuccession, int nurseId = -1, double cost = DBL_MAX, double dualCost = DBL_MAX) :
        Pattern(nurseId, firstDay, shiftSuccession.size()),
        id_(s_count), cost_(cost),
        consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
        dualCost_(dualCost), timeDuration_(shiftSuccession.size())
    {
      ++s_count;
      for(int k=0; k<length_; k++)
        shifts_[firstDay+k] = shiftSuccession[k];
    }

    Rotation(const std::vector<double>& compactPattern) :
        Pattern(compactPattern),
        id_(s_count), cost_(DBL_MAX),
        consShiftsCost_(0), consDaysWorkedCost_(0), completeWeekendCost_(0), preferenceCost_(0), initRestCost_(0),
        dualCost_(DBL_MAX), timeDuration_((int)compactPattern.back())
    {
      ++s_count;
      for(int k=0; k<length_; k++)
        shifts_[firstDay_+k] = (int)compactPattern[k+3];
    }

    Rotation(const Rotation& rotation, int nurseId) :
        Pattern(nurseId, rotation.firstDay_, rotation.length_),
        shifts_(rotation.shifts_), id_(rotation.id_), cost_(rotation.cost_),
        consShiftsCost_(rotation.consShiftsCost_), consDaysWorkedCost_(rotation.consDaysWorkedCost_),
        completeWeekendCost_(rotation.completeWeekendCost_), preferenceCost_(rotation.preferenceCost_), initRestCost_(rotation.initRestCost_),
        dualCost_(rotation.dualCost_), timeDuration_(rotation.timeDuration_)
    {
      if(rotation.nurseId_ != nurseId_){
        cost_ = DBL_MAX;
        dualCost_ = DBL_MAX;
      }
    }

    ~Rotation(){}

    int getShift(int day) const override {
      return shifts_.at(day);
    }

    //count rotations
    //
    static unsigned int s_count;

    // Shifts to be performed
    //
    std::map<int,int> shifts_;

    //Id of the rotation
    //
    long id_;

    // Cost
    //
    double cost_;
    double consShiftsCost_ , consDaysWorkedCost_, completeWeekendCost_, preferenceCost_, initRestCost_ ;

    // Dual cost as found in the subproblem
    //
    double dualCost_;

    // Level of the branch and bound tree where the rotation has been generated
    //
    int treeLevel_=0;

    // Time duration (in a certain unit: day, hours, half-hours, ...)

    int timeDuration_;

    //compact the rotation in a vector
    std::vector<double> getCompactPattern() const override {
      std::vector<double> pattern = Pattern::getCompactPattern();
      for(std::pair<int,int> p: shifts_) pattern.push_back(p.second);
      pattern.push_back(timeDuration_);
      return pattern;
    }

    //Compute the cost of a rotation
    //
    void computeCost(Scenario* pScenario, Preferences* pPreferences, const std::vector<LiveNurse*>& liveNurses, int horizon);

    //Compute the dual cost of a rotation
    //
    void checkDualCost(DualCosts& costs);

    // calcule le nombre d'heures d'une rotation
    void computeTimeDuration(Scenario* pScenario) {
      timeDuration_ = 0;
      for (std::pair<int, int> p: shifts_) {
        timeDuration_ += pScenario->timeDurationToWork_[p.second];
      }
    }

    std::string toString(int nbDays = -1, std::vector<int> shiftIDToShiftTypeID={}) const override;

    //Compare rotations on index
    //
    static bool compareId(const Rotation& rot1, const Rotation& rot2);

    //Compare rotations on cost
    //
    static bool compareCost(const Rotation& rot1, const Rotation& rot2);

    //Compare rotations on dual cost
    //
    static bool compareDualCost(const Rotation& rot1, const Rotation& rot2);
};

//-----------------------------------------------------------------------------
//
//  C l a s s   M a s t e r P r o b l e m
//
// Build and solve the master problem of the column generation scheme
//
//-----------------------------------------------------------------------------
class RotationMP: public MasterProblem {
  public:
    RotationMP(Scenario* pScenario, Demand* pDemand, Preferences* pPreferences, std::vector<State> *pInitState,
               MySolverType solver);
    virtual ~RotationMP();

    PPattern getPattern(const std::vector<double>& pattern) const override;

    MyVar* addColumn(int nurseId, const RCSolution& solution) override;

    // STAB
    // Update all the upper bounds of the stabilization variables by multiplying
    // them by an input factor
    void stabUpdateBound(OsiSolverInterface* solver, double factor) override;

    // STAB
    // Update all the costs of the stabilization variables to the values
    // corresponding dual variables with a small margin in input
    void stabUpdateCost(OsiSolverInterface* solver, double margin) override;

    // STAB
    // Check the stopping criterion of the relaxation solution specific to the
    // the stabilization
    // The point is that current solution can be infeasible if  stabilization
    // variables are non zero
    bool stabCheckStoppingCriterion() const override ;

    // STAB: compute the lagrangian bound
    //
    double computeLagrangianBound(double objVal,double sumRedCost) const override ;

    // STAB: reset the costs and bounds of the stabilization variables
    //
    void stabResetBoundAndCost(OsiSolverInterface* solver,
                               const SolverParam& parameters) override ;

    //get a reference to the restsPerDay_ for a Nurse
    inline const std::vector<MyVar*>& getRestVarsPerDay(LiveNurse* pNurse, int day) const override {
      return restsPerDay_[pNurse->id_][day];
    }

    const std::vector<MyVar*>& getMinWorkedDaysVars() const {return minWorkedDaysVars_;}
    const std::vector<MyVar*>& getMaxWorkedDaysVars() const {return maxWorkedDaysVars_;}
    const std::vector<MyVar*>& getMaxWorkedWeekendVars() const {return maxWorkedWeekendVars_;}

  protected:
    // Main method to build the rostering problem for a given input
    void build(const SolverParam& parameters) override ;

    // Provide an initial solution to the solver. If empty, add artificial columns
    void initializeSolution(const std::vector<Roster>& solution) override ;

    //Create a new rotation variable
    //add the correct constraints and coefficients for the nurse i working on a rotation
    //if s=-1, the nurse works on all shifts
    MyVar* addRotation(const Rotation& rotation, const char* baseName, bool coreVar = false);

    //compute and add the last rotation finishing on the day just before the first one
    Rotation computeInitStateRotation(LiveNurse* pNurse);

    /* Build each set of constraints - Add also the coefficient of a column for each set */
    void buildRotationCons(const SolverParam& parameters);
    int addRotationConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs,
                             int i, int k, bool firstDay, bool lastDay);

    void buildMinMaxCons(const SolverParam& parameters);
    int addMinMaxConsToCol(std::vector<MyCons*>& cons, std::vector<double>& coeffs, int i, int nbDays, int nbWeekends);

    // return the costs of all active columns associated to the type
    double getColumnsCost(CostType costType, bool historicalCosts) const override ;
    double getColumnsCost(CostType costType, const std::vector<MyVar*>& vars) const;

    double getMinDaysCost() const override;
    double getMaxDaysCost() const override;
    double getMaxWeekendCost() const override;

    /* retrieve the dual values */
    vector2D<double> getShiftsDualValues(LiveNurse*  pNurse) const override;
    std::vector<double> getStartWorkDualValues(LiveNurse* pNurse) const override ;
    std::vector<double> getEndWorkDualValues(LiveNurse* pNurse) const override ;
    double getWorkedWeekendDualValue(LiveNurse* pNurse) const override;

    /*
    * Variables
    */
    vector3D<MyVar*> restsPerDay_; //stores all the arcs that are resting on a day for each nurse
    vector2D<MyVar*> restingVars_; //binary variables for the resting arcs in the rotation network
    vector3D<MyVar*> longRestingVars_; //binary variables for the resting arcs in the rotation network
    std::vector<MyVar*> initialStateVars_; //stores all the initial rotations finishing on the first day

    std::vector<MyVar*> minWorkedDaysVars_; //count the number of missing worked days per nurse
    std::vector<MyVar*> maxWorkedDaysVars_; //count the number of exceeding worked days per nurse
    std::vector<MyVar*> maxWorkedWeekendVars_; //count the number of exceeding worked weekends per nurse

    std::vector<MyVar*> minWorkedDaysAvgVars_; //count the number of missing worked days from average per nurse
    std::vector<MyVar*> maxWorkedDaysAvgVars_; // count the number of exceeding worked days from average per nurse
    std::vector<MyVar*> maxWorkedWeekendAvgVars_; //count the number of exceeding worked weekends from average per nurse

    std::vector<MyVar*> minWorkedDaysContractAvgVars_; //count the number of missing worked days from average per contract
    std::vector<MyVar*> maxWorkedDaysContractAvgVars_; // count the number of exceeding worked days from average per contract
    std::vector<MyVar*> maxWorkedWeekendContractAvgVars_; //count the number of exceeding worked weekends from average per contract


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

    // vectors of booleans indicating whether some above constraints are present
    // in the model
    std::vector<bool> isMinWorkedDaysAvgCons_,
        isMaxWorkedDaysAvgCons_,
        isMaxWorkedWeekendAvgCons_,
        isMinWorkedDaysContractAvgCons_,
        isMaxWorkedDaysContractAvgCons_,
        isMaxWorkedWeekendContractAvgCons_;
};


#endif //NURSESCHEDULER_ROTATIONMP_H
