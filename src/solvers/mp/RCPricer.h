/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#ifndef SRC_SOLVERS_MP_RCPRICER_H_
#define SRC_SOLVERS_MP_RCPRICER_H_

#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "tools/MyTools.h"
#include "MasterProblem.h"
#include "solvers/mp/sp/SubProblem.h"
#include "solvers/mp/modeler/Modeler.h"

//---------------------------------------------------------------------------
//
// C l a s s   R o t a t i o n P r i c e r
//
// Contains the pricer that generates new columns.
//
//---------------------------------------------------------------------------
class RCPricer : public MyPricer {
 public:
  RCPricer(MasterProblem *master, const char *name, const SolverParam &param);

  virtual ~RCPricer();

  /* perform pricing */
  // at root node, after_fathom and backtracked are true
  std::vector<MyVar *> pricing(double bound = 0,
                               bool before_fathom = false,
                               bool after_fathom = false,
                               bool backtracked = false) override;

  const std::vector<double> &getLastMinOptimalReducedCost() const override {
    return minOptimalReducedCosts_;
  }

  bool isLastRunOptimal() const override {
    return optimal_;
  }

 protected:
  // DATA - instance-related data
  MasterProblem *pMaster_;
  PScenario pScenario_;
  int nbDays_;
  Modeler *pModel_;
  std::vector<PLiveNurse> nursesToSolve_;

  // One subproblem list per contract because the consecutive same shift
  // constraints vary by contract.
  std::map<PConstContract, std::list<SubProblem *>> subProblems_;

  // DATA - Solutions, rotations, etc.
  std::vector<MyVar *> allNewColumns_;

  // True if all subproblems have been solved to optimality
  bool optimal_;

  // Stats on the number of subproblems solved and successfully solved
  int nbSPTried_;
  int nbSPSolvedWithSuccess_;

  // SETTINGS - Options for forbidden shifts, nurses, starting days, etc.
  std::set<std::pair<int, int> > forbiddenShifts_;
  std::set<int> forbiddenNursesIds_;
  std::set<int> forbiddenStartingDays_;
  std::set<int> forbiddenEndingDays_;

  // SETTINGS - Options for the neighborhood. need of an original to reset at
  // the end of each node when optimality has been reached.
  // Here, we could have all the parameters as fields
  // but they would be too many.
  std::vector<int> currentSubproblemStrategy_;  // by nurse

  // store the min reduced cost find for each subproblem solved
  std::vector<double> minOptimalReducedCosts_;
  double minReducedCost_;

  // mutex for concurrency (can be locked several times by the same thread)
  std::recursive_mutex m_subproblem_;

 public:
  // METHODS - Solutions, rotations, etc.
  void resetSolutions() {
    allNewColumns_.clear();
    forbiddenShifts_.clear();
    optimal_ = true;  // will be set to false whenever possible
    nbSPSolvedWithSuccess_ = 0;
    nbSPTried_ = 0;
    Tools::initVector<double>(&minOptimalReducedCosts_,
        pMaster_->getNbNurses(),
        -DBL_MAX);
    minReducedCost_ = 0;
  }

  // Retrieve the right subproblem
  SubProblem *retriveSubproblem(PLiveNurse pNurse);

  // release the subproblem
  void releaseSubproblem(PLiveNurse pNurse, SubProblem *subProblem);

  // METHODS - Forbidden shifts, nurses, starting days, etc.
  //
  // !!! WARNING !!! : SOME METHODS ARE NOT YET IMPLEMENTED IN THE SUBPROBLEM
  // (ALTHOUGH THE NECESSARY STRUCTURES MAY ALREADY BE THERE !!!
  //
  // update current nurse strategy and reduced costs based on the solutions
  // found.
  // return true if the strategy has been updated
  bool updateCurrentStrategyAndRedCost(PLiveNurse pNurse,
                                       const std::vector<RCSolution> &solutions,
                                       bool disjointForbidden);

  // add nurses to nursesToSolve_ in the reverse order
  template<typename T>
  void reversePushBackNurses(T *array) {
    std::reverse(array->begin(), array->end());  // reverse array
    nursesToSolve_.insert(nursesToSolve_.end(), array->begin(), array->end());
  }

  // Shifts
  void forbidShift(int k, int s) override {
    forbiddenShifts_.insert(std::pair<int,
                                      int>(k, s));
  }
  void forbidShifts(const std::set<std::pair<int, int> > &shifts) {
    for (auto s : shifts)forbidShift(s.first, s.second);
  }
  void authorizeShift(int k, int s) {
    forbiddenShifts_.erase(std::pair<int,
                                     int>(k,
                                          s));
  }
  void clearForbiddenShifts() { forbiddenShifts_.clear(); }
  // Nurses
  void forbidNurse(int nurseNum) override {
    forbiddenNursesIds_.insert(nurseNum);
  }
  void forbidNurses(const std::set<int> &nurses) {
    for (auto n : nurses)forbidNurse(n);
  }
  void authorizeNurse(int nurseNum) override {
    forbiddenNursesIds_.erase(nurseNum);
  }
  void clearForbiddenNurses() override { forbiddenNursesIds_.clear(); }
  // Starting days
  void forbidStartingDay(int k) override { forbiddenStartingDays_.insert(k); }
  void forbidStartingDays(const std::set<int> &days) {
    for (int d : days)forbidStartingDay(d);
  }
  void authorizeStartingDay(int k) override { forbiddenStartingDays_.erase(k); }
  void clearForbiddenStartingDays() override { forbiddenStartingDays_.clear(); }
  // Ending days
  void forbidEndingDay(int k) override { forbiddenEndingDays_.insert(k); }
  void forbidEndingDays(const std::set<int> &days) {
    for (int d : days)forbidEndingDay(d);
  }
  void authorizeEndingDay(int k) { forbiddenEndingDays_.erase(k); }
  void clearForbiddenEndingDays() { forbiddenEndingDays_.clear(); }

  // Test functions
  bool isShiftForbidden(int k, int n) {
    return (forbiddenShifts_.find(std::pair<int, int>(k, n))
        != forbiddenShifts_.end());
  }
  bool isNurseForbidden(int n) {
    return (forbiddenNursesIds_.find(n) != forbiddenNursesIds_.end());
  }
  bool isStartingDayForbidden(int k) {
    return (forbiddenStartingDays_.find(k) != forbiddenStartingDays_.end());
  }
  bool isEndingDayForbidden(int k) {
    return (forbiddenEndingDays_.find(k) != forbiddenEndingDays_.end());
  }

 protected:
  /*
   * Methods
   */

  // get the duals values per day and per shift for a nurse
  vector2D<double> getShiftsDualValues(PLiveNurse pNurse);
  std::vector<double> getStartWorkDualValues(PLiveNurse pNurse);
  std::vector<double> getEndWorkDualValues(PLiveNurse pNurse);
  double getWorkedWeekendDualValue(PLiveNurse pNurse);

  // compute some forbidden shifts from the lasts solutions and add them to
  // forbidden shifts.
  // only the shift from the best solution with a negative dual costs are
  // added to the forbidden shifts.
  bool addForbiddenShifts(const std::vector<RCSolution> &solutions,
                          const DualCosts &duals);

  // Add the rotations to the master problem
  int addColumnsToMaster(int nurseNum, std::vector<RCSolution> *solutions);

  // Sort the rotations that just were generated for a nurse.
  // Default option is sort by increasing reduced cost but we
  // could try something else (involving disjoint columns for ex.)
  void sortGeneratedSolutions(std::vector<RCSolution> *solutions) const;

  // Set the subproblem options depending on the parameters
//  void setSubproblemOptions(vector<SolveOption> &options,
//                            int &maxRotationLengthForSubproblem,
//                            pLiveNurse pNurse);

  int nb_int_solutions_;

  // Print functions
  //
  void printStatSPSolutions();

  // ----------------------------------------
  //
  // DBG - DEBUG & STATS DATA AND FUNCTIONS
  //
  // ----------------------------------------
  // DBG / Stats
  double timeInExSubproblems_ = 0;
  double timeForS_ = 0;
  double timeForN_ = 0;
  double timeForNL_ = 0;
  int nbSubproblems_ = 0;
  int nbExSubproblems_ = 0;
  int nbS_ = 0;
  int nbN_ = 0;
  int nbNL_ = 0;
  std::minstd_rand rand_;

  // DBG functions
  void generateRandomForbiddenStartingDays();
};

#endif  // SRC_SOLVERS_MP_RCPRICER_H_
