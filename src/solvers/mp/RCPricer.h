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

#include "tools/Tools.h"
#include "MasterProblem.h"
#include "solvers/mp/sp/SubProblem.h"
#include "solvers/mp/modeler/Modeler.h"

//---------------------------------------------------------------------------
//
// C l a s s   R C P r i c e r
//
// Contains the pricer that generates new columns.
//
//---------------------------------------------------------------------------
class RCPricer : public MyPricer {
 public:
  RCPricer(MasterProblem *pMaster, const char *name, const SolverParam &param);

  virtual ~RCPricer();

  /* perform pricing */
  // at root node, after_fathom and backtracked are true
  std::vector<MyVar *> pricing(double bound = 0,
                               bool before_fathom = false,
                               bool after_fathom = false,
                               bool backtracked = false) override;

  void updateParameters(bool useMoreTime) override;

  double getLastMinReducedCost() const override {
    return minReducedCost_;
  }

  const std::vector<double> &getLastReducedCostLBs() const override {
    return reducedCostLBs_;
  }

  bool isLastRunOptimal() const override {
    return optimal_;
  }

  bool isLastRunLowerBounded() const override {
    return lowerBounded_;
  }

  void initNursesAvailabilities() override;

 protected:
  // DATA - instance-related data
  MasterProblem *pMaster_;
  PScenario pScenario_;
  int nbDays_;
  Modeler *pModel_;
  std::vector<PLiveNurse> nursesToSolve_;

  // One subproblem list per contract because the consecutive same shift
  // constraints vary by contract.
  std::map<PLiveNurse, SubProblem*> subProblems_;

  // DATA - Solutions, rotations, etc.
  std::vector<MyVar *> allNewColumns_;

  // True if all subproblems have been solved to optimality
  bool optimal_;

  // True if all subproblems have a lower bound
  bool lowerBounded_;

  // Stats on the number of subproblems solved and successfully solved
  int nbSPTried_;
  int nSPSolvedWithSuccess_;

  // SETTINGS - Options for forbidden shifts, nurses, starting days, etc.
  std::vector<std::set<std::pair<int, int> >> nursesForbiddenShifts_;
  std::set<std::pair<int, int> > forbiddenShifts_;
  std::set<int> forbiddenNursesIds_;

  // store the min reduced cost find for each subproblem solved
  std::vector<double> reducedCostLBs_;
  double minReducedCost_;

  // mutex for concurrency (can be locked several times by the same thread)
  std::recursive_mutex m_subproblem_;

 public:
  MasterProblem* pMaster() const {
    return pMaster_;
  }

  // METHODS - Solutions, rotations, etc.
  void resetSolutions() {
    allNewColumns_.clear();
    forbiddenShifts_.clear();
    // will be set to false whenever possible
    optimal_ = true;
    lowerBounded_ = false;
    nSPSolvedWithSuccess_ = 0;
    nbSPTried_ = 0;
    Tools::initVector(&reducedCostLBs_, pMaster_->nNurses(), -DBL_MAX);
    minReducedCost_ = 0;
  }

  // Retrieve the right subproblem
  SubProblem *retrieveSubproblem(const PLiveNurse &pNurse,
                                 const SubProblemParam &spParam);

  SubProblem *buildSubproblem(const PLiveNurse &pNurse,
                              const SubProblemParam &spParam) const;

  void computeCost(Column *col);

  // METHODS - Forbidden shifts, nurses, starting days, etc.
  //
  // !!! WARNING !!! : SOME METHODS ARE NOT YET IMPLEMENTED IN THE SUBPROBLEM
  // (ALTHOUGH THE NECESSARY STRUCTURES MAY ALREADY BE THERE !!!
  //
  // update reduced costs based on the solutions found.
  void updateRedCost(SubProblem *pSP,
                     const std::vector<RCSolution> &solutions,
                     bool disjointForbidden);

  // add nurses to nursesToSolve_ in the reverse order
  template<typename T>
  void reversePushBackNurses(vector<T> *vec) {
    std::reverse(vec->begin(), vec->end());  // reverse array
    nursesToSolve_ = Tools::appendVectors(nursesToSolve_, *vec);
  }

  template<typename T>
  void reverseOrderPushBackNurses(vector<T> *vec,
                                  const std::map<T, int> &order) {
    std::stable_sort(vec->begin(), vec->end(),
                     [&order](const T &v1, const T &v2) {
                       return order.at(v1) > order.at(v2);
                     });
    nursesToSolve_ = Tools::appendVectors(nursesToSolve_, *vec);
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

  // Test functions
  bool isShiftForbidden(int k, int n) {
    return (forbiddenShifts_.find(std::pair<int, int>(k, n))
        != forbiddenShifts_.end());
  }
  bool isNurseForbidden(int n) {
    return (forbiddenNursesIds_.find(n) != forbiddenNursesIds_.end());
  }

 protected:
  /*
   * Methods
   */

  // compute some forbidden shifts from the lasts solutions and add them to
  // forbidden shifts.
  // only the shift from the best solution with a negative dual costs are
  // added to the forbidden shifts.
  bool addForbiddenShifts(const std::vector<RCSolution> &solutions);

  // Add the rotations to the master problem
  int addColumnsToMaster(int nurseNum, std::vector<RCSolution> *solutions);

  // Sort the rotations that just were generated for a nurse.
  // Default option is sort by increasing reduced cost but we
  // could try something else (involving disjoint columns for ex.)
  void sortGeneratedSolutions(std::vector<RCSolution> *solutions);

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

 private:
  std::minstd_rand rdm_;
};

#endif  // SRC_SOLVERS_MP_RCPRICER_H_
