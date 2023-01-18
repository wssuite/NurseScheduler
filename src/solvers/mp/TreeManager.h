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

#ifndef SRC_SOLVERS_MP_TREEMANAGER_H_
#define SRC_SOLVERS_MP_TREEMANAGER_H_

#include <algorithm>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "tools/Tools.h"
#include "solvers/mp/modeler/Modeler.h"
#include "solvers/mp/MasterProblem.h"


struct Score {
  explicit Score(double sc = -DBL_MAX) : score(sc) {}
  double score;
  PLiveNurse pNurse = nullptr;
  int k = -1, s = -1, sk = -1;
  vector<int> forbiddenShifts;
  MyVar *var;
};

// Node that correspond to branching on a set of columns (diving)
//
struct ColumnsNode : public MyNode {
  ColumnsNode(MyPNode pParent,
              std::vector<PColumn> Columns) :
      MyNode(std::move(pParent)), columns_(std::move(Columns)) {}

  std::string write() const {
    std::stringstream out;
    out << "ColumnsNode: (" << getInfo()
        << ", #Cols=" << columns_.size() << ")";
    return out.str();
  }

  // vector of the columns on which we have branched.
  const std::vector<PColumn> columns_;
};

// Node that correspond to branching on a resting day
//
struct RestNode : public MyNode {
  RestNode(MyPNode pParent, PLiveNurse pNurse, int day, bool rest) :
      MyNode(std::move(pParent)), pNurse_(std::move(pNurse)),
      day_(day), rest_(rest) {}

  std::string write() const {
    std::stringstream out;
    out << "RestNode: (" << getInfo()
        << ", Nurse=" << pNurse_->num_ << ", Day=" << day_ << ", ";
    if (rest_) out << "rest";
    else
      out << "work";
    out << ")";
    return out.str();
  }

  // nurse and day on which we have branched for rest or work. pNurse_ can be 0
  const PLiveNurse pNurse_;
  const int day_;
  const bool rest_;
};

// Node that correspond to branching on working shifts
//
struct ShiftNode : public MyNode {
  ShiftNode(MyPNode pParent,
            PLiveNurse pNurse,
            int day,
            bool work,
            const std::vector<int> &forbiddenShifts) :
      MyNode(std::move(pParent)),
      pNurse_(pNurse),
      day_(day),
      work_(work),
      forbiddenShifts_(forbiddenShifts) {}

  std::string write() const {
    std::stringstream out;
    out << "ShiftNode: (" << getInfo()
        << ", Nurse=" << pNurse_->num_ << ", Day=" << day_ << ", ";
    if (work_) out << "work";
    else
      out << "N.A.";
    out << "). Forbidden shifts:";
    for (int s : forbiddenShifts_)
      out << " " << s;
    return out.str();
  }

  // nurse and day on which we have branched for rest or work. pNurse_ can be 0
  const PLiveNurse pNurse_;
  const int day_;
  const bool work_;
  std::vector<int> forbiddenShifts_;
};

// Node that corresponds to branching on a penalized original variable
// The variable can relate to cover constraints, total shifts or total weekends
//
struct PenaltyNode : public MyNode {
  PenaltyNode(MyPNode pParent,
              PLiveNurse pNurse,
              int day,
              bool work,
              std::vector<int> forbiddenShifts) :
      MyNode(std::move(pParent)),
      pNurse_(std::move(pNurse)),
      day_(day),
      work_(work),
      forbiddenShifts_(std::move(forbiddenShifts)) {}

  std::string write() const {
    std::stringstream out;
    out << "PenaltyNode: (" << getInfo()
        << ", Nurse=" << pNurse_->num_ << ", Day=" << day_ << ", ";
    if (work_) out << "work";
    else
      out << "N.A.";
    out << "). Forbidden shifts:";
    for (int s : forbiddenShifts_)
      out << " " << s;
    return out.str();
  }

  // nurse and day on which we have branched for rest or work. pNurse_ can be 0
  const PLiveNurse pNurse_;
  const int day_;
  const bool work_;
  std::vector<int> forbiddenShifts_;
};

struct VarNode : public MyNode {
  VarNode(MyPNode pParent, MyVar *var, double lb, double ub) :
      MyNode(std::move(pParent)),
      var_(var),
      lb_(lb),
      ub_(ub) {}

  std::string write() const {
    std::stringstream out;
    out << "VarNode: (" << getInfo()
        << ", Var=" << var_->name_ << ", LB=" << lb_ << ", UB=" << ub_ << ")";
    return out.str();
  }

  // number of nurse on which we have branched. pNumberOfNurses_ can be 0
  const MyVar *var_;
  const double lb_, ub_;
};

struct CoverageNode : public MyNode {
  CoverageNode(MyPNode pParent,
      string cutName, double lhs, double rhs) :
      MyNode(std::move(pParent)),
    cutName_(std::move(cutName)),
    lhs_(lhs),
    rhs_(rhs) {}

  std::string write() const {
    std::stringstream out;
    out << "CoverageNode: (" << getInfo()
        << ", Var=" << cutName_ << ", LHS=" << lhs_ << ", RHS=" << rhs_ << ")";
    return out.str();
  }

  // number of nurse on which we have branched. pNumberOfNurses_ can be 0
  const string cutName_;
  const double lhs_, rhs_;
};


class RCBranchingCandidate : public MyBranchingCandidate {
 public:
  explicit RCBranchingCandidate(MyPNode currentNode):
      MyBranchingCandidate(), currentNode_(std::move(currentNode)) {}

  void addVarNode(int index, MyVar *var, double lb, double ub) {
    children_[index]->addPNode(
        std::make_shared<VarNode>(currentNode_, var, lb, ub));
  }

  void addDemandNode(int index, const char * cutName, double lhs, double rhs) {
    children_[index]->addPNode(
        std::make_shared<CoverageNode>(currentNode_, cutName, lhs, rhs));
  }

  void addColumnsNode(int index, const std::vector<PColumn> &Columns) {
    children_[index]->addPNode(
        std::make_shared<ColumnsNode>(currentNode_, Columns));
  }

  void addRestNode(int index, PLiveNurse pNurse, int day, bool rest) {
    children_[index]->addPNode(
        std::make_shared<RestNode>(currentNode_, pNurse, day, rest));
  }

  void addShiftNode(int index,
                    PLiveNurse pNurse,
                    int day,
                    bool work,
                    const std::vector<int> &shifts) {
    children_[index]->addPNode(
        std::make_shared<ShiftNode>(currentNode_, pNurse, day, work, shifts));
  }

 protected:
  MyPNode currentNode_;
};


struct RestTree : public MyTree {
  RestTree(PScenario pScenario, PDemand pDemand,
           double epsilon, bool printCurrentNode);

  void addForbiddenShifts(PLiveNurse pNurse,
                          std::set<std::pair<int, int> > *forbidenShifts);

  bool isColumnsNode() const {
    return dynamic_cast<ColumnsNode*>(currentNode_.get()) != nullptr;
  }

  bool continueDiving() const {
    if (isColumnsNode()) return true;
    return (diveDepth_ <= std::max(minDepth_ + 10, diveLength_));
  }

  void reset();

  void updateStats(MyNode *node);

  std::string writeBranchStats() const;

  std::string writeOneStat(
      std::string name,
      const std::vector<std::pair<int, double> > &stats) const;

 protected:
  PScenario pScenario_;
  PDemand pDemand_;
  // Tree statistics
  // each pair represents first the number of times of the branching and
  // second the increase of the LB
  std::vector<std::pair<int, double>> statsRestByDay_, statsWorkByDay_,
      statsRestByNurse_, statsWorkByNurse_;
  std::pair<int, double> statsCols_;
};

struct ColumnsComparator {
  virtual bool is_disjoint(PColumn col1, PColumn col2) const {
    return true;
  }
};

struct DayDisjointComparator : public ColumnsComparator {
  bool is_disjoint(PColumn col1, PColumn col2) const {
    return col1->isDisjointWith(col2, true);
  }
};

struct ShiftDisjointComparator : public ColumnsComparator {
  bool is_disjoint(PColumn col1, PColumn col2) const {
    return col1->nurseNum() != col2->nurseNum()
        && col1->isShiftDisjointWith(col2, true);
  }
};

class DiveBranchingRule;

struct ScoreVar {
  explicit ScoreVar(const DiveBranchingRule *pRule) : pRule(pRule) {}
  virtual ~ScoreVar() {}

  virtual double score(PLiveNurse pNurse,
                       int day,
                       const std::vector<int> &shifts,
                       const std::vector<double> & values,
                       double baseScore = 0) const = 0;

  const DiveBranchingRule *pRule;
};

struct ScoreVarCloseHalf : ScoreVar {
  explicit ScoreVarCloseHalf(const DiveBranchingRule *pRule,
                             double advantage = .1) :
      ScoreVar(pRule), weekendAdvantage_(advantage) {}

  double score(PLiveNurse pNurse,
               int day,
               const std::vector<int> &shifts,
               const std::vector<double> & values,
               double baseScore = 0) const override;

  double weekendAdvantage_;
};

struct ScoreVarCloseHalfWeekendDecrement : ScoreVar {
  explicit ScoreVarCloseHalfWeekendDecrement(const DiveBranchingRule *pRule) :
      ScoreVar(pRule) {}

  double score(PLiveNurse pNurse,
               int day,
               const std::vector<int> &shifts,
               const std::vector<double> & values,
               double baseScore = 0) const override;
};

class DiveBranchingRule : public MyBranchingRule {
 public:
  DiveBranchingRule(MasterProblem *master,
                    RestTree *tree,
                    const char *name,
                    bool randomSwapOfChilfren = false);
  virtual ~DiveBranchingRule() {}

  /* compute branching decisions */
  bool branching_candidates(
      int nCandidates,
      std::vector<MyPBranchingCandidate> *candidates) override;

  /* branch on a number of nurses working on a shift */
  void branchOnNumberNurses(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates);

  /* Branch on opt var */
  void branchOnOptDemand(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates);

  /* branch on a set of resting arcs */
  void branchOnRestDay(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates);

  /* branch on a set of shifts */
  void branchOnShifts(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates);

  /* compute fixing decisions */
  bool column_node(std::vector<MyPBranchingCandidate> *candidates) override;

  /* Choose columns */
  std::vector<MyVar *> chooseColumns(
      const std::vector<std::pair<MyVar *, double>> &candidates,
      std::vector<PColumn> *columns,
      double *maxValue,
      const ColumnsComparator &comparator);

  /* compare columns */
  static bool compareColumnCloseToInt(const std::pair<MyVar *, double> &obj1,
                                      const std::pair<MyVar *, double> &obj2);

  static bool compareColumnCloseTo5(const std::pair<MyVar *, double> &obj1,
                                    const std::pair<MyVar *, double> &obj2);

  MasterProblem *getMaster() const;
  Modeler *getModel() const;

 protected:
  // Pointer to the master problem to link the master and the sub problems
  MasterProblem *pMaster_;

  // Pointer to the branching tree
  RestTree *pTree_;

  // Pointers to the data
  Modeler *pModel_;

  // random swap
  // if true randomly swap children before inserting them in the tree
  bool randomSwapOfChilfren_;

  std::unique_ptr<ScoreVar> scoreFunc_;

  // branching decisions to test
  typedef void (DiveBranchingRule::*BranchFunc)(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates);
  std::vector<BranchFunc> branchFunctions_ =
      {&DiveBranchingRule::branchOnRestDay,
       &DiveBranchingRule::branchOnShifts};

  virtual void buildRestNodesCut(const MyPBranchingCandidate &candidate,
                                 PLiveNurse pNurse,
                                 int day,
                                 bool forceRest,
                                 const MyPBranchingNode &restNode,
                                 const MyPBranchingNode &workNode) const;

  void deactivateColumns(const MyPBranchingCandidate &candidate,
                         int nurseNum,
                         int day,
                         std::vector<int> forbiddenShifts,
                         const MyPBranchingNode &forbiddenNode,
                         const MyPBranchingNode &complementaryNode) const;


  void randomSwapLastChildrenIfEnable(const MyPBranchingCandidate &candidate) {
    if (!randomSwapOfChilfren_) return;
    // Here : random choice to decide the order of the children
    if (Tools::randomInt(0, 1)) {
      candidate->swapLastChildren();
      pTree_->swapLastNodes();
    }
  }

  // compute a base score for each nurse in order to advantage certain nurses
  // in the branching selection process
  std::vector<double> computeNurseBaseScore(double coeff);
};

class RosterBranchingRule : public DiveBranchingRule {
 public:
  RosterBranchingRule(MasterProblem *master, RestTree *tree, const char *name) :
      DiveBranchingRule(master, tree, name) {}

  virtual ~RosterBranchingRule() {}

 protected:
  void buildRestNodesCut(const MyPBranchingCandidate &candidate,
                         PLiveNurse pNurse,
                         int day,
                         bool forceRest,
                         const MyPBranchingNode &restNode,
                         const MyPBranchingNode &workNode) const override {}
};

#endif  // SRC_SOLVERS_MP_TREEMANAGER_H_
