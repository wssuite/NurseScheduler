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
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "tools/Tools.h"
#include "solvers/mp/modeler/Modeler.h"
#include "solvers/mp/MasterProblem.h"

// Node that correspond to branching on a set of columns (diving)
//
struct ColumnsNode : public MyNode {
  ColumnsNode(int index,
              MyNode *pParent,
              const std::vector<PPattern> &patterns) :
      MyNode(index, pParent), patterns_(patterns) {}

  std::string write() const {
    std::stringstream out;
    out << "ColumnsNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",#Cols="
        << patterns_.size() << ")";
    return out.str();
  }

  // vector of the columns on which we have branched.
  const std::vector<PPattern> patterns_;
};

// Node that correspond to branching on a resting day
//
struct RestNode : public MyNode {
  RestNode(int index, MyNode *pParent, PLiveNurse pNurse, int day, bool rest) :
      MyNode(index, pParent), pNurse_(pNurse), day_(day), rest_(rest) {}

  std::string write() const {
    std::stringstream out;
    out << "RestNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse="
        << pNurse_->num_ << ",Day=" << day_ << ",";
    if (rest_) std::cout << "rest";
    else
      std::cout << "work";
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
  ShiftNode(int index,
            MyNode *pParent,
            PLiveNurse pNurse,
            int day,
            bool work,
            const std::vector<int> &forbiddenShifts) :
      MyNode(index, pParent),
      pNurse_(pNurse),
      day_(day),
      work_(work),
      forbiddenShifts_(forbiddenShifts) {}

  std::string write() const {
    std::stringstream out;
    out << "ShiftNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse="
        << pNurse_->num_ << ",Day=" << day_ << ",";
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
// The variable can relate to cover constraints, total shifts or total week-ends
//
struct PenaltyNode : public MyNode {
  PenaltyNode(int index,
              MyNode *pParent,
              PLiveNurse pNurse,
              int day,
              bool work,
              const std::vector<int> &forbiddenShifts) :
      MyNode(index, pParent),
      pNurse_(pNurse),
      day_(day),
      work_(work),
      forbiddenShifts_(forbiddenShifts) {}

  std::string write() const {
    std::stringstream out;
    out << "PenaltyNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse="
        << pNurse_->num_ << ",Day=" << day_ << ",";
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
  VarNode(int index, MyNode *pParent, MyVar *var, double lb, double ub) :
      MyNode(index, pParent),
      var_(var),
      lb_(lb),
      ub_(ub) {}

  std::string write() const {
    std::stringstream out;
    out << "VarNode: (depth=" << depth_ << ",LB=" << bestLB_;
    out << ",Var=" << var_->name_ << ",LB=" << lb_ << ",UB=" << ub_ << ")";
    return out.str();
  }

  // number of nurse on which we have branched. pNumberOfNurses_ can be 0
  const MyVar *var_;
  const double lb_, ub_;
};

struct CoverageNode : public MyNode {
  CoverageNode(int index, MyNode *pParent,
      const char * cutName, double lhs, double rhs) :
      MyNode(index, pParent),
    cutName_(cutName),
    lhs_(lhs),
    rhs_(rhs) {}

  std::string write() const {
    std::stringstream out;
    out << "CoverageNode: (depth=" << depth_ << ",LB=" << bestLB_;
    out << ",Var=" << cutName_ << ",LHS=" << lhs_ << ",RHS="
        << rhs_ << ")";
    return out.str();
  }

  // number of nurse on which we have branched. pNumberOfNurses_ can be 0
  const char * cutName_;
  const double lhs_, rhs_;
};


struct RestTree : public MyTree {
  RestTree(PScenario pScenario, PDemand pDemand, double epsilon);

  void addForbiddenShifts(PLiveNurse pNurse,
                          std::set<std::pair<int, int> > *forbidenShifts);

  void pushBackNewVarNode(MyVar *var, double lb, double ub) {
    VarNode *node = new VarNode(tree_.size(), currentNode_, var, lb, ub);
    pushBackNode(node);
  }

  void pushBackNewDemandNode(const char * cutName, double lhs, double rhs) {
    CoverageNode *node =
        new CoverageNode(tree_.size(), currentNode_, cutName, lhs, rhs);
    pushBackNode(node);
  }

  void pushBackNewColumnsNode(const std::vector<PPattern> &patterns) {
    ColumnsNode *node = new ColumnsNode(tree_.size(), currentNode_, patterns);
    pushBackNode(node);
  }

  void pushBackNewRestNode(PLiveNurse pNurse, int day, bool rest) {
    RestNode
        *node = new RestNode(tree_.size(), currentNode_, pNurse, day, rest);
    pushBackNode(node);
  }

  void pushBackNewShiftNode(PLiveNurse pNurse,
                            int day,
                            bool work,
                            const std::vector<int> &shifts) {
    ShiftNode *node =
        new ShiftNode(tree_.size(), currentNode_, pNurse, day, work, shifts);
    pushBackNode(node);
  }

  bool is_columns_node() const {
    return dynamic_cast<ColumnsNode *>(currentNode_) != nullptr;
  }

  bool continueDiving() const {
    if (is_columns_node()) return true;
    return (diveDepth_ <= std::max(min_depth_ + 10, diveLength_));
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
  virtual bool is_disjoint(PPattern col1, PPattern col2) const = 0;
};

struct DayDisjointComparator : public ColumnsComparator {
  bool is_disjoint(PPattern col1, PPattern col2) const {
    return col1->isDisjointWith(col2);
  }
};

struct ShiftDisjointComparator : public ColumnsComparator {
  bool is_disjoint(PPattern col1, PPattern col2) const {
    return col1->nurseNum() != col2->nurseNum()
        && col1->isShiftDisjointWith(col2);
  }
};

class DiveBranchingRule;

struct ScoreVar {
  explicit ScoreVar(const DiveBranchingRule *pRule) : pRule(pRule) {}
  virtual ~ScoreVar() {}

  virtual double score(PLiveNurse pNurse,
                       int day,
                       const std::vector<int> &shifts,
                       const std::vector<double> & values) const = 0;

  const DiveBranchingRule *pRule;
};

struct ScoreVarCloseHalf : ScoreVar {
  explicit ScoreVarCloseHalf(const DiveBranchingRule *pRule,
                             double advantage = .1) :
      ScoreVar(pRule), advantage(advantage) {}

  double score(PLiveNurse pNurse,
               int day,
               const std::vector<int> &shifts,
               const std::vector<double> & values) const override;

  double advantage;
};

struct ScoreVarBestExpectedLBImprovement : ScoreVar {
  explicit ScoreVarBestExpectedLBImprovement(const DiveBranchingRule *pRule) :
      ScoreVar(pRule) {}

  double score(PLiveNurse pNurse,
               int day,
               const std::vector<int> &shifts,
               const std::vector<double> & values) const override;
};

class DiveBranchingRule : public MyBranchingRule {
 public:
  DiveBranchingRule(MasterProblem *master, RestTree *tree, const char *name);
  virtual ~DiveBranchingRule() {}

  /* compute branching decisions */
  bool branching_candidates(MyBranchingCandidate *candidate) override;

  /* branch on a number of nurses working on a shift */
  void branchOnNumberNurses(MyBranchingCandidate *candidate);

  /* Branch on opt var */
  void branchOnOptDemand(MyBranchingCandidate *candidate);

  /* branch on a set of resting arcs */
  void branchOnRestDay(MyBranchingCandidate *candidate);

  /* branch on a set of shifts */
  void branchOnShifts(MyBranchingCandidate *candidate);

  /* compute fixing decisions */
  bool column_candidates(MyBranchingCandidate *candidate) override;

  /* Choose columns */
  std::vector<MyVar *> chooseColumns(
      const std::vector<std::pair<MyVar *, double>> &candidates,
      std::vector<PPattern> *patterns,
      double *maxValue,
      const ColumnsComparator &comparator);

  /* compare columns */
  static bool compareColumnCloseToInt(const std::pair<MyVar *, double> &obj1,
                                      const std::pair<MyVar *, double> &obj2);

  static bool compareColumnCloseTo5(const std::pair<MyVar *, double> &obj1,
                                    const std::pair<MyVar *, double> &obj2);

  MasterProblem *getMaster() const;
  Modeler *getModel() const;
  const std::vector<PDualCosts> &getPDualCosts() const;

 protected:
  // Pointer to the master problem to link the master and the sub problems
  MasterProblem *pMaster_;

  // Pointer to the rest tree
  RestTree *tree_;

  // Pointers to the data
  Modeler *pModel_;

  // store dual costs
  std::vector<PDualCosts> pDualCosts_;

  std::unique_ptr<ScoreVar> scoreFunc_;

  // branching decisions to test
  typedef void (DiveBranchingRule::*BranchFunc)(MyBranchingCandidate *);
  std::vector<BranchFunc> branchFunctions_ =
      {&DiveBranchingRule::branchOnRestDay,
       &DiveBranchingRule::branchOnShifts};

  virtual void buildRestNodesCut(MyBranchingCandidate *candidate,
                                 PLiveNurse pNurse,
                                 int day,
                                 bool forceRest,
                                 MyBranchingNode *restNode,
                                 MyBranchingNode *workNode) const;

  void deactivateColumns(MyBranchingCandidate *candidate,
                         int nurseNum,
                         int day,
                         std::vector<int> forbiddenShifts,
                         MyBranchingNode *forbiddenNode,
                         MyBranchingNode *complementaryNode) const;


  void randomSwapLastChildren(MyBranchingCandidate *candidate) {
    // Here : random choice to decide the order of the children
    if (Tools::randomInt(0, 1)) {
      candidate->swapLastChildren();
      tree_->swapLastNodes();
    }
  }
};

class RosterBranchingRule : public DiveBranchingRule {
 public:
  RosterBranchingRule(MasterProblem *master, RestTree *tree, const char *name) :
      DiveBranchingRule(master, tree, name) {}

  virtual ~RosterBranchingRule() {}

 protected:
  void buildRestNodesCut(MyBranchingCandidate *candidate,
                         PLiveNurse pNurse,
                         int day,
                         bool forceRest,
                         MyBranchingNode *restNode,
                         MyBranchingNode *workNode) const override {}
};

#endif  // SRC_SOLVERS_MP_TREEMANAGER_H_
