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

#ifndef SRC_SOLVERS_MP_MODELER_MODELER_H_
#define SRC_SOLVERS_MP_MODELER_MODELER_H_

#include <stdio.h>

#include <cassert>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include <set>
#include <list>
#include <map>
#include <memory>
#include <utility>

#include "solvers/Solver.h"
#include "solvers/mp/modeler/Stabilization.h"
#include "tools/Tools.h"

/*
 * Var types
 */
enum VarType { VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY };

/*
 * Rule Search Strategy
 */

enum SearchStrategy {
  BestFirstSearch,
  BreadthFirstSearch,
  DepthFirstSearch,
  HighestGapFirst
};

/*
 * My Modeling objects
 * If the object is added to the vector objects_ of the Modeler,
 * the modeler will also delete it at the end.
 */
struct MyObject {
  explicit MyObject(const char *name) {
    char *tName = new char[255];
    strncpy(tName, name, 255);
    name_ = tName;
  }

  MyObject(const MyObject &myObject) : num_(myObject.num_) {
    char *name2 = new char[255];
    strncpy(name2, myObject.name_, 255);
    name_ = name2;
  }

  virtual ~MyObject() {
    delete[] name_;
  }

  // for the map of objects
  int operator<(const MyObject &m) const { return num_ < m.num_; }

  int num() const { return num_; }

  const char *name_;

 private:
  int num_ = -1;

  void num(int num) {
    if (num < -1)
      Tools::throwError("Negative num for MyObject: %d", num);
    num_ = num;
  }

  friend class Modeler;
};

struct MyVar : public MyObject {
  MyVar(const char *name,
        int index,
        double cost,
        VarType type,
        double lb,
        double ub,
        std::vector<double> column = {}) :
      MyObject(name),
      index_(index),
      current_epoch_(-1),
      current_index_(-1),
      type_(type),
      cost_(cost),
      lb_(lb),
      ub_(ub),
      compactColumn_(std::move(column)),
      iteration_creation_(-1),
      active_count_(0),
      last_active_(0) {}

  MyVar(const MyVar &var) :
      MyObject(var),
      index_(var.index_),
      current_epoch_(var.current_epoch_),
      current_index_(var.current_index_),
      type_(var.type_),
      cost_(var.cost_),
      lb_(var.lb_),
      ub_(var.ub_),
      compactColumn_(var.compactColumn_),
      iteration_creation_(var.iteration_creation_),
      active_count_(var.active_count_),
      last_active_(var.last_active_) {}

  ~MyVar() override {}

  int getIndex() const { return index_; }

  void setCurrentIndex(int epoch, int c_index) {
    current_epoch_ = epoch;
    current_index_ = c_index;
  }

  int getCurrentIndex() const { return current_index_; }
  int getCurrentIndex(int epoch) const {
    // epoch == 0 is only used for the core variables
    if (current_epoch_ == 0 || current_epoch_ == epoch)
      return current_index_;
    return -1;
  }

  virtual VarType getVarType() const { return type_; }

  double getCost() const { return cost_; }

  virtual double getLB() const { return lb_; }

  virtual double getUB() const { return ub_; }

  void setCost(double cost) { cost_ = cost; }

  virtual void setLB(double lb) { lb_ = lb; }

  virtual void setUB(double ub) { ub_ = ub; }

  virtual void setVarType(VarType newType) { type_ = newType; }

  bool is_integer() const { return type_ != VARTYPE_CONTINUOUS; }

  const std::vector<double> &getCompactColumn() const { return compactColumn_; }

  int getIterationCreation() const { return iteration_creation_; }

  int getActiveCount() const { return active_count_; }

  int getLastActive() const { return last_active_; }

  int getNbConsInactiveIteration(int currentLPIteration) const {
    return currentLPIteration - last_active_;
  }
  double getActivityRate(int currentLPIteration) const {
    return active_count_ * 1.0 / (currentLPIteration - iteration_creation_);
  }

  void addActiveIteration(int iteration);

  void initActiveCounters(int iteration_creation);

  // get the first day of the rotation corresponding to the column
  virtual int getFirstDay() const {
    return static_cast<int>(compactColumn_[1]);
  }

  // get the id of the nurse in charge of the rotation
  virtual int getNurseNum() const {
    return static_cast<int>(compactColumn_[0]);
  }

 protected:
  int index_;  // count var
  int current_epoch_;  // epoch when current_index_ has been set
  int current_index_;  // position of the variables in the current epoch
  VarType type_;  // type of the variable
  double cost_;  // cost of the variable
  double lb_;  // lower bound
  double ub_;  // upper bound
  const std::vector<double> compactColumn_;  // column for a column
  // save the iteration number at the creation of the variable
  int iteration_creation_;
  // count the number of times where the variable is present in the solution
  int active_count_;
  // save the iteration number of the last activity
  int last_active_;

  virtual void setIndex(int index) {
    index_ = index;
  }

  friend class Modeler;
};

struct MyCons : public MyObject {
  MyCons(const char *name, int index, double lhs, double rhs) :
      MyObject(name),
      index_(index),
      lhs_(lhs),
      rhs_(rhs) {}

  MyCons(const MyCons &cons) :
      MyObject(cons.name_),
      index_(cons.index_),
      lhs_(cons.lhs_),
      rhs_(cons.rhs_) {}

  ~MyCons() override {}

  /*
   * Getters
   */

  double getLhs() const { return lhs_; }

  double getRhs() const { return rhs_; }

  void setLhs(double lhs) { lhs_ = lhs; }

  void setRhs(double rhs) { rhs_ = rhs; }

  int getIndex() const { return index_; }

 protected:
  int index_;
  double lhs_;  // left hand side == lower bound
  double rhs_;  // rihgt hand side == upper bound
};

/*
 * My pricer
 */
struct MyPricer {
  explicit MyPricer(const char *name) : name_(name) {}
  virtual ~MyPricer() {}

  // name of the pricer handler
  //
  const char *name_;

  /* perform pricing */
  // return true if optimal
  virtual std::vector<MyVar *> pricing(double bound = 0,
                                       bool before_fathom = false,
                                       bool after_fathom = false,
                                       bool backtracked = false) = 0;

  virtual void updateParameters(bool useMoreTime) = 0;

  virtual double getLastMinReducedCost() const = 0;
  virtual const std::vector<double> &getLastReducedCostLBs() const = 0;

  // return true if all subproblems have been solved to optimality
  virtual bool isLastRunOptimal() const = 0;

  // return true if all subproblems have a lower bound
  virtual bool isLastRunLowerBounded() const = 0;

  // METHODS - Forbidden shifts, nurses, starting days, etc.
  //
  // !!! WARNING !!! : SOME METHODS ARE NOT YET IMPLEMENTED IN THE SUBPROBLEM
  // (ALTHOUGH THE NECESSARY STRUCTURES MAY ALREADY BE THERE !!!
  //
  // Shifts
  virtual void forbidShift(int k, int s) = 0;
//   virtual void forbidShifts(set<pair<int,int> > shifts)=0;
//   virtual void authorizeShift(int k, int s)=0;
//   virtual void clearForbiddenShifts()=0;
  // Nurses
  virtual void forbidNurse(int nurseNum) = 0;
//   virtual void forbidNurses(set<int,int> nurses)=0;
  virtual void authorizeNurse(int nurseNum) = 0;
  virtual void clearForbiddenNurses() = 0;

  virtual void initNursesAvailabilities() = 0;
};

/* Represents the nodes of the branching tree */
struct MyNode;
typedef std::shared_ptr<MyNode> MyPNode;

struct MyNode {
  MyNode() : index_(0),
             pParent_(nullptr),
             bestLB_(-XLARGE_SCORE),
             processed_(false),
             gap_(XLARGE_SCORE),
             smallestGap_(XLARGE_SCORE),
             depth_(0),
             bestLagLB_(-XLARGE_SCORE),
             lastLagLB_(-XLARGE_SCORE),
             presolvedUB_(XLARGE_SCORE) {}
  explicit MyNode(const MyPNode &pParent) :
      index_(-1),
      pParent_(pParent.get()),
      bestLB_(pParent_->bestLB_),
      processed_(false),
      gap_(XLARGE_SCORE),
      smallestGap_(XLARGE_SCORE),
      depth_(pParent_->depth_ + 1),
      bestLagLB_(pParent_->bestLB_),
      lastLagLB_(-XLARGE_SCORE),
      presolvedUB_(XLARGE_SCORE) {}
  virtual ~MyNode() = default;

  // parent
  MyNode *pParent_;

  void pushBackChild(const MyPNode &child) {
    children_.push_back(child);
  }

  const std::vector<MyPNode>& getChildren() const {
    return children_;
  }

  bool updateLB(double newLB);

  void computeGap(double treeBestLB) {
    gap_ = (bestLB_ - treeBestLB) / treeBestLB;
  }

  double getLPGap() const {
    return gap_;
  }

  double getHighestGap() const {
    // if root, it is the best
    if (!pParent_)
      return XLARGE_SCORE;

    // otherwise compare the current gap
    return pParent_->smallestGap_;
  }

  // The quality of a node is used to sort the candidate list
  // (siblings of the branching tree)
  // They are sorted in ascending order
  double getQuality() const {
    // if root, does not apply
    if (!pParent_)
      return XLARGE_SCORE;
    // otherwise return parent's best LB
    return pParent_->getBestLB() + 1e-6 * depth_;
  }

  double getBestLB() const { return bestLB_; }

  // STAB
  double getBestLagLB() const { return bestLagLB_; }
  void setBestLagLB(double bestLagLB) { bestLagLB_ = bestLagLB; }

  // LAGLB
  double getLastLagLB() const { return lastLagLB_; }
  void setLastLagLB(double lb) { lastLagLB_ = lb; }

  // PRESOLVED
  double getPresolvedUB() const { return presolvedUB_; }
  void setPresolvedUB(double presolvedUB) { presolvedUB_ = presolvedUB; }

  void setDepth(int depth) { depth_ = depth; }

  int getDepth() const { return depth_; }

  bool isProcessed() const {
    return processed_;
  }

  virtual std::string write() const {
    return "MyNode: ("+getInfo()+")";
  }

  virtual std::string writeInheritance() const {
    std::stringstream out;
    const MyNode *pNode = this;
    while (pNode != nullptr) {
      out <<  pNode->write() << std::endl;
      pNode = pNode->pParent_;
    }
    return out.str();
  }

  int getIndex() const { return index_; }

 protected:
  int index_;
  double bestLB_;
  bool processed_;
  // gap: between bestLB_ and current best lb
  // highest gap: between the bestLB_ and the computed bestLB_ of the children
  double gap_, smallestGap_;
  int depth_;

  // STAB : best lagrangian bound obtained at this node
  double bestLagLB_;

  // LAGLB: last Lagrangian lower bound in the column generation
  double lastLagLB_;

  // PRESOLVED : UB obtained when performing strong banching
  double presolvedUB_;

  // children
  std::vector<MyPNode> children_;

  friend class MyTree;
  void setIndex(int index) { index_ = index; }

  std::string getInfo() const {
    std::stringstream out;
    out << "depth=" << depth_ << ", LB=" << bestLB_;
    if (presolvedUB_ < XLARGE_SCORE - 1)
      out << ", presolved=" << presolvedUB_;
    return out.str();
  }
};

struct MyTree {
  explicit MyTree(double epsilon, bool printCurrentNode = true)
      : epsilon_(epsilon),
        treeSize_(0),
        nbNodesProcessed_(0),
        nbNodesLastIncumbent_(0),
        diveDepth_(0),
        diveLength_(LARGE_SCORE),
        minDepth_(0),
        nbNodesSinceDive_(0),
        currentNode_(nullptr),
        printCurrentNode_(printCurrentNode),
        bestLbInRoot(-XLARGE_SCORE),
        bestLb(-XLARGE_SCORE),
        maxBestLb(XLARGE_SCORE),
        bestUb(XLARGE_SCORE),
        bestLbMinTreeLevel_(0) {}

  virtual ~MyTree() {}

  void setRootLB(double bestLBRoot) { bestLbInRoot = bestLBRoot; }

  double getRootLB() const { return bestLbInRoot; }

  void setCurrentNode(const MyPNode &currentNode);

  /* clear tree */
  void clear() {
    tree_.clear();
    activeTreeMapping_.clear();
  }

  std::string writeCurrentNode() {
    if (currentNode_)
      return currentNode_->write();
    else
      return "";
  }

  void keepFirstChild() {
    // set current node to first child
    setCurrentNode(currentNode_->getChildren()[0]);
  }

  void eraseCurrentSiblings();

  const MyPNode &getNode(const int nodeIndex) const {
    return tree_[nodeIndex];
  }

  double getBestLB() const { return bestLb; }

  int getBestLBMinTreeLevel() const { return bestLbMinTreeLevel_; }

  void setBestUB(double ub) {
    /* reinitialize nbNodesLastIncumbent_ */
    if (ub + 1 < bestUb) nbNodesLastIncumbent_ = 0;
    if (ub < bestUb) bestUb = ub;
  }

  double getBestUB() const { return bestUb; }

  double getCurrentLB() const { return currentNode_->getBestLB(); }

  // fix maximum value of the best LB. In general, it's due to a branch that has
  // been cut, but should have not. It can only happen when not solving
  // subproblems to optimality.
  void fixMaxBestLb(double lb) {
    if (lb < maxBestLb) maxBestLb = lb;
  }

  void swapNodes(int index1, int index2) {
    std::swap(tree_[index1], tree_[index2]);
  }

  void swapLastNodes() {
    int index1 = tree_.size() - 2,
        index2 = tree_.size() - 1;
    swapNodes(index1, index2);
  }

  // Reset and clear solving parameters
  virtual void reset() {
    clear();
    bestUb = XLARGE_SCORE;
    currentNode_ = nullptr;
    treeSize_ = 0;
    nbNodesProcessed_ = 0;
    nbNodesLastIncumbent_ = 0;
    nbNodesSinceDive_ = 0;
    diveDepth_ = 0;
    diveLength_ = LARGE_SCORE;
    bestLbInRoot = -XLARGE_SCORE;
    bestLb = -XLARGE_SCORE;
    maxBestLb = XLARGE_SCORE;
  }

  int getTreeSize() const { return treeSize_; }

  int getNbNodesProcessed() const { return nbNodesProcessed_; }

  int getNbNodesSinceLastIncumbent() const { return nbNodesLastIncumbent_; }

  void resetNbNodesSinceLastDive() { nbNodesSinceDive_ = 0; }

  int getDiveLength() const { return diveLength_; }

  int getNbDives() const {
    return nbNodesSinceDive_ / diveLength_;
  }

  double updateNodeLB(double lb);

  // compute best lb by visiting all unprocessed nodes
  // return the best lb and its level in tree
  // (in case of tie, choose the highest one in the tree)
  std::pair<double, int> computeCurrentBestLB();

  // STAB: getter and setter for lagrangian bound
  double getNodeBestLagLB() const { return currentNode_->getBestLagLB(); }

  double getNodeLastLagLB() const { return currentNode_->getLastLagLB(); }

  bool updateNodeLagLB(double lb) {
    currentNode_->setLastLagLB(lb);
    if (lb > currentNode_->getBestLagLB() + epsilon_) {
      currentNode_->setBestLagLB(lb);
      return true;
    }
    return false;
  }

  double getObjective() const { return bestUb; }

  double getRelaxedObjective() const { return bestLbInRoot; }

  void createRootNode() {
    tree_.push_back(std::make_shared<MyNode>());
    ++treeSize_;
  }

  void pushBackNode(const MyPNode &node) {
    node->setIndex(tree_.size());
    tree_.push_back(node);
    treeSize_++;
    currentNode_->pushBackChild(node);
  }

  const MyPNode &getCurrentNode() const { return currentNode_; }

  double getIntegralityGap() const {
    return (getBestUB() - getBestLB()) / getBestLB();
  }

  virtual bool isColumnsNode() const { return false; }

  virtual void addForbiddenShifts(
      PLiveNurse pNurse,
      std::set<std::pair<int, int> > *forbidenShifts) {}

  /*
   * Stats
   */

  virtual void updateStats(const MyPNode &node) {}

  virtual std::string writeBranchStats() const { return ""; }

  void printStats() const { std::cout << writeBranchStats(); }

 protected:
  double epsilon_;
  // mapping between the Siblings and MyNode*
  // a siblings contains a list of all its leaves MyNode
  std::map<MyNode*, std::vector<MyPNode>> activeTreeMapping_;
  // branching tree
  std::vector<MyPNode> tree_;
  // tree size, number of nodes since last incumbent,
  // depth of the current dive, length of a dive
  int treeSize_, nbNodesProcessed_, nbNodesLastIncumbent_, diveDepth_,
      diveLength_, minDepth_, nbNodesSinceDive_, bestLbMinTreeLevel_;
  // current node
  MyPNode currentNode_;
  // print current node every time it has been changed
  bool printCurrentNode_;

  // best lb in root and current best lb
  double bestLbInRoot, bestLb, maxBestLb;
  // best current upper bound found
  double bestUb;
};

/* Structures to store a branching candidate and its potential children */
struct MyBranchingNode {
  friend struct MyBranchingCandidate;
  MyBranchingNode() = default;

  const std::vector<double> &getLb() const { return LB; }

  void setLb(int index, double lb) { LB[index] = lb; }

  const std::vector<double> &getLhs() const { return lhs; }

  void setLhs(int index, double lhs) { this->lhs[index] = lhs; }

  const std::vector<double> &getRhs() const { return rhs; }

  void setRhs(int index, double rhs) { this->rhs[index] = rhs; }

  const std::vector<double> &getUb() const { return UB; }

  void setUb(int index, double ub) {  UB[index] = ub; }

  void addPNode(MyPNode pNode) { pNode_ = std::move(pNode); }

  const MyPNode &pNode() const  {  return pNode_; }

  virtual void removeVariables(const std::vector<int> &indexToRemove);

 protected:
  MyPNode pNode_ = nullptr;
  std::vector<double> LB, UB;
  std::vector<double> rhs, lhs;
};
typedef std::shared_ptr<MyBranchingNode> MyPBranchingNode;

struct MyBranchingCandidate {
  MyBranchingCandidate() = default;
  ~MyBranchingCandidate() {
    for (MyVar *pV : newBranchingVars_) delete pV;
    for (MyCons *pC : newBranchingCons_) delete pC;
  }

  virtual int createNewChild() {
    children_.push_back(std::make_shared<MyBranchingNode>());
    initialize(children_.back());
    return children_.size() - 1;
  }

  virtual void createNewChild(int position) {
    children_.insert(children_.begin() + position,
                     std::make_shared<MyBranchingNode>());
    initialize(children_[position]);
  }

  const MyPBranchingNode &getChild(int index) {
    return children_[index];
  }

  void swapChildren(int index1, int index2) {
    std::swap(children_[index1], children_[index2]);
  }

  void swapLastChildren() {
    int index1 = children_.size() - 2,
        index2 = children_.size() - 1;
    swapChildren(index1, index2);
  }

  void initialize(const MyPBranchingNode &pNode) {
    for (MyVar *var : branchingVars_) {
      pNode->LB.push_back(var->getLB());
      pNode->UB.push_back(var->getUB());
    }

    for (MyCons *cons : branchingCons_) {
      pNode->lhs.push_back(cons->getLhs());
      pNode->rhs.push_back(cons->getRhs());
    }
  }

  int addBranchingVar(MyVar *var) {
    branchingVars_.push_back(var);
    for (const MyPBranchingNode &pChild : children_) {
      pChild->LB.push_back(var->getLB());
      pChild->UB.push_back(var->getUB());
    }
    return branchingVars_.size() - 1;
  }

  int addNewBranchingVar(MyVar *var) {
    newBranchingVars_.push_back(var);
    return addBranchingVar(var);
  }

  int addBranchingCons(MyCons *cons) {
    branchingCons_.push_back(cons);
    for (const MyPBranchingNode &pChild : children_) {
      pChild->lhs.push_back(cons->getLhs());
      pChild->rhs.push_back(cons->getRhs());
    }
    return branchingCons_.size() - 1;
  }

  int addNewBranchingCons(MyCons *cons) {
    newBranchingCons_.push_back(cons);
    return addBranchingCons(cons);
  }

  const std::vector<MyCons *> &getBranchingCons() const {
    return branchingCons_;
  }

  const std::vector<MyVar *> &getBranchingVars() const {
    return branchingVars_;
  }

  const std::vector<MyPBranchingNode> &getChildren() const {
    return children_;
  }

  const std::vector<MyCons *> &getNewBranchingCons() const {
    return newBranchingCons_;
  }

  const std::vector<MyVar *> &getNewBranchingVars() const {
    return newBranchingVars_;
  }

  virtual void removeDeactivatedColumns(
      Modeler* pModel, const std::function<bool(MyVar*)> &shouldRemoveCol);

 protected:
  std::vector<MyVar *> branchingVars_, newBranchingVars_;
  std::vector<MyCons *> branchingCons_, newBranchingCons_;
  std::vector<MyPBranchingNode> children_;
};
typedef std::shared_ptr<MyBranchingCandidate> MyPBranchingCandidate;

/*
 * My branching rule
 */
struct MyBranchingRule {
  explicit MyBranchingRule(const char *name)
      : name_(name), searchStrategy_(BestFirstSearch) {}
  virtual ~MyBranchingRule() {}

  // name of the branching rule handler
  const char *name_;

  /* Compute logical fixing decisions */
  // Add the new children to the candidate (just column node).
  // Return true a child if created, false otherwise
  virtual bool column_node(std::vector<MyPBranchingCandidate> *candidates) = 0;

  /* Compute branching decisions */
  // Add the new children to the candidate (other nodes).
  // Return true if a child is created
  virtual bool branching_candidates(
      int nCandidates,
      std::vector<MyPBranchingCandidate> *candidates) = 0;

  void set_search_strategy(SearchStrategy searchStrategy) {
    searchStrategy_ = searchStrategy;
  }

 protected:
  SearchStrategy searchStrategy_;
};


//-----------------------------------------------------------------------------
//
// C l a s s   M o d e l e r
//
// Generic class of models
//
//
//-----------------------------------------------------------------------------


class Modeler {
 public:
  Modeler():
      pPricer_(nullptr),
      pBranchingRule_(nullptr),
      pTree_(nullptr),
      stab_(this) {}

  virtual ~Modeler() {
    for (MyVar *var : initialColumnVars_)
      delete var;
    initialColumnVars_.clear();
    for (MyObject *object : objects_)
      delete object;
    objects_.clear();
  }

  // solve the model
  virtual int solve(bool relaxation = false) = 0;

  // Reset and clear solving parameters and objects
  virtual void reset();

  // Add a pricer
  virtual int addPricer(MyPricer *pPricer) {
    delete pPricer_;
    pPricer_ = pPricer;
    return 1;
  }

  MyPricer *pPricer() const {
    return pPricer_;
  }

  // Add a branching rule
  virtual int addBranchingRule(MyBranchingRule *pBranchingRule) {
    delete pBranchingRule_;
    pBranchingRule_ = pBranchingRule;
    pBranchingRule_->set_search_strategy(searchStrategy_);
    return 1;
  }

  MyBranchingRule *pBranchingRule() const {
    return pBranchingRule_;
  }

  // Add a tree
  virtual int addTree(MyTree *pTree) {
    delete pTree_;
    pTree_ = pTree;
    return 1;
  }

  MyTree *pTree() const {
    return pTree_;
  }

  virtual void addForbiddenShifts(
      PLiveNurse pNurse,
      std::set<std::pair<int, int> > *forbidenShifts) {
    pTree_->addForbiddenShifts(pNurse, forbidenShifts);
  }


  /*
   * Class methods for pricer and branching rule
   */

  // return true if optimal
  std::vector<MyVar *> pricing(double bound = 0,
                               bool before_fathom = false,
                               bool after_fathom = false,
                               bool backtracked = false) {
    if (pPricer_)
      return pPricer_->pricing(bound, before_fathom, after_fathom, backtracked);
    return {};
  }

  bool branching_candidates(
      int nCandidates, std::vector<MyPBranchingCandidate> *candidates) {
    if (pBranchingRule_)
      return pBranchingRule_->branching_candidates(nCandidates, candidates);
    return false;
  }

  // remove all bad candidates from fixingCandidates
  bool branch_on_column(std::vector<MyPBranchingCandidate> *candidates) {
    if (pBranchingRule_)
      return pBranchingRule_->column_node(candidates);
    return false;
  }

  // Set search strategy
  void set_search_strategy(SearchStrategy searchStrategy) {
    if (pBranchingRule_)
      pBranchingRule_->set_search_strategy(searchStrategy);
  }

  double epsilon() const {
    return parameters_.epsilon_;
  }

  /*
   * Store and set num for objects objects
   */
 protected:
  std::recursive_mutex mutex_;

  void registerObject(MyObject *o) {
    std::lock_guard<std::recursive_mutex> lock(mutex_);
    if (o->num() < 0)
      o->num(objectsCount++);
  }

  void addObject(MyObject *o) {
    std::lock_guard<std::recursive_mutex> lock(mutex_);
    registerObject(o);
    objects_.push_back(o);
  }

  /*
   * Create variable:
   *    var is a pointer to the pointer of the variable
   *    var_name is the name of the variable
   *    lhs, rhs are the lower and upper bound of the variable
   *    vartype is the type of the variable:
   *    VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
   */

 protected:
  virtual int createVar(MyVar **var,
                        const char *var_name,
                        int index,
                        double objCoeff,
                        double lb,
                        double ub,
                        VarType vartype,
                        const std::vector<double> &column,
                        double score) = 0;

 public:
  void createVar(MyVar **var,
                 const char *var_name,
                 double objCoeff,
                 double lb,
                 double ub,
                 VarType vartype,
                 const std::vector<double> &column,
                 double score);

  void createPositiveVar(MyVar **var,
                         const char *var_name,
                         double objCoeff,
                         const std::vector<double> &column = {},
                         double score = 0,
                         double ub = DBL_MAX);

  void createIntVar(MyVar **var,
                    const char *var_name,
                    double objCoeff,
                    const std::vector<double> &column = {},
                    double score = 0,
                    double ub = DBL_MAX);

  void createBinaryVar(MyVar **var,
                       const char *var_name,
                       double objCoeff,
                       const std::vector<double> &column = {},
                       double score = 0);

  void createPositiveFeasibilityVar(
      MyVar **var,
      const char *var_name,
      const std::vector<double> &column = {},
      double score = 0,
      double ub = DBL_MAX);

 protected:
  virtual int createColumnVar(MyVar **var,
                              const char *var_name,
                              int index,
                              double objCoeff,
                              const std::vector<double> &column,
                              double dualObj,
                              double lb,
                              double ub,
                              VarType vartype,
                              double score) = 0;

 public:
  void createColumnVar(MyVar **var,
                       const char *var_name,
                       double objCoeff,
                       const std::vector<double> &column,
                       double dualObj,
                       double lb,
                       double ub,
                       VarType vartype,
                       double score);

  void createPositiveColumnVar(MyVar **var,
                               const char *var_name,
                               double objCoeff,
                               const std::vector<double> &column,
                               double dualObj = LARGE_SCORE,
                               double score = 0,
                               double ub = DBL_MAX);

  void createIntColumnVar(MyVar **var,
                          const char *var_name,
                          double objCoeff,
                          const std::vector<double> &column,
                          double dualObj = LARGE_SCORE,
                          double score = 0,
                          double ub = DBL_MAX);

  void createBinaryColumnVar(MyVar **var,
                             const char *var_name,
                             double objCoeff,
                             const std::vector<double> &column,
                             double dualObj = LARGE_SCORE,
                             double score = 0);

  /*
   * Create linear constraint present at the beginning) and cut (generated):
   *    con is a pointer to the pointer of the constraint
   *    con_name is the name of the constraint
   *    lhs, rhs are the lower and upper bound of the constraint
   *    nonZeroVars is the number of non-zero coefficients to add to the constraint
   *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
   *    coeffs is the array of coefficient to add to the constraints
   */

  // Add linear constraints
 protected:
  virtual int createConsLinear(MyCons **cons,
                               const char *con_name,
                               int index,
                               double lhs,
                               double rhs,
                               std::vector<MyVar *> vars = {},
                               std::vector<double> coeffs = {}) = 0;

 public:
  void createConsLinear(MyCons **cons,
                        const char *con_name,
                        double lhs,
                        double rhs,
                        std::vector<MyVar *> vars = {},
                        std::vector<double> coeffs = {});

  void createLEConsLinear(MyCons **cons,
                          const char *con_name,
                          double rhs,
                          std::vector<MyVar *> vars = {},
                          std::vector<double> coeffs = {});

  void createGEConsLinear(MyCons **cons,
                          const char *con_name,
                          double lhs,
                          std::vector<MyVar *> vars = {},
                          std::vector<double> coeffs = {});

  void createEQConsLinear(MyCons **cons,
                          const char *con_name,
                          double eq,
                          std::vector<MyVar *> vars = {},
                          std::vector<double> coeffs = {});

 protected:
  virtual int createCutLinear(MyCons **cons,
                              const char *con_name,
                              int index,
                              double lhs,
                              double rhs,
                              std::vector<MyVar *> vars = {},
                              std::vector<double> coeffs = {}) = 0;

 public:
  // set manually lhs and rhs
  void createCutLinear(MyCons **cons,
                       const char *con_name,
                       double lhs,
                       double rhs,
                       std::vector<MyVar *> vars = {},
                       std::vector<double> coeffs = {});

  // Add a lower or equal constraint
  void createLECutLinear(MyCons **cons,
                         const char *con_name,
                         double rhs,
                         std::vector<MyVar *> vars = {},
                         std::vector<double> coeffs = {});

  // Add a greater or equal constraint
  void createGECutLinear(MyCons **cons,
                         const char *con_name,
                         double lhs,
                         std::vector<MyVar *> vars = {},
                         std::vector<double> coeffs = {});

  // Add an equality constraint
  void createEQCutLinear(MyCons **cons,
                         const char *con_name,
                         double eq,
                         std::vector<MyVar *> vars = {},
                         std::vector<double> coeffs = {});

  /*
   * Add variables to constraints
   */

  virtual int addCoefLinear(MyCons *cons,
                            MyVar *var,
                            double coeff,
                            bool transformed = false) = 0;

  /*
   * Add new Column to the problem
   */

  void createColumn(MyVar **var,
                    const char *var_name,
                    double objCoeff,
                    const std::vector<double> &column,
                    double dualObj,
                    VarType vartype,
                    std::vector<MyCons *> cons = {},
                    std::vector<double> coeffs = {},
                    bool transformed = false,
                    double score = 0);

  void createPositiveColumn(MyVar **var,
                            const char *var_name,
                            double objCoeff,
                            const std::vector<double> &column,
                            double dualObj,
                            std::vector<MyCons *> cons = {},
                            std::vector<double> coeffs = {},
                            bool transformed = false,
                            double score = 0);

  void createBinaryColumn(MyVar **var,
                          const char *var_name,
                          double objCoeff,
                          const std::vector<double> &column,
                          double dualObj,
                          std::vector<MyCons *> cons = {},
                          std::vector<double> coeffs = {},
                          bool transformed = false,
                          double score = 0);

  void createIntColumn(MyVar **var,
                       const char *var_name,
                       double objCoeff,
                       const std::vector<double> &column,
                       double dualObj,
                       std::vector<MyCons *> cons = {},
                       std::vector<double> coeffs = {},
                       bool transformed = false,
                       double score = 0);

  /*
   * get the primal values
   */

  virtual bool isInteger(MyVar *var) const {
    return isInteger(getVarValue(var));
  }

  virtual bool isInteger(double v) const {
    return std::fabs(round(v) - v) < epsilon();
  }

  virtual double getVarValue(MyVar *var) const = 0;

  // compute the total cost of a multiple vectors of MyObject* in the solution
  double getVarValue(const std::vector<MyVar *> &vector) const {
    double value = 0;
    for (MyVar *v : vector)
      if (v) value += getVarValue(v);
    return value;
  }

  template<typename V>
  double getVarValue(const std::vector<V> &vector) const {
    double value = 0;
    for (const V &vect : vector)
      value += getVarValue(vect);
    return value;
  }

  std::vector<double> getVarValues(const std::vector<MyVar *> &vars) const {
    std::vector<double> values(vars.size(), 0);
    for (unsigned int i = 0; i < vars.size(); ++i)
      if (vars[i]) values[i] = getVarValue(vars[i]);
    return values;
  }

  /*
   * Get the dual variables
   */

  virtual double getDual(MyCons *cons, bool transformed = false) const = 0;

  double getDual(const std::vector<MyCons *> &vector,
                 bool transformed = false) const {
    double value = 0;
    for (MyCons *c : vector)
      if (c) value += getDual(c);
    return value;
  }

  template<typename T, typename V>
  double getDual(const std::pair<T, V> &p,
                 bool transformed = false) const {
    double value =
        getDual(p.first, transformed) + getDual(p.second, transformed);
    return value;
  }

  std::vector<double> getDuals(const std::vector<MyCons *> &cons,
                               bool transformed = false) const {
    std::vector<double> dualValues(cons.size(), 0);
    for (unsigned int i = 0; i < cons.size(); ++i)
      if (cons[i]) dualValues[i] = getDual(cons[i], transformed);
    return dualValues;
  }

  template<typename T, typename V>
  std::vector<double> getDuals(const std::map<T, V> &cons,
                               bool transformed = false) const {
    std::vector<double> dualValues(cons.size(), 0);
    for (unsigned int i = 0; i < cons.size(); ++i)
      if (cons[i]) dualValues[i] = getDual(cons[i], transformed);
    return dualValues;
  }

  /*
   * Get the reduced cost
   */

  virtual double getReducedCost(MyVar *var) const = 0;

  /*
   * Getters and setters
   */

  // compute the total cost of MyObject* in the solution
  virtual double getTotalCost(MyVar *var, bool print = false) const {
    if (var == nullptr) return 0;
    double value = getVarValue(var);
    if (print && value > epsilon())
      std::cout << var->name_ << ": " << value << "*" << var->getCost()
                << std::endl;
    return value * var->getCost();
  }

  // compute the total cost of a vector of MyObject* in the solution
  template<typename T>
  double getTotalCost(const std::map<MyVar *, T> &map0,
                      bool print = false) const {
    double value = 0;
    for (const std::pair<MyVar *, T> &var : map0)
      value += getTotalCost(var.first, print);
    return value;
  }

  template<typename T, typename V>
  double getTotalCost(const std::map<T, V> &m, bool print = false) const {
    double value = 0;
    for (const auto &p : m)
      value += getTotalCost(p.second, print);
    return value;
  }

  // compute the total cost of a multiple vectors of MyObject* in the solution
  template<typename V>
  double getTotalCost(const std::vector<V> &vector, bool print = false) const {
    double value = 0;
    for (const V &vect : vector)
      value += getTotalCost(vect, print);
    return value;
  }

  template<typename T, typename V>
  double getTotalCost(const std::pair<T, V> &p, bool print = false) const {
    double value =
        getTotalCost(p.first, print) + getTotalCost(p.second, print);
    return value;
  }

  /**************
   * Parameters *
   *************/
  virtual int setVerbosity(int v) = 0;

  /**************
   * Outputs *
   *************/

  virtual int printStats() const {
    pTree_->printStats();
    return 1;
  }

  virtual int printBestSol();

  virtual int writeProblem(std::string fileName) const = 0;

  virtual int writeLP(std::string fileName) const = 0;

  virtual int writeMPS(std::string fileName) const = 0;

  virtual void toString(MyObject *obj) const {
    std::cout << obj->name_ << std::endl;
  }

  /**************
   * Getters *
   *************/

  template<typename M>
  M getModel() const {
    std::string error = "This template has not been implemented.";
    Tools::throwError(error.c_str());
  }

  const std::vector<MyVar *> &getBinaryCoreVars() const {
    return binaryCoreVars_;
  }

  const std::vector<MyVar *> &getIntegerCoreVars() const {
    return integerCoreVars_;
  }

  const std::vector<MyVar *> &getPositiveCoreVars() const {
    return positiveCoreVars_;
  }

  const std::vector<MyVar *> &getFeasibilityCoreVars() const {
    return feasibilityCoreVars_;
  }

  bool isFeasible() const {
    for (MyVar *var : feasibilityCoreVars_) {
      double v = getVarValue(var);
      if (v > epsilon() || v < -epsilon())
        return false;
    }
    return true;
  }

  bool isInfeasible() const { return !isFeasible(); }

  void printInfeasibleVars() const {
    for (MyVar *var : feasibilityCoreVars_) {
      double v = getVarValue(var);
      if (v > epsilon() || v < -epsilon()) {
        std::cout << var->name_ << ": "
                  << std::setprecision(3) << v << std::endl;
      }
    }
  }

  int getCurrentIndex(MyVar *pVar) const {
    return pVar->getCurrentIndex(current_epoch_);
  }

  int getVerbosity() const { return verbosity_; }

  virtual void setBestUB(double ub) { pTree_->setBestUB(ub); }

  virtual double getObjective() const { return pTree_->getBestUB(); }

  virtual double getBestUB() const { return pTree_->getBestUB(); }

  virtual double getObjective(int index) const { return XLARGE_SCORE; }

  virtual double getRelaxedObjective() const { return pTree_->getRootLB(); }

  double getCurrentLB() const { return pTree_->getCurrentLB(); }

  virtual double updateNodeLB(double lb) { return pTree_->updateNodeLB(lb); }

  const MyPNode &getNode(const int nodeIndex) const {
    return pTree_->getNode(nodeIndex);
  }

  // STAB
  double getNodeBestLagLB() { return pTree_->getNodeBestLagLB(); }

  double getNodeLastLagLB() { return pTree_->getNodeLastLagLB(); }

  bool updateNodeLagLB(double lb) { return pTree_->updateNodeLagLB(lb); }

  double getRootLB() const { return pTree_->getRootLB(); }

  double getBestLB() const { return pTree_->getBestLB(); }

  int getBestLBMinTreeLevel() const { return pTree_->getBestLBMinTreeLevel(); }

  int getDiveLength() const { return pTree_->getDiveLength(); }

  double getIntegralityGap() const { return pTree_->getIntegralityGap(); }

  int getTreeSize() const { return pTree_->getTreeSize(); }

  int getNbNodesProcessed() const { return pTree_->getNbNodesProcessed(); }

  std::string writeCurrentNode() const { return pTree_->writeCurrentNode(); }

  bool isColumnsNode() const { return pTree_->isColumnsNode(); }

  void keepFirstChild() { return pTree_->keepFirstChild(); }

  int getNbDives() const { return pTree_->getNbDives(); }

  void resetNbNodesSinceLastDive() { pTree_->resetNbNodesSinceLastDive(); }

  virtual int nbSolutions() const { return 0; }

  void setSearchStrategy(SearchStrategy searchStrategy) {
    searchStrategy_ = searchStrategy;
    set_search_strategy(searchStrategy);
  }

  SearchStrategy getSearchStrategy() const { return searchStrategy_; }

  virtual void setParameters(const SolverParam &parameters,
                             PrintSolution *func) {
    parameters_ = parameters;
    parameters_.saveFunction_ = func;
    setVerbosity(parameters_.verbose_);
    logfile_ = parameters.logfile_;
  }

  std::string logfile() const { return logfile_; }

  const SolverParam &getParameters() const { return parameters_; }

  SolverParam* getpParameters() { return &parameters_; }

  void setLogFile(std::string fileName) { logfile_ = fileName; }

  void setInfinity(double inf) { infinity_ = inf; }
  double getInfinity() const { return infinity_; }

  // get the variables that are always present in the model
  const std::vector<MyVar *> &getCoreVars() const { return coreVars_; }

  const std::vector<MyCons *> &getCoreCons() const { return coreCons_; }

  // get the variables that are generating during the resolution (columns)
  // and currently active (i.e. have a positive value in the current solution)
  const std::vector<MyVar *> &getActiveColumns() const {
    return activeColumnVars_;
  }

  std::pair<int, int> getFractionalAndPositiveColumns() const;

  bool hasFractionalColumns() const;

  // get the variables that are generating during the resolution (columns)
  // and currently active (i.e. have a positive value in the current solution)
  const std::vector<MyVar *> &getInitialColumns() const {
    return initialColumnVars_;
  }

  // add the column and mark it owned by the modeler
  void addInitialColumn(MyVar *var) {
    initialColumnVars_.push_back(var);
  }

  void setInitialColumns(const std::vector<MyVar *> &vars) {
    initialColumnVars_ = vars;
  }

  virtual void clear() {
    // reset stabilization
    stab_ = Stabilization(this);
    // clear and delete columns
    clearActiveColumns();
    // clear and delete core vars and cons
    coreVars_.clear();
    binaryCoreVars_.clear();
    integerCoreVars_.clear();
    positiveCoreVars_.clear();
    feasibilityCoreVars_.clear();
    coreCons_.clear();
    for (MyObject *obj : objects_)
      delete obj;
    objects_.clear();
  }

  virtual void clearActiveColumns(size_t new_size = 0) {
    activeColumnVars_ = std::vector<MyVar *>();
    if (new_size > 0) activeColumnVars_.reserve(new_size);
    current_epoch_++;  // go to the next epoch
  }

  virtual void clearInitialColumns() {
    initialColumnVars_.clear();
  }

  virtual void resetInitialColumns() {
    // reset the columns index in the Modeler after a reset
    for (MyVar *v : initialColumnVars_)
      v->setIndex(var_count++);
  }

  Stabilization &stab() {
    return stab_;
  }

  virtual void copyActiveToInitialColumns() {}

  // Clear the active column, set the active columns with those
  // in the best solution, and set the primal values accordingly
  virtual bool loadBestSol(bool integer) { return false; }
  virtual bool loadBestSol() { return loadBestSol(true); }

  virtual bool isSolutionInteger() const { return false; }

  // Get the current level in the branch and bound tree
  virtual int getCurrentTreeLevel() const { return 0; }

  virtual void addNode(
      const MyPBranchingNode &pNode, double presolvedUB = DBL_MAX) {
    if (presolvedUB < DBL_MAX - 1) pNode->pNode()->setPresolvedUB(presolvedUB);
    pTree_->pushBackNode(pNode->pNode());
  }

  // get the current number of objects
  int getNbObjectsInMemory() const { return objects_.size(); }

 protected:
  virtual void addActiveColumn(MyVar *var, int c_index) {
    activeColumnVars_.push_back(var);
    var->setCurrentIndex(current_epoch_, c_index);
  }

  // store the variables and the constraints that belongs to the core of the
  // model (the one that are not generated).
  std::vector<MyVar *> coreVars_;
  std::vector<MyVar *> binaryCoreVars_;
  std::vector<MyVar *> integerCoreVars_;
  std::vector<MyVar *> positiveCoreVars_;
  // vector of variables that are used to soften some hard constraints to
  // ensure feasibility
  std::vector<MyVar *> feasibilityCoreVars_;
  std::vector<MyCons *> coreCons_;
  size_t var_count = 0, cons_count = 0;
  // as soon as columns are added, core variables cannot be added anymore
  // as soon as cut are added, core constraints cannot be added anymore
  // The reason why we implement such a mechanism is that core variables and
  // constraints need to be indexed in a sequence from 0 to N.
  // If columns or cuts are generated in the middle, it breaks the sequence.
  // Even more, it goes against the logic that core variables and constraints
  // are instantiated at the beginning of the problem.
  bool column_added = false, cut_added = false;

  // store the active columns (the one that has a positive value)
  std::vector<MyVar *> activeColumnVars_;

  // When starting a new solve, the columns present in this vector will be
  // added back to the model
  // the vector will also be cleared
  std::vector<MyVar *> initialColumnVars_;

  MyPricer *pPricer_;  // prices the rotations
  MyTree *pTree_;  // store the tree information
  // choose the variables on which we should branch
  MyBranchingRule *pBranchingRule_;
  Stabilization stab_;

  int verbosity_ = 0;
  SearchStrategy searchStrategy_ = BestFirstSearch;

  SolverParam parameters_;

  // log file where outputs must be written
  std::string logfile_ = "";

  // Coin data. It is very important to set this value big enough (>1e27)
  // as otherwise the model will interpret it as real bound instead of infinity
  double infinity_ = INFINITY;

 protected:
  // store all MyObject* (objects owned by modeler)
  std::vector<MyObject *> objects_;
  // count the objects created by the modeler
  // (could be different of objects_.size())
  int objectsCount = 0;
  // store the current epoch for the active solution stored
  // current index of the variables are only valid for this epoch
  int current_epoch_ = 0;
};

#endif  // SRC_SOLVERS_MP_MODELER_MODELER_H_
