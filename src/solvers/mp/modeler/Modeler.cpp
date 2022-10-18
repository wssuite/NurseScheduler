/*
 * Copyright (C) 2022 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "solvers/mp/modeler/Modeler.h"

#include <cassert>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <list>
#include <memory>
#include <utility>

#include "tools/Tools.h"


void MyVar::addActiveIteration(int iteration) {
  if (last_active_ != iteration) {
    last_active_ = iteration;
    ++active_count_;
  }
}

void MyVar::initActiveCounters(int iteration_creation) {
  iteration_creation_ = iteration_creation;
  active_count_ = 0;
  last_active_ = iteration_creation;
}


bool MyNode::updateLB(double newLB) {
  // if already processed, send a warning
  if (processed_)
    std::cerr << "The node LB has already been updated." << std::endl;
  // if not root and the LB is decreasing -> error.
  // BUG: need to be found. For the moment, return false
  bool increased = true;
  if (pParent_ && pParent_->bestLB_ > newLB + 1e-2) {
    std::cerr << "The node lower bound (" << newLB
              << ") is smaller than its parent one ("
              << pParent_->bestLB_ << ")." << std::endl;
    increased = false;
    newLB = pParent_->bestLB_;
  }
  bestLB_ = newLB;
  processed_ = true;
  // if not root
  if (pParent_) {
    double parentGap = (bestLB_ - pParent_->bestLB_) / pParent_->bestLB_;
    if (parentGap < pParent_->smallestGap_)
      pParent_->smallestGap_ = parentGap;
  }
  return increased;
}

void MyTree::setCurrentNode(const MyPNode &currentNode) {
  // if dive length has not been updated and we are not diving
  // this would be called once at the end of the first dive
  bool diving = currentNode->pParent_ &&
      currentNode->pParent_->children_.size() <= 1;
  if (!diving && diveDepth_ > 0 && diveLength_ == LARGE_SCORE)
    diveLength_ = 1 + diveDepth_;

  // update these parameters only if current node exists
  // (i.e., currentNode is not the root node)
  if (currentNode_) {
    currentNode_->processed_ = true;
    ++nb_nodes_processed_;
    // one more node without new incumbent
    ++nb_nodes_last_incumbent_;
    ++nb_nodes_since_dive_;
    // set dive depth in case of diving
    diveDepth_ = currentNode_->getDepth();
    // add the children if any
    if (!currentNode_->getChildren().empty())
      activeTreeMapping_[currentNode_.get()] = currentNode_->getChildren();
  }
  // update current node
  currentNode_ = currentNode;
  --tree_size_;
  // remove the parent from the mapping if not the root node and
  // if the last one to be processed
  if (currentNode_->pParent_ &&
      activeTreeMapping_.at(currentNode_->pParent_).back() == currentNode_)
    eraseCurrentSiblings();
//  #ifdef DBG
  if (currentNode_ && printCurrentNode_)
    std::cout << currentNode_->write() << std::endl;
//  #endif
}

void MyTree::eraseCurrentSiblings() {
  // should remain only one node not process (the current one)
  for (const MyPNode &pNode : activeTreeMapping_.at(currentNode_->pParent_))
    if (pNode != currentNode_ && !pNode->isProcessed())
      Tools::throwException("Erase a sibling that has not yet "
                            "been processed: "+pNode->write());
  size_t n = activeTreeMapping_.erase(currentNode_->pParent_);
  if (n != 1)
    Tools::throwException("The active mapping wasn't containing "
                          "the parent node (" + currentNode_->pParent_->write()
                          + ") of: " + currentNode_->write());
  // update min_depth_
  min_depth_ = LARGE_SCORE;
  for (const auto &p : activeTreeMapping_)
    if (p.first->getDepth() < min_depth_) min_depth_ = p.first->getDepth();
}

double MyTree::updateNodeLB(double lb) {
  // set lb and gap
  currentNode_->updateLB(lb);
  // if root -> set root and best lb
  double parentLB;
  if (!currentNode_->pParent_) {
    best_lb_in_root = lb;
    best_lb = lb;
    best_lb_min_tree_level_ = currentNode_->getDepth();
    parentLB = 0;
  } else {
    // update best_lb
    auto bestLBAndMinLevel = computeCurrentBestLB();
    if (bestLBAndMinLevel.first + 1e-5 < best_lb)
      std::cerr << "Best LB is decreasing from " << std::setprecision(10)
                << best_lb << " to " << bestLBAndMinLevel.first << std::endl;
    best_lb = bestLBAndMinLevel.first;
    best_lb_min_tree_level_ = bestLBAndMinLevel.second;
    parentLB = currentNode_->pParent_->getBestLB();
  }

  // compute gap
  currentNode_->computeGap(best_lb);
  updateStats(currentNode_);

  return lb - parentLB;
}

// compute best lb by visiting all unprocessed nodes
// return the best lb and its level in tree
// (in case of tie, choose the highest one in the tree)
std::pair<double, int> MyTree::computeCurrentBestLB() {
  double new_best_lb = currentNode_->getBestLB();
  int min_level = currentNode_->getDepth();

  auto processLB = [&new_best_lb, &min_level](MyNode* node) {
    if (node && !node->isProcessed()) {
      if (new_best_lb >= node->getBestLB() + 1e-6) {
        new_best_lb = node->getBestLB();
        min_level = node->getDepth();
      } else if (new_best_lb >= node->getBestLB() - 1e-6 &&
                 min_level > node->getDepth()) {
        min_level = node->getDepth();
      }
    }
  };

  for (const auto &p : activeTreeMapping_) {
    processLB(p.first);
    for (const MyPNode &node : p.second)
      processLB(node.get());
  }

  return {new_best_lb, min_level};
}

void MyBranchingCandidate::removeDeactivatedColumns(
    Modeler* pModel, const std::function<bool(MyVar*)> &doRemoveCol) {
  int index = -1;
  vector<int> indexToRemove;
  for (MyVar* pVar : branchingVars_) {
    index++;
    // if in use, do nothing
    if (pModel->getVarValue(pVar) > pModel->epsilon()) continue;
    // check if any of the LB > 1
    bool active = false;
    for (const MyPBranchingNode& pNode : children_)
      if (pNode->getLb()[index] > pModel->epsilon()) {
        active = true;
        break;
      }
    if (active) continue;
    // check if one of the UB = 0
    active = true;
    for (const MyPBranchingNode& pNode : children_)
      if (pNode->getUb()[index] < pModel->epsilon()) {
        active = false;
        break;
      }
    if (active) continue;
    // if all LB = 0 and one UB = 0 -> mark to be removed if wished
    if (doRemoveCol(pVar))
      indexToRemove.push_back(index);
  }

  // update variables and children to remove these indices
  vector<MyVar*> newBranchingVars;
  auto it = indexToRemove.begin();
  for (int i = 0; i < branchingVars_.size(); i++) {
    if (it == indexToRemove.end() || i != *it) {
      newBranchingVars.push_back(branchingVars_[i]);
    } else if (i == *it) {
      ++it;
    }
  }
  branchingVars_ = newBranchingVars;

  for (const MyPBranchingNode& pNode : children_)
    pNode->removeVariables(indexToRemove);
}

void MyBranchingNode::removeVariables(const std::vector<int> &indexToRemove) {
  std::vector<double> newLB, newUB;
  auto it = indexToRemove.begin();
  for (int i = 0; i < LB.size(); i++) {
    if (it == indexToRemove.end() || i != *it) {
      newLB.push_back(LB[i]);
      newUB.push_back(UB[i]);
    } else if (i == *it) {
      ++it;
    }
  }
  // set new bounds
  LB = newLB;
  UB = newUB;
}

/*
  * Create variable:
  *    var is a pointer to the pointer of the variable
  *    var_name is the name of the variable
  *    lhs, rhs are the lower and upper bound of the variable
  *    vartype is the type of the variable:
  *    VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
  */
void Modeler::createVar(
    MyVar **var, const char *var_name, double objCoeff,
    double lb, double ub, VarType vartype,
    const std::vector<double> &column, double score) {
  // check flag column_added
  if (column_added)
    Tools::throwError("Cannot create variable %s, as some columns have "
                      "already been generated.", var_name);
  assert(lb <= ub);
  createVar(var, var_name, var_count++, objCoeff, lb, ub,
            vartype, column, score);
  coreVars_.push_back(*var);
  addObject(*var);
}

void Modeler::createPositiveVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double score, double ub) {
  ub = (ub == DBL_MAX) ? infinity_ : ub;
  createVar(var, var_name, objCoeff, 0.0, ub,
            VARTYPE_CONTINUOUS, column, score);
  positiveCoreVars_.push_back(*var);
}

void Modeler::createIntVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double score, double ub) {
  ub = (ub == DBL_MAX) ? infinity_ : ub;
  createVar(var, var_name, objCoeff, 0, ub, VARTYPE_INTEGER, column, score);
  integerCoreVars_.push_back(*var);
}

void Modeler::createBinaryVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double score) {
  createVar(var, var_name, objCoeff, 0, 1, VARTYPE_BINARY, column, score);
  binaryCoreVars_.push_back(*var);
}

void Modeler::createPositiveFeasibilityVar(
    MyVar **var, const char *var_name,
    const std::vector<double> &column, double score, double ub) {
  createPositiveVar(var, var_name, LARGE_SCORE, column, score, ub);
  feasibilityCoreVars_.push_back(*var);
}

void Modeler::createColumnVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double dualObj,
    double lb, double ub, VarType vartype, double score) {
  // set flag column_added to true
  column_added = true;
  assert(lb <= ub);
  createColumnVar(var, var_name, var_count++, objCoeff, column, dualObj,
                  lb, ub, vartype, score);
  registerObject(*var);
}

void Modeler::createPositiveColumnVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column,
    double dualObj, double score, double ub) {
  ub = (ub == DBL_MAX) ? infinity_ : ub;
  createColumnVar(var, var_name, objCoeff, column, dualObj,
                  0, ub, VARTYPE_CONTINUOUS, score);
}

void Modeler::createIntColumnVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column,
    double dualObj, double score, double ub) {
  ub = (ub == DBL_MAX) ? infinity_ : ub;
  createColumnVar(var, var_name, objCoeff, column, dualObj,
                  0, ub, VARTYPE_INTEGER, score);
}

void Modeler::createBinaryColumnVar(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column,
    double dualObj, double score) {
  createColumnVar(var, var_name, objCoeff, column, dualObj,
                  0, 1, VARTYPE_BINARY, score);
}

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
void Modeler::createConsLinear(
    MyCons **cons, const char *con_name, double lhs, double rhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  // check flag cut_added
  if (cut_added)
    Tools::throwError("Cannot create constraint %s, as some cuts "
                      "have already been generated.", con_name);
  assert(lhs <= rhs);
  createConsLinear(cons, con_name, cons_count++, lhs, rhs,
                   std::move(vars), std::move(coeffs));
  coreCons_.push_back(*cons);
  addObject(*cons);
}

void Modeler::createLEConsLinear(
    MyCons **cons, const char *con_name, double rhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createConsLinear(cons, con_name, -infinity_, rhs,
                   std::move(vars), std::move(coeffs));
}

void Modeler::createGEConsLinear(
    MyCons **cons, const char *con_name, double lhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createConsLinear(cons, con_name, lhs, infinity_,
                   std::move(vars), std::move(coeffs));
}

void Modeler::createEQConsLinear(
    MyCons **cons, const char *con_name, double eq,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createConsLinear(cons, con_name, eq, eq, std::move(vars), std::move(coeffs));
}

// set manually lhs and rhs
void Modeler::createCutLinear(
    MyCons **cons, const char *con_name, double lhs, double rhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  // set flag cut_added to true
  cut_added = true;
  assert(lhs <= rhs);
  createCutLinear(cons, con_name, cons_count++, lhs, rhs,
                  std::move(vars), std::move(coeffs));
  registerObject(*cons);
}

// Add a lower or equal constraint
void Modeler::createLECutLinear(
    MyCons **cons, const char *con_name, double rhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createCutLinear(cons, con_name, -infinity_, rhs,
                  std::move(vars), std::move(coeffs));
}

// Add a greater or equal constraint
void Modeler::createGECutLinear(
    MyCons **cons, const char *con_name, double lhs,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createCutLinear(cons, con_name, lhs, infinity_,
                  std::move(vars), std::move(coeffs));
}

// Add an equality constraint
void Modeler::createEQCutLinear(
    MyCons **cons, const char *con_name, double eq,
    std::vector<MyVar *> vars, std::vector<double> coeffs) {
  createCutLinear(cons, con_name, eq, eq,
                  std::move(vars), std::move(coeffs));
}

/*
 * Add new Column to the problem
 */
void Modeler::createColumn(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double dualObj, VarType vartype,
    std::vector<MyCons *> cons, std::vector<double> coeffs,
    bool transformed, double score) {
  switch (vartype) {
    case VARTYPE_BINARY:
      createBinaryColumnVar(var, var_name, objCoeff, column, dualObj, score);
      break;
    case VARTYPE_INTEGER:
      createIntColumnVar(var, var_name, objCoeff, column, dualObj, score);
      break;
    default:
      createPositiveColumnVar(var, var_name, objCoeff, column, dualObj, score);
      break;
  }

  for (unsigned int i = 0; i < cons.size(); i++)
    addCoefLinear(cons[i], *var, coeffs[i], transformed);
}

void Modeler::createPositiveColumn(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double dualObj,
    std::vector<MyCons *> cons, std::vector<double> coeffs,
    bool transformed, double score) {
  createColumn(var, var_name, objCoeff, column, dualObj, VARTYPE_CONTINUOUS,
               std::move(cons), std::move(coeffs), transformed, score);
}

void Modeler::createBinaryColumn(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double dualObj,
    std::vector<MyCons *> cons, std::vector<double> coeffs,
    bool transformed, double score) {
  createColumn(var, var_name, objCoeff, column, dualObj, VARTYPE_BINARY,
               std::move(cons), std::move(coeffs), transformed, score);
}

void Modeler::createIntColumn(
    MyVar **var, const char *var_name, double objCoeff,
    const std::vector<double> &column, double dualObj,
    std::vector<MyCons *> cons, std::vector<double> coeffs,
    bool transformed, double score) {
  createColumn(var, var_name, objCoeff, column, dualObj, VARTYPE_INTEGER,
               std::move(cons), std::move(coeffs), transformed, score);
}

// Reset and clear solving parameters and objects
void Modeler::reset() {
  // reset tree
  pTree_->reset();
  // clear active columns
  clearActiveColumns();
  // reset counters and flags
  var_count = coreVars_.size();
  column_added = false;
  cons_count = coreCons_.size();
  cut_added = false;
  // check it's an ordered sequence
  if (coreVars_.back()->getIndex() != (var_count - 1))
    Tools::throwError("coreVars_ is not an ordered sequence: the last "
                      "variable does not have the right index.");
  if (coreCons_.back()->getIndex() != (cons_count - 1))
    Tools::throwError("coreCons_ is not an ordered sequence: the last "
                      "constraint does not have the right index.");
  // reset initial columns to ensure the right index
  resetInitialColumns();
}

int Modeler::printBestSol() {
  Tools::LogOutput log(logfile_, true);

  // print the value of the relaxation
  log.printnl("Relaxation: %4.2f", getRelaxedObjective());

  if (!loadBestSol())
    return 0;

  // print the objective value
  log.printnl("Objective: %4.2f; %d solutions found.",
              Modeler::getObjective(), nbSolutions());

  if (verbosity_ >= 2) {
    // display all objective solution found
    log.printnl("%-30s", "All Objective:");
    for (int n = 0; n < nbSolutions(); ++n)
      log.printnl("Solution number: %-5d: %4.2f", n, getObjective(n));

    // print the value of the positive variables
    log.printnl("%-30s", "Variables:");
    double tolerance = pow(.1, DECIMALS);
    // iterate on core variables
    for (MyVar *var : coreVars_) {
      double value = getVarValue(var);
      if (value > tolerance)
        log.printnl("%-30s %4.2f (%6.0f)",
                    var->name_,
                    value,
                    var->getCost());
    }
    // iterate on column variables
    for (MyVar *var : activeColumnVars_) {
      double value = getVarValue(var);
      if (value > tolerance)
        log.printnl("%-30s %4.2f (%6.0f)",
                    var->name_,
                    value,
                    var->getCost());
    }
    log.printnl("");
  }

  return 1;
}

std::pair<int, int> Modeler::getFractionalAndPositiveColumns() const {
  int frac = 0;
  int non_zero = 0;
  for (MyVar *var : getActiveColumns()) {
    double value = getVarValue(var);
    if (value < epsilon())
      continue;
    non_zero++;
    if (value < 1 - epsilon())
      frac++;
  }
  return {frac, non_zero};
}

bool Modeler::hasFractionalColumns() const {
  for (MyVar *var : getActiveColumns()) {
    double value = getVarValue(var);
    if (value >= epsilon() && value < 1 - epsilon())
      return true;
  }
  return false;
}
