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

#include "MyRCSPP.h"

#include <algorithm>
#include <list>
#include <memory>
#include <vector>


PRCLabel LabelPool::getNewLabel() {
  PRCLabel pL;
  // if at the end of the vector
  if (endLabel_ == pLabels_.end()) {
    // add it to the vector
    pLabels_.emplace_back(pFactory_->makePRCLabel());
    // set pL
    pL = pLabels_.back();
    // update end: to be done at the end as vector can be moved
    endLabel_ = pLabels_.end();
  } else {
//    // recreate labels all the time ...
//    pL = pFactory_->makePRCLabel();
//    *endLabel_  = pL;
    // reuse an already initialized label
    pL = *endLabel_;
    // update the num
    pFactory_->updateNumlabel(pL);
    // update the endLabel
    endLabel_++;
  }
  return pL;
}

void LabelPool::sort() {
  // Sorting the vector of labels by increasing cost
  std::sort(begin(), end(), LabelCostIncreasing());
}

MyRCSPPSolver::MyRCSPPSolver(MyRCGraph *pRcGraph) :
    pRcGraph_(pRcGraph),
    pFactory_(std::make_shared<RCLabelFactory>()),
    labelPool_(pFactory_, 100),
    nParetoLabels_(0),
    number_of_infeasible_deleted_labels_(0),
    total_number_of_generated_labels_(0),
    total_number_of_dominations_(0),
    verbose_(false),
    maxReducedCostBound_(0),
    epsilon_(1.0e-5),
    strategy_(SP_BREADTH_FIRST),
    nb_max_paths_(1),
    sortLabelsOption_(false),
    minimumCostToSinksOption_(false),
    worstCaseCostOption_(false) {}

void MyRCSPPSolver::reset() {
  verbose_ = false;
  maxReducedCostBound_ = 0;
  epsilon_ = 1.0e-5;
  strategy_ = SP_BREADTH_FIRST;
  nb_max_paths_ = 1;
  pLabelsPerNode_.clear();
  pLabelsPerNode_.resize(pRcGraph_->nNodes());
}

std::vector<RCSolution> MyRCSPPSolver::solve(double maxReducedCostBound,
                                             int verbose,
                                             double epsilon,
                                             SPSearchStrategy strategy,
                                             int nb_max_paths,
                                             PScenario pS,
                                             vector<int> maxConsShifts,
                                             vector<int> minConsShifts,
                                             bool enumeratedSubPathOption,
                                             bool sortLabelsOption,
                                             bool minimumCostToSinksOption,
                                             bool worstCaseCostOption) {
  maxReducedCostBound_ = maxReducedCostBound;
  verbose_ = verbose;
  epsilon_ = epsilon;
  strategy_ = strategy;
  nb_max_paths_ = nb_max_paths;
  sortLabelsOption_ = sortLabelsOption;
  minimumCostToSinksOption_ = minimumCostToSinksOption;
  worstCaseCostOption_ = worstCaseCostOption;
  enumeratedSubPathOption_ = enumeratedSubPathOption;

  // Displaying of the activated options for the algorithm
  if (verbose_ >= 3) {
    std::cout << "\nOPTION(S): " << std::endl;
    if (sortLabelsOption_)
      std::cout << " - sortLabelsOption is activated " << std::endl;
    if (minimumCostToSinksOption_)
      std::cout << " - minimumCostToSinksOption is activated " << std::endl;
    if (worstCaseCostOption_)
      std::cout << " - worstCaseCostOption is activated " << std::endl;
    if (enumeratedSubPathOption_)
      std::cout << " - enumeratedSubPathOption is activated " << std::endl;

    std::cout << "\nGRAPH : " << std::endl;
    std::cout << " - number of nodes: " << pRcGraph_->nNodes() << std::endl;
    std::cout << " - number of arcs: " << pRcGraph_->nArcs() << std::endl;
  }

  timer_.start();

  // sort nodes as the graph should be acyclic
  // Allow to visit each node once
  std::vector<PRCNode> sortedNodes = pRcGraph_->sortNodes();
  for (const auto &pN : sortedNodes) {
    // expand labels from all predecessors
    for (const auto &pArc : pN->inArcs) {
      if (minimumCostToSinksOption_) {
        // If we computed the minimum costs to the sinks, we expand only the
        // labels from which a negative path to the sinks can be built
        for (const auto &pL : pLabelsPerNode_.at(pArc->origin->id)) {
          // check if possible to obtain a negative cost for the current label
          // when enable
          if (hasPotentialNegativePathToSinks(pL)) {
            // Expansion on the current arc
            if (expand(pL, pArc, labelPool_.getNewLabel()))
              total_number_of_generated_labels_++;
            else
              number_of_infeasible_deleted_labels_++;
          } else {
            number_of_infeasible_deleted_labels_++;
          }
        }
      } else {
        // Expand all the labels
        for (const auto &pL : pLabelsPerNode_.at(pArc->origin->id)) {
          // Expansion on the current arc
          if (expand(pL, pArc, labelPool_.getNewLabel())) {
            // has been able to expand
            total_number_of_generated_labels_++;
          } else {
            // has not been able to expand: release label
            labelPool_.releaseLastLabel();
            number_of_infeasible_deleted_labels_++;
          }
        }
      }
    }

    // Sorting by increasing cost of labels just being expanded
    if (sortLabelsOption_)
      labelPool_.sort();

    // Domination of all the labels of the current node to get a
    // pareto-optimal set of labels
    checkAllDominations(labelPool_.begin(), labelPool_.end());

    // clear the pool (allow to reuse labels for next node)
    labelPool_.clear();
  }

  if (verbose_ >= 3) {
    std::cout << "\nLABELS STATISTICS: " << std::endl;
    // We add 1 to nParetoLabels in order to take account the first labels of
    // the source node
    std::cout << " - Total number of not dominated labels: "
                 "" << nParetoLabels_ + 1 << std::endl;
    std::cout << " - Number of infeasible labels deleted: "
                 "" << number_of_infeasible_deleted_labels_ << std::endl;
    std::cout << " - Total number of generated labels: " <<
              total_number_of_generated_labels_ << std::endl;
    std::cout << " - Total number of dominations: " <<
              total_number_of_dominations_ << std::endl << std::endl;
  }

  // Collecting of all the labels of the sink nodes
  vector<PRCLabel> pLabelsSinks;
  for (const auto &pN : pRcGraph_->pSinks()) {
    for (const auto &pl : pLabelsPerNode_.at(pN->id))
      pLabelsSinks.push_back(pl);
  }

  // Sorting by increasing cost all the labels of the sink nodes
  std::sort(pLabelsSinks.begin(),
            pLabelsSinks.end(),
            LabelCostIncreasing());

  // Creation of the vector containing the solution(s) of the problem
  vector<RCSolution> finalSolutions;
  // TODO(AL): I do not understand
//  int max_path = nb_max_paths_ > 0 ? nb_max_paths_ : pLabelsSinks.size();
  int max_path = nb_max_paths_ > 0 && nb_max_paths_ < pLabelsSinks.size() ?
      nb_max_paths_ : pLabelsSinks.size();
  for (int n = 0; n < max_path; ++n) {
    PRCLabel currentLabel = pLabelsSinks[n];
    RCSolution solution = createSolution(currentLabel);
    finalSolutions.push_back(solution);
  }

  timer_.stop();



  return finalSolutions;
}

RCSolution MyRCSPPSolver::createSolution(const PRCLabel& finalLabel) {
  // Backtracking
  PRCLabel currentLabel = finalLabel;
  vector<int> shifts;
  while (currentLabel->getPreviousLabel() != nullptr) {
    for (const auto &pS : currentLabel->getInArc()->stretch.pShifts())
      shifts.insert(shifts.begin(), pS->type);
    currentLabel = currentLabel->getPreviousLabel();
  }
  return RCSolution(0, shifts, finalLabel->cost());
}



PRCNode MyRCSPPSolver::getTargetNode(const PRCLabel &pl) const {
  return pl->getInArc()->target;
}



void MyRCSPPSolver::checkAllDominations(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end) {
  // if no labels, returns
  if (begin == end)
    return;

  // compute active resources
  auto pN = (*begin)->getNode();
  vector<PResource> activeResources;
  activeResources.reserve(pRcGraph_->nResources());
  for (const auto &pR : pRcGraph_->pResources())
    if (pR->isActive(pN->day, *pN->pShift))
      activeResources.push_back(pR);

  // compare labels pairwise for domination
  size_t nLabels = std::distance(begin, end);
  vector<bool> isDominated(nLabels, false);
  int nNonDominatedLabels = nLabels;
  auto itL1 = begin;
  for (auto itD1 = isDominated.begin();
       itD1 != isDominated.end();
       ++itD1, ++itL1) {
    if (*itD1) continue;

    auto itL2 = itL1 + 1;
    for (auto itD2 = itD1 + 1; itD2 != isDominated.end(); ++itD2, ++itL2) {
      if (*itD2) continue;

      if (dominate(*itL1, *itL2, activeResources)) {
        *itD2 = true;
        nNonDominatedLabels--;
      }
    }
    itL2 = itL1 + 1;
    if (sortLabelsOption_) {
      for (auto itD2 = itD1 + 1; itD2 != isDominated.end(); ++itD2, ++itL2) {
        if (*itD2) continue;

        if (dominate(*itL2, *itL1, activeResources)) {
          *itD1 = true;
          nNonDominatedLabels--;
          break;
        }
        // when the labels are sorted by increasing cost, domination of L2
        // over L1 can stopped being checked as soon as we find L2 with
        // larger cost than L1
        if ((*itL2)->cost() > (*itL1)->cost())
          break;
      }
    } else {
      for (auto itD2 = itD1 + 1; itD2 != isDominated.end(); ++itD2, ++itL2) {
        if (*itD2) continue;

        if (dominate(*itL2, *itL1, activeResources)) {
          *itD1 = true;
          nNonDominatedLabels--;
          break;
        }
      }
    }
  }

  // set the list of non-dominated labels for current node
  if (nNonDominatedLabels > 0) {
    std::vector<PRCLabel> &nodeLabels = pLabelsPerNode_[pN->id];
    nodeLabels.reserve(nNonDominatedLabels);
    itL1 = begin;
    for (bool d : isDominated) {
      if (!d) {
        // create a new pointer from the current label
        nodeLabels.push_back(pFactory_->makePRCLabel(**itL1));
        nParetoLabels_++;
      }
      ++itL1;
    }
  }
}

bool MyRCSPPSolver::dominate(const PRCLabel &pL1,
                             const PRCLabel &pL2,
                             const vector<PResource>& resources) {
  total_number_of_dominations_++;

  if (worstCaseCostOption_) {
    // Improved domination
    double cost1 = pL1->cost(), cost2 = pL2->cost();
    for (const auto &r : resources) {
      // if at anytime cost1 > cost2, return false
      if (cost1 > cost2) return false;
      double ubDiff =
          pL1->getWorstUbCost(r->id()) - pL2->getWorstUbCost(r->id());
      double lbDiff =
          pL1->getWorstLbCost(r->id()) - pL2->getWorstLbCost(r->id());
      double maxDiff = std::max(.0, std::max(ubDiff, lbDiff));
      cost1 += maxDiff;
    }
    return cost1 <= cost2;
  } else {
    // Basic domination
    if (pL1->cost() > pL2->cost()) return false;
    for (const auto &r : resources) {
      if (!r->dominates(pL1->getConsumption(r->id()),
                        pL2->getConsumption(r->id())))
        return false;
    }
    return true;
  }
}


bool MyRCSPPSolver::expand(
    const PRCLabel &pLParent, const PRCArc &pArc, const PRCLabel &pLChild) {
  pLChild->copy(*pLParent);
  pLChild->setParentLabel(pLParent);
  pLChild->setInArc(pArc);
  // Expansion of each resource
  for (auto &e : pArc->expanders) {
    // get child ResourceValues (copy of the parent ResourceValues)
    ResourceValues *v = &pLChild->getResourceValues(e->resourceId);
    if (!e->expand(pLParent->getResourceValues(e->resourceId), pLChild, v))
      return false;  // has not been able to expand
  }
  pLChild->addBaseCost(pArc->baseCost);
  pLChild->addDualCost(pArc->dualCost);
  return true;  // has been able to expand
}




void MyRCSPPSolver::displayLabel(const PRCLabel &pL) const {
  std::cout << pL->getNum() << " ";
  std::cout << "Cost : " << pL->cost() << " ";
  std::cout << "Resources : ( ";
  for (int idRes = 0; idRes < pL->nResLabels(); ++idRes) {
    std::cout << pL->getConsumption(idRes) << " ";
    std::cout << "[" << pL->getWorstLbCost(idRes) << ","
              << pL->getWorstUbCost(idRes)
              << "] ";
  }
  std::cout << ")" << std::endl;
}




void MyRCSPPSolver::backTracking(const PRCLabel &pL) {
  PRCLabel currentLabel = pL;
  std::cout << "Total shift cost : "
            << currentLabel->totalShiftCost() << std::endl;
  std::cout << "Total week-end cost : "
            << currentLabel->totalWeekendCost() << std::endl;
  std::cout << "Total cons shift cost : "
            << currentLabel->consShiftCost() << std::endl;

  while (currentLabel->getPreviousLabel() != nullptr) {
    std::cout << "================" << std::endl;
    std::cout << Tools::intToDay(getTargetNode(currentLabel)->day)
              << " Node " << getTargetNode(currentLabel)->id << std::endl;
    getTargetNode(currentLabel)->print();
    std::cout << std::endl;
    currentLabel->getInArc()->print();
    std::cout << std::endl;
    displayLabel(currentLabel);
    currentLabel = currentLabel->getPreviousLabel();
  }

  std::cout << "================" << std::endl;
  pRcGraph_->pSource()->print();
  std::cout << std::endl;
  displayLabel(currentLabel);
}



vector<double> MyRCSPPSolver::shortestPathAcyclic(const PRCNode &sink) {
  // Initialization of all the costs of the shortest paths at infinity
  double inf = std::numeric_limits<double>::infinity();
  vector<double> shortestPathToSink(pRcGraph_->nNodes(), inf);

  // Only the cost of the shortest path to the sink node is set to 0
  shortestPathToSink.at(sink->id) = 0.0;
  // Iteration through all the nodes of the acyclic graph in the reverse
  // order of the topological order
  std::list<PRCNode> nodesToProcess = {sink};
  int i = 0;
  while (!nodesToProcess.empty()) {
    PRCNode pN = nodesToProcess.front();
    nodesToProcess.pop_front();
    double sp = shortestPathToSink[pN->id];
    for (const auto &pArc : pN->inArcs) {
      int predecessor = pArc->origin->id;
      // update shortest path if needed
      if (shortestPathToSink[predecessor] > sp + pArc->cost) {
        shortestPathToSink[predecessor] = sp + pArc->cost;
        nodesToProcess.push_back(pArc->origin);
      }
    }
    if (i++ > pRcGraph_->nArcs())
      Tools::throwError("The RC graph is cyclic.");
  }
  return shortestPathToSink;
}



void MyRCSPPSolver::computeMinimumCostToSinks() {
  minimumCostToSinks_.clear();
  minimumCostToSinks_.resize(pRcGraph_->nNodes(), LARGE_SCORE);
  for (const auto &pS : pRcGraph_->pSinks()) {
    vector<double> sp = shortestPathAcyclic(pS);
    for (int i = 0; i < sp.size(); ++i)
      if (sp[i] < minimumCostToSinks_[i])
        minimumCostToSinks_[i] = sp[i];
  }
}



bool MyRCSPPSolver::hasPotentialNegativePathToSinks(
    const PRCLabel &pl) const {
  if (pl->getNode()->type == SOURCE_NODE)
    return true;
  PRCNode pN  = getTargetNode(pl);
  if (pN->type == SINK_NODE)
    return true;

  // If a path to a sink node with a negative cost exists, this label can be
  // expanded
  return pl->cost() + minimumCostToSinks_[pN->id] < 0;
}
