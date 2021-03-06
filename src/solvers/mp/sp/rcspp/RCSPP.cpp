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

#include "solvers/mp/sp/rcspp/RCSPP.h"

#include <algorithm>
#include <memory>
#include <list>
#include <vector>

#include "solvers/mp/RCPricer.h"

// TODO(JO): works only for  soft resources. AL commented
int dominate(const PRCLabel &pL1,
             const PRCLabel &pL2,
             const vector<PResource> &resources) {
  // Improved domination
  double cost1 = pL1->cost(), cost2 = pL2->cost();
  for (const auto &r : resources) {
    // if at anytime cost1 > cost2, return false
    if (cost1 > cost2) return 0;
    // if pL1 doesn't dominate pL2 on resource r, stop
    if (!r->dominates(pL1, pL2, &cost1)) return 0;
  }
  return cost1 <= cost2;
}

// Variants of the label setting functions when using the decrease state-space
// relaxation where lower bounds are first ignored for the constraints on the
// total number of shifts and the number of shifts per rotation
// returned statuses:
// - 0 if no domination when LBs are ignored
// - 1 if domination with every constraint
// - 2 if domination when LBs are ignored but not with every constraint
int dominateNoLB(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 const vector<PResource> &resources)  {
  // We implement only improved domination for DSSR
  double cost1noLB = pL1->cost(), costOfLBs = 0;
  double cost2 = pL2->cost();
  for (const auto &r : resources) {
    // if at anytime cost1 > cost2, return false
    if (cost1noLB > cost2) return 0;
    // if hard resource, just check both bounds
    if (r->isHard()) {
      if (!r->dominates(pL1, pL2)) return 0;
      continue;
    }
    // if soft resource, evaluate worst case for UB
    double ubDiff =
        pL1->getWorstUbCost(r->id()) - pL2->getWorstUbCost(r->id());
    double lbDiff =
        pL1->getWorstLbCost(r->id()) - pL2->getWorstLbCost(r->id());
#ifdef DBG
    if (((lbDiff > .0) && (ubDiff > .0)) || ((lbDiff < .0) && (ubDiff < .0))) {
      std::cout << "ubDiff = " << ubDiff << "; lbDiff = " << lbDiff <<
                std::endl;
      Tools::throwError("ubDiff and lbDiff should never have the same sign");
    }
#endif
    double maxDiff;
    if (r->isAnyWorkShiftResource()) {
      // ignore the lower bounds for any resource on any working shift
      if (lbDiff > .0) {
        costOfLBs += lbDiff;
        maxDiff = .0;
      } else {
        maxDiff = std::max(.0, ubDiff);
      }
    } else {
      // for other resources, do as usual
      maxDiff = std::max(ubDiff, lbDiff);
    }
    cost1noLB += maxDiff;
  }
  return cost1noLB > cost2 ? 0 : ((cost1noLB + costOfLBs <= cost2) ? 1 : 2);
}

vector<double> MyRCSPPSolver::shortestPathToSinksAcyclic(
    const RCGraph *pRCGraph) {
  // Initialization of all the costs of the shortest paths at infinity
  double inf = std::numeric_limits<double>::infinity();
  vector<double> shortestPathToSinks(pRCGraph->nNodes(), inf);
  // Only the cost of the shortest path to the sink node is set to 0
  for (const auto &pN : pRCGraph->pSinks())
    shortestPathToSinks.at(pN->id) = 0.0;

  // Iteration through all the nodes of the acyclic graph in the reverse
  // order of the topological order
  vector<PRCNode> sortedNodes = pRCGraph->sortNodes();
  auto itN = sortedNodes.rbegin();
  for (; itN != sortedNodes.rend(); itN++) {
    double sp = shortestPathToSinks[(*itN)->id];
    for (const auto &pArc : (*itN)->inArcs) {
      if (pArc->forbidden) continue;  // do not use forbidden arcs
      int predecessor = pArc->origin->id;
      // update shortest path if needed
      if (shortestPathToSinks[predecessor] > sp + pArc->cost) {
        shortestPathToSinks[predecessor] = sp + pArc->cost;
      }
    }
  }
  return shortestPathToSinks;
}

void MyRCSPPSolver::computeMinimumCostToSinks(RCGraph *pRCGraph) {
  minimumCostToSinks_.clear();
  minimumCostToSinks_ = shortestPathToSinksAcyclic(pRCGraph);
}

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

MyRCSPPSolver::MyRCSPPSolver(RCGraph *pRcGraph,
                             const SubproblemParam &param) :
    pRcGraph_(pRcGraph),
    pFactory_(std::make_shared<RCLabelFactory>()),
    labelPool_(pFactory_, 100),
    nParetoLabels_(0),
    bestPrimalBound_(0.0),
    number_of_infeasible_deleted_labels_(0),
    total_number_of_generated_labels_(0),
    total_number_of_dominations_(0),
    param_(param),
    verbose_(false),
    maxReducedCostBound_(0),
    epsilon_(1.0e-5),
    strategy_(SP_BREADTH_FIRST),
    nb_max_paths_(1) {}

void MyRCSPPSolver::reset() {
  verbose_ = false;
  maxReducedCostBound_ = 0;
  epsilon_ = 1.0e-5;
  strategy_ = SP_BREADTH_FIRST;
  nb_max_paths_ = 1;
  this->pExpandedLabelsPerNode_.clear();
  this->pLabelsToExpandPerNode_.clear();
  this->pLabelsNoExpandPerNode_.clear();
}

void MyRCSPPSolver::displaySolveInputInfo() {
  if (verbose_ >= 3) {
    std::cout << "\nOPTION(S): " << std::endl;
    if (param_.rcsppSortLabels_)
      std::cout << " - sortLabels option is activated " << std::endl;
    if (param_.rcsppMinCostToSinks_)
      std::cout << " - minimumCostToSinks option is activated " << std::endl;
    if (param_.rcsppImprovedDomination_)
      std::cout << " - worstCaseCost option is activated " << std::endl;
    if (param_.rcsppEnumSubpaths_)
      std::cout << " - enumeratedSubPath option is activated " << std::endl;
    if (param_.rcsppDssr_)
      std::cout << " - decrementalSearchSpace option is activated " <<
                std::endl;

    std::cout << "\nGRAPH : " << std::endl;
    std::cout << " - number of nodes: " << pRcGraph_->nNodes() << std::endl;
    std::cout << " - number of arcs: " << pRcGraph_->nArcs() << std::endl;
  }
}

void MyRCSPPSolver::displaySolveStatistics() {
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
}

std::vector<RCSolution> MyRCSPPSolver::solve(double maxReducedCostBound,
                                             int verbose,
                                             double epsilon,
                                             SPSearchStrategy strategy,
                                             int nb_max_paths) {
  maxReducedCostBound_ = maxReducedCostBound;
  verbose_ = verbose;
  epsilon_ = epsilon;
  strategy_ = strategy;
  nb_max_paths_ = nb_max_paths;

  // Displaying of the activated options for the algorithm
  displaySolveInputInfo();

  bestPrimalBound_ = DBL_MAX;
  vector<PRCLabel> pLabelsSinks;
  // sort nodes as the graph should be acyclic
  vector<PRCNode> sortedNodes = pRcGraph_->sortNodes();

  while (true) {
    vector<PRCLabel> destinationLabels = forwardLabelSetting(sortedNodes, 0);
    pLabelsSinks.insert(pLabelsSinks.end(),
                        destinationLabels.begin(),
                        destinationLabels.end());
    // Sorting by increasing cost all the labels of the sink nodes
    int nbNegativeLabels = getNbNegativeCostLabels(pLabelsSinks,
                                                   &bestPrimalBound_);
    if (verbose_ >= 3) {
      std::cout << "Current best primal bound = " << bestPrimalBound_ <<
                std::endl;
      std::cout << "Number of negative cost labels = " << nbNegativeLabels <<
                std::endl;
    }

    if (param_.rcsppToOptimality_) {
      // if no heuristic option is activated, we can break, because
      // optimality must be reached
      // otherwise, we prepare for a new iteration where we aim at optimality
      if (!param_.rcsppDssr_ && !param_.rcsppNbToExpand_)
        break;
      else
        prepareForNextExecution(pRcGraph_->pNodes());
    } else {
      // if any heuristic option is activated and we got enough negative cost
      // rosters, we can break, because we got what we were looking for
      // otherwise, we deactivate the most restrictive heuristic and get
      // prepared for a new solution
      if (nbNegativeLabels <= param_.rcsppMinNegativeLabels_)
        break;
      else
        prepareForNextExecution(pRcGraph_->pNodes());
    }
  }

  // Creation of the vector containing the negative cost rosters
  // first, sort by increasing cost all the labels of the sink nodes
  std::sort(pLabelsSinks.begin(),
            pLabelsSinks.end(),
            LabelCostIncreasing());

  vector<RCSolution> finalSolutions;
  for (const auto& pL : pLabelsSinks) {
    if (pL->cost() + epsilon_ < maxReducedCostBound_) {
      RCSolution solution = createSolution(pL);
      finalSolutions.push_back(solution);
    } else {
      break;
    }
  }
  displaySolveStatistics();

  return finalSolutions;
}

std::vector<PRCLabel> MyRCSPPSolver::forwardLabelSetting(const std::vector<
    PRCNode> &sortedNodes,
                                                         int finalDay) {
  vector<PRCLabel> destinationLabels;
  // Allow to visit each node once
  for (const auto &pN : sortedNodes) {
    // Expand all the labels of the predecessors of the node through the
    // arcs entering the node
    pullLabelsFromPredecessors(pN);

    // information of domination of newly generated labels
    size_t nLabels = labelPool_.nLabels();
    vector<int> domStatusNew(nLabels, 0);

    // information of domination of former labels that have not been
    // expanded yet
    std::vector<PRCLabel> &labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
    size_t nLabelsNoExpand = labelsNoExpand.size();
    vector<int> domStatusNewNoExpand(nLabelsNoExpand, 0);

    // if expand all the labels at the node, dominate all labels first
    if (param_.rcsppNbToExpand_ == 0) {
      // Domination of all the labels of the current node to get a
      // pareto-optimal set of labels
      checkAllDominations(labelPool_.begin(),
                          labelPool_.end(),
                          !param_.rcsppDssr_ ? dominate : dominateNoLB,
                          &domStatusNew,
                          &domStatusNewNoExpand);
    }

    // For heuristic search, select only a subset of labels for expansion
    selectLabelsToExpand(labelPool_.begin(),
                         labelPool_.end(),
                         domStatusNew,
                         domStatusNewNoExpand, pN);

    // clear the pool (allow to reuse labels for next node)
    labelPool_.clear();
  }

  // Collect all the labels of the sink nodes
  for (const auto &pN : pRcGraph_->pSinks()) {
    for (const auto &pl : pLabelsToExpandPerNode_.at(pN->id))
      destinationLabels.push_back(pl);
  }
  return destinationLabels;
}

void MyRCSPPSolver::pullLabelsFromPredecessors(const PRCNode& pN) {
  // Define default expansion function
  std::function<bool(const PRCLabel &, const PRCArc &)> func  =
      [&](const PRCLabel &pL, const PRCArc &pArc) {
        return  expand(pL, pArc, labelPool_.getNewLabel());
      };
  // If we computed the minimum costs to the sinks, we expand only the
  // labels from which a negative path to the sinks can be built
  if (param_.rcsppMinCostToSinks_) {
    func = [&](const PRCLabel &pL, const PRCArc &pArc) {
      // check if possible to obtain a negative cost for the current label
      // when enable, then, expansion on the current arc
      return hasPotentialImprovingPathToSinks(
          pL, pArc->origin, bestPrimalBound_) &&
          expand(pL, pArc, labelPool_.getNewLabel());
    };
  }

  // expand labels from all predecessors
  for (const auto &pArc : pN->inArcs) {
    if (pArc->forbidden) continue;  // do not use forbidden arcs
    for (const auto &pL : pLabelsToExpandPerNode_.at(pArc->origin->id)) {
      // Expand the labels
      if (func(pL, pArc)) total_number_of_generated_labels_++;
      else
        number_of_infeasible_deleted_labels_++;
    }
  }
}

bool MyRCSPPSolver::expand(
    const PRCLabel &pLParent, const PRCArc &pArc, const PRCLabel &pLChild) {
  pLChild->setAsNext(pLParent, pArc);

  // Expansion of each resource
  for (auto &e : pArc->expanders) {
    // get child ResourceValues (copy of the parent ResourceValues)
    ResourceValues &v = pLChild->getResourceValues(e->resourceId);
    if (!e->expand(pLChild, &v)) {
      labelPool_.releaseLastLabel();
      return false;
    }
  }
  pLChild->addCost(pArc->baseCost + pArc->dualCost);
#ifdef DBG
  pLChild->addBaseCost(pArc->baseCost);
  pLChild->addDualCost(pArc->dualCost);
#endif
  return true;
}


void MyRCSPPSolver::checkDominationsWithPreviousLabels(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    const vector<PResource> &resources,
    int (*domFunction)(const PRCLabel &,
                       const PRCLabel &,
                       const vector<PResource> &),
    vector<int> *domStatus,
    vector<int> *domStatusNoExpand) {


  auto pN = (*begin)->getNode();
  std::vector<PRCLabel> &expandedLabels = pExpandedLabelsPerNode_[pN->id];
  std::vector<PRCLabel> &labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];

  // if no label is present on the node, no comparison is necessary
  if (expandedLabels.empty() && labelsNoExpand.empty()) return;

  // check for domination between the newly expanded labels and those that
  // are already stored at the node
  // it is useless to check if the stored labels that have already been
  // expanded are dominated
  // also, the stored labels have already been compared pairwise for domination
  auto itL1 = begin;
  int status;
  for (auto itD1 = domStatus->begin();
       itD1 != domStatus->end();
       ++itD1, ++itL1) {
    if (*itD1 == 1) continue;

    // check if the expanded labels dominate the newly generated ones
    for (const auto& pL2 : expandedLabels) {
      total_number_of_dominations_++;
      status = domFunction(pL2, *itL1, resources);
      *itD1 = (status == 1) ? 1 : std::max(*itD1, status);
      if (*itD1 == 1) break;
    }
    if (*itD1 == 1) continue;

    // check if the non-expanded labels dominate the newly generated ones
    for (const auto& pL2 : labelsNoExpand) {
      total_number_of_dominations_++;
      status = domFunction(pL2, *itL1, resources);
      *itD1 = (status == 1) ? 1 : std::max(*itD1, status);
      if (*itD1 == 1) break;
    }
    if (*itD1 == 1) continue;

    // check if the non-dominated newly generated labels dominate the stored
    // labels that have not been expanded yet
    auto itL2 = labelsNoExpand.begin();
    for (auto itD2 = domStatusNoExpand->begin();
         itD2 != domStatusNoExpand->end();
         ++itD2, ++itL2) {
      if (*itD2 == 1) continue;
      total_number_of_dominations_++;
      status = domFunction(*itL1, *itL2, resources);
      *itD2 = (status == 1) ? 1 : std::max(*itD2, status);
    }
  }
}

void MyRCSPPSolver::checkDominationsPairwise
    (const vector<PRCLabel>::iterator &begin,
     const vector<PRCLabel>::iterator &end,
     const vector<PResource> &resources,
     int (*domFunction)(const PRCLabel &,
                        const PRCLabel &,
                        const vector<PResource> &),
     vector<int> *domStatus) {
  // compare newly generated labels pairwise for domination
  auto itL1 = begin;
  for (auto itD1 = domStatus->begin();
       itD1 != domStatus->end();
       ++itD1, ++itL1) {
    if (*itD1 == 1) continue;

    auto itL2 = itL1 + 1;
    int status;
    for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
      if (*itD2 == 1) continue;

      total_number_of_dominations_++;
      status = domFunction(*itL1, *itL2, resources);
      *itD2 = (status == 1) ? 1 : std::max(*itD2, status);
    }
    itL2 = itL1 + 1;
    if (param_.rcsppSortLabels_) {
      for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
        if (*itD2 == 1) continue;

        total_number_of_dominations_++;
        status = domFunction(*itL2, *itL1, resources);
        *itD1 = (status == 1) ? 1 : std::max(*itD1, status);
        if (*itD1 == 1) break;

        // when the labels are sorted by increasing cost, domination of L2
        // over L1 can be stopped being checked as soon as we find L2 with
        // larger cost than L1
        if ((*itL2)->cost() > (*itL1)->cost()) break;
      }
    } else {
      for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
        if (*itD2 == 1) continue;

        total_number_of_dominations_++;
        status = domFunction(*itL2, *itL1, resources);
        *itD1 = (status == 1) ? 1 : std::max(*itD1, status);
        if (*itD1 == 1) break;
      }
    }
  }
}


void MyRCSPPSolver::checkAllDominations(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    int (*domFunction)(const PRCLabel &,
                       const PRCLabel &,
                       const vector<PResource> &),
    vector<int> *domStatus,
    vector<int> *domStatusNoExpand) {
  // if no labels, returns
  if (begin == end) return;

  // Sort the labels by increasing order before checking domination
  if (param_.rcsppSortLabels_)
    std::sort(begin, end, LabelCostIncreasing());


  // compute active resources
  auto pN = (*begin)->getNode();
  vector<PResource> activeResources;
  activeResources.reserve(pRcGraph_->nResources());
  for (const auto &pR : pRcGraph_->pResources()) {
    if (param_.rcsppImprovedDomination_) pR->useAltenativeDomination();
    else
      pR->useDefaultDomination();
    if (pR->isActive(pN->day, *pN->pAShift))
      activeResources.push_back(pR);
  }

  // check for domination between the newly generated labels and those that
  // are already stored at the node
  checkDominationsWithPreviousLabels(begin, end,
                                     activeResources,
                                     domFunction,
                                     domStatus, domStatusNoExpand);

  // compare newly generated labels pairwise for domination
  checkDominationsPairwise(begin, end, activeResources, domFunction, domStatus);
}


void MyRCSPPSolver::selectLabelsToExpand(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    const vector<int>& domStatus,
    const vector<int>& domStatusNoExpand,
    const PRCNode &pN) {

  // if no labels, returns
  size_t nLabels = std::distance(begin, end);
  vector<PRCLabel>& labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
  size_t nLabelsNoExpand = labelsNoExpand.size();
  if (nLabels + nLabelsNoExpand == 0) return;

  // get the vector of labels to expand
  vector<PRCLabel>& labelsToExpand = pLabelsToExpandPerNode_[pN->id];

  // set the list of non-dominated labels for current node
  if (param_.rcsppNbToExpand_ >= 1) {
    std::sort(begin, end, SortForFewExpansions());
    int cnt = 0;
    for (auto itL = begin; itL != end; ++itL) {
      if (cnt >= param_.rcsppNbToExpand_) break;
      cnt++;
      // expand the "best" nbExpanded labels
      labelsToExpand.push_back(pFactory_->makePRCLabel(**itL));
    }
  } else if (param_.rcsppDssr_) {
    auto itL = begin;
    // add the newly generated non-dominated labels
    for (int d : domStatus) {
      if (d == 0) {
        // expand the label if not dominated with relaxed constraints
        labelsToExpand.push_back(pFactory_->makePRCLabel(**itL));
        nParetoLabels_++;
      } else if (d == 2) {
        // store the label without expanding it if not dominated with every
        // constraint
        labelsNoExpand.push_back(pFactory_->makePRCLabel(**itL));
      }
      ++itL;
    }
  } else {
    auto itL = labelsNoExpand.begin();
    // add the non-expanded non-dominated labels
    for (int d : domStatusNoExpand) {
      if (d == 0) {
        // expand the label if not dominated with every constraint
        labelsToExpand.push_back(*itL);
        nParetoLabels_++;
      }
      ++itL;
    }
    labelsNoExpand.clear();
    itL = begin;
    // add the newly generated non-dominated labels
    for (int d : domStatus) {
      if (d == 0) {
        // expand the label if not dominated with every constraint
        labelsToExpand.push_back(pFactory_->makePRCLabel(**itL));
        nParetoLabels_++;
      }
      ++itL;
    }
  }
}

void MyRCSPPSolver::prepareForNextExecution(const vector<PRCNode> &nodes) {
  // go over activated heuristic options from the most to the least agressive
  if (param_.rcsppNbToExpand_ >= 1) {
    // if we have expanded only a  small subset of labels at each node, we did
    // not perform any domination, so we just start the algorithm over
    param_.rcsppNbToExpand_ = 0;
    for (const auto &pN : nodes) {
      if (pN->type != SOURCE_NODE) {
        pLabelsToExpandPerNode_[pN->id].clear();
      }
    }
  } else if (param_.rcsppDssr_) {
    // if we used a decremental state-space relaxation, we wish to keep the
    // expansion and domination information we got during the execution and
    // warm-start the new execution of the label-setting algorithm
    param_.rcsppDssr_ = false;
    // move the labels that were expanded to the list of expanded labels
    for (const auto &pN : nodes) {
      std::vector<PRCLabel>
          &expandedLabels = pExpandedLabelsPerNode_[pN->id];
      std::vector<PRCLabel>
          &labelsToExpand = pLabelsToExpandPerNode_[pN->id];
      expandedLabels.reserve(expandedLabels.size() + labelsToExpand.size());
      for (const auto &pL : labelsToExpand) {
        expandedLabels.push_back(pL);
      }
      labelsToExpand.clear();
    }
  }
}

RCSolution MyRCSPPSolver::createSolution(const PRCLabel &finalLabel) {
  // Backtrack from the label
  PRCLabel pL = finalLabel;
  int firstDay;
  vector<int> shifts;
  while (pL->getPreviousLabel() != nullptr) {
    for (const auto &pS : pL->getInArc()->stretch.pShifts())
      shifts.insert(shifts.begin(), pS->type);
    if (pL->getNode()->day >= 0)
      firstDay = pL->getNode()->day;
    pL = pL->getPreviousLabel();
  }
  // "Track forward" from the label if some bakcward expansion
  pL = finalLabel;
  while (pL->getNextLabel() != nullptr) {
    for (const auto &pS : pL->getOutArc()->stretch.pShifts())
      shifts.push_back(pS->type);
    pL = pL->getNextLabel();
  }

  return RCSolution(firstDay, shifts, finalLabel->cost());
}



bool MyRCSPPSolver::hasPotentialImprovingPathToSinks(
    const PRCLabel &pl, const PRCNode& pN, double primalBound) const {
  // If a path to a sink node with a negative cost exists, this label can be
  // expanded
  return pl->cost() + minimumCostToSinks_[pN->id] < primalBound;
}

int MyRCSPPSolver::getNbNegativeCostLabels(const vector<PRCLabel> &labels,
                                           double *pMinCost) {
  *pMinCost = 0;
  int nbNegativeCosts = 0;
  for (const auto &pL : labels) {
    if (pL->cost() <= -epsilon_) {
      nbNegativeCosts++;
      *pMinCost = std::min(*pMinCost, pL->cost());
    }
  }
  return nbNegativeCosts;
}

std::vector<PRCLabel> MyRCSPPSolver::backwardLabelSetting(
    const std::vector<PRCNode> &sortedNodes, int finalDay) {
  vector<PRCLabel> destinationLabels;
  // Allow to visit each node once
  auto itN = sortedNodes.rbegin();
  for (; itN != sortedNodes.rend(); ++itN) {
    PRCNode pN = *itN;
    if (pN->type == SOURCE_NODE) break;
    if (pN->day < finalDay) break;

    // Expand all the labels of the predecessors of the node through the
    // arcs entering the node
    pushLabelsToPredecessors(pN);

    // information of domination of newly generated labels
    size_t nLabels = labelPool_.nLabels();
    vector<int> domStatusNew(nLabels, 0);

    // information of domination of former labels that have not been
    // expanded yet
    std::vector<PRCLabel> &labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
    size_t nLabelsNoExpand = labelsNoExpand.size();
    vector<int> domStatusNewNoExpand(nLabelsNoExpand, 0);

    // if expand all the labels at the node, dominate all labels first
    if (param_.rcsppNbToExpand_ == 0) {
      // Domination of all the labels of the current node to get a
      // pareto-optimal set of labels
      checkAllDominations(labelPool_.begin(),
                          labelPool_.end(),
                          !param_.rcsppDssr_ ? dominate : dominateNoLB,
                          &domStatusNew,
                          &domStatusNewNoExpand);
    }

    // For heuristic search, select only a subset of labels for expansion
    selectLabelsToExpand(labelPool_.begin(),
                         labelPool_.end(),
                         domStatusNew,
                         domStatusNewNoExpand, pN);

    // clear the pool (allow to reuse labels for next node)
    labelPool_.clear();
  }

  // Collect all the labels of the sink nodes
  for (const auto &pN : sortedNodes) {
    if (pN->day > finalDay) break;
    if (pN->day == finalDay) {
      for (const auto &pl : pLabelsToExpandPerNode_.at(pN->id))
        destinationLabels.push_back(pl);
    }
  }
  return destinationLabels;
}

void MyRCSPPSolver::pushLabelsToPredecessors(const PRCNode &pN) {
  // expand labels to all predecessors
  for (const auto &pArc : pN->inArcs) {
    if (pArc->forbidden) continue;  // do not use forbidden arcs
    // Expand all the labels
    for (const auto &pL : pLabelsToExpandPerNode_.at(pN->id)) {
      // Expansion on the current arc
      if (expandBack(pL, pArc, labelPool_.getNewLabel()))
        total_number_of_generated_labels_++;
      else
        number_of_infeasible_deleted_labels_++;
    }
  }
}
bool MyRCSPPSolver::expandBack(
    const PRCLabel &pLNext, const PRCArc &pArc, const PRCLabel &pLPrevious) {
  pLPrevious->setAsPrevious(pLNext, pArc);

  // Expansion of each resource
  for (auto &e : pArc->expanders) {
    // get child ResourceValues (copy of the parent ResourceValues)
    ResourceValues &v = pLPrevious->getResourceValues(e->resourceId);
    if (!e->expandBack(pLPrevious, &v)) {
      labelPool_.releaseLastLabel();
      return false;
    }
  }
  pLPrevious->addCost(pArc->baseCost + pArc->dualCost);
  return true;
}
