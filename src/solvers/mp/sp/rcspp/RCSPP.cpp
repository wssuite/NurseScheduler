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
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

// return domination status of pL2
DominationStatus dominate(RCLabel *pL1, RCLabel *pL2) {
  // Improved domination
  double cost1 = pL1->cost(), cost2 = pL2->cost();
  for (const auto &pR : pL1->getNode()->activeResources) {
    // if at anytime cost1 > cost2, return false
    if (cost1 >= cost2 + 1e-6) return NOT_DOMINATED;
    // if pL1 doesn't dominate pL2 on resource pR, stop
    // when using default domination, only DOMINATED counts
    // i.e.: UB_DOMINATED is like a NOT_DOMINATED when LB counts
    if (pR->isDefaultDomination()) {
      if (pR->dominates(pL1, pL2) != DOMINATED)
        return NOT_DOMINATED;
    } else {
      // will be just paying for LB as only active for soft resources
      if (pR->dominates(pL1, pL2, &cost1) == NOT_DOMINATED)
        return NOT_DOMINATED;
    }
  }
  return cost1 <= cost2 + 1e-6 ? DOMINATED : NOT_DOMINATED;
}

// Variants of the label setting functions when using the decrease state-space
// relaxation where lower bounds are first ignored for the constraints on the
// total number of shifts and the number of shifts per rotation
// returned statuses for pL2:
// - NOT_DOMINATED if no domination in any case
// - DOMINATED if domination with every constraint
// - UB_DOMINATED if domination when LBs are ignored
DominationStatus dominateNoLB(RCLabel *pL1, RCLabel *pL2)  {
  // We implement only improved domination for DSSR
  double cost1noLB = pL1->cost(), costOfLBs = 0;
  double cost2 = pL2->cost();
  bool dominated = true;
  // if cost1 > cost2, return false
  if (cost1noLB >= cost2  + 1e-6) return NOT_DOMINATED;
  for (const auto &pR : pL1->getNode()->activeResources) {
    // if hard resource, just check both bounds
    if (pR->isDefaultDomination()) {
      DominationStatus status = pR->dominates(pL1, pL2);
      if (status == NOT_DOMINATED) return NOT_DOMINATED;
      if (status == UB_DOMINATED) dominated = false;
      continue;
    }
    // evaluate resource
    double c = 0;
    DominationStatus status = pR->dominates(pL1, pL2, &c);
    if (status == NOT_DOMINATED) return NOT_DOMINATED;
    // ignore the lower bounds for any resource
    if (status == UB_DOMINATED)
      costOfLBs += c;
    // if at anytime cost1noLB > cost2, return false
    if (cost1noLB >= cost2  + 1e-6) return NOT_DOMINATED;
  }
  if (dominated && cost1noLB + costOfLBs <= cost2 + 1e-6)
    return DOMINATED;
  return UB_DOMINATED;
}

PRCLabel LabelPool::getNewLabel() {
  PRCLabel pL;
  // if at the end of the vector
  if (endLabel_ == pLabels_.end()) {
    // add it to the vector
    pLabels_.push_back(pFactory_->makePRCLabel());
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

const std::map<std::string, RCSPPSolver::SearchLevel>
    RCSPPSolver::searchLevelByName_ = {
    {"NB_TO_EXPAND", NB_TO_EXPAND_},
    {"DSSR_NO_LB", DSSR_NO_LB_},
    {"DSSR_LB", DSSR_LB_},
    {"INCR_DSSR_NO_LB", INCR_DSSR_NO_LB_},
    {"INCR_DSSR_LB", INCR_DSSR_LB_},
    {"OPTIMAL", OPTIMAL_}
};

const std::map<RCSPPSolver::SearchLevel, std::string>
    RCSPPSolver::namesBySearchLevel_ =
        Tools::buildNamesByType(RCSPPSolver::searchLevelByName_);

std::recursive_mutex RCSPPSolver::gMutex_;
bool RCSPPSolver::gCanComputeLB_ = true;
bool RCSPPSolver::gWarningHasBeenPrinted_ = false;
RCSPPSolver::SearchLevel RCSPPSolver::gMaxSearchLevel_ = OPTIMAL_;

RCSPPSolver::RCSPPSolver(PRCGraph pRcGraph,
                         SubProblemParam param,
                         int seed,
                         RCSPPSolver* pSolver) :
    pRcGraph_(std::move(pRcGraph)),
    pFactory_(std::make_shared<RCLabelFactory>()),
    labelPool_(pFactory_, 100),
    bestPrimalBound_(0.0),
    dssrLvl_(0),
    param_(param),
    maxReducedCostBound_(0),
    redCostLB_(-DBL_MAX),
    useDominateNoLb_(true),
    total_number_of_nondominated_labels_(0),
    number_of_infeasible_deleted_labels_(0),
    total_number_of_generated_labels_(0),
    total_number_of_dominations_(0),
    rdm_(Tools::getANewRandomGenerator(seed)),
    timer_("RCSPP solver"),
    maxSearchLevel_(pSolver ? pSolver->maxSearchLevel_ : OPTIMAL_),
    computingRelaxation_(false) {
  // if we are not solving to optimality and rcsppMinNegativeLabels_ has is
  // default value, override it based on the number of columns to generate
  if (param_.rcsppToOptimality_ && !param.rcsppRandomStartDay_) {
    param_.rcsppMaxNegativeLabels_ = XLARGE_SCORE;
  } else {
    if (param_.rcsppMinNegativeLabels_ == -1)
      param_.rcsppMinNegativeLabels_ = static_cast<int>(round(
          param_.spMinColumnsRatioForIncrease_ * param_.nbMaxColumnsToAdd_));
    if (param_.rcsppMaxNegativeLabels_ == -1)
      param_.rcsppMaxNegativeLabels_ = static_cast<int>(round(
          param_.spColumnsRatioForNumberPaths_ * param_.nbMaxColumnsToAdd_));
  }
}

void RCSPPSolver::resetLabels() {
  pExpandedLabelsPerNode_.clear();
  pLabelsToExpandPerNode_.clear();
  pLabelsNoExpandPerNode_.clear();
  int nNodes = pRcGraph_->nNodes();
  pExpandedLabelsPerNode_.resize(nNodes);
  pLabelsToExpandPerNode_.resize(nNodes);
  pLabelsNoExpandPerNode_.resize(nNodes);
}

std::vector<RCSolution> RCSPPSolver::solve(
    double maxReducedCostBound, bool relaxation) {
  maxReducedCostBound_ = maxReducedCostBound;
  maxSolvingTime_ = param_.spMaxSolvingTimeSeconds_;

  // set domination rules
  for (const auto &pR : pRcGraph_->pResources()) {
    if (param_.rcsppImprovedDomination_) pR->useAltenativeDomination();
    else
      pR->useDefaultDomination();
  }

  // compute a LB if possible (i.e., if it has enough time to compute an LB)
  if (relaxation) return solveRelaxation();

  // initialize the resources used for dominance
  // if dssrLvl_ == 0 (not set or to optimality)
  if (dssrLvl_ == 0) {
    dssrLvl_ = param_.rcsppDssr_ ? pRcGraph_->nDays() : 0;
    pRcGraph_->initializeDominance(dssrLvl_);
  }
  if (!param_.rcsppDssr_) useDominateNoLb_ = false;

  // normal solve
  return solve();
}

std::vector<RCSolution> RCSPPSolver::solve() {
  redCostLB_ = -DBL_MAX;
  timeRanOut_ = false;
  timer_.start();

  // Displaying of the activated options for the algorithm
  displaySolveInputInfo();

  bestPrimalBound_ = DBL_MAX;
  vector<PRCLabel> pLabelsSinks;
  vector<RCSolution> finalSolutions;
  // sort nodes as the graph should be acyclic
  vector<PRCNode> sortedNodes = pRcGraph_->sortNodes();
  while (true) {
    // reset counters for infeasible labels
    number_of_infeasible_deleted_labels_per_resource_.clear();
    number_of_infeasible_deleted_labels_per_resource_.resize(
        pRcGraph_->pResources().size());

    // initialize the search level based on the current parameters
    initSearchLevel();

    // should normally be solved until optimality if no heuristic activated
    isOptimal_ = param_.rcsppToOptimality_;
    vector<PRCLabel> destinationLabels;
    // use bidirectional, only if rcsppNbToExpand_ is deactivated
    if (param_.rcsppBidirectional_ && param_.rcsppNbToExpand_ == 0) {
      destinationLabels = bidirectionalLabelSetting(sortedNodes);
    } else {
      destinationLabels = forwardLabelSetting(sortedNodes,
                                              pRcGraph_->nDays() - 1);
    }
    pLabelsSinks = Tools::appendVectors(pLabelsSinks, destinationLabels);

    // Sorting by increasing cost all the labels of the sink nodes
    int nLabels = getNbLabelsBelowMaxBound(pLabelsSinks, &bestPrimalBound_);
    if (param_.verbose_ >= 3) {
      std::cout << "Current best primal bound = "
                << bestPrimalBound_ << std::endl;
      std::cout << "Number of negative cost labels = " << nLabels
                << std::endl;
    }

    // check if enough time left
    if (!enoughTimeLeft()) {
      isOptimal_ = false;
      break;
    }

    if (param_.rcsppToOptimality_) {
      // if no heuristic option is activated, we can break, because
      // optimality must be reached
      if (!param_.rcsppDssr_ && !param_.rcsppNbToExpand_)
        break;
      else if (!prepareForNextExecution(pRcGraph_->pNodes()))
        break;  // should never be reached
    } else {
      // if any heuristic option is activated, and we got enough negative cost
      // rosters, we can break, because we got what we were looking for
      // otherwise, we deactivate the most restrictive heuristic and get
      // prepared for a new solution
      if (nLabels >= param_.rcsppMinNegativeLabels_)
        break;
      else if (!prepareForNextExecution(pRcGraph_->pNodes()))
        break;
    }
  }

  // stop timer
  timer_.stop();

  // Creation of the vector containing the negative cost rosters
  // first, sort by increasing cost all the labels of the sink nodes
  std::sort(pLabelsSinks.begin(), pLabelsSinks.end(), LabelCostIncreasing());

  for (const auto &pL : pLabelsSinks) {
    RCSolution solution = createSolution(pL);
#ifdef NS_DEBUG
    //    if (finalSolutions.empty()) {
    //      std::cout << "Label bound: " << maxReducedCostBound << std::endl;
    //      createSolution(pL, pRcGraph_);
    //    }
#endif
    if (pL->cost() + param_.epsilon_ < maxReducedCostBound_) {
      finalSolutions.push_back(solution);
    } else {
      break;
    }
  }
  displaySolveStatistics();

  if (isOptimal_) {
    redCostLB_ = maxReducedCostBound_;
    for (const RCSolution &sol : finalSolutions)
      if (sol.reducedCost() < redCostLB_) redCostLB_ = sol.reducedCost();
  }

  return finalSolutions;
}

std::vector<RCSolution> RCSPPSolver::solveRelaxation() {
  // if wants a LB, check if was able to solve with a search level greater
  // than DSSR_NO_LB_
  if (gCanComputeLB_) {
    maxSolvingTime_ =
        maxSolvingTime_ * param_.spMaxSolvingTimeRatioForRelaxation_;
    bool relaxed = maxSearchLevel_ <= DSSR_NO_LB_;
    auto sols = computeValidLB(relaxed);
    // stop now if asked
    if (job_.shouldStop()) return sols;
    // if wasn't able to compute an LB with enforced LBs, try without LBs
    if (timeRanOut_ && !relaxed)
      sols = computeValidLB(true);
    // if time ran out, stop computing an LB
    if (timeRanOut_) {
      std::lock_guard<std::recursive_mutex> l(gMutex_);
      gCanComputeLB_ = false;
    }
    // clean the graph as infeasible labels could now exist
    resetLabels();
    return sols;
  }
  return {};
}

// alert other RCSPPSolver that time ran out for this solver
void RCSPPSolver::timeRanOut() {
  // if already ran out of time, nothing to do
  if (timeRanOut_) return;
  timeRanOut_ = true;
  std::lock_guard<std::recursive_mutex> l(gMutex_);
  // if computing LB, nothing else
  if (computingRelaxation_)
    return;
  // if already at minimum level, nothing to do
  if (maxSearchLevel_ == NB_TO_EXPAND_) return;
  SearchLevel newLevel = SearchLevel(searchLevel_ - 1);
  if (newLevel == NB_TO_EXPAND_) {
    newLevel = DSSR_NO_LB_;
    if (!gWarningHasBeenPrinted_)
      std::cout << "Max search was not decrease to "
                << namesBySearchLevel_.at(NB_TO_EXPAND_)
                << ", but " << namesBySearchLevel_.at(DSSR_NO_LB_)
                << " instead. You should consider giving more computational "
                   "time to the subproblems (currently "
                << param_.spMaxSolvingTimeSeconds_
                << " seconds)." << std::endl;
    gWarningHasBeenPrinted_ = true;
  }
  if (maxSearchLevel_ > newLevel) {
    if (param_.verbose_ >= 2)
      std::cout << "RCSPPSolver: Decrease max search to level "
                << namesBySearchLevel_.at(newLevel) << std::endl;
    maxSearchLevel_ = newLevel;
  }
  if (gMaxSearchLevel_ > newLevel) {
    if (param_.verbose_ >= 1)
      std::cout << "RCSPPSolver: Decrease global max search to level "
                << namesBySearchLevel_.at(newLevel) << std::endl;
    gMaxSearchLevel_ = newLevel;
  }
}

// check if there is enough time left to continue the resolution
// if not, update timeRanOut
bool RCSPPSolver::enoughTimeLeft() {
  // always perform at least NB_TO_EXPAND
  if (searchLevel_ == NB_TO_EXPAND_) return true;
  // if exceed time limit
  if (timer_.dSinceStart() > maxSolvingTime_)
    timeRanOut();
  else if (computingRelaxation_)
    return gCanComputeLB_ && !job_.shouldStop();
  if (param_.rcsppUseGlobalMaxLevel_ && gMaxSearchLevel_ < maxSearchLevel_)
    maxSearchLevel_ = gMaxSearchLevel_;
  if (searchLevel_ <= maxSearchLevel_) return true;
  return false;
}

void RCSPPSolver::initSearchLevel() {
  if (param_.rcsppNbToExpand_ > 0) {
    searchLevel_ = NB_TO_EXPAND_;
  } else if (param_.rcsppToOptimality_) {
    searchLevel_ = OPTIMAL_;
  } else if (param_.rcsppDssr_) {
    if (dssrLvl_ == pRcGraph_->nDays())
      searchLevel_ = useDominateNoLb_ ? DSSR_NO_LB_ : DSSR_LB_;
    else
      searchLevel_ = useDominateNoLb_ ? INCR_DSSR_NO_LB_ : INCR_DSSR_LB_;
  } else {
    Tools::throwError("There is no valid RCSPP search level for those "
                      "current parameters");
  }
}

// compute a LB by using relaxing many resources
std::vector<RCSolution> RCSPPSolver::computeValidLB(bool relaxedLBs) {
  // use only necessary resources
  std::set<int> alwaysActiveResources;
  std::vector<PResource> otherResources;
  for (const auto &pR : pRcGraph_->pResources()) {
    if (pR->isActive(LARGE_SCORE))
      alwaysActiveResources.insert(pR->id());
    else
      otherResources.push_back(pR);
  }

  // deactivate the others expanders
  vector2D<PExpander> expanders(pRcGraph_->nArcs());
  for (const PRCArc &pA : pRcGraph_->pArcs()) {
    expanders[pA->id] = pA->expanders;
    std::vector<PExpander> activeExpanders;
    for (const PExpander &pE : pA->expanders)
      if (alwaysActiveResources.find(pE->indResource) !=
          alwaysActiveResources.end())
        activeExpanders.push_back(pE);
    pA->expanders = activeExpanders;
  }

  // clean the graph as a looser dominance would be used
  resetLabels();

  // solve the graph (and update redCostLB_ in the solve function)
  SubProblemParam spParam = param_;
  param_.rcsppToOptimality_ = true;
  param_.rcsppNbToExpand_ = 0;
  param_.rcsppDssr_ = false;
  param_.rcsppMaxNegativeLabels_ = XLARGE_SCORE;
  dssrLvl_ = 0;
  pRcGraph_->initializeDominance(0);
  useDominateNoLb_ = relaxedLBs;
  computingRelaxation_ = true;
  vector<RCSolution> sols = solve();
  computingRelaxation_ = false;
  param_ = spParam;

  if (param_.verbose_ >= 2)
    std::cout << "LB found in " << timer_.dSinceStart() << "seconds"
              << (timeRanOut_ ? " (time ran out)" : "")
              << (job_.shouldStop() ? " (has been stopped)" : "") << ": "
              << redCostLB_ << " for " << (relaxedLBs ? "relaxed" : "enforced")
              << " LBs." << std::endl;

  // reactivate the others expanders
  for (const PRCArc &pA : pRcGraph_->pArcs())
    pA->expanders = expanders[pA->id];

  // if time ran out or should stop, nothing has been produced
  if (timeRanOut_ || job_.shouldStop())
    return {};

  // check feasibility of the paths
  vector<RCSolution> finalSolutions;
  double bestSolRedCost = maxReducedCostBound_;
  for (const RCSolution &sol : sols) {
    bool feasible = true;
    // all the other resources must be feasible
    const PRCNode &pO = pRcGraph_->pSource(sol.firstDayId()),
        &pT = pRcGraph_->pNodesPerDayId(sol.lastDayId()).at(
        sol.pShifts().back()->type);
    PRCArc pA = std::make_shared<RCArc>(0, pO, pT, sol, 0, TO_SINK);
    PRCLabel pL = std::make_shared<RCLabel>(*pLSources_[sol.firstDayId()]);
    try {
      for (const PResource &pR : otherResources) {
        auto pE = pR->initialize(*pO->pAShift, sol, pA, pR->id());
        ResourceValues &v = pL->getResourceValues(pR->id());
        if (pE && !pE->expand(pL, &v)) {
          feasible = false;
          break;
        }
      }
      if (feasible) {
        finalSolutions.push_back(sol);
        if (sol.reducedCost() < bestSolRedCost)
          bestSolRedCost = sol.reducedCost();
      }
    } catch (...) {}
  }

  // update optimality
  isOptimal_ = (redCostLB_ + param_.epsilon_ > bestSolRedCost);

  return finalSolutions;
}

std::vector<PRCLabel> RCSPPSolver::forwardLabelSetting(
    const std::vector<PRCNode> &sortedNodes, int finalDay) {
  // vector of non-dominated labels at the nodes of day=finalDay
  vector<PRCLabel> destinationLabels;

  // add the initial labels to the expansion list
  for (const PRCLabel &pL : pLSources_) addLabelToExpand(pL);

  // Allow to visit each node once
  int nLabelsFound = 0;
  int k = param_.rcsppRandomStartDay_ ? rdm_() % pLSources_.size() : 0;
  const PRCNode &pFirstSource = pRcGraph_->pSources()[k];
  auto it = find(sortedNodes.begin(), sortedNodes.end(), pFirstSource);
  bool hasExpandedLabels = false;
  int currentDay = (*it)->dayId;
  int nProcessedNodes = 0;
  // loop while has done one whole tour or
  // has finished to treat nodes of the current day or
  // at least one label has been expanded on current day
  while (nProcessedNodes < sortedNodes.size() ||
         currentDay == (*it)->dayId || hasExpandedLabels) {
    // check if enough time left
    if (!enoughTimeLeft())
      break;

    const PRCNode &pN = *it;

    // update day and hasExpandedLabels
    if (currentDay != pN->dayId) {
      hasExpandedLabels = false;
      currentDay = pN->dayId;
    }

    // Expand all the labels of the predecessors of the node through the
    // arcs entering the node if starts early enough
    pullLabelsFromPredecessors(pN, finalDay - 1);

    // information of domination of newly generated labels
    size_t nLabels = labelPool_.nLabels();
    vector<DominationStatus> domStatusNew(nLabels, NOT_DOMINATED);

    // information of domination of former labels that have not been
    // expanded yet
    std::vector<PRCLabel> &labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
    size_t nLabelsNoExpand = labelsNoExpand.size();
    vector<DominationStatus> domStatusNewNoExpand(
        nLabelsNoExpand, NOT_DOMINATED);

    // if expand all the labels at the node, dominate all labels first
    if (param_.rcsppNbToExpand_ == 0) {
      // copy any labelsToExpand in expanded (could be there from previous pass)
      if (pN->type != SOURCE_NODE &&
          !pLabelsToExpandPerNode_[pN->id].empty()) {
        pExpandedLabelsPerNode_[pN->id].insert(
            pExpandedLabelsPerNode_[pN->id].end(),
            pLabelsToExpandPerNode_[pN->id].begin(),
            pLabelsToExpandPerNode_[pN->id].end());
        pLabelsToExpandPerNode_[pN->id].clear();
      }

      // Domination of all the labels of the current node to get a
      // pareto-optimal set of labels
      checkAllDominations(labelPool_.begin(),
                          labelPool_.end(),
                          useDominateNoLb_ ? dominateNoLB : dominate,
                          &domStatusNew,
                          &domStatusNewNoExpand);
    }

    // check if any not_dominated new labels
    for (auto s : domStatusNew)
      if (s == NOT_DOMINATED) {
        hasExpandedLabels = true;
        break;
      }

    // For heuristic search, select only a subset of labels for expansion
    selectLabelsToExpand(labelPool_.begin(),
                         labelPool_.end(),
                         domStatusNew,
                         domStatusNewNoExpand, pN);

    if (param_.verbose_ >= 4)
      std::cout << "RCSSP: process node " << pN->toString()
                << " (active resources=" << pN->activeResources.size()
                << ", expanded labels=" << nLabels
                << ", non dominated=" <<  pLabelsToExpandPerNode_[pN->id].size()
                << ")" << std::endl;

    // clear the pool (allow to reuse labels for next node)
    labelPool_.clear();

    // Collect all the labels of the final day nodes
    // we test on the dayId for the bidirectional label-setting
    // we test also on sink for the rotations
    if (pN->dayId >= finalDay || pN->type == SINK_NODE) {
      const auto &pLabels = pLabelsToExpandPerNode_.at(pN->id);
      for (const auto &pL : pLabels)
        if (pL->getOutArc() == nullptr)
          destinationLabels.push_back(pL);
      // count labels that have reached a sink
      if (pN->type == SINK_NODE) {
        double minCost;
        nLabelsFound += getNbLabelsBelowMaxBound(pLabels, &minCost);
        if (nLabelsFound >= param_.rcsppMaxNegativeLabels_) {
          isOptimal_ = false;
          break;
        }
      }
    }

    // next node
    if (++it == sortedNodes.end())
      it = sortedNodes.begin();
    ++nProcessedNodes;
  }

  return destinationLabels;
}

void RCSPPSolver::pullLabelsFromPredecessors(
    const PRCNode& pN, int predecessorMaxDay) {
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
      return hasPotentialImprovingPathToSinks(pL, pArc->origin->id, -1e-6) &&
          expand(pL, pArc, labelPool_.getNewLabel());
    };
  }

  // expand labels from all predecessors
  for (const auto &pArc : pN->inArcs) {
    if (pArc->forbidden) continue;  // do not use forbidden arcs
    // arcs that do not start early enough
    if (pArc->origin->dayId > predecessorMaxDay) continue;
    for (const auto &pL : pLabelsToExpandPerNode_.at(pArc->origin->id)) {
      // Expand the labels
      if (func(pL, pArc)) total_number_of_generated_labels_++;
      else
        number_of_infeasible_deleted_labels_++;
    }
  }
}

bool RCSPPSolver::expand(
    const PRCLabel &pLParent, const PRCArc &pArc, const PRCLabel &pLChild) {
  pLChild->setAsNext(pLParent, pArc);

  // Expansion of each resource
  for (auto &e : pArc->expanders) {
    // get child ResourceValues (copy of the parent ResourceValues)
    ResourceValues &v = pLChild->getResourceValues(e->indResource);
    if (!e->expand(pLChild, &v)) {
      labelPool_.releaseLastLabel();
      number_of_infeasible_deleted_labels_per_resource_[e->indResource]++;
      return false;
    }
  }

  pLChild->addBaseCost(pArc->baseCost);
  pLChild->addCost(-pArc->dualCost);
#ifdef NS_DEBUG
  pLChild->addDualCost(pArc->dualCost);
  pLChild->addPreferencesCost(pArc->baseCost);
#endif
  return true;
}

void RCSPPSolver::checkDominationsWithPreviousLabels(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    DominationStatus (*domFunction)(RCLabel *, RCLabel *),
    vector<DominationStatus> *domStatus,
    vector<DominationStatus> *domStatusNoExpand) {
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
    if (*itD1 == DOMINATED) continue;

    // check if the expanded labels dominate the newly generated ones
    for (const auto& pL2 : expandedLabels) {
      total_number_of_dominations_++;
      status = domFunction(pL2.get(), itL1->get());
      if (status == DOMINATED) {
        *itD1 = DOMINATED;
        break;
      }
      if (status == UB_DOMINATED) *itD1 = UB_DOMINATED;
    }
    if (*itD1 == DOMINATED) continue;

    // check if the non-expanded labels dominate the newly generated ones
    for (const auto& pL2 : labelsNoExpand) {
      total_number_of_dominations_++;
      status = domFunction(pL2.get(), itL1->get());
      if (status == DOMINATED) {
        *itD1 = DOMINATED;
        break;
      }
      if (status == UB_DOMINATED) *itD1 = UB_DOMINATED;
    }
    if (*itD1 == DOMINATED) continue;

    // check if the non-dominated newly generated labels dominate the stored
    // labels that have not been expanded yet
    auto itL2 = labelsNoExpand.begin();
    for (auto itD2 = domStatusNoExpand->begin();
         itD2 != domStatusNoExpand->end();
         ++itD2, ++itL2) {
      if (*itD2 == DOMINATED) continue;
      total_number_of_dominations_++;
      status = domFunction(itL1->get(), itL2->get());
      if (status == DOMINATED) *itD2 = DOMINATED;
      else if (status == UB_DOMINATED) *itD2 = UB_DOMINATED;
    }
  }
}

void RCSPPSolver::checkDominationsPairwise(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    DominationStatus (*domFunction)(RCLabel *, RCLabel *),
    vector<DominationStatus> *domStatus) {
  // compare newly generated labels pairwise for domination
  auto itL1 = begin;
  for (auto itD1 = domStatus->begin();
       itD1 != domStatus->end();
       ++itD1, ++itL1) {
    if (*itD1 == DOMINATED) continue;

    auto itL2 = itL1 + 1;
    int status;
    for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
      if (*itD2 == DOMINATED) continue;

      total_number_of_dominations_++;
      status = domFunction(itL1->get(), itL2->get());
      if (status == DOMINATED) *itD2 = DOMINATED;
      else if (status == UB_DOMINATED) *itD2 = UB_DOMINATED;
    }
    itL2 = itL1 + 1;
    if (param_.rcsppSortLabels_) {
      for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
        if (*itD2 == DOMINATED) continue;

        total_number_of_dominations_++;
        status = domFunction(itL2->get(), itL1->get());
        if (status == DOMINATED) *itD1 = DOMINATED;
        else if (status == UB_DOMINATED) *itD1 = UB_DOMINATED;
        if (*itD1 == DOMINATED) break;

        // when the labels are sorted by increasing cost, domination of L2
        // over L1 can be stopped being checked as soon as we find L2 with
        // larger cost than L1
        if ((*itL2)->cost() > (*itL1)->cost()) break;
      }
    } else {
      for (auto itD2 = itD1 + 1; itD2 != domStatus->end(); ++itD2, ++itL2) {
        if (*itD2 == DOMINATED) continue;

        total_number_of_dominations_++;
        status = domFunction(itL2->get(), itL1->get());
        if (status == DOMINATED) *itD1 = DOMINATED;
        else if (status == UB_DOMINATED) *itD1 = UB_DOMINATED;
        if (*itD1 == DOMINATED) break;
      }
    }
  }
}

void RCSPPSolver::checkAllDominations(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    DominationStatus (*domFunction)(RCLabel *, RCLabel *),
    vector<DominationStatus> *domStatus,
    vector<DominationStatus> *domStatusNoExpand) {
  // if no labels, returns
  if (begin == end) return;

  // Sort the labels by increasing order before checking domination
  if (param_.rcsppSortLabels_)
    std::sort(begin, end, LabelCostIncreasing());

  // check for domination between the newly generated labels and those that
  // are already stored at the node
  checkDominationsWithPreviousLabels(begin, end,
                                     domFunction,
                                     domStatus, domStatusNoExpand);

  // compare newly generated labels pairwise for domination
  checkDominationsPairwise(begin, end, domFunction, domStatus);
}


void RCSPPSolver::selectLabelsToExpand(
    const vector<PRCLabel>::iterator &begin,
    const vector<PRCLabel>::iterator &end,
    const vector<DominationStatus>& domStatus,
    const vector<DominationStatus>& domStatusNoExpand,
    const PRCNode &pN) {
  // if no labels, returns
  size_t nLabels = std::distance(begin, end);
  vector<PRCLabel>& labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
  size_t nLabelsNoExpand = labelsNoExpand.size();
  if (nLabels + nLabelsNoExpand == 0) return;

  // get the vector of labels to expand
  vector<PRCLabel>& labelsToExpand = pLabelsToExpandPerNode_[pN->id];

  // set the list of non-dominated labels for current node
  if (param_.rcsppNbToExpand_ > 0) {
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
    for (DominationStatus d : domStatus) {
      if (d == NOT_DOMINATED) {
        // expand the label if not dominated with relaxed constraints
        labelsToExpand.push_back(pFactory_->makePRCLabel(**itL));
        total_number_of_nondominated_labels_++;
      } else if (d == UB_DOMINATED) {
        // store the label without expanding it if not dominated with every
        // constraint
        labelsNoExpand.push_back(pFactory_->makePRCLabel(**itL));
      }
      ++itL;
    }
  } else {
    auto itL = labelsNoExpand.begin();
    // add the non-expanded non-dominated labels
    for (DominationStatus d : domStatusNoExpand) {
      if (d == NOT_DOMINATED) {
        // expand the label if not dominated with every constraint
        labelsToExpand.push_back(*itL);
        total_number_of_nondominated_labels_++;
      }
      ++itL;
    }
    labelsNoExpand.clear();
    itL = begin;
    // add the newly generated non-dominated labels
    for (int d : domStatus) {
      if (d == NOT_DOMINATED) {
        // expand the label if not dominated with every constraint
        labelsToExpand.push_back(pFactory_->makePRCLabel(**itL));
        total_number_of_nondominated_labels_++;
      }
      ++itL;
    }
  }
}

bool RCSPPSolver::prepareForNextExecution(const vector<PRCNode> &nodes) {
  // set optimality to false, as a next execution is needed
  isOptimal_ = false;
  // go over activated heuristic options from the most to the least aggressive
  if (param_.rcsppNbToExpand_ > 0) {
    // if we have expanded only a  small subset of labels at each node, we did
    // not perform any domination, so we just start the algorithm over
    if (param_.verbose_ >= 3) {
      std::cout << "Did not find negative cost roster by expanding "
                << param_.rcsppNbToExpand_ << " labels, switch to "
                << (param_.rcsppDssr_ ? "dssr" : "optimality") << std::endl;
    }
    param_.rcsppNbToExpand_ = 0;
    if (!param_.rcsppDssr_)
      param_.rcsppToOptimality_ = true;

    for (const auto &pN : nodes)
      if (pN->type != SOURCE_NODE)
        pLabelsToExpandPerNode_[pN->id].clear();
    // stop current execution (return false) if enable and
    // max search level is DSSR_NO_LB_
    return !param_.rcsspWaitBeforeStartingNextExecution_
        || maxSearchLevel_ != DSSR_NO_LB_;
  } else if (param_.rcsppDssr_) {
    // if we used a decremental state-space relaxation, we may wish to keep the
    // expansion and domination information we got during the execution and
    // warm-start the new execution of the label-setting algorithm
    // move the labels that were expanded to the list of expanded labels
    for (const auto &pN : nodes) {
      std::vector<PRCLabel>
          &expandedLabels = pExpandedLabelsPerNode_[pN->id];
      std::vector<PRCLabel>
          &labelsToExpand = pLabelsToExpandPerNode_[pN->id];
      expandedLabels = Tools::appendVectors(expandedLabels, labelsToExpand);
      labelsToExpand.clear();
    }

    // DSSR will deactivate some hard constraints
    // compute always active constraints
    bool anyHardConstraints = false, hasLB = false;
    std::set<int> alwaysActiveResources;
    for (const auto &pR : pRcGraph_->pResources()) {
      if (pR->isActive(LARGE_SCORE))
        alwaysActiveResources.insert(pR->id());
      if (pR->isHard())
        anyHardConstraints = true;
      auto pRB = std::dynamic_pointer_cast<BoundedResource>(pR);
      if (pRB && pRB->getLb() > 0) hasLB = true;
    }

    // if all resources are always active: go to optimality (dssr lvl = 0)
    bool sameSearchLevel = false;
    if (pRcGraph_->pResources().size() == alwaysActiveResources.size())
      dssrLvl_ = 0;

    // if not incremental: go to no lb or optimality (dssr lvl = 0)
    if (!param_.rcsppIncrementalDssr_) {
      // DSSR_NO_LB -> DSSR_LB
      if (searchLevel_ == DSSR_NO_LB_ && hasLB) {
        useDominateNoLb_ = false;
        if (param_.verbose_ >= 3)
          std::cout << "Solving with dssr and LBs." << std::endl;
        return !param_.rcsspWaitBeforeStartingNextExecution_;
      }
      // optimality
      dssrLvl_ = 0;
    } else {
      // Incremental dssr
      // First, find active resources
      vector<bool> activeResources(pRcGraph_->pResources().size(), false);
      for (const auto &pN : nodes)
        for (const PResource &pR : pN->activeResources)
          activeResources[pR->id()] = true;

      // Second: DSSR_NO_LB -> INCR_DSSR_NO_LB (ignore LBs and adding resources)
      // or already INCR_DSSR_NO_LB
      bool activated = false;
      if (useDominateNoLb_) {
        if (searchLevel_ == DSSR_NO_LB_ && param_.verbose_ >= 3)
          std::cout << "Solving with incremental dssr and without LB."
                    << std::endl;
        useDominateNoLb_ = true;
        sameSearchLevel = (searchLevel_ == INCR_DSSR_NO_LB_);
        // check if any inactive resource for domination has eliminated
        // some infeasible labels. If the case, add it to the domination
        for (const auto &pR : pRcGraph_->pResources()) {
          if (!pR->isActive(dssrLvl_) &&
              number_of_infeasible_deleted_labels_per_resource_[pR->id()] > 0) {
            for (const auto &pN : nodes)
              if (pN->activateResourceIfNecessary(pR))
                activated = true;
            if (activated) {
              if (param_.verbose_ >= 3)
                std::cout << "Activating resource " << pR->id() << " ("
                          << pR->name << ") for dssr no LB" << std::endl;
              return true;
            }
          }
        }
        // Then, decrease dssrLvl_ until one new resource is activated
        int dssrLvl = dssrLvl_;
        while (dssrLvl_ > 0) {
          for (const auto &pR : pRcGraph_->pResources()) {
            if (pR->isActive(dssrLvl_) && !activeResources[pR->id()]) {
              for (const auto &pN : nodes)
                pN->activateResourceIfNecessary(pR);
              if (param_.verbose_ >= 3)
                std::cout << "Activating resource " << pR->id() << " ("
                          << pR->name << ") for dssr no LB." << std::endl;
              activated = true;
              break;
            }
          }
          if (activated) break;
          dssrLvl_--;
        }
        // if reach 0 (i.e., not activated), reset if for DSSR_LB
        if (dssrLvl_ == 0) {
          // DSSR_LB if any resources with LB
          if (hasLB) {
            useDominateNoLb_ = false;
            dssrLvl_ = pRcGraph_->nDays();
            // make all resources inactive except the always present ones
            pRcGraph_->initializeDominance(dssrLvl_);
            if (param_.verbose_ >= 3)
              std::cout << "Solving with dssr and LBs." << std::endl;
            return !param_.rcsspWaitBeforeStartingNextExecution_;
          }
          // optimality
          dssrLvl_ = 0;
        } else if (dssrLvl_ < dssrLvl) {
          if (param_.verbose_ >= 3)
            std::cout << "Increase dssr lvl to " << dssrLvl_ << std::endl;
        }
      }
      // Third, INCR_DSSR_LB: Activate next resource, deactivate current one
      // We must have more than one resource to activate,
      // otherwise it's the same than solving to optimality
      if (!useDominateNoLb_ &&
          pRcGraph_->pResources().size() > alwaysActiveResources.size() + 1) {
        // find current activated resource if any
        int maxR = 0;
        for (const auto &pN : nodes)
          for (const PResource &pR : pN->activeResources)
            if (alwaysActiveResources.find(pR->id()) ==
                alwaysActiveResources.end() && pR->id() > maxR)
              maxR = pR->id();
        // initialize dominance to remove all the unnecessary resources
        pRcGraph_->initializeDominance(pRcGraph_->nDays());
        // find next resource
        dssrLvl_ = 1;
        for (const auto &pR : pRcGraph_->pResources()) {
          // select an inactive resource bigger than current activated one
          if (activeResources[pR->id()] || !pR->isActive(dssrLvl_)
              || pR->id() < maxR)
            continue;
          // if no always active constraints has any LB,
          // ensure the new resource has one
          if (!hasLB) {
            auto pRB = std::dynamic_pointer_cast<BoundedResource>(pR);
            if (!pRB || pRB->getLb() == 0) continue;
          }
          // activate resource at every node
          for (const auto &pN : nodes)
            pN->activateResourceIfNecessary(pR);
          if (param_.verbose_ >= 3)
            std::cout << "Activating resource " << pR->id() << " ("
                      << pR->name << ") for dssr with LBs." << std::endl;
          activated = true;
          sameSearchLevel = (searchLevel_ == INCR_DSSR_LB_);
          break;
        }
      }

      // if new resources are taken into account in the dominance,
      // labels must be reset
      if (activated)
        resetLabels();
      else
        dssrLvl_ = 0;  // solve to optimality
    }

    // try other subproblems when going from dssr to incremental dssr
    if (dssrLvl_ > 0)
      return !param_.rcsspWaitBeforeStartingNextExecution_ || sameSearchLevel;

    // set to optimality (dssrLvl = 0)
    useDominateNoLb_ = false;
    param_.rcsppDssr_ = false;
    dssrLvl_ = 0;
    param_.rcsppToOptimality_ = true;
    pRcGraph_->initializeDominance(0);
    // Has LB were totally ignored, domination issues could happen
    // -> reset
    if (anyHardConstraints)
      resetLabels();
    if (param_.verbose_ >= 2)
      std::cout << "Solving to optimality." << std::endl;
    // stop process here if set this way, and wait to being called again
    return !param_.rcsspWaitBeforeStartingNextExecution_;
  }
  // it has been solved to optimality even if it wasn't a parameters option
  // should never be reached
  Tools::throwError("RCSPP parameters are not set correctly, as no heuristics "
                    "are enable and, at the same time, the optimality was not "
                    "set to true.");
  param_.rcsppToOptimality_ = true;
  return false;
}

RCSolution RCSPPSolver::createSolution(
    const PRCLabel &finalLabel, const PRCGraph pRcGraph) {
  // Backtrack from the label
  PRCLabel pL = finalLabel;
  Stretch stretch(finalLabel->getNode()->dayId);
  if (pRcGraph)
    std::cout << "===========================================" << std::endl;
  if (pRcGraph)
    std::cout << pL->toString(pRcGraph->pResources()) << std::endl;
  while (pL->getPreviousLabel() != nullptr) {
    stretch.pushFront(pL->getInArc()->stretch);
    pL = pL->getPreviousLabel();
    if (pRcGraph)
      std::cout << pL->toString(pRcGraph->pResources()) << std::endl;
  }
  if (pRcGraph)
    std::cout << "===========================================" << std::endl;
  // "Track forward" from the label if some backward expansion
  pL = finalLabel;
  while (pL->getNextLabel() != nullptr) {
    stretch.pushBack(pL->getOutArc()->stretch);
    pL = pL->getNextLabel();
  }
  return RCSolution(stretch, finalLabel);
}

bool RCSPPSolver::hasPotentialImprovingPathToSinks(
    const PRCLabel &pl, int nodeId, double primalBound) const {
  // If a path to a sink node with a negative cost exists, this label can be
  // expanded
  return pl->cost() + minimumCostToSinks_[nodeId] < primalBound;
}

int RCSPPSolver::getNbLabelsBelowMaxBound(const vector<PRCLabel> &labels,
                                          double *pMinCost) const {
  *pMinCost = 0;
  int nbNegativeCosts = 0;
  for (const auto &pL : labels) {
    if (pL->cost() <= maxReducedCostBound_ - param_.epsilon_) {
      nbNegativeCosts++;
      *pMinCost = std::min(*pMinCost, pL->cost());
    }
  }
  return nbNegativeCosts;
}


std::vector<PRCLabel> RCSPPSolver::bidirectionalLabelSetting(
    const std::vector<PRCNode> &sortedNodes) {
  int middleDay = this->pRcGraph_->nDays()/2;
  this->initializeLabels();
  vector<PRCLabel> forwardLabels = forwardLabelSetting(sortedNodes, middleDay);
  // store nodes where a merge should happen
  std::vector<PRCNode> mergingNodes;
  for (const PRCNode &pN : sortedNodes)
    if (pN->dayId >= middleDay && !pLabelsToExpandPerNode_[pN->id].empty())
      mergingNodes.push_back(pN);
//  resetLabels();
  vector<PRCLabel> backwardLabels =
      backwardLabelSetting(sortedNodes, mergingNodes, middleDay);

  // check if enough time left
  if (!enoughTimeLeft())
    return {};

  // merge the labels built by forward and backward propagation
  if (param_.verbose_ >= 4)
    std::cout << "RCSSP: merging ..." << std::flush;
  mergeLabels(&forwardLabels, &backwardLabels);
  if (param_.verbose_ >= 4)
    std::cout << " Done" << std::endl;

  // store the resulting negative cost labels
  vector<PRCLabel> negativeLabels;
  auto itL = labelPool_.begin();
  for (; itL != labelPool_.end(); ++itL)
    negativeLabels.push_back(pFactory_->makePRCLabel(**itL));
  labelPool_.clear();

  return negativeLabels;
}

void RCSPPSolver::mergeLabels(
    vector<PRCLabel> *pForwardLabels,
    vector<PRCLabel> *pBackwardLabels
) {
  // sort the labels by increasing costs to stop the merging process
  // prematurely when no negative cost roster can happen from merging
  std::sort(pForwardLabels->begin(), pForwardLabels->end(),
            LabelCostIncreasing());
  std::sort(pBackwardLabels->begin(), pBackwardLabels->end(),
            LabelCostIncreasing());

  for (const auto& pLForward : *pForwardLabels) {
    for (const auto& pLBackward : *pBackwardLabels) {
      // break if no future merge can yield a better negative cost roster
      if (pLForward->cost() + pLBackward->cost() >
          std::min(bestPrimalBound_, 0.0) - param_.epsilon_) break;
      // merge only the labels that are hosted by the same node
      if (pLForward->getNode() != pLBackward->getNode()) continue;
      // merge the two labels
      merge(pLForward, pLBackward);
    }
  }
}

bool RCSPPSolver::merge(const PRCLabel& pLForward,
                        const PRCLabel& pLBackward) {
  // get a new label from the pool for the result of merging
  PRCLabel pLMerged = labelPool_.getNewLabel();
  if (pLMerged->num() == 19931)
    int b = 0;
  pLMerged->setAsMerged(pLForward, pLBackward);
  // merge every resource consumption
  for (const auto &pR : pRcGraph_->pResources()) {
    const ResourceValues &vForward = pLForward->getResourceValues(pR->id());
    const ResourceValues &vBackward = pLBackward->getResourceValues(pR->id());
    ResourceValues *vMerged = &pLMerged->getResourceValues(pR->id());
    if (!pR->merge(vForward, vBackward, vMerged, pLMerged)) {
      labelPool_.releaseLastLabel();
      number_of_infeasible_deleted_labels_per_resource_[pR->id()]++;
      return false;
    }
  }
  // keep only the labels with negative costs
  if (pLMerged->cost() > bestPrimalBound_ -param_.epsilon_) {
    labelPool_.releaseLastLabel();
    return false;
  } else {
    bestPrimalBound_ = pLMerged->cost();
  }
  return true;
}

std::vector<PRCLabel> RCSPPSolver::backwardLabelSetting(
    const std::vector<PRCNode> &sortedNodes,
    const std::vector<PRCNode> &mergingNodes, int finalDay) {
  // vector of non-dominated labels at the nodes of day=finalDay
  vector<PRCLabel> destinationLabels;

  // add zero labels to the sink nodes and to the expansion list
  for (const auto &pN : pRcGraph_->pSinks()) {
    auto pL = std::make_shared<RCLabel>(pRcGraph_->nResources());
    pL->setNode(pN);
    pLabelsToExpandPerNode_[pN->id].push_back(pL);
  }

  // Allow to visit each node once
  auto itN = sortedNodes.rbegin();
  for (; itN != sortedNodes.rend(); ++itN) {
    // check if enough time left
    if (!enoughTimeLeft())
      break;

    PRCNode pN = *itN;
    if (pN->dayId < finalDay) break;

    // Expand all the labels of the predecessors of the node through the
    // arcs entering the node
    pullLabelsFromSuccessors(pN);

    // information of domination of newly generated labels
    size_t nLabels = labelPool_.nLabels();
    vector<DominationStatus> domStatusNew(nLabels, NOT_DOMINATED);

    // information of domination of former labels that have not been
    // expanded yet
    std::vector<PRCLabel> &labelsNoExpand = pLabelsNoExpandPerNode_[pN->id];
    size_t nLabelsNoExpand = labelsNoExpand.size();
    vector<DominationStatus> domStatusNewNoExpand(
        nLabelsNoExpand, NOT_DOMINATED);

    if (param_.verbose_ >= 4)
      std::cout << "RCSSP: process node " << pN->toString()
                << " (active resources=" << pN->activeResources.size()
                << ", expanded labels=" << nLabels
                << ", non dominated=" <<  pLabelsToExpandPerNode_[pN->id].size()
                << ")" << std::endl;

    // if expand all the labels at the node, dominate all labels first
    if (param_.rcsppNbToExpand_ == 0) {
      // Domination of all the labels of the current node to get a
      // pareto-optimal set of labels
      checkAllDominations(labelPool_.begin(),
                          labelPool_.end(),
                          useDominateNoLb_ ? dominateNoLB : dominate,
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

  // Collect all the labels of the final day nodes
  for (const auto &pN : mergingNodes)
    for (const auto &pL : pLabelsToExpandPerNode_.at(pN->id))
      if (pL->getInArc() == nullptr)
        destinationLabels.push_back(pL);

  return destinationLabels;
}

void RCSPPSolver::pullLabelsFromSuccessors(const PRCNode &pN) {
  // expand labels to all predecessors
  for (const auto &pArc : pN->outArcs) {
    if (pArc->forbidden) continue;  // do not use forbidden arcs
    // expand all the labels
    for (const auto &pL : pLabelsToExpandPerNode_.at(pArc->target->id)) {
      // expansion on the current arc
      if (expandBack(pL, pArc, labelPool_.getNewLabel()))
        total_number_of_generated_labels_++;
      else
        number_of_infeasible_deleted_labels_++;
    }
  }
}

bool RCSPPSolver::expandBack(
    const PRCLabel &pLNext, const PRCArc &pArc, const PRCLabel &pLPrevious) {
  pLPrevious->setAsPrevious(pLNext, pArc);

  // Expansion of each resource
  for (auto &e : pArc->expanders) {
    // get child ResourceValues (copy of the parent ResourceValues)
    ResourceValues &v = pLPrevious->getResourceValues(e->indResource);
    if (!e->expandBack(pLPrevious, &v)) {
      labelPool_.releaseLastLabel();
      number_of_infeasible_deleted_labels_per_resource_[e->indResource]++;
      return false;
    }
  }

  pLPrevious->addBaseCost(pArc->baseCost);
  pLPrevious->addCost(-pArc->dualCost);
#ifdef NS_DEBUG
  pLPrevious->addDualCost(pArc->dualCost);
#endif
  return true;
}


void RCSPPSolver::displaySolveInputInfo() {
  if (param_.verbose_ >= 3) {
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

void RCSPPSolver::displaySolveStatistics() const {
  if (param_.verbose_ >= 3) {
    std::cout << "\nLABELS STATISTICS: " << std::endl;
    // We add 1 to nParetoLabels in order to take account the first labels of
    // the source node
    std::cout << " - Total number of not dominated labels: "
                 "" << total_number_of_nondominated_labels_ + 1 << std::endl;
    std::cout << " - Number of infeasible labels deleted: "
                 "" << number_of_infeasible_deleted_labels_ << std::endl;
    std::cout << " - Total number of generated labels: " <<
              total_number_of_generated_labels_ << std::endl;
    std::cout << " - Total number of dominations: " <<
              total_number_of_dominations_ << std::endl << std::endl;
  }
}
