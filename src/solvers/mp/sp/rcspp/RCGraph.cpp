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

#include "solvers/mp/sp/rcspp/RCGraph.h"

#include <list>
#include <utility>

void RCGraph::initializeExpanders() {
  vector<PResource> softResources;
  softResources.reserve(nResources());
  // initialize hard resources first
  for (const auto &pR : pNonEnumResources_) {
    if (!pR->isHard()) {
      softResources.push_back(pR);
      continue;
    }
    for (const auto &pA : pArcs_)
      pR->initialize(*pA->origin->pAShift, pA->stretch, pA);
  }

  // then, initialize soft resources
  for (const auto &pR : softResources)
    for (const auto &pA : pArcs_)
      pR->initialize(*pA->origin->pAShift, pA->stretch, pA);
}

void RCGraph::initializeDominance() {
  for (const auto &pN : pNodes_) {
    pN->activeResources.clear();
    for (const auto &pR : pNonEnumResources_)
      if (pR->isActive(pN->day, *pN->pAShift))
        pN->activeResources.push_back(pR);
  }
}

std::vector<PRCNode> RCGraph::sortNodes() const {
  // compute the depth of each node from the sinks
  vector<int> depths(nNodes(), 0);
  std::list<PRCNode> nodesToProcess(pSinks_.begin(), pSinks_.end());
  int maxDepth = -1;
  while (!nodesToProcess.empty()) {
    PRCNode pN = nodesToProcess.front();  // retrieve first node
    nodesToProcess.pop_front();  // remove first node
    // update predecessors
    int depth = depths[pN->id] + 1;
    for (const PRCArc &pArc : pN->inArcs) {
      // if depth is updated -> reprocess node to update predecessors
      if (depths[pArc->origin->id] < depth) {
        depths[pArc->origin->id] = depth;
        nodesToProcess.push_back(pArc->origin);
      }
    }
    // update max depth
    if (depth > maxDepth)
      maxDepth = depth;

    if (depth > nNodes())
      Tools::throwError("The RC graph is cyclic.");
  }

  // order the nodes starting from the deepest nodee
  std::vector<PRCNode> sortedNodes;
  int depth = maxDepth;
  for (; depth >= 0; depth--) {
    // find all nodes of the right depth
    for (int i = 0; i < depths.size(); i++)
      if (depths[i] == depth)
        sortedNodes.push_back(pNode(i));
  }
  return sortedNodes;
}

void RCGraph::addResource(const PResource& pR) {
  pR->setId(pResources_.size());
  pResources_.push_back(pR);
  pNonEnumResources_.push_back(pR);
}

PRCNode RCGraph::addSingleNode(
    NodeType type, int day, const PAbstractShift& pAS) {
  pNodes_.emplace_back(
      std::make_shared<RCNode>(pNodes_.size(), type, day, pAS));
  if (type == SINK_NODE)  // if a sink
    pSinks_.push_back(pNodes_.back());
  else if (type == SOURCE_NODE)  // if a source
    pSources_.push_back(pNodes_.back());
  return pNodes_.back();
}

PRCArc RCGraph::addSingleArc(PRCNode o,
                             PRCNode d,
                             const Stretch &s,
                             double cost) {
  pArcs_.emplace_back(std::make_shared<RCArc>(
      pArcs_.size(), std::move(o), std::move(d), s, cost));
  // Adding the id of the last arc created in the vector of the incident
  // arcs ids of the target node
  PRCArc pArc = pArcs_.back();
  pArc->target->inArcs.push_back(pArcs_.back());
  pArc->origin->outArcs.push_back(pArcs_.back());
  // add the arc to pArcsPerDayShift_ for each day/shift of the stretch
  const Stretch &st = pArc->stretch;
  for (int k = st.firstDay(); k <= st.lastDay(); k++)
    pArcsPerDayShift_[k][st.pShift(k)->id].push_back(pArc);
  return pArc;
}

void RCGraph::printSummaryOfGraph() const {
  std::cout <<  "# SUMMARY OF THE CONFLICT GRAPH" << std::endl;
  printAllNodes();
  printAllArcs();
}

void RCGraph::printAllNodes() const {
  std::cout << "List of nodes: " << std::endl;
  for (const auto& pN : pNodes_) {
    std::cout << "\t";
    pN->print();
    std::cout << std::endl;
  }
}
void RCGraph::printAllArcs() const {
  std::cout << "List of arcs: " << std::endl;
  for (const auto& pA : pArcs_) {
    std::cout << "\t";
    pA->print();
    std::cout << std::endl;
  }
}

const vector<PResource> &RCGraph::pResources() const {
  return pResources_;
}

const vector<PResource> &RCGraph::pNonEnumResources() const {
  return pNonEnumResources_;
}

void RCGraph::pNonEnumResources(const vector<PResource> &pNonEnumResources) {
  pNonEnumResources_ = pNonEnumResources;
}

PRCArc RCGraph::getArc(const PRCNode& origin, const PRCNode& target) const {
  for (auto inArc : target->inArcs) {
    if (inArc->origin->id == origin->id)
      return inArc;
  }
  Tools::throwException("Arc %d not found.", origin->id);
  return nullptr;
}

// Forbid a node / arc
void RCGraph::forbidDayShift(int k, int s) {
  for (const PRCArc &pA : pArcsPerDayShift_[k][s])
    forbidArc(pA);
}

// Authorize a node / arc
void RCGraph::authorizeDayShift(int k, int s) {
  for (const PRCArc &pA : pArcsPerDayShift_[k][s])
    authorizeArc(pA);
}

void RCGraph::forbidArc(const PRCArc &pA) {
  pA->forbidden = true;
  pForbiddenArcs_.insert(pA);
}

// Authorize a node / arc
void RCGraph::authorizeArc(const PRCArc &pA) {
  pA->forbidden = false;
  pForbiddenArcs_.erase(pA);
}

void RCGraph::resetAuthorizationsArcs() {
  for (const auto &pArc : pForbiddenArcs_)
    pArc->forbidden = false;
  pForbiddenArcs_.clear();
}
