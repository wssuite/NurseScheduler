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

void RCGraph::initializeExpanders() {
  // initialize hard resources first
  int ind = -1;
  for (const auto &pR : pResources_) {
    ind++;
    if (!pR->isHard()) continue;
    for (const auto &pA : pArcs_) {
      pR->initialize(*pA->origin->pAShift, pA->stretch, pA, ind);
    }
  }

  // then, initialize soft resources
  ind = -1;
  for (const auto &pR : pResources_) {
    ind++;
    if (pR->isHard()) continue;
    for (const auto &pA : pArcs_) {
      pR->initialize(*pA->origin->pAShift, pA->stretch, pA, ind);
    }
  }
}

void RCGraph::initializeDominance() {
  for (const auto &pN : pNodes_) {
    pN->activeResources.clear();
    for (const auto &pR : pResources_)
      if (pR->isActive(pN->dayId, pN->pAShift))
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
        nodesToProcess.push_back(pNode(pArc->origin->id));
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
}

PRCNode RCGraph::addSingleNode(
    NodeType type, const PDay& pDay, const PAbstractShift& pAS) {
  pNodes_.push_back(
      std::make_shared<RCNode>(pNodes_.size(), type, pDay, pAS));
  if (type == SINK_NODE) {
    // if a sink
    pSinks_.push_back(pNodes_.back());
    pNodesPerDay_[pDay->id - firstDayId_].push_back(pNodes_.back());
  } else if (type == SOURCE_NODE) {
    // if a source
    pSources_.push_back(pNodes_.back());
  } else {
    pNodesPerDay_[pDay->id - firstDayId_].push_back(pNodes_.back());
  }

  return pNodes_.back();
}

PRCArc RCGraph::findArc(const PRCNode& o, const Stretch &s) {
  for (auto pA : o->outArcs) {
    bool isEqual = true;
    if (pA->stretch.nDays() == s.nDays()) {
      auto itS = s.pShifts().begin();
      for (const auto &pS : pA->stretch.pShifts()) {
        if (!pS->equals(**itS)) {
          isEqual = false;
          break;
        }
        itS++;
      }
      if (isEqual) return pA;
    }
  }
  return nullptr;
}

PRCArc RCGraph::addSingleArc(const PRCNode& o,
                             const PRCNode& d,
                             const Stretch &s,
                             double cost) {
#ifdef DBG
  if (o->id == d->id) Tools::throwError(
        "It is forbidden to create loop arc for an acyclic graph. "
        "Creating an arc with the same origin and destination: %s",
        o->toString().c_str());
#endif

  // we first check that the arc is not already present in the graph, there
  // will be errors if the same arc appears twice in the graph
  PRCArc pArc = findArc(o, s);
  if (pArc != nullptr) return pArc;

  pArcs_.push_back(std::make_shared<RCArc>(
      pArcs_.size(), o, d, s, cost));
  // Adding the id of the last arc created in the vector of the incident
  // arcs ids of the target node
  pArc = pArcs_.back();
  pArc->target->inArcs.push_back(pArc);
  pArc->origin->outArcs.push_back(pArc);
  // add the arc to pArcsPerDayShift_ for each day/shift of the stretch
  const Stretch &st = pArc->stretch;
  for (int k = st.firstDayId(); k <= st.lastDayId(); k++)
    pArcsPerDayShift_[k][st.pShift(k)->id].push_back(pArc);
  return pArc;
}

// delete an arc from every vector where it appears
void RCGraph::deleteArc(const PRCArc& pA) {
  RCNode *pO = pA->origin, *pT = pA->target;
  Tools::erase(&pArcs_, pA, true);
  Tools::erase(&pO->outArcs, pA, true);
  Tools::erase(&pT->inArcs, pA, true);

  const Stretch &st = pA->stretch;
  for (int k = st.firstDayId(); k <= st.lastDayId(); k++)
    Tools::erase(&pArcsPerDayShift_[k][st.pShift(k)->id],
                 pA, true);
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
