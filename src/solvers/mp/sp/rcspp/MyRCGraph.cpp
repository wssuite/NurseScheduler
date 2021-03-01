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

#include "solvers/mp/sp/rcspp/MyRCGraph.h"

#include <list>
#include <utility>

void MyRCGraph::initializeExpanders() {
  vector<PResource> softResources;
  softResources.reserve(nResources());
  // initialize hard resources first
  for (const auto &pR : pResources_) {
    if (!pR->isHard()) {
      softResources.push_back(pR);
      continue;
    }
    for (const auto &pA : pArcs_)
      pR->initialize(*pA->origin->pShift, pA->stretch, pA);
  }

  // then, initialize soft resources
  for (auto &pR : softResources)
    for (auto &pA : pArcs_)
      pR->initialize(*pA->origin->pShift, pA->stretch, pA);
}

void MyRCGraph::initializeDominance() {
  for (const auto &pN : pNodes_)
    for (const auto &pR : pResources_)
      if (pR->isActive(pN->day, *pN->pShift))
        pN->indActiveResources.push_back(pR->id());
}

std::vector<PRCNode> MyRCGraph::sortNodes() const {
  // compute the depth of each node from the sinks
  vector<int> depths(nNodes(), 0);
  std::list<PRCNode> nodesToProcess(pSinks_.begin(), pSinks_.end());
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
    if (depth > nNodes())
      Tools::throwError("The RC graph is cyclic.");
  }

  // order the nodes starting from the source
  std::vector<PRCNode> sortedNodes;
  int depth = depths[pSource_->id];
  for (; depth >= 0; depth--) {
    // find all nodes of the right depth
    for (int i = 0; i < depths.size(); i++)
      if (depths[i] == depth)
        sortedNodes.push_back(pNode(i));
  }
  return sortedNodes;
}

void MyRCGraph::addResource(PResource pR) {
  pR->setId(pResources_.size());
  pResources_.push_back(pR);
}

PRCNode MyRCGraph::addSingleNode(NodeType type, int day, PShift pS) {
  pNodes_.emplace_back(std::make_shared<RCNode>(pNodes_.size(), type, day, pS));

  if (type == SINK_NODE) {
    pSinks_.push_back(pNodes_.back());
  } else if (type == SOURCE_NODE) {
    // TODO(AL): use a vector of sources instead of one source
    if (pSource_ != nullptr) Tools::throwError("A source is already defined");
    pSource_ = pNodes_.back();
  }

  return pNodes_.back();
}

PRCArc MyRCGraph::addSingleArc(PRCNode o,
                             PRCNode d,
                             const Stretch &s,
                             double cost) {
  pArcs_.emplace_back(std::make_shared<RCArc>(
      pArcs_.size(), std::move(o), std::move(d), s, cost));
  // Adding the id of the last arc created in the vector of the incident
  // arcs ids of the target node
  pArcs_.back()->target->inArcs.push_back(pArcs_.back());
  return pArcs_.back();
}

void MyRCGraph::printSummaryOfGraph() const {
  std::cout <<  "# SUMMARY OF THE CONFLICT GRAPH" << std::endl;
  printAllNodes();
  printAllArcs();
}

void MyRCGraph::printAllNodes() const {
  std::cout << "List of nodes: " << std::endl;
  for (auto pN : pNodes_) {
    std::cout << "\t";
    pN->print();
    std::cout << std::endl;
  }
}
void MyRCGraph::printAllArcs() const {
  std::cout << "List of arcs: " << std::endl;
  for (auto pA : pArcs_) {
    std::cout << "\t";
    pA->print();
    std::cout << std::endl;
  }
}
const vector<PResource> &MyRCGraph::pResources() const {
  return pResources_;
}

PRCArc MyRCGraph::getArc(PRCNode origin, PRCNode target) const {
  for (auto inArc : target->inArcs) {
    if (inArc->origin->id == origin->id)
      return inArc;
  }
  Tools::throwException("Arc %d not found.", origin->id);
  return nullptr;
}
