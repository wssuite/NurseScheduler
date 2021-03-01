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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_MYRCGRAPH_H_
#define SRC_SOLVERS_MP_SP_RCSPP_MYRCGRAPH_H_

#include <utility>
#include <set>
#include <memory>
#include <vector>

#include "solvers/mp/sp/rcspp/MyRCLabel.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

#include <boost/graph/adjacency_list.hpp>

using std::vector;
using std::shared_ptr;
using std::unique_ptr;

/**
 * Description of a node of the RCGraph
 *
 */

struct RCArc;
typedef shared_ptr<RCArc> PRCArc;

struct RCNode {
  RCNode(int id, NodeType t, int d, PShift pS) : id(id), type(t), day(d),
                                                 pShift(pS) {
  }

  bool isWorkNode() {return pShift->isWork();}

  void print() const {
    std::cout << "Id: " << id << ", type: " << type << ", day: " <<
    day << ", shift: " << pShift->name;
  }

  const int id;
  NodeType type;
  const int day;
  PShift pShift;
  vector<PRCArc> inArcs;  // Ids of the arcs entering in this node
  vector<int> indActiveResources;
};

typedef shared_ptr<RCNode> PRCNode;

/**
 * Description of an arc of the RCGraph
 *
 */
struct RCArc {
  // Constructor
  //
  explicit RCArc(int n,
                 const PRCNode o,
                 const PRCNode t,
                 Stretch s,
                 double c,
                 ArcType type = NONE_ARC) :
      id(n),
      origin(o),
      target(t),
      stretch(std::move(s)),
      type(type),
      cost(c),
      baseCost(c),
      dualCost(0),
      forbidden(false) {}

  void print() const {
    std::cout << id <<" : "
              << "(" << origin->id << "," << target->id
              << "): base cost = " << baseCost
              << ", dual cost = " << dualCost
              << ", first day = " << Tools::intToDay(stretch.firstDay())
              << ", origin shift = " << origin->pShift->name
              << ", stretch =";
    for (auto pS : stretch.pShifts()) {
      std::cout << " " << pS->name;
    }
  }

  // Set the dual costs on the arc, this should be done when solving the RCSPP
  void resetDualCost() {
    cost = baseCost;
    dualCost = 0;
  }

  void addDualCost(double c) {
    dualCost += c;
    cost += c;
  }

  void addBaseCost(double c) {
    baseCost += c;
    cost += c;
  }

  // Return true if the arc was added during the enumeration process
  bool isEnumeratedArc() {
    return stretch.nDays() > 1;
  }

  int id;
  const PRCNode origin, target;   // end vertices of the arc

  // the arc symbolizes the affectation of the following stretch
  const Stretch stretch;

  // global information on the type of arc, see definition of enumerated type
  ArcType type;

  // traversal costs of the arc
  double cost;  // total traversal cost
  double baseCost;  // base cost due to soft constraints
  double dualCost;  // cost due to the dual values of the master problem

  // vector of resources expanded along the arc
  vector<PExpander> expanders;

  // true if the traversal of the arc is forbidden (by branching rule or LNS)
  bool forbidden;
};

class MyRCGraph {
 public:
  explicit MyRCGraph(int nDays) : nDays_(nDays), pSource_(nullptr) {}

  // Getters of private attributes
  int nNodes() const {return pNodes_.size();}
  const vector<PRCNode>& pNodes() const {return pNodes_;}
  const PRCNode& pNode(int id) const {return pNodes_[id];}
  PRCNode pSource() const {return pSource_;}
  const vector<PRCNode> &pSinks() const { return pSinks_; }
  int nArcs() const {return pArcs_.size();}
  int nDays() const {return nDays_;}
  const vector<PRCArc> pArcs() const {return pArcs_;}
  const PRCArc& pArc(int id) const {return pArcs_[id];}
  int nResources() const {return pResources_.size();}
  const vector<PResource> &pResources() const;
  PRCArc getArc(PRCNode origin, PRCNode target) const;

  void reset() {
    clearResources();
    clearNodes();
    clearArcs();
  }

  void clearNodes() {
    pNodes_.clear();
  }

  void clearArcs() {
    pArcs_.clear();
  }

  void clearResources() {
    pResources_.clear();
  }

  const vector<int> & getActiveResourcesNode(int idNode) {
    return pNode(idNode)->indActiveResources;
  }

  void addResource(PResource pR);
  PRCNode addSingleNode(NodeType type, int day, PShift pS);

  PRCArc addSingleArc(PRCNode o, PRCNode d, const Stretch &s, double d1);

  // preprocess the expansion of every resource on every arc to save cpu at
  // the actual expansion of labels
  void initializeExpanders();

  // preprocess dominance of every resource on every node to save cpu at the
  // actual dominance of labels
  void initializeDominance();

  // return a sorted vector of the nodes from sources to sinks to ensure
  // that a given node cannont reach previous nodes in the graph
  std::vector<PRCNode> sortNodes() const;

  // print utilities
  void printSummaryOfGraph() const;
  void printAllNodes() const;
  void printAllArcs() const;

 protected:
  // THE GRAPH
  int nDays_;
  vector<PRCNode> pNodes_;
  PRCNode pSource_;  // source node
  vector<PRCNode> pSinks_;   // sink nodes

  vector<PRCArc> pArcs_;

 protected:
  std::set<int> forbiddenNodes_;
  std::set<int> forbiddenArcs_;

  // RESOURCES
  vector<PResource> pResources_;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_MYRCGRAPH_H_
