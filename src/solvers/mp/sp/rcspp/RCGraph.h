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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RCGRAPH_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RCGRAPH_H_

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "solvers/mp/sp/rcspp/RCLabel.h"

using std::vector;
using std::shared_ptr;
using std::unique_ptr;

/**
 * Description of a node of the RCGraph
 *
 */

// Different node types and their names
enum NodeType {
  SOURCE_NODE, PRINCIPAL_NETWORK, SINK_NODE, NONE_NODE
};

static const std::vector<std::string> nodeTypeName = {
    "SOURCE_NODE", "PPL_NETWORK", "SINK_NODE  ", "NONE       "};

// Different arc types and their names
enum ArcType {
  SOURCE_TO_PRINCIPAL,
  SHIFT_TO_NEWSHIFT,
  SHIFT_TO_SAMESHIFT,
  TO_SINK,
  NONE_ARC
};

static const std::vector<std::string> arcTypeName = {
    "SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH",
    "TO_SINK        ", "NONE           "};

struct RCArc;
typedef shared_ptr<RCArc> PRCArc;

struct RCNode {
  RCNode(int id, NodeType t, int d, PAbstractShift pAS) :
      id(id), type(t), day(d), pAShift(pAS) {}

  bool isWorkNode() {return pAShift->isWork();}

  void print() const {
    std::cout << toString();
  }

  std::string toString() const {
    std::stringstream rep;
    rep << "Id: " << id << ", type: " << type << ", day: " <<
        day << ", shift: " << pAShift->name;
    return rep.str();
  }

  const int id;
  NodeType type;
  const int day;
  PAbstractShift pAShift;
  vector<PRCArc> inArcs;  // pointers to the arcs entering in this node
  vector<PRCArc> outArcs;  // pointers to the arcs exiting in this node
  vector<PResource> activeResources;
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
                 PRCNode  o,
                 PRCNode  t,
                 Stretch s,
                 double c,
                 ArcType type = NONE_ARC) :
      id(n),
      origin(std::move(o)),
      target(std::move(t)),
      stretch(std::move(s)),
      type(type),
      cost(c),
      baseCost(c),
      dualCost(0),
      forbidden(false) {}

  std::string toString() const {
    std::stringstream rep;
    rep << id <<" : "
        << "(" << origin->id << "," << target->id
        << "): base cost = " << baseCost
        << ", dual cost = " << dualCost
        << ", first day = " << stretch.firstDay()
        << "(" << Tools::intToDay(stretch.firstDay()) << ")"
        << ", origin shift = " << origin->pAShift->name
        << ", stretch =";
    for (const auto& pS : stretch.pShifts())
      rep << " " << pS->name;
    return rep.str();
  }

  void print() const {
    std::cout << toString();
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

class RCGraph {
 public:
  explicit RCGraph(int nDays, int nShifts) :
      nDays_(nDays), nShifts_(nShifts) {
    Tools::initVector3D(&pArcsPerDayShift_, nDays, nShifts, 0);
  }

  void copy(const RCGraph &g) {
    nDays_ = g.nDays_;
    nShifts_ = g.nShifts_;
    for (const auto &pN : g.pNodes_) {
      pNodes_.emplace_back(std::make_shared<RCNode>(*pN));
      if (pN->type == SOURCE_NODE) {
        pSources_.push_back(pNodes_.back());
      }
      if (pN->type == SINK_NODE) {
        pSinks_.push_back(pNodes_.back());
      }
    }
    for (const auto &pA : g.pArcs_) {
      pArcs_.emplace_back(std::make_shared<RCArc>(*pA));
    }
    Tools::initVector3D(&pArcsPerDayShift_, nDays_, nShifts_, 0);
    for (int k=0; k < nDays_; k++)
      for (int s=0; s < nShifts_; s++)
        for (const PRCArc &pA : g.pArcsPerDayShift_[k][s])
          pArcsPerDayShift_[k][s].push_back(pArc(pA->id));

    for (const PRCNode &pN : g.pForbiddenNodes_)
      pForbiddenNodes_.insert(pNode(pN->id));
    for (const PRCArc &pA : g.pForbiddenArcs_)
      pForbiddenArcs_.insert(pArc(pA->id));
  }

  // Getters of private attributes
  int nNodes() const {return pNodes_.size();}
  const vector<PRCNode>& pNodes() const {return pNodes_;}
  const PRCNode& pNode(int id) const {return pNodes_[id];}
  const vector<PRCNode> &pSources() const { return pSources_; }
  PRCNode pSource(int d) const {return pSources_[d];}
  const vector<PRCNode> &pSinks() const { return pSinks_; }
  PRCNode pSink(int d) const {return pSinks_[d];}
  int nArcs() const {return pArcs_.size();}
  int nDays() const {return nDays_;}
  int nShifts() const {return nShifts_;}
  const vector<PRCArc> pArcs() const {return pArcs_;}
  const PRCArc& pArc(int id) const {return pArcs_[id];}
  const vector<PRCArc>& pArcs(int day, int shift) const {
    return pArcsPerDayShift_[day][shift];
  }
  int nResources() const {return pResources_.size();}
  const vector<PResource> &pResources() const;
  const vector<PResource> &pNonEnumResources() const;
  void pNonEnumResources(const vector<PResource> &pNonEnumResources);
  PRCArc getArc(const PRCNode& origin, const PRCNode& target) const;
  const std::set<PRCArc>& pForbiddenArcs() const { return pForbiddenArcs_; }

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
    pNonEnumResources_.clear();
  }

  void addResource(const PResource& pR);
  PRCNode addSingleNode(NodeType type, int day, const PAbstractShift& pAS);

  PRCArc addSingleArc(PRCNode o, PRCNode d, const Stretch &s, double d1);

  // preprocess the expansion of every resource on every arc to save cpu at
  // the actual expansion of labels
  void initializeExpanders();

  // preprocess dominance of every resource on every node to save cpu at the
  // actual dominance of labels
  void initializeDominance();

  // return a sorted vector of the nodes from sources to sinks to ensure
  // that a given node cannot reach previous nodes in the graph
  std::vector<PRCNode> sortNodes() const;

  // Forbid all node / arc using day k and shift s
  void forbidDayShift(int k, int s);

  // Authorize a node / arc using day k and shift s
  void authorizeDayShift(int k, int s);

  // Forbid arc
  void forbidArc(const PRCArc &pA);

  // Authorize a arc
  void authorizeArc(const PRCArc &pA);

  // remove all forbidden arcs
  void resetAuthorizationsArcs();

  // print utilities
  void printSummaryOfGraph() const;
  void printAllNodes() const;
  void printAllArcs() const;

 protected:
  // THE GRAPH
  int nDays_, nShifts_;
  vector<PRCNode> pNodes_;
  // source and sink nodes
  vector<PRCNode> pSources_, pSinks_;

  vector<PRCArc> pArcs_;
  // arcs stored by day and shifts contains within the arc stretch
  vector3D<PRCArc> pArcsPerDayShift_;

 protected:
  std::set<PRCNode> pForbiddenNodes_;
  std::set<PRCArc> pForbiddenArcs_;

  // RESOURCES
  // Non enum resources contains all the resources that has not been enumerated
  // in the RC graph when enable
  vector<PResource> pResources_, pNonEnumResources_;
};

typedef shared_ptr<RCGraph> PRCGraph;

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCGRAPH_H_
