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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCGRAPH_H_
#define SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCGRAPH_H_

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include "boost/config.hpp"

#include "tools/Tools.h"
#include "data/Shift.h"

namespace boostRCSPP {

enum LABEL {
  CONS_DAYS = 0,
  DAYS = 1,
  WEEKEND = 2
};

// true if the dominance is done with a descending order (lower the better),
// false otherwise
static const std::vector<std::string> labelsName = {
    "CONS_DAYS", "DAYS     ", "WEEKEND  "};

// Penalties structure
class Penalties {
 public:
  Penalties() = default;

  Penalties(std::vector<int> labelsMinLevel,
            std::vector<int> labelsMaxLevel,
            std::vector<double> labelsWeight) :
      labelsMinLevel_(std::move(labelsMinLevel)),
      labelsMaxLevel_(std::move(labelsMaxLevel)),
      labelsWeight_(std::move(labelsWeight)) {}

  void addLabel(LABEL label, const Penalties &penalties) {
    labelsMinLevel_.push_back(penalties.minLevel(label));
    labelsMaxLevel_.push_back(penalties.maxLevel(label));
    labelsWeight_.push_back(penalties.weight(label));
  }

  void addLabel(int minLvl, int maxLvl, double weight) {
    labelsMinLevel_.push_back(minLvl);
    labelsMaxLevel_.push_back(maxLvl);
    labelsWeight_.push_back(weight);
  }

  int minLevel(int label) const {
    return labelsMinLevel_[label];
  }

  void minLevel(int label, int lvl) {
    labelsMinLevel_[label] = lvl;
  }

  int maxLevel(int label) const {
    return labelsMaxLevel_[label];
  }

  void maxLevel(int label, int lvl) {
    labelsMaxLevel_[label] = lvl;
  }

  double weight(int label) const {
    return labelsWeight_[label];
  }

  double penalty(int label, int labelLvl) const {
    if (labelLvl < minLevel(label)) return minPenalty(label, labelLvl);
    if (labelLvl > maxLevel(label)) return maxPenalty(label, labelLvl);
    return 0;
  }

  double minPenalty(int label, int labelLvl) const {
    return std::max(0, minLevel(label) - labelLvl) * weight(label);
  }

  double maxPenalty(int label, int labelLvl) const {
    return std::max(0, labelLvl - maxLevel(label)) * weight(label);
  }

 private:
  // levels for which a label will have a penalty if it is respectively
  // below the min or above the max
  std::vector<int> labelsMinLevel_, labelsMaxLevel_;
  // Weight associated to a given label
  std::vector<double> labelsWeight_;
};

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
  SHIFT_TO_ENDSEQUENCE,
  REPEATSHIFT,
  TO_SINK,
  NONE_ARC
};

static const std::vector<std::string> arcTypeName = {
    "SOURCE_TO_PPL  ", "SHIFT_TO_NEWSH ", "SHIFT_TO_SAMESH", "SHIFT_TO_ENDSEQ",
    "REPEATSHIFT    ", "TO_SINK        ", "NONE           "};

//////////////////////////////////////////////////////////////////////////
//
// The following structures are the properties of the nodes and arcs:
//   For the nodes: id, earliest arrival time, and latest arrival time
//   For the arcs : id, cost and time (time is the resource here)
//   For the rcspp: usual stuff + nodes and and special properties
//
//////////////////////////////////////////////////////////////////////////

// Nodes specific properties for RC
struct Vertex_Properties {
  // Constructor
  explicit Vertex_Properties(int n = 0,
                    NodeType t = NONE_NODE,
                    std::vector<int> lbs = {},
                    std::vector<int> ubs = {},
                    bool hard_lbs = false,
                    bool forbidden = false) :
      num(n),
      type(t),
      lbs(lbs),
      ubs(ubs),
      hard_lbs_(hard_lbs),
      forbidden(forbidden) {}

  // id
  int num;

  // type
  NodeType type;

  // label lower bounds
  std::vector<int> lbs;

  // label upper bounds
  std::vector<int> ubs;

  // hard lb
  bool hard_lbs_;

  // forbidden
  bool forbidden;

  int lb(LABEL l) const {
    return lbs.at(l);
  }

  int ub(LABEL l) const {
    return ubs.at(l);
  }

  int size() const {
    return lbs.size();
  }
};

// Arcs specific properties for RC
struct Arc_Properties {
  // Constructor
  Arc_Properties(int n = 0,
                 int origin = -1,
                 int destination = -1,
                 ArcType ty = NONE_ARC,
                 double c = 0,
                 std::vector<int> consumptions = {},
                 int day = -1,
                 std::vector<PShift> pShifts = {},
                 bool forbidden = false) :
      num(n), origin(origin), destination(destination),
      type(ty), cost(c), initialCost(c), consumptions(consumptions),
      day(day), pShifts(pShifts), forbidden(forbidden) {}

  // id
  int num;

  // VERTICES
  int origin, destination;

  // type
  ArcType type;

  // traversal cost
  double cost;
  double initialCost;

  // label consumption
  std::vector<int> consumptions;

  int day;  // day
  std::vector<PShift> pShifts;  // shifts id
  bool forbidden;

  int consumption(LABEL l) const {
    return consumptions.at(l);
  }

  int size() const {
    return consumptions.size();
  }

  std::set<LABEL> labelsToPrice;
  Penalties penalties;

  bool findLabelToPrice(LABEL l) const {
    return labelsToPrice.find(l) != labelsToPrice.end();
  }
};

// Graph with RC generic structure
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
Vertex_Properties, Arc_Properties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex;
typedef boost::graph_traits<Graph>::edge_descriptor edge;

//---------------------------------------------------------------------------
//
// C l a s s e s   G r a p h s
//
// Contains a interface for subgraphs and the main rcspp
//
//---------------------------------------------------------------------------
class SubGraph {
 public:
  SubGraph() {}
  virtual ~SubGraph() {}

  virtual int entrance(int day = -1) const = 0;
  virtual int exit(int day = -1) const = 0;

  // link two sub graphs together:
  // 1. create an arc from the exit of inSubGraph to the current entrance
  // 2. the entrance method should now returns the entrance of the  inSubGraph
  virtual void linkInSubGraph(SubGraph *inSubGraph, int day = -1) = 0;

  // link two sub graphs together:
  // 1. create an arc from the current exit to the entrance of outSubGraph
  // 2. the exit method should now returns the exit of the outSubGraph
  virtual void linkOutSubGraph(SubGraph *outSubGraph, int day = -1) = 0;
};

class RCGraph {
 public:
  explicit RCGraph(int nDays = 0);

  virtual ~RCGraph();

  // remove every forbidden arc before solving the rcspp; there is a need to
  // add the arcs backs after solving
  std::map<int, Arc_Properties> removeForbiddenArcsFromBoost();
  void restoreForbiddenArcsToBoost(std::map<int, Arc_Properties> arcs_removed);

  int nDays() const { return nDays_; }

  ////////////////////
  //  NODES
  ////////////////////
  Graph &g() { return g_; }

  // Basic function for adding a node
  int addSingleNode(NodeType type,
                    std::vector<int> lbs,
                    std::vector<int> ubs,
                    bool hard_lbs = false);

  void setSource(int v) { source_ = v; }

  int source() const { return source_; }

  void addSink(int v) { sinks_.push_back(v); }

  int sink(int k = 0) const { return sinks_.at(k); }

  int lastSink() const { return sinks_.back(); }

  const std::vector<vertex> &sinks() const { return sinks_; }

  // Get info from the node ID
  int nodesSize() const { return nNodes_; }

  const Vertex_Properties &node(int v) const {
    if (v == -1) Tools::throwError("Cannot retrieve a node for index -1;");
    return get(boost::vertex_bundle, g_)[v];
  }

  NodeType nodeType(int v) const {
    return get(&Vertex_Properties::type, g_)[v];
  }

  const std::vector<int> &nodeLBs(int v) const {
    return get(&Vertex_Properties::lbs, g_)[v];
  }

  const std::vector<int> &nodeUBs(int v) const {
    return get(&Vertex_Properties::ubs, g_)[v];
  }

  bool nodeForbidden(int v) const {
    return forbiddenNodes_.find(v) != forbiddenNodes_.end();
  }

  void updateUBs(int v, const std::vector<int> &ubs) {
    boost::put(&Vertex_Properties::ubs, g_, v, ubs);
  }

  void forbidNode(int v) {
    boost::put(&Vertex_Properties::forbidden, g_, v, true);
    forbiddenNodes_.insert(v);
  }

  void authorizeNode(int v) {
    boost::put(&Vertex_Properties::forbidden, g_, v, false);
    forbiddenNodes_.erase(v);
  }

  ////////////////////
  //  ARCS
  ////////////////////
  // Basic function for adding an arc
  int addSingleArc(int origin,
                   int destination,
                   double baseCost,
                   std::vector<int> consumptions,
                   ArcType type,
                   int day = -1,
                   std::vector<PShift> shifts = {});

  int addPricingArc(int origin,
                    int destination,
                    double baseCost,
                    std::vector<int> consumptions,
                    ArcType type,
                    int day,
                    std::set<LABEL> labelsToPrice,
                    const Penalties &penalties);

  // Get info with the arc ID
  int arcsSize() const {
    return nArcs_;
  }

  const Arc_Properties &arc(int a) const {
    if (a == -1) Tools::throwError("Cannot retrieve an arc for index -1;");
    return get(boost::edge_bundle, g_)[arcsDescriptors_[a]];
  }

  ArcType arcType(int a) const {
    return get(&Arc_Properties::type, g_, arcsDescriptors_[a]);
  }

  int arcOrigin(int a) const {
    return boost::source(arcsDescriptors_[a], g_);
  }

  int arcDestination(int a) const {
    return target(arcsDescriptors_[a], g_);
  }

  const std::vector<int> &arcConsumptions(int a) const {
    return get(&Arc_Properties::consumptions, g_, arcsDescriptors_[a]);
  }

  double arcCost(int a) const {
    return get(&Arc_Properties::cost, g_, arcsDescriptors_[a]);
  }

  double arcInitialCost(int a) const {
    return get(&Arc_Properties::initialCost, g_, arcsDescriptors_[a]);
  }

  const std::vector<PShift> &arcShifts(int a) const {
    return get(&Arc_Properties::pShifts, g_, arcsDescriptors_[a]);
  }

  int arcDay(int a) const {
    return get(&Arc_Properties::day, g_, arcsDescriptors_[a]);
  }

  bool arcForbidden(int a) const {
    return get(&Arc_Properties::forbidden, g_, arcsDescriptors_[a]);
  }

  const std::set<LABEL> &arcLabelsToPrice(int a) {
    return get(&Arc_Properties::labelsToPrice, g_, arcsDescriptors_[a]);
  }

  void updateConsumptions(int a, const std::vector<int> &consumptions) {
    boost::put(&Arc_Properties::consumptions,
               g_,
               arcsDescriptors_[a],
               consumptions);
  }

  void updateShifts(int a, const std::vector<PShift> &pShifts) {
    boost::put(&Arc_Properties::pShifts, g_, arcsDescriptors_[a], pShifts);
  }

  void updateCost(int a, double cost) {
    boost::put(&Arc_Properties::cost, g_, arcsDescriptors_[a], cost);
  }

  void forbidArc(int a) {
    boost::put(&Arc_Properties::forbidden, g_, arcsDescriptors_[a], true);
    forbiddenArcs_.insert(a);
  }

  void authorizeArc(int a) {
    boost::put(&Arc_Properties::forbidden, g_, arcsDescriptors_[a], false);
    forbiddenArcs_.erase(a);
  }

  void resetAuthorizations();

  // Print functions.
  void printGraph(int nLabel = -1, int nShiftsToDisplay = 1) const;

  std::string printNode(int v, int nLabel = -1) const;

  void printAllNodes(int nLabel = -1) const;

  std::string printArc(int a, int nLabel = -1, int nShiftsToDisplay = 1) const;
  std::string printArc(const Arc_Properties &arc_prop,
                       int nLabel = -1,
                       int nShiftsToDisplay = 1) const;

  void printAllArcs(int nLabel = -1, int nShiftsToDisplay = 1) const;

  std::string shortNameNode(int v) const;

  std::string printSummaryOfGraph() const;

 protected:
  // THE GRAPH
  Graph g_;
  int nDays_;
  int nNodes_;  // Total number of nodes in the rcspp
  // Source
  vertex source_;
  // Sink Node
  std::vector<vertex> sinks_;

  int nArcs_;  // Total number of arcs in the rcspp
  std::vector<edge> arcsDescriptors_;

  std::set<int> forbiddenNodes_;
  std::set<int> forbiddenArcs_;
};

}  // namespace boostRCSPP

#endif  // SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCGRAPH_H_
