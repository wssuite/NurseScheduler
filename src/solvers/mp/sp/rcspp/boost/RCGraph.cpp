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

#include "RCGraph.h"

#include <iostream>
#include <map>
#include <utility>

#include "RCSPP.h"

namespace boostRCSPP {

/*
 * Main definitions of the RCGraph
 */

RCGraph::RCGraph(int nDays) : nDays_(nDays), nNodes_(0), nArcs_(0) {}
RCGraph::~RCGraph() {}

std::map<int, Arc_Properties> RCGraph::removeForbiddenArcsFromBoost() {
  // get the boost edges of all forbidden arcs
  std::map<int, Arc_Properties> arcs_removed;
  std::vector<edge> edges_to_remove;
  for (int a = 0; a < nArcs_; ++a) {
    const Arc_Properties &arc_prop = arc(a);
    if (arc_prop.forbidden) {
      edges_to_remove.push_back(arcsDescriptors_[a]);
      arcs_removed[a] = arc_prop;
    }
  }
  // remove the edges
  for (auto &e : edges_to_remove)
    boost::remove_edge(e, g_);

  return arcs_removed;
}

void RCGraph::restoreForbiddenArcsToBoost(
    std::map<int, Arc_Properties> arcs_removed) {
  // restore the boost edges of all forbidden arcs
  for (auto &p : arcs_removed) {
    edge e =
        add_edge(p.second.origin, p.second.destination, p.second, g_).first;
    arcsDescriptors_[p.first] = e;
  }
}

// Addition of a single node
int RCGraph::addSingleNode(NodeType type,
                           std::vector<int> lbs,
                           std::vector<int> ubs,
                           bool hard_lbs) {
  add_vertex(Vertex_Properties(nNodes_, type, lbs, ubs, hard_lbs), g_);
  return nNodes_++;
}

// Adds a single arc (origin, destination, cost, travel time, type)
int RCGraph::addSingleArc(int o,
                          int d,
                          double baseCost,
                          std::vector<int> consumptions,
                          ArcType type,
                          int day,
                          std::vector<PShift> shifts) {
  Arc_Properties a(nArcs_, o, d, type, baseCost, consumptions, day, shifts);
  edge e = add_edge(o, d, a, g_).first;
  arcsDescriptors_.push_back(e);
  return nArcs_++;
}

int RCGraph::addPricingArc(int o,
                           int d,
                           double baseCost,
                           std::vector<int> consumptions,
                           ArcType type,
                           int day,
                           std::set<LABEL> labelsToPrice,
                           const Penalties &penalties) {
  Arc_Properties a(nArcs_, o, d, type, baseCost, consumptions, day, {});
  a.labelsToPrice = labelsToPrice;
  a.penalties = penalties;
  edge e = add_edge(o, d, a, g_).first;
  arcsDescriptors_.push_back(e);
  return nArcs_++;
}

void RCGraph::resetAuthorizations() {
  for (int v : forbiddenNodes_)
    boost::put(&Vertex_Properties::forbidden, g_, v, false);
  forbiddenNodes_.clear();
  for (int a : forbiddenArcs_)
    boost::put(&Arc_Properties::forbidden, g_, arcsDescriptors_[a], false);
  forbiddenArcs_.clear();
}

// Print the rcspp
void RCGraph::printGraph(int nLabel, int nShiftsToDisplay) const {
  // TITLE
  std::cout << "# " << std::endl;
  std::cout << "# GRAPH OF THE SUBPROBLEM " << std::endl;
  std::cout << "# " << std::endl;

  // THE NODES
  printAllNodes(nLabel);

  // THE ARCS
  printAllArcs(nLabel, nShiftsToDisplay);

  // SUMMARY
  std::cout << printSummaryOfGraph();
}

// Prints the line of a node
std::string RCGraph::printNode(int v, int nLabel) const {
  std::stringstream rep;
  const Vertex_Properties &vert_prop = node(v);
  char buff[255];
  snprintf(buff, sizeof(buff), "# NODE   %5d  %15s",
           v, nodeTypeName[vert_prop.type].c_str());
  rep << buff;

  if (nLabel == -1) nLabel = vert_prop.size();
  for (int l = 0; l < nLabel; ++l) {
    snprintf(buff, sizeof(buff),
             "  %13s=%s%3d,%3d]",
             labelsName[l].c_str(),
             (vert_prop.hard_lbs_ ? "[" : "("),
             vert_prop.lbs[l],
             vert_prop.ubs[l]);
    rep << buff;
  }
  rep << "  " << shortNameNode(v);
  return rep.str();
}

// Prints all nodes
void RCGraph::printAllNodes(int nLabel) const {
  std::cout << "#   NODES (" << nNodes_ << ")" << std::endl;
  for (int v = 0; v < nNodes_; v++)
    std::cout << printNode(v, nLabel) << std::endl;
  std::cout << "# " << std::endl;
}

// Prints the line of an arc
std::string RCGraph::printArc(int a, int nLabel, int nShiftsToDisplay) const {
  const Arc_Properties &arc_prop = arc(a);
  return printArc(arc_prop, nLabel, nShiftsToDisplay);
}

std::string RCGraph::printArc(
    const Arc_Properties &arc_prop, int nLabel, int nShiftsToDisplay) const {
  std::stringstream rep;
  char buff[255];
  snprintf(buff, sizeof(buff),
           "# ARC   %5d  %15s  (%5d,%5d)  c=%12.2f  ",
           arc_prop.num,
           arcTypeName[arc_prop.type].c_str(),
           arc_prop.origin,
           arc_prop.destination,
           arc_prop.cost);
  rep << buff;
  snprintf(buff, sizeof(buff),
           "%10s  %2d:",
           (arc_prop.forbidden ? "forbidden" : "authorized"),
           arc_prop.day);
  rep << buff;
  for (unsigned int s = 0; static_cast<int>(s) < nShiftsToDisplay; ++s) {
    if (s < arc_prop.pShifts.size())
      snprintf(buff, sizeof(buff), " %3d", arc_prop.pShifts[s]->id);
    else
      snprintf(buff, sizeof(buff), " %3s", "");
    rep << buff;
  }
  rep << (nShiftsToDisplay < static_cast<int>(arc_prop.pShifts.size()) ? " .."
                                                                        "."
                                                                      : "    ");
  snprintf(buff, sizeof(buff), "  [%20s] -> [%20s]",
           shortNameNode(arc_prop.origin).c_str(),
           shortNameNode(arc_prop.destination).c_str());
  rep << buff;

  if (nLabel == -1) nLabel = arc_prop.size();
  for (int l = 0; l < nLabel; ++l) {
    if (l < arc_prop.size())
      snprintf(buff, sizeof(buff),
               "  %13s=%3d",
               labelsName[l].c_str(),
               arc_prop.consumptions[l]);
    else
      snprintf(buff, sizeof(buff), "  %13s=%3s", labelsName[l].c_str(), "");
    rep << buff;
  }

  rep << "   Price:";
  for (int l : arc_prop.labelsToPrice)
    rep << " " << labelsName[l].c_str();

  return rep.str();
}

// Prints all arcs
void RCGraph::printAllArcs(int nLabel, int nShiftsToDisplay) const {
  std::cout << "#   ARCS (" << nArcs_ << "]" << std::endl;
  for (int a = 0; a < nArcs_; a++)
    std::cout << printArc(a, nLabel, nShiftsToDisplay) << std::endl;
  std::cout << "# " << std::endl;
}

// Short name for a node
std::string RCGraph::shortNameNode(int v) const {
  std::stringstream rep;
  NodeType type_v = get(&Vertex_Properties::type, g_)[v];

  if (type_v == SOURCE_NODE) {
    rep << "SOURCE";
  } else if (type_v == PRINCIPAL_NETWORK) {
    rep << "PRINCIPAL_NETWORK";
  } else if (type_v == SINK_NODE) {
    rep << "SINK";
  } else {
    rep << "NONE";
  }

  return rep.str();
}

// Summary of the rcspp
std::string RCGraph::printSummaryOfGraph() const {
  std::stringstream rep;
  std::map<NodeType, int> nNodesPerType;
  std::map<ArcType, int> nArcsPerType;
  rep << "# +------------------+" << std::endl;
  rep << "# | SUBPROBLEM GRAPH |" << std::endl;
  rep << "# +------------------+" << std::endl;
  rep << "# " << std::endl;
  rep << "#     " << nDays_ << " days" << std::endl;
  rep << "# " << std::endl;
  // COUNT THE NODES
  for (int t = SOURCE_NODE; t != NONE_NODE; t++) {
    NodeType ty = static_cast<NodeType>(t);
    nNodesPerType.insert(std::pair<NodeType, int>(ty, 0));
  }
  for (int v = 0; v < nNodes_; v++) nNodesPerType.at(nodeType(v))++;
  // DISPLAY NODES
  rep << "#     -------------------------" << std::endl;
  rep << "#   > NODES                    " << std::endl;
  rep << "#     -------------------------" << std::endl;
  for (int t = SOURCE_NODE; t != NONE_NODE; t++) {
    NodeType ty = static_cast<NodeType>(t);
    rep << "#        " << nodeTypeName[ty] << "      " << nNodesPerType.at(ty)
        << std::endl;
  }
  rep << "#     -------------------------" << std::endl;
  rep << "#        TOTAL            " << nNodes_ << std::endl;
  rep << "#     -------------------------" << std::endl;
  rep << "# " << std::endl;
  rep << "# " << std::endl;

  // COUNT THE ARCS
  for (int a = 0; a < nArcs_; a++) nArcsPerType[arcType(a)]++;
  // DISPLAY ARCS
  rep << "#     -------------------------" << std::endl;
  rep << "#   > ARCS                     " << std::endl;
  rep << "#     -------------------------" << std::endl;
  for (auto p : nArcsPerType)
    rep << "#        " << arcTypeName[p.first] << "  " << p.second << std::endl;
  rep << "#     -------------------------" << std::endl;
  rep << "#        TOTAL            " << nArcs_ << std::endl;
  rep << "#     -------------------------" << std::endl;

  return rep.str();
}

}  // namespace boostRCSPP
