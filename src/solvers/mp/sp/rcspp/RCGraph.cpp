//
// Created by antoine legrain on 2020-04-10.
//

#include "RCGraph.h"
#include <iostream>


// print a solution
std::string RCSolution::toString(std::vector<int> shiftIDToShiftTypeID) const {
  std::stringstream buff;
  buff << "RC solution of cost " << cost
       << " starting on day " << firstDay << ":" << std::endl;
  for(int k=0; k<firstDay; k++) buff << "|     ";

  char buffer[500];
  if(shiftIDToShiftTypeID.empty())
    for(int s: shifts) {
      std::sprintf(buffer,  "| -:%2d", s);
      buff << buffer;
    }
  else
    for(int s: shifts) {
      std::sprintf(buffer, "|%2d:%2d", shiftIDToShiftTypeID[s], s);
      buff << buffer;
    }
  buff << "|" << std::endl;
  return buff.str();
}

/*
 * Main definitions of the RCGraph
 */

RCGraph::RCGraph(int nDays): nDays_(nDays), nNodes_(0), nArcs_(0) {}
RCGraph::~RCGraph() {}

std::vector<RCSolution> RCGraph::solve(RCSPPSolver *rcspp, int nLabels,
                                       const std::vector<int>& labelsMinLevel,
                                       std::vector<vertex> sinks) {
  // 0 - Remove all fordidden edges
  //
  std::map<int, Arc_Properties> arcs_to_remove;
  std::vector<edge> edges_to_remove;
  for (int a = 0; a < nArcs_; ++a) {
    const Arc_Properties &arc_prop = arc(a);
    if (arc_prop.forbidden) {
      arcs_to_remove[a] = arc_prop;
      edges_to_remove.push_back(arcsDescriptors_[a]);
    }
  }
  for (auto &e: edges_to_remove)
    boost::remove_edge(e, g_);

  // 1 - solve the resource constraints shortest path problem
  //
  if (sinks.empty()) sinks = sinks_;
  std::vector<RCSolution> rc_solutions = rcspp->solve(nLabels, labelsMinLevel, sinks);

  // 2 - Add back all fordidden edges
  //
  for (auto &p: arcs_to_remove) {
    edge e = (add_edge(p.second.origin, p.second.destination, p.second, g_)).first;
    arcsDescriptors_[p.first] = e;
  }

  return rc_solutions;
}

// Addition of a single node
int RCGraph::addSingleNode(NodeType type, std::vector<int> lbs, std::vector<int> ubs, bool hard_lbs){
  add_vertex( Vertex_Properties( nNodes_, type, lbs, ubs, hard_lbs ), g_ );
  return nNodes_++;
}

// Adds a single arc (origin, destination, cost, travel time, type)
int RCGraph::addSingleArc(int o, int d, double baseCost, std::vector<int> consumptions,
    ArcType type, int day, std::vector<int> shifts){
  edge e = (add_edge( o, d, Arc_Properties( nArcs_, o, d, type, baseCost, consumptions, day, shifts ), g_ )).first;
  arcsDescriptors_.push_back(e);
  return nArcs_++;
}

void RCGraph::resetAuthorizations() {
  for(int v: forbiddenNodes_)
    boost::put( &Vertex_Properties::forbidden, g_, v, false);
  forbiddenNodes_.clear();
  for(int a: forbiddenArcs_)
    boost::put( &Arc_Properties::forbidden, g_, arcsDescriptors_[a], false);
  forbiddenArcs_.clear();
}

// Print the rcspp
void RCGraph::printGraph(int nLabel, int nShiftsToDisplay) const {

  // TITLE
  std::cout << "# " << std::endl;
  std::cout << "# GRAPH OF THE SUBPROBLEM " << std::endl;
  std::cout << "# " << std::endl;

  // THE NODES
  //
  printAllNodes(nLabel);

  // THE ARCS
  //
  printAllArcs(nLabel, nShiftsToDisplay);

  // SUMMARY
  //
  std::cout << printSummaryOfGraph();



}

// Prints the line of a node
std::string RCGraph::printNode(int v, int nLabel) const {
  std::stringstream rep;
  const Vertex_Properties& vert_prop = node(v);
  char buff[255];
  sprintf(buff, "# NODE   %5d  %15s", v, nodeTypeName[vert_prop.type].c_str());
  rep << buff;

  if(nLabel==-1) nLabel = vert_prop.size();
  for(int l=0; l<nLabel; ++l) {
    sprintf(buff, "  %13s=%s%3d,%3d]", labelName[l].c_str(), (vert_prop.hard_lbs_ ? "[": "("),
        vert_prop.lb(l), vert_prop.ub(l));
    rep << buff;
  }
  rep << "  " << shortNameNode(v);
  return rep.str();
}

// Prints all nodes
void RCGraph::printAllNodes(int nLabel) const {
  std::cout << "#   NODES (" << nNodes_ << ")" << std::endl;
  for(int v=0; v<nNodes_; v++) std::cout << printNode(v, nLabel) << std::endl;
  std::cout << "# " << std::endl;

}
// Prints the line of an arc

std::string RCGraph::printArc(int a, int nLabel, int nShiftsToDisplay) const {
  std::stringstream rep;
  const Arc_Properties& arc_prop = arc(a);
  char buff[255];
  sprintf(buff, "# ARC   %5d  %15s  (%5d,%5d)  c=%12.2f  ", a, arcTypeName[arc_prop.type].c_str(),
      arcOrigin(a),  arcDestination(a), arc_prop.cost);
  rep << buff;
  sprintf(buff, "%10s  %2d:", (arc_prop.forbidden ? "forbidden" : "authorized"), arc_prop.day);
  rep << buff;
  for(unsigned int s=0; s<nShiftsToDisplay; ++s){
    if(s < arc_prop.shifts.size())
      sprintf(buff, " %3d", arc_prop.shifts[s]);
    else sprintf(buff, " %3s", "");
    rep << buff;
  }
  rep << (nShiftsToDisplay < (int) arc_prop.shifts.size() ? " ..." : "    ");
  sprintf(buff, "  [%20s] -> [%20s]", shortNameNode(arcOrigin(a)).c_str(),
      shortNameNode(arcDestination(a)).c_str());
  rep << buff;

  if(nLabel==-1) nLabel = arc_prop.size();
  for(int l=0; l<nLabel; ++l) {
    if(l<arc_prop.size())
      sprintf(buff, "  %13s=%3d", labelName[l].c_str(), arc_prop.consumptions[l]);
    else sprintf(buff, "  %13s=%3s", labelName[l].c_str(), "");
    rep << buff;
  }
  return rep.str();

}

// Prints all arcs
void RCGraph::printAllArcs(int nLabel, int nShiftsToDisplay) const {
  std::cout << "#   ARCS (" << nArcs_ << "]" << std::endl;
  for(int a=0; a<nArcs_; a++) std::cout << printArc(a, nLabel, nShiftsToDisplay) << std::endl;
  std::cout << "# " << std::endl;
}

// Short name for a node
std::string RCGraph::shortNameNode(int v) const {

  std::stringstream rep;
  NodeType type_v = get( &Vertex_Properties::type, g_)[v];

  if(type_v == SOURCE_NODE){
    rep << "SOURCE";
  }

  else if (type_v == PRINCIPAL_NETWORK){
//    int k = get( &Vertex_Properties::day, g_)[v];
//    int cons = principalToCons_.at(v);
//    rep << (pScenario_->intToShiftType_[principalToShift_.at(v)])[0] << "-" << k << "-" << cons;
    rep << "PRINCIPAL_NETWORK";
  }

  else if (type_v == PRICE_LABEL_ENTRANCE){
    rep << "LEN_IN";
  }

  else if (type_v == PRICE_LABEL){
    rep << "LEN_<=";
    for(int ub: nodeUBs(v)) rep << " " << ub;
  }

  else if (type_v == PRICE_LABEL_EXIT){
    rep << "ROTSIZE OUT";
  }

  else if (type_v == SINK_NODE){
    rep << "SINK";
  }

  else{
    rep << "NONE";
  }

  return rep.str();
}

// Summary of the rcspp
std::string RCGraph::printSummaryOfGraph() const {
  std::stringstream rep;
  std::map<NodeType,int> nNodesPerType;
  std::map<ArcType,int> nArcsPerType;
  rep << "# +------------------+" << std::endl;
  rep << "# | SUBPROBLEM GRAPH |" << std::endl;
  rep << "# +------------------+" << std::endl;
  rep << "# " << std::endl;
  rep << "#     " << nDays_ << " days" << std::endl;
  rep << "# " << std::endl;
  // COUNT THE NODES
  for(int t = SOURCE_NODE; t!=NONE_NODE; t++){
    NodeType ty = static_cast<NodeType>(t);
    nNodesPerType.insert(std::pair<NodeType,int>(ty,0));
  }
  for(int v=0; v<nNodes_; v++) nNodesPerType.at(nodeType(v))++;
  // DISPLAY NODES
  rep << "#     -------------------------" << std::endl;
  rep << "#   > NODES                    " << std::endl;
  rep << "#     -------------------------" << std::endl;
  for(int t = SOURCE_NODE; t!=NONE_NODE; t++){
    NodeType ty = static_cast<NodeType>(t);
    rep << "#        " << nodeTypeName[ty] << "      " <<nNodesPerType.at(ty) << std::endl;
  }
  rep << "#     -------------------------" << std::endl;
  rep << "#        TOTAL            " << nNodes_ << std::endl;
  rep << "#     -------------------------" << std::endl;
  rep << "# " << std::endl;
  rep << "# " << std::endl;

  // COUNT THE ARCS
  for(int a=0; a<nArcs_; a++) nArcsPerType[arcType(a)]++;
  // DISPLAY ARCS
  rep << "#     -------------------------" << std::endl;
  rep << "#   > ARCS                     " << std::endl;
  rep << "#     -------------------------" << std::endl;
  for(auto p: nArcsPerType){
    rep << "#        " << arcTypeName[p.first] << "  " << p.second << std::endl;
  }
  rep << "#     -------------------------" << std::endl;
  rep << "#        TOTAL            " << nArcs_ << std::endl;
  rep << "#     -------------------------" << std::endl;

  return rep.str();
}