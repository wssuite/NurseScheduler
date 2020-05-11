//
// Created by antoine legrain on 2020-04-10.
//

#include "solvers/mp/rcspp/RCGraph.h"
#include <iostream>


#ifdef DBG
void spp_res_cont::print() const {
  RCSolution sol(first_day, shifts_, cost);
  std::cout << sol.toString();
  for(int l=0; l<size(); l++) {
    std::cout << labelName[l].c_str() << "=" << label_value(l) << "  ";
  }
  std::cout << std::endl;
  std::cout << "Arc taken: " << pred_arc << std::endl;
}
#endif

// Resources extension model (arc has cost + label consumptions)
    bool ref_spp::operator()( const Graph& g,
                            spp_res_cont& new_cont,
                            const spp_res_cont& old_cont,
                            boost::graph_traits<Graph>::edge_descriptor ed ) const {
  const Arc_Properties &arc_prop = get(boost::edge_bundle, g)[ed];
  const Vertex_Properties &vert_prop = get(boost::vertex_bundle, g)[target(ed, g)];
#ifdef DBG
  if (arc_prop.forbidden)
    return false;
  if (vert_prop.forbidden)
    return false;
#endif
  new_cont.cost = old_cont.cost + arc_prop.cost;
  if(new_cont.first_day == -1) new_cont.first_day = arc_prop.day;
  if(arc_prop.day != -1) new_cont.day = arc_prop.day;

#ifdef DBG
  new_cont.pred_arc = arc_prop.num;
  if(!arc_prop.shifts.empty()) {
    new_cont.shifts_.insert(new_cont.shifts_.end(),
                            arc_prop.shifts.begin(),
                            arc_prop.shifts.end());
  }
#endif

#ifdef DBG
//  if(new_cont.pred_arc == 368 || new_cont.pred_arc == 68 || new_cont.pred_arc == 524)
    assert(old_cont.size() == new_cont.size());
#endif

  for (int l = 0; l < old_cont.size(); ++l) {
    int lv = old_cont.label_value(l) + arc_prop.consumption(l);
    int lb = vert_prop.lb(l);
    if(lv < lb) {
      if(vert_prop.hard_lbs_)
        return false;
      lv = lb;
    }
    if (lv > vert_prop.ub(l))
      return false;
    new_cont.label_values[l] = lv;
  }

  return true;
}

// Dominance function model
bool dominance_spp::operator()( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ) const {
  // must be "<=" here!!!
  // must NOT be "<"!!!
  if (res_cont_1.cost > res_cont_2.cost + EPSILON) return false;

  /* Dominance:
   * if label is increasing -> can dominate after a certain level:
   *    if the label do not reach the minimum level, there is no additional cost
   * if label is decreasing -> cannot dominate (min level should be an ub).
   *    the label could continue to decrease in the future and not imply any additional cost.
   * So, res_cont_1 dominates res_cont_2 in one of these 3 situations:
   * a- cost1 < cost2 at epsilon and all label1 <= label2
   * b- cost1 == cost2 at epsilon and all label1 <= label2 and
   *    it exists one label1 < label2 with label2 > minLevel (significant)
   * c- res_cont_1 == res_cont_2
   */

  bool dominate = (res_cont_1.cost < res_cont_2.cost - EPSILON),
      biggerThanMinLevel = false,
      equal = !dominate;
  for (int l = 0; l < res_cont_1.size(); ++l) {
    int label1 = res_cont_1.label_value(l),
        label2 = res_cont_2.label_value(l),
        minLevel = labelsMinLevel_.at(l);
    // label1 > label2 -> cannot be dominated in any case
    if (label1 > label2) return false;

    // b or not c
    if (!dominate && label1 < label2) {
      equal = false;
      // b: check if label is bigger than the min level
      if (label2 > minLevel)
        biggerThanMinLevel = true;
    }
  }

#ifdef DBG
//  if(dominate || biggerThanMinLevel) {
//    std::cout << "**************** DOMINATED ****************" << std::endl;
//    res_cont_2.print();
//    std::cout << "******************* BY ********************" << std::endl;
//    res_cont_1.print();
//    std::cout << "*******************************************" << std::endl;
//  }
#endif

  // a
  if(dominate) return true;
  // b
  if(biggerThanMinLevel) return true;
  // c
  return equal;

  // this is not a contradiction to the documentation
  // the documentation says:
  // "A label $l_1$ dominates a label $l_2$ if and only if both are resident
  // at the same vertex, and if, for each resource, the resource consumption
  // of $l_1$ is less than or equal to the resource consumption of $l_2$,
  // and if there is at least one resource where $l_1$ has a lower resource
  // consumption than $l_2$."
  // one can think of a new label with a resource consumption equal to that
  // of an old label as being dominated by that old label, because the new
  // one will have a higher number and is created at a later point in time,
  // so one can implicitly use the number or the creation time as a resource
  // for tie-breaking
}


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

// comparator to order the processing of the nodes in boost rc spp
typedef ks_smart_pointer<boost::r_c_shortest_paths_label<Graph, spp_res_cont> > Spplabel;
bool SpplabelComparator::operator()(const Spplabel& splabel1, const Spplabel& splabel2) const {
  if(splabel1->cumulated_resource_consumption.day > splabel2->cumulated_resource_consumption.day)
    return true;
  if(splabel1->cumulated_resource_consumption.day < splabel2->cumulated_resource_consumption.day)
    return false;
  return splabel1->cumulated_resource_consumption.cost > splabel2->cumulated_resource_consumption.cost;
}

/*
 * Main definitions of the RCGraph
 */

RCGraph::RCGraph(int nDays): nDays_(nDays), nNodes_(0), nArcs_(0) {}
RCGraph::~RCGraph() {}

std::vector<RCSolution> RCGraph::solve(int nLabels, double maxReducedCostBound,
                                       const std::vector<int>& labelsMinLevel,
                                       std::vector<vertex> sinks,
                                       std::function<void (spp_res_cont&)> post_process_rc) {
  // 0 - Remove all fordidden edges
  //
  std::map<int,Arc_Properties> arcs_to_remove;
  std::vector<edge> edges_to_remove;
  for(int a=0; a<nArcs_; ++a) {
    const Arc_Properties& arc_prop = arc(a);
    if(arc_prop.forbidden) {
      arcs_to_remove[a] = arc_prop;
      edges_to_remove.push_back(arcsDescriptors_[a]);
    }
  }
  for(auto& e: edges_to_remove)
    boost::remove_edge(e, g_);

  // 1 - solve the resource constraints shortest path problem
  //
  if(sinks.empty()) sinks = sinks_;
  dominance_spp dominance(labelsOrder, labelsMinLevel);
  vector2D<edge> opt_solutions_spp;
  std::vector<spp_res_cont> pareto_opt_rcs_spp;
  SpplabelComparator comp;
  r_c_shortest_paths_dispatch( g_,
      get( &Vertex_Properties::num, g_ ),
      source_,
      sinks,
      opt_solutions_spp,
      pareto_opt_rcs_spp,
      spp_res_cont (0, std::vector<int>(nLabels)),
      ref_spp(),
      dominance,
      std::allocator< boost::r_c_shortest_paths_label< Graph, spp_res_cont> >(),
      boost::default_r_c_shortest_paths_visitor(),
      comp); // std::greater<Spplabel>()

  // process paths if needed (process_solution is not used)
  std::vector<RCSolution> rc_solutions;
  // For each path of the list, record the corresponding rotation (if negativeOnly=true, do it only if the dualCost < 0)
  for(unsigned int p=0; p < opt_solutions_spp.size(); ++p) {
    spp_res_cont &rc = pareto_opt_rcs_spp[p];
    if (processPath(opt_solutions_spp[p], rc, post_process_rc))
      if (rc.cost < maxReducedCostBound)
        rc_solutions.push_back(solution(opt_solutions_spp[p], rc));
  }

  // 2 - Add back all fordidden edges
  //
  for(auto& p: arcs_to_remove) {
    edge e = (add_edge( p.second.origin, p.second.destination, p.second, g_ )).first;
    arcsDescriptors_[p.first] = e;
  }

  return rc_solutions;
}

bool RCGraph::processPath(std::vector<edge>& path, spp_res_cont& rc,
    std::function<void (spp_res_cont&)> post_process_rc) const {
  // a. Check if it is valid
  bool b_is_a_path_at_all = false;
  bool b_feasible = false;
  bool b_correctly_extended = false;
  spp_res_cont actual_final_resource_levels(0, std::vector<int>(rc.size()));
  boost::graph_traits<Graph>::edge_descriptor ed_last_extended_arc;
  check_r_c_path(g_,
                 path,
                 spp_res_cont(0, std::vector<int>(rc.size())),
                 true,
                 rc,
                 actual_final_resource_levels,
                 ref_spp(),
                 b_is_a_path_at_all,
                 b_feasible,
                 b_correctly_extended,
                 ed_last_extended_arc);
  // b. if feasible, add the solution
  if (b_is_a_path_at_all && b_feasible && b_correctly_extended) {
    if (post_process_rc) post_process_rc(rc);
    return true;
  }
    // c. print a warning as it shouldn't be the case
  else {
    if (!b_is_a_path_at_all)
      std::cerr << "Not a path." << std::endl;
    if (!b_feasible)
      std::cerr << "Not a feasible path." << std::endl;
    if (!b_correctly_extended)
      std::cerr << "Not correctly extended." << std::endl;
    printPath(std::cerr, path, rc);
    return false;
  }
}

RCSolution RCGraph::solution(const std::vector<edge>& path, const spp_res_cont& resource) const {
#ifdef DBG
//  printPath(std::cout, path, resource);
#endif

  RCSolution sol(resource.cost);

  // All arcs are consecutively considered
  //
  for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
    int a = boost::get(&Arc_Properties::num, g_, path[j]);
    int day = arcDay(a);
    const std::vector<int>& shifts = arcShifts(a);

    if(day != -1 && !shifts.empty()) {
      // if it's the first day
      if(sol.firstDay == -1) sol.firstDay = day;
      // append the shifts
      for(int s: shifts) sol.shifts.push_back(s);
    }
  }

  return sol;
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
void RCGraph::printGraph(int nLabel, int nShifts) const {

  // TITLE
  std::cout << "# " << std::endl;
  std::cout << "# GRAPH OF THE SUBPROBLEM " << std::endl;
  std::cout << "# " << std::endl;

  // THE NODES
  //
  printAllNodes(nLabel);

  // THE ARCS
  //
  printAllArcs(nLabel, nShifts);

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

std::string RCGraph::printArc(int a, int nLabel, int nShifts) const {
  std::stringstream rep;
  const Arc_Properties& arc_prop = arc(a);
  char buff[255];
  sprintf(buff, "# ARC   %5d  %15s  (%5d,%5d)  c=%12.2f  ", a, arcTypeName[arc_prop.type].c_str(),
      arcOrigin(a),  arcDestination(a), arc_prop.cost);
  rep << buff;
  sprintf(buff, "%10s  %2d:", (arc_prop.forbidden ? "forbidden" : "authorized"), arc_prop.day);
  rep << buff;
  for(int s=0; s<nShifts; ++s){
    if(s < arc_prop.shifts.size())
      sprintf(buff, " %3d", arc_prop.shifts[s]);
    else sprintf(buff, " %3s", "");
    rep << buff;
  }
  rep << (nShifts < arc_prop.shifts.size() ? " ..." : "    ");
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
void RCGraph::printAllArcs(int nLabel, int nShifts) const {
  std::cout << "#   ARCS (" << nArcs_ << "]" << std::endl;
  for(int a=0; a<nArcs_; a++) std::cout << printArc(a, nLabel, nShifts) << std::endl;
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

// Print the path (arcs, nodes, cost of each arc in the current network, etc.)
//
void RCGraph::printPath(std::ostream& out, std::vector<edge> path, spp_res_cont resource) const {

  // The successive nodes, and corresponding arc costs / time
  //
//  out << "# " << std::endl;
//  for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
//    int a = boost::get(&Arc_Properties::num, g_, path[j]);
//    out << "# \t| [ " << shortNameNode(boost::source( path[j], g_ )) << " ]";
//    out << "\t\tCost:  " << arcCost(a);
//    const std::vector<int> & consumptions = arcConsumptions(a);
//    for(int l=0; l<consumptions.size(); ++l)
//      out << "\t\t" << labelName[l] << ":" << consumptions[l];
//    out << "\t\t[" << (arcForbidden(a) ? "forbidden" : " allowed ") << "]" << std::endl;
//  }
//  out << "# " << std::endl;
  for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
    int a = boost::get(&Arc_Properties::num, g_, path[j]);
    out << printArc(a, resource.size()) << std::endl;
  }

  // Last node and total
  //
  out << "# \t| [";
  for(int s: sinks_) out << shortNameNode(s) << " ";
  out << "]" << std::endl;
  out << "# \t| ~TOTAL~   \t\tCost:   " << resource.cost;
  for(int l=0; l<resource.size(); ++l)
    out << "\t\t" << labelName[l] << ":" << resource.label_value(l);
  out << std::endl << "# \t| " << std::endl;
  out << "# \t| RC Solution: |";

  // Print it
  //
  int k=0;
  int firstDay = -1;
  for( int j = static_cast<int>( path.size() ) - 1; j >= 0;	--j){
    if(firstDay == -1 && boost::source( path[j], g_ ) == source_) {
      firstDay = boost::get(&Arc_Properties::day, g_, path[j]);
      while (k < firstDay) {
        out << " |";
        k++;
      }
    }
    for(int s: get( &Arc_Properties::shifts, g_, path[j])) {
      out << s << "|";
      k++;
    }
  }
  while(k < nDays_){
    out << " |";
    k++;
  }
  out << std::endl;
  out << std::endl;
}


//
// Comparison override for 1 and 2 resource paths
//
/////////////////////////////////////////////////

// 1 resource comparisons (== and <)
bool operator==( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ) {
  if ( res_cont_1.cost != res_cont_2.cost )
    return false;
  for (int l=0; l<res_cont_1.size(); ++l)
    if ( res_cont_1.label_value(l) != res_cont_2.label_value(l) )
      return false;
  return true;
}

bool operator<( const spp_res_cont& res_cont_1, const spp_res_cont& res_cont_2 ){
  if ( res_cont_1.cost < res_cont_2.cost )
    return true;
  if ( res_cont_1.cost > res_cont_2.cost )
    return false;
  for (int l=0; l<res_cont_1.size(); ++l) {
    int v1 = res_cont_1.label_value(l), v2 = res_cont_2.label_value(l);
    if (v1 < v2)
      return labelsOrder[l]; // if order is descending -> true
    if (v1 > v2)
      return !labelsOrder[l]; // if order is descending -> !true
  }
  // are equal -> false
  return false;
}


// modified boost::r_c_shortest_paths_dispatch function (body/implementation)
template<class Graph,
    class VertexIndexMap,
    class Resource_Container,
    class Resource_Extension_Function,
    class Dominance_Function,
    class Label_Allocator,
    class Visitor,
    class Splabel_Comparator>
void r_c_shortest_paths_dispatch(const Graph& g,
                                 const VertexIndexMap& vertex_index_map,
                                 vertex s,
                                 std::vector<vertex> t,
                                 vector2D<edge>& opt_solutions_spp,
                                 std::vector<spp_res_cont>& pareto_opt_rcs_spp,
                                 const Resource_Container& rc,
                                 const Resource_Extension_Function& ref,
                                 const Dominance_Function& dominance,
                                 Label_Allocator /*la*/,
                                 Visitor vis,
                                 Splabel_Comparator /*comp*/) {

  size_t i_label_num = 0;
  typedef
  typename
  Label_Allocator::template rebind
      <boost::r_c_shortest_paths_label
          <Graph, Resource_Container> >::other LAlloc;
  LAlloc l_alloc;
  typedef
  ks_smart_pointer
      <boost::r_c_shortest_paths_label<Graph, Resource_Container> > Splabel;
  std::priority_queue<Splabel, std::vector<Splabel>, Splabel_Comparator > unprocessed_labels;

  int nbCreatedLabels = 0;
  int nbDeletedLabels = 0;
  bool b_feasible = true;
  boost::r_c_shortest_paths_label<Graph, Resource_Container> *first_label =
      l_alloc.allocate(1);
  l_alloc.construct
      (first_label,
       boost::r_c_shortest_paths_label
           <Graph, Resource_Container>(i_label_num++,
                                       rc,
                                       0,
                                       typename boost::graph_traits<Graph>::
                                       edge_descriptor(),
                                       s));
  nbCreatedLabels++;

  Splabel splabel_first_label = Splabel(first_label);
  unprocessed_labels.push(splabel_first_label);
  std::vector<std::list<Splabel> > vec_vertex_labels_data(num_vertices(g));
  boost::iterator_property_map<typename std::vector<std::list<Splabel> >::iterator,
      VertexIndexMap>
      vec_vertex_labels(vec_vertex_labels_data.begin(), vertex_index_map);
  vec_vertex_labels[s].push_back(splabel_first_label);
  typedef
  std::vector<typename std::list<Splabel>::iterator>
      vec_last_valid_positions_for_dominance_data_type;
  vec_last_valid_positions_for_dominance_data_type
      vec_last_valid_positions_for_dominance_data(num_vertices(g));
  boost::iterator_property_map<
      typename vec_last_valid_positions_for_dominance_data_type::iterator,
      VertexIndexMap>
      vec_last_valid_positions_for_dominance
      (vec_last_valid_positions_for_dominance_data.begin(),
       vertex_index_map);
  BGL_FORALL_VERTICES_T(v, g, Graph) {
      put(vec_last_valid_positions_for_dominance, v, vec_vertex_labels[v].begin());
    }
  std::vector<size_t> vec_last_valid_index_for_dominance_data(num_vertices(g), 0);
  boost::iterator_property_map<std::vector<size_t>::iterator, VertexIndexMap>
      vec_last_valid_index_for_dominance
      (vec_last_valid_index_for_dominance_data.begin(), vertex_index_map);
  std::vector<bool>
      b_vec_vertex_already_checked_for_dominance_data(num_vertices(g), false);
  boost::iterator_property_map<std::vector<bool>::iterator, VertexIndexMap>
      b_vec_vertex_already_checked_for_dominance
      (b_vec_vertex_already_checked_for_dominance_data.begin(),
       vertex_index_map);

  std::vector<Splabel> label_trash;
#ifdef DBG
//  int iter = 0;
#endif
  while (!unprocessed_labels.empty() && vis.on_enter_loop(unprocessed_labels, g)) {
#ifdef DBG
//    if(iter++ % 1000 == 0)
//      std::cout << "RC SPP iteration " << iter << ": number of labels to process = " << unprocessed_labels.size() << std::endl;
#endif
    Splabel cur_label = unprocessed_labels.top();
    assert (cur_label->b_is_valid);
    unprocessed_labels.pop();
    vis.on_label_popped(*cur_label, g);
    // an Splabel object in unprocessed_labels and the respective Splabel
    // object in the respective list<Splabel> of vec_vertex_labels share their
    // embedded r_c_shortest_paths_label object
    // to avoid memory leaks, dominated
    // r_c_shortest_paths_label objects are marked and deleted when popped
    // from unprocessed_labels, as they can no longer be deleted at the end of
    // the function; only the Splabel object in unprocessed_labels still
    // references the r_c_shortest_paths_label object
    // this is also for efficiency, because the else branch is executed only
    // if there is a chance that extending the
    // label leads to new undominated labels, which in turn is possible only
    // if the label to be extended is undominated
    assert (cur_label->b_is_valid);
    if (!cur_label->b_is_dominated) {
      typename boost::graph_traits<Graph>::vertex_descriptor
          i_cur_resident_vertex = cur_label->resident_vertex;
      std::list<Splabel> &list_labels_cur_vertex =
          get(vec_vertex_labels, i_cur_resident_vertex);
      if (list_labels_cur_vertex.size() >= 2
          && vec_last_valid_index_for_dominance[i_cur_resident_vertex]
             < list_labels_cur_vertex.size()) {
        typename std::list<Splabel>::iterator outer_iter =
            list_labels_cur_vertex.begin();
        bool b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = false;
        while (outer_iter != list_labels_cur_vertex.end()) {
          Splabel cur_outer_splabel = *outer_iter;
          assert (cur_outer_splabel->b_is_valid);
          typename std::list<Splabel>::iterator inner_iter = outer_iter;
          if (!b_outer_iter_at_or_beyond_last_valid_pos_for_dominance
              && outer_iter ==
                 get(vec_last_valid_positions_for_dominance,
                     i_cur_resident_vertex))
            b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = true;
          if (!get(b_vec_vertex_already_checked_for_dominance, i_cur_resident_vertex)
              || b_outer_iter_at_or_beyond_last_valid_pos_for_dominance) {
            ++inner_iter;
          } else {
            inner_iter =
                get(vec_last_valid_positions_for_dominance,
                    i_cur_resident_vertex);
            ++inner_iter;
          }
          bool b_outer_iter_erased = false;
          while (inner_iter != list_labels_cur_vertex.end()) {
            Splabel cur_inner_splabel = *inner_iter;
            assert (cur_inner_splabel->b_is_valid);
            if (dominance(cur_outer_splabel->
                              cumulated_resource_consumption,
                          cur_inner_splabel->
                              cumulated_resource_consumption)) {
              typename std::list<Splabel>::iterator buf = inner_iter;
              ++inner_iter;
              list_labels_cur_vertex.erase(buf);
              if (cur_inner_splabel->b_is_processed) {
                cur_inner_splabel->b_is_valid = false;
                l_alloc.destroy(cur_inner_splabel.get());
                l_alloc.deallocate(cur_inner_splabel.get(), 1);
                nbDeletedLabels++;
              } else
                cur_inner_splabel->b_is_dominated = true;
              continue;
            } else
              ++inner_iter;
            if (dominance(cur_inner_splabel->
                              cumulated_resource_consumption,
                          cur_outer_splabel->
                              cumulated_resource_consumption)) {
              typename std::list<Splabel>::iterator buf = outer_iter;
              ++outer_iter;
              list_labels_cur_vertex.erase(buf);
              b_outer_iter_erased = true;
              assert (cur_outer_splabel->b_is_valid);
              if (cur_outer_splabel->b_is_processed) {
//                label_trash.push_back(cur_outer_splabel);
                cur_outer_splabel->b_is_valid = false;
                l_alloc.destroy(cur_outer_splabel.get());
                l_alloc.deallocate(cur_outer_splabel.get(), 1);
                nbDeletedLabels++;
              } else
                cur_outer_splabel->b_is_dominated = true;
              break;
            }
          }
          if (!b_outer_iter_erased)
            ++outer_iter;
        }
        if (list_labels_cur_vertex.size() > 1)
          put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex,
              (--(list_labels_cur_vertex.end())));
        else
          put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex,
              list_labels_cur_vertex.begin());
        put(b_vec_vertex_already_checked_for_dominance,
            i_cur_resident_vertex, true);
        put(vec_last_valid_index_for_dominance, i_cur_resident_vertex,
            list_labels_cur_vertex.size() - 1);
      }
    }
    assert (cur_label->b_is_valid);

    // ------------------------------------------------------------------------- START MODIFICATION
//    if( !b_all_pareto_optimal_solutions && cur_label->resident_vertex == t ) {
//      // the devil don't sleep
//      if (cur_label->b_is_dominated) {
//        cur_label->b_is_valid = false;
//        l_alloc.destroy(cur_label.get());
//        l_alloc.deallocate(cur_label.get(), 1);
//        nbDeletedLabels++;
//      }
//      while (unprocessed_labels.size()) {
//        Splabel l = unprocessed_labels.top();
//        assert (l->b_is_valid);
//        unprocessed_labels.pop();
//        // delete only dominated labels, because nondominated labels are
//        // deleted at the end of the function
//        if (l->b_is_dominated) {
//          l->b_is_valid = false;
//          l_alloc.destroy(l.get());
//          l_alloc.deallocate(l.get(), 1);
//          nbDeletedLabels++;
//        }
//      }
//    }
    // ------------------------------------------------------------------------- END MODIFICATION
    if (!cur_label->b_is_dominated) {
      cur_label->b_is_processed = true;
      vis.on_label_not_dominated(*cur_label, g);
      typename boost::graph_traits<Graph>::vertex_descriptor cur_vertex =
          cur_label->resident_vertex;
      typename boost::graph_traits<Graph>::out_edge_iterator oei, oei_end;
      for (boost::tie(oei, oei_end) = out_edges(cur_vertex, g);
           oei != oei_end;
           ++oei) {
        b_feasible = true;
        boost::r_c_shortest_paths_label<Graph, Resource_Container> *new_label =
            l_alloc.allocate(1);
        l_alloc.construct(new_label,
                          boost::r_c_shortest_paths_label
                              <Graph, Resource_Container>
                              (i_label_num++,
                               cur_label->cumulated_resource_consumption,
                               cur_label.get(),
                               *oei,
                               target(*oei, g)));
        nbCreatedLabels++;

        b_feasible =
            ref(g,
                new_label->cumulated_resource_consumption,
                new_label->p_pred_label->cumulated_resource_consumption,
                new_label->pred_edge);

        if (!b_feasible) {
          vis.on_label_not_feasible(*new_label, g);
          new_label->b_is_valid = false;
          l_alloc.destroy(new_label);
          l_alloc.deallocate(new_label, 1);
          nbDeletedLabels++;
        } else {
          const boost::r_c_shortest_paths_label<Graph, Resource_Container> &
              ref_new_label = *new_label;
          vis.on_label_feasible(ref_new_label, g);
          Splabel new_sp_label(new_label);
          vec_vertex_labels[new_sp_label->resident_vertex].
              push_back(new_sp_label);
          unprocessed_labels.push(new_sp_label);
        }
      }
    } else {
      assert (cur_label->b_is_valid);
      vis.on_label_dominated(*cur_label, g);
      cur_label->b_is_valid = false;
      l_alloc.destroy(cur_label.get());
      l_alloc.deallocate(cur_label.get(), 1);
      nbDeletedLabels++;
    }
  }

  typename std::list<Splabel>::const_iterator csi;
  typename std::list<Splabel>::const_iterator csi_end;

  // ------------------------------------------------------------------------- START MODIFICATION
  for (int sink = 0; sink < t.size(); sink++) {
    std::list<Splabel> dsplabels = get(vec_vertex_labels, t[sink]);
    csi = dsplabels.begin();
    csi_end = dsplabels.end();
    // if d could be reached from o
    if (!dsplabels.empty()) {
      for (; csi != csi_end; ++csi) {
        std::vector<typename boost::graph_traits<Graph>::edge_descriptor>
            cur_pareto_optimal_path;
        const boost::r_c_shortest_paths_label<Graph, Resource_Container> *p_cur_label =
            (*csi).get();
        assert (p_cur_label->b_is_valid);
        pareto_opt_rcs_spp.
            push_back(p_cur_label->cumulated_resource_consumption);
        while (p_cur_label->num != 0) {
          cur_pareto_optimal_path.push_back(p_cur_label->pred_edge);
          p_cur_label = p_cur_label->p_pred_label;
          assert (p_cur_label->b_is_valid);
        }
        opt_solutions_spp.push_back(cur_pareto_optimal_path);
      }
    }
  }
  // ------------------------------------------------------------------------- END MODIFICATION

  BGL_FORALL_VERTICES_T(i, g, Graph) {
      const std::list<Splabel> &list_labels_cur_vertex = vec_vertex_labels[i];
      csi_end = list_labels_cur_vertex.end();
      for (csi = list_labels_cur_vertex.begin(); csi != csi_end; ++csi) {
        assert ((*csi)->b_is_valid);
        (*csi)->b_is_valid = false;
        l_alloc.destroy((*csi).get());
        l_alloc.deallocate((*csi).get(), 1);
        nbDeletedLabels++;
      }
    }
//  for (Splabel label: label_trash) {
//    assert(label->b_is_valid);
//    label->b_is_valid = false;
//    l_alloc.destroy(label.get());
//    l_alloc.deallocate(label.get(), 1);
//  }
  // std::cout << "Nb created labels " << nbCreatedLabels << std::endl;
  // std::cout << "Nb deleted labels " << nbDeletedLabels << std::endl;
} // r_c_shortest_paths_dispatch


///////////////////////////////////////////
