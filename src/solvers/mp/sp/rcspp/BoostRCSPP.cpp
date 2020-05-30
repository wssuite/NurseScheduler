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

#include "BoostRCSPP.h"

#include <memory>

void spp_res_cont::print(std::ostream &out) const {
  out << "Cost: " << cost << std::endl;
  for (int l = 0; l < size(); l++) {
    out << labelName[l].c_str() << "=" << label_value(l) << "  ";
  }
  out << std::endl;
#ifdef DBG
  RCSolution sol(first_day, shifts_, cost);
  out << sol.toString();
  out << "Arc taken: " << pred_arc << std::endl;
#endif
}

void spp_res_cont::dominate(const spp_res_cont &res) {
  // if first dominance -> set to parent lvl
  // the goal is to ensure that this label continue to be processed first
  // while they dominate other labels.
  if (dominanceLvl < parentDominanceLvl)
    dominanceLvl = parentDominanceLvl;
  if (dominanceLvl <= res.parentDominanceLvl)
    dominanceLvl = res.parentDominanceLvl + 1;
}

//
// Comparison override for 1 and 2 resource paths
//
/////////////////////////////////////////////////

// 1 resource comparisons (== and <)
bool operator==(const spp_res_cont &res_cont_1,
                const spp_res_cont &res_cont_2) {
  if (res_cont_1.cost != res_cont_2.cost)
    return false;
  for (int l = 0; l < res_cont_1.size(); ++l)
    if (res_cont_1.label_value(l) != res_cont_2.label_value(l))
      return false;
  return true;
}

bool operator!=(const spp_res_cont &res_cont_1,
                const spp_res_cont &res_cont_2) {
  return !(res_cont_1 == res_cont_2);
}

bool operator<(const spp_res_cont &res_cont_1, const spp_res_cont &res_cont_2) {
  if (res_cont_1.cost < res_cont_2.cost)
    return true;
  if (res_cont_1.cost > res_cont_2.cost)
    return false;
  for (int l = 0; l < res_cont_1.size(); ++l) {
    int v1 = res_cont_1.label_value(l), v2 = res_cont_2.label_value(l);
    if (v1 < v2)
      return labelsOrder[l];  // if order is descending -> true
    if (v1 > v2)
      return !labelsOrder[l];  // if order is descending -> !true
  }
  // are equal -> false
  return false;
}

// Resources extension model (arc has cost + label consumptions)
bool ref_spp::operator()(const Graph &g,
                         spp_res_cont *new_cont,
                         const spp_res_cont &old_cont,
                         const edge &ed) const {
  const Arc_Properties &arc_prop = get(boost::edge_bundle, g)[ed];
  const Vertex_Properties
      &vert_prop = get(boost::vertex_bundle, g)[target(ed, g)];
#ifdef DBG
  if (arc_prop.forbidden)
    return false;
  if (vert_prop.forbidden)
    return false;
#endif
  new_cont->cost = old_cont.cost + arc_prop.cost;
  if (new_cont->first_day == -1) new_cont->first_day = arc_prop.day;
  if (arc_prop.day != -1) new_cont->day = arc_prop.day;

#ifdef DBG
  new_cont->pred_arc = arc_prop.num;
  if (!arc_prop.shifts.empty()) {
    new_cont->shifts_.insert(new_cont->shifts_.end(),
                             arc_prop.shifts.begin(),
                             arc_prop.shifts.end());
  }
#endif

  assert(old_cont.size() == new_cont->size());

  for (int l = 0; l < old_cont.size(); ++l) {
    int lv = old_cont.label_value(l) + arc_prop.consumption(l);
    int lb = vert_prop.lb(l);
    if (lv < lb) {
      if (vert_prop.hard_lbs_)
        return false;
      lv = lb;
    }
    if (lv > vert_prop.ub(l))
      return false;
    new_cont->label_values[l] = lv;
  }

  // set parent dominance lvl and reset dominance lvl
  // When his children will be created, the  dominance lvl will be higher than
  // all dominated, and thus they will be processed before the children of
  // the dominated path.
  new_cont->dominanceLvl = 0;
  new_cont->parentDominanceLvl = old_cont.dominanceLvl;

  return true;
}

// Dominance function model
bool dominance_spp::operator()(spp_res_cont *res_cont_1,
                               spp_res_cont *res_cont_2) const {
  // must be "<=" here!!!
  // must NOT be "<"!!!
  if (res_cont_1->cost > res_cont_2->cost + epsilon_) return false;

  /* Dominance:
   * if label is increasing -> can dominate after a certain level:
   *    if the label do not reach the minimum level, there is no additional cost
   * if label is decreasing -> cannot dominate (min level should be an ub).
   *    the label could continue to decrease in the future and not imply
   *    any additional cost.
   * So, res_cont_1 dominates res_cont_2 in one of these 3 situations:
   * a- cost1 < cost2 at epsilon and all label1 <= label2
   * b- cost1 == cost2 at epsilon and all label1 <= label2 and
   *    it exists one label1 < label2 with label2 > minLevel (significant)
   * c- res_cont_1 == res_cont_2
   */

  bool dominate = (res_cont_1->cost < res_cont_2->cost - epsilon_),
      biggerThanMinLevel = false,
      equal = !dominate;
  for (int l = 0; l < res_cont_1->size(); ++l) {
    int label1 = res_cont_1->label_value(l),
        label2 = res_cont_2->label_value(l),
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
//    res_cont_2->print();
//    std::cout << "******************* BY ********************" << std::endl;
//    res_cont_1->print();
//    std::cout << "*******************************************" << std::endl;
//  }
#endif

  // a, b, c
  if (dominate || biggerThanMinLevel || equal) {
    res_cont_1->dominate(*res_cont_2);
    return true;
  }
  return false;
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

// Comparator for the priority queue used by boost to process the labels.
// Note that the Compare parameter of a priority queue is defined such that
// it returns true if its first argument comes before its second argument in
// a weak ordering.
// But because the priority queue outputs largest elements first,
// the elements that "come before" are actually output last.
// That is, the front of the queue contains the "last" element according to the
// weak ordering imposed by Compare
// -> return True if element should be processed last (so first in the queue)

// Breath First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the smallest day (return true if day1 > day2)
// 2. break tie with cost (return true if cost1 > cost2)
bool SpplabelBreadthFirstComparator::operator()(
    const Spplabel &splabel1,
    const Spplabel &splabel2) const {
  if (splabel1->cumulated_resource_consumption.day
      > splabel2->cumulated_resource_consumption.day)
    return true;
  if (splabel1->cumulated_resource_consumption.day
      < splabel2->cumulated_resource_consumption.day)
    return false;
  return splabel1->cumulated_resource_consumption.cost
      > splabel2->cumulated_resource_consumption.cost + 1e-5;
}

// Depth First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the biggest day (day1 < day2)
// 2. break tie with cost (return true if cost1 > cost2)
bool SpplabelDepthFirstComparator::operator()(
    const Spplabel &splabel1,
    const Spplabel &splabel2) const {
  if (splabel1->cumulated_resource_consumption.day
      < splabel2->cumulated_resource_consumption.day)
    return true;
  if (splabel1->cumulated_resource_consumption.day
      > splabel2->cumulated_resource_consumption.day)
    return false;
  return splabel1->cumulated_resource_consumption.cost
      > splabel2->cumulated_resource_consumption.cost + 1e-5;
}

// Best First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the smallest cost (return true if cost1 > cost2)
bool SpplabelBestFirstComparator::operator()(
    const Spplabel &splabel1,
    const Spplabel &splabel2) const {
  return splabel1->cumulated_resource_consumption.cost
      > splabel2->cumulated_resource_consumption.cost + 1e-5;
}

// Dominant First Comparator to order the processing of the nodes
// in boost rc spp:
// 1. starts with the ones with the most parent dominant level (order the label
// by dominance)
// 2. starts with the ones with the smallest cost (return true if cost1 > cost2)
bool SpplabelDominantFirstComparator::operator()(
    const Spplabel &splabel1,
    const Spplabel &splabel2) const {
  if (splabel1->cumulated_resource_consumption.parentDominanceLvl >
      splabel2->cumulated_resource_consumption.parentDominanceLvl)
    return false;  // should be processed first
  if (splabel1->cumulated_resource_consumption.parentDominanceLvl <
      splabel2->cumulated_resource_consumption.parentDominanceLvl)
    return true;  // should be processed last
  return splabel1->cumulated_resource_consumption.cost
      > splabel2->cumulated_resource_consumption.cost + 1e-5;
}

BoostRCSPPSolver::BoostRCSPPSolver(
    RCGraph *rcg,
    double maxReducedCostBound,
    double epsilon,
    SPSearchStrategy strategy,
    int nb_max_paths,
    std::function<void(spp_res_cont *)> post_process_rc) :
    rcg_(rcg),
    maxReducedCostBound_(maxReducedCostBound),
    epsilon_(epsilon),
    strategy_(strategy),
    nb_max_paths_(nb_max_paths),
    post_process_rc_(post_process_rc) {}

std::vector<RCSolution> BoostRCSPPSolver::solve(
    int nLabels,
    const std::vector<int> &labelsMinLevel,
    std::vector<vertex> sinks) {
  // 1 - solve the resource constraints shortest path problem
  dominance_spp dominance(labelsOrder, labelsMinLevel, epsilon_);
  vector2D<edge> opt_solutions_spp;
  std::vector<spp_res_cont> pareto_opt_rcs_spp;
  rc_spp_visitor
      vis(nb_max_paths_, sinks, post_process_rc_, maxReducedCostBound_);
  r_c_shortest_paths_solve(
      rcg_->g(),
      rcg_->source(),
      sinks,
      &opt_solutions_spp,
      &pareto_opt_rcs_spp,
      spp_res_cont(0, std::vector<int>(nLabels)),
      ref_spp(),
      dominance,
      std::allocator<boost::r_c_shortest_paths_label<Graph, spp_res_cont> >(),
      &vis,  // boost::default_r_c_shortest_paths_visitor(),
      strategy_);

  // 2 - Post process the solutions
  // process paths if needed (process_solution is not used)
  std::vector<RCSolution> rc_solutions;
  // For each path of the list, store the solutions with a negative cost
  for (unsigned int p = 0; p < opt_solutions_spp.size(); ++p) {
    spp_res_cont &rc = pareto_opt_rcs_spp[p];
    if (processPath(&opt_solutions_spp[p], &rc, post_process_rc_))
      if (rc.cost < maxReducedCostBound_)
        rc_solutions.push_back(solution(opt_solutions_spp[p], rc));
  }

  return rc_solutions;
}

bool BoostRCSPPSolver::processPath(
    std::vector<edge> *path,
    spp_res_cont *rc,
    std::function<void(spp_res_cont *)> post_process_rc,
    bool printBadPath) const {
  // a. Check if it is valid
  bool b_is_a_path_at_all = false;
  bool b_feasible = false;
  bool b_correctly_extended = false;
  spp_res_cont actual_final_resource_levels(0, std::vector<int>(rc->size()));
  boost::graph_traits<Graph>::edge_descriptor ed_last_extended_arc;
  ref_spp ref;
  check_r_c_path(rcg_->g(),
                 *path,
                 spp_res_cont(0, std::vector<int>(rc->size())),
                 true,
                 *rc,
                 actual_final_resource_levels,
                 [&ref](const Graph &g,
                        spp_res_cont &new_cont,  // NOLINT
                        const spp_res_cont &old_cont,
                        const edge &ed) {
                   return ref(g, &new_cont, old_cont, ed);
                 },
                 b_is_a_path_at_all,
                 b_feasible,
                 b_correctly_extended,
                 ed_last_extended_arc);
  // b. if feasible, postprocess the solution
  if (b_is_a_path_at_all && b_feasible && b_correctly_extended) {
    if (post_process_rc) post_process_rc(rc);
    return true;
  } else {
    // c. print a warning as it shouldn't be the case
    if (printBadPath) {
      if (!b_is_a_path_at_all)
        std::cerr << "Not a path." << std::endl;
      if (!b_feasible)
        std::cerr << "Not a feasible path." << std::endl;
      if (!b_correctly_extended)
        std::cerr << "Not correctly extended." << std::endl;
      printPath(std::cerr, *path, *rc);
      std::cerr << "vs" << std::endl;
      actual_final_resource_levels.print(std::cerr);
    }
    return false;
  }
}

RCSolution BoostRCSPPSolver::solution(const std::vector<edge> &path,
                                      const spp_res_cont &resource) const {
#ifdef DBG
  //  printPath(std::cout, path, resource);
#endif

  RCSolution sol(resource.cost);

  // All arcs are consecutively considered
  for (int j = path.size() - 1; j >= 0; --j) {
    int a = boost::get(&Arc_Properties::num, rcg_->g(), path[j]);
    int day = rcg_->arcDay(a);
    const std::vector<int> &shifts = rcg_->arcShifts(a);

    if (day != -1 && !shifts.empty()) {
      // if it's the first day
      if (sol.firstDay == -1) sol.firstDay = day;
      // append the shifts
      for (int s : shifts) sol.shifts.push_back(s);
    }
  }
  return sol;
}

// Print the path (arcs, nodes, cost of each arc in the current network, etc.)
void BoostRCSPPSolver::printPath(std::ostream &out,
                                 std::vector<edge> path,
                                 spp_res_cont resource) const {
  // The successive nodes, and corresponding arc costs / time
  for (int j = path.size() - 1; j >= 0; --j) {
    int a = boost::get(&Arc_Properties::num, rcg_->g(), path[j]);
    out << rcg_->printArc(a, resource.size()) << std::endl;
  }

  // Last node and total
  out << "# \t| ~TOTAL~   \t\tCost:   " << resource.cost;
  for (int l = 0; l < resource.size(); ++l)
    out << "\t\t" << labelName[l] << ":" << resource.label_value(l);
  out << std::endl << "# \t| " << std::endl;
  out << "# \t| RC Solution: |";

  // Print it
  int k = 0;
  int firstDay = -1;
  for (int j = path.size() - 1; j >= 0; --j) {
    if (firstDay == -1 && boost::source(path[j], rcg_->g()) == rcg_->source()) {
      firstDay = boost::get(&Arc_Properties::day, rcg_->g(), path[j]);
      while (k < firstDay) {
        out << " |";
        k++;
      }
    }
    for (int s : get(&Arc_Properties::shifts, rcg_->g(), path[j])) {
      out << s << "|";
      k++;
    }
  }
  while (k < rcg_->nDays()) {
    out << " |";
    k++;
  }
  out << std::endl;
  out << std::endl;
}

// rc_spp_visitor struct: derived from boost::default_r_c_shortest_paths_visitor
rc_spp_visitor::rc_spp_visitor(
    int nPaths,
    const std::vector<vertex> &sinks,
    std::function<void(spp_res_cont *)> post_process_rc,
    double maxReducedCostBound) :
    nMaxPaths_(nPaths),
    sinks_(sinks),
    postProcessRc_(post_process_rc),
    maxReducedCostBound_(maxReducedCostBound) {}

void rc_spp_visitor::on_label_popped(Label &spplabel, const Graph &) {}

void rc_spp_visitor::on_label_feasible(Label &spplabel, const Graph &) {}

void rc_spp_visitor::on_label_not_feasible(Label &, const Graph &) {}

void rc_spp_visitor::on_label_dominated(Label &spplabel, const Graph &) {}

void rc_spp_visitor::on_label_not_dominated(Label &spplabel, const Graph &) {
  // if seek all paths, do nothing
  if (nMaxPaths_ == -1) return;
  // check if we found a non-dominated path from a source to a sink
  if (find(sinks_.begin(), sinks_.end(), spplabel.resident_vertex)
      != sinks_.end()) {
    paths_[spplabel.num] = spplabel.cumulated_resource_consumption;
  }
}

template<typename Queue>
bool rc_spp_visitor::on_enter_loop(const Queue &queue, const Graph &graph) {
  // if go until optimality
  if (nMaxPaths_ == -1) return true;
  // if exceed the number of paths searched
  if (paths_.size() >= nMaxPaths_)
    return false;
  // stop if have explored too many feasible nodes (avoid stalling)
  return nLoop_++ < nMaxPaths_ * num_vertices(graph);
}
