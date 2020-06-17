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

#include <algorithm>
#include <memory>

void spp_res_cont::print(std::ostream &out) const {
  out << "Cost: " << cost << "\t\tPostprocess cost: "
      << postprocessCost << std::endl;
  for (int l = 0; l < size(); l++) {
    out << "\t\t" << labelsName[l].c_str() << "=" << label_value(l);
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
    if (v1 < v2) return true;
    if (v1 > v2) return false;
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
    LABEL label = labels_[l];
    int lv = old_cont.label_value(l);
    // check if should price label
    if (arc_prop.findLabelToPrice(label))
      new_cont->cost += arc_prop.penalties.penalty(label, lv);
    lv += arc_prop.consumption(label);
    int lb = vert_prop.lb(label);
    if (lv < lb) {
      if (vert_prop.hard_lbs_)
        return false;
      lv = lb;
    }
    if (lv > vert_prop.ub(label))
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
                               const spp_res_cont *res_cont_2) const {
  // must be "<=" here at epsilon!!!
  // must NOT be "<"!!!
  double worst_penalty = worstPenalty(*res_cont_1, *res_cont_2);
  if (res_cont_1->cost + worst_penalty <= res_cont_2->cost + epsilon_) {
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

double dominance_spp::worstPenalty(
    const spp_res_cont &res_cont_1,
    const spp_res_cont &res_cont_2) const {
  double penalty = 0;
  for (int i = 0; i < res_cont_1.size(); i++)
    penalty += worstPenalty(i,
                            res_cont_1.label_value(i),
                            res_cont_2.label_value(i));
  return penalty;
}

double dominance_spp::worstPenalty(int l, int level1, int level2) const {
  // i. if level 1 < level2
  if (level1 < level2) {
    int minLvl = penalties_.minLevel(l), maxLvl = penalties_.maxLevel(l);
    // a. compute penalty with respect to the min level.
    // worst case: stop working immediately (as the penalty is linear)
    // -> pay the difference between penalties
    double min_penalty = penalties_.weight(l) *
        (std::max(0, minLvl - level1) - std::max(0, minLvl - level2));
    // b. compute penalty with respect to the max level.
    // worst case: stop working immediately (as the penalty is linear)
    // -> pay the difference between penalties
    double max_penalty = penalties_.weight(l) *
        (std::max(0, level1 - maxLvl) - std::max(0, level2 - maxLvl));
    // return sum of both
    // if both are different of 0:
    //   - lvl1 < minLvl: min_penalty > 0
    //   - lvl2 > maxLvl: max_penalty < 0
    // the sum remains the worst case as there is the same weight penalty
    // associated to the min and the max (otherwise, it would be a little
    // more complex).
    assert(min_penalty == 0 || max_penalty == 0 ||
        (min_penalty > 0 && max_penalty < 0));
    return min_penalty + max_penalty;
  } else {
    // ii. if level 1 >= level2
    // a. compute penalty with respect to the min level.
    // worst case: continue to work until level2 >= minLvl
    // -> no penalty
    // b. compute penalty with respect to the max level.
    // worst case: continue working until level2 reaches maxLvl
    // -> look at the difference between level (as the penalty is linear)
    return penalties_.weight(l) * (level1 - level2);
  }
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
    int verbose,
    double epsilon,
    SPSearchStrategy strategy,
    int nb_max_paths,
    std::function<void(spp_res_cont *)> post_process_rc) :
    rcg_(rcg),
    verbose_(verbose),
    maxReducedCostBound_(maxReducedCostBound),
    epsilon_(epsilon),
    strategy_(strategy),
    nb_max_paths_(nb_max_paths),
    post_process_rc_(post_process_rc) {}

std::vector<RCSolution> BoostRCSPPSolver::solve(
    std::vector<LABEL> labels,
    const Penalties &penalties,
    std::vector<vertex> sinks) {
  timer_.start();
  // 1 - solve the resource constraints shortest path problem
  ref_spp ref(labels);
  dominance_spp dominance(labels, penalties, epsilon_);
  vector2D<edge> opt_solutions_spp;
  std::vector<spp_res_cont> pareto_opt_rcs_spp;
  rc_spp_visitor vis(nb_max_paths_,
                     sinks,
                     post_process_rc_,
                     maxReducedCostBound_);
  r_c_shortest_paths_solve(
      rcg_->g(),
      rcg_->source(),
      sinks,
      &opt_solutions_spp,
      &pareto_opt_rcs_spp,
      spp_res_cont(0, std::vector<int>(labels.size())),
      ref,
      dominance,
      std::allocator<boost::r_c_shortest_paths_label<Graph, spp_res_cont> >(),
      &vis,  // boost::default_r_c_shortest_paths_visitor(),
      strategy_);

  // 2 - Post process the solutions
  // process paths if needed (process_solution is not used)
  std::vector<RCSolution> rc_solutions;
  // For each path of the list, store the solutions with a negative cost
  double bestCost = maxReducedCostBound_;
  for (unsigned int p = 0; p < opt_solutions_spp.size(); ++p) {
    spp_res_cont &rc = pareto_opt_rcs_spp[p];
    if (processPath(&opt_solutions_spp[p], &rc, ref, post_process_rc_)) {
      if (rc.cost < maxReducedCostBound_) {
        if (verbose_ >= 4)
          printPath(std::cout, opt_solutions_spp[p], rc);
        rc_solutions.push_back(solution(opt_solutions_spp[p], rc));
        if (rc.cost < bestCost) bestCost = rc.cost;
      }
    }
  }

  timer_.stop();
  if (verbose_ >= 3) {
    vis.printStats(std::cout);
    std::cout << "Time spent in the RCSPP: " << timer_.dSinceStart()
              << " and best solution cost: " << bestCost << std::endl;
  }

  return rc_solutions;
}

bool BoostRCSPPSolver::processPath(
    std::vector<edge> *path,
    spp_res_cont *rc,
    const ref_spp &ref,
    std::function<void(spp_res_cont *)> post_process_rc,
    bool printBadPath) const {
  // a. Check if it is valid
  spp_res_cont checked_final_resource_levels(0, std::vector<int>(rc->size()));
  bool valid = check_r_c_path(*path,
                              *rc,
                              &checked_final_resource_levels,
                              ref);
  // b. if valid, postprocess the solution
  if (valid) {
    if (post_process_rc) post_process_rc(rc);
    return true;
  } else {
    // c. print a warning as it shouldn't be the case
    if (printBadPath) {
      printPath(std::cerr, *path, *rc);
      std::cerr << "vs" << std::endl;
      checked_final_resource_levels.print(std::cerr);
    }
    return false;
  }
}

RCSolution BoostRCSPPSolver::solution(const std::vector<edge> &path,
                                      const spp_res_cont &resource) const {
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
                                 const std::vector<edge> &path,
                                 const spp_res_cont &resource) const {
  // The successive nodes, and corresponding arc costs / time
  for (int j = path.size() - 1; j >= 0; --j) {
    int a = boost::get(&Arc_Properties::num, rcg_->g(), path[j]);
    out << rcg_->printArc(a, resource.size()) << std::endl;
  }

  // Last node and total
  out << "# \t| ~TOTAL~   \t\t";
  resource.print(out);
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

// modified boost::check_r_c_path function
bool BoostRCSPPSolver::check_r_c_path(
    std::vector<edge> ed_vec_path,
    // if true, computed accumulated final resource levels must
    // be equal to desired_final_resource_levels
    // if false, computed accumulated final resource levels must
    // be less than or equal to desired_final_resource_levels
    const spp_res_cont &found_final_resource_levels,
    spp_res_cont *checked_final_resource_levels,
    const ref_spp &ref) const {
  // if empty, return checked <= found
  if (ed_vec_path.empty())
    return *checked_final_resource_levels < found_final_resource_levels
        || *checked_final_resource_levels == found_final_resource_levels;

  // if size == 1, check if target(0) == source(1)
  // return checked <= found
  if (ed_vec_path.size() == 1) {
    // if target(0) == source(1)
    if (target(ed_vec_path[0], rcg_->g()) ==
        source(ed_vec_path[1], rcg_->g()))
      return *checked_final_resource_levels < found_final_resource_levels
          || *checked_final_resource_levels == found_final_resource_levels;
    // otherwise, return false
    std::cerr << "Target of arc 0 is different of source of arc 1."
              << std::endl;
    Arc_Properties arc_prop = get(
        boost::edge_bundle, rcg_->g())[ed_vec_path[0]];
    std::cerr << "Arc 0: " << rcg_->printArc(
        arc_prop, found_final_resource_levels.size()) << std::endl;
    arc_prop = get(boost::edge_bundle, rcg_->g())[ed_vec_path[1]];
    std::cerr << "Arc 1: " << rcg_->printArc(
        arc_prop, found_final_resource_levels.size()) << std::endl;
    return false;
  }

  // reverse the path
  std::reverse(ed_vec_path.begin(), ed_vec_path.end());

  // check if target(i) == source(i+1)
  for (int i = 0; i < ed_vec_path.size() - 1; ++i)
    if (target(ed_vec_path[i], rcg_->g()) !=
        source(ed_vec_path[i + 1], rcg_->g())) {
      std::cerr << "Target of arc " << i << "is different of source of arc "
                << i + 1 << "." << std::endl;
      Arc_Properties arc_prop = get(
          boost::edge_bundle, rcg_->g())[ed_vec_path[i]];
      std::cerr << "Arc " << i << ": " << rcg_->printArc(
          arc_prop, found_final_resource_levels.size()) << std::endl;
      arc_prop = get(boost::edge_bundle, rcg_->g())[ed_vec_path[i + 1]];
      std::cerr << "Arc " << i + 1 << ": " << rcg_->printArc(
          arc_prop, found_final_resource_levels.size()) << std::endl;
      return false;
    }

  // extension of path
  for (int i = 0; i < ed_vec_path.size(); ++i) {
    spp_res_cont current_resource_levels = *checked_final_resource_levels;
    // if extension infeasible
    if (!ref(rcg_->g(),
             checked_final_resource_levels,
             current_resource_levels,
             ed_vec_path[i])) {
      std::cerr << "Extension of arc " << i << " is infeasible:" << std::endl;
      Arc_Properties arc_prop = get(
          boost::edge_bundle, rcg_->g())[ed_vec_path[i]];
      std::cerr << rcg_->printArc(
          arc_prop, found_final_resource_levels.size()) << std::endl;
      std::cerr << "For resource: ";
      current_resource_levels.print(std::cerr);
      return false;
    }
  }
  // return checked <= final
  return *checked_final_resource_levels < found_final_resource_levels
      || *checked_final_resource_levels == found_final_resource_levels;
}  // check_path

// build a path and the corresponding label by backtracking an original one
void BoostRCSPPSolver::backtrack(
    const Graph &g,
    const std::vector<std::list<Spplabel> > &vec_vertex_labels,
    const Label *p_original_label,
    const Label *p_cur_label,
    const ref_spp &ref,
    const dominance_spp &dominance,
    vector2D<edge> *opt_solutions_spp,
    std::vector<spp_res_cont> *pareto_opt_rcs_spp,
    std::vector<edge> path) const {
  assert(p_cur_label->b_is_valid);
  spp_res_cont original_cont = p_cur_label->cumulated_resource_consumption;

  // 1 - if source, store path and build resource container
  if (p_cur_label->num == 0) {
    spp_res_cont cur_cont = p_cur_label->cumulated_resource_consumption;
    for (int i = path.size() - 1; i >= 0; --i) {
      spp_res_cont old_cont = cur_cont;
      // extend path. If not feasible, return (could happen when solving
      // heuristically and the bounds are tighter)
      if (!ref(g, &cur_cont, old_cont, path[i])) return;
    }
#ifdef DBG
    assert(dominance(&cur_cont,
                     &(p_original_label->cumulated_resource_consumption)));
#endif
    opt_solutions_spp->push_back(path);
    pareto_opt_rcs_spp->push_back(cur_cont);
    return;
  }

  // 2 - find dominating labels when backtracking
  // a - add the edge to the path
  path.push_back(p_cur_label->pred_edge);
  // b - verify target vertex
  assert(
      target(p_cur_label->pred_edge, g) == p_cur_label->resident_vertex);
  // c - find a dominating label at current vertex
  vertex cur_vertex = source(p_cur_label->pred_edge, g);
  spp_res_cont cur_cont = p_cur_label->cumulated_resource_consumption;
  bool found_label = false, infeasible_once = false;
  for (Spplabel label : vec_vertex_labels[cur_vertex]) {
    spp_res_cont pred_cont = label->cumulated_resource_consumption,
        new_cont = pred_cont;
    // d - extend label: if feasible, check if new label dominate the
    // current one and then backtrack
    if (ref(g, &new_cont, pred_cont, p_cur_label->pred_edge)) {
      if (dominance(&new_cont, &cur_cont)) {
        assert(label->b_is_valid);
        assert(cur_vertex == label->resident_vertex);
        // backtrack
        backtrack(g,
                  vec_vertex_labels,
                  p_original_label,
                  label,
                  ref,
                  dominance,
                  opt_solutions_spp,
                  pareto_opt_rcs_spp,
                  path);
        found_label = true;
      }
    } else {
      infeasible_once = true;
    }
  }
  // e - verify that a new label has been found and is valid
  // The main reason a label could not be found is that the predecessor
  // label has been dominated and deleted, and the label residing at the
  // current vertex could not be expanded in a feasible way.
  // That could happen if the LB of a label is hard
  if (!found_label) {
    // if one expanded label was infeasible (UB reached or hard LB),
    // discard the path
    if (infeasible_once) {
//      std::cerr << "Discard: ";
//      cur_cont.print(std::cerr);
      return;
    } else {
//      Tools::throwError("Predecessor label could not be found when "
//                        "reconstructing solution path in the RCGraph.");
      std::cerr << "Predecessor label could not be found when "
                   "reconstructing solution path in the RCGraph." << std::endl;
    }
  }
}

// rc_spp_visitor struct: derived from boost::default_r_c_shortest_paths_visitor
rc_spp_visitor::rc_spp_visitor(
    int nMax,
    const std::vector<vertex> &sinks,
    std::function<void(spp_res_cont *)> post_process_rc,
    double maxReducedCostBound) :
    nMax_(nMax),
    sinks_(sinks),
    postProcessRc_(post_process_rc),
    maxReducedCostBound_(maxReducedCostBound) {}

void rc_spp_visitor::on_label_popped(Label &spplabel, const Graph &) {
  ++nPoppedLabels_;
}

void rc_spp_visitor::on_label_feasible(Label &spplabel, const Graph &) {
  ++nFeasibleLabels_;
}

void rc_spp_visitor::on_label_not_feasible(Label &, const Graph &) {
  ++nInfeasibleLabels_;
}

void rc_spp_visitor::on_label_dominated(Label &spplabel, const Graph &) {
  ++nDominatedLabels_;
}

void rc_spp_visitor::on_label_not_dominated(Label &spplabel, const Graph &) {
  ++nNotDominatedLabels_;
  // if seek all paths, do nothing
  if (nMax_ == -1) return;
  // check if we found a non-dominated path from a source to a sink
  if (find(sinks_.begin(), sinks_.end(), spplabel.resident_vertex)
      != sinks_.end()) {
    paths_[spplabel.num] = spplabel.cumulated_resource_consumption;
  }
}

template<typename Queue>
bool rc_spp_visitor::on_enter_loop(const Queue &queue, const Graph &graph) {
  // if go until optimality
  if (nMax_ == -1) return true;
  // if exceed the number of paths searched
  if (paths_.size() >= nMax_)
    return false;
  // stop if have explored too many nodes (avoid stalling)
  return nPoppedLabels_ < nMax_ * num_vertices(graph);
}

void rc_spp_visitor::printStats(std::ostream &out) const {
  out << "Labels statistics:\tpopped=" << nPoppedLabels_
      << "\tfeasible=" << nFeasibleLabels_
      << "\tinfeasible=" << nInfeasibleLabels_
      << "\tdominated=" << nDominatedLabels_
      << "\tnot dominated=" << nNotDominatedLabels_ << std::endl;
}
