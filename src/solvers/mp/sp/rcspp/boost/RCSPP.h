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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCSPP_H_
#define SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCSPP_H_

#include <map>
#include <memory>
#include <queue>
#include <vector>
#include <list>

#include "RCGraph.h"
#include "solvers/Solver.h"
#include "solvers/mp/MasterProblem.h"

#include "boost/graph/r_c_shortest_paths.hpp"

namespace boostRCSPP {

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// The following functions / data structures are used for the time resource:
//   Resource container ("resource" = cost + time)
//   Comparison override --> in SubProblem.cpp
//   Cost extension function
//   Dominance function
//
/////////////////////////////////////////////////////////////////////////////

// data structures for shortest path problem with labels
// ResourceContainer model
struct spp_res_cont {
  // Constructor
  spp_res_cont() : cost(XLARGE_SCORE), postprocessCost(0) {}

  spp_res_cont(double c, const std::vector<int> &label_values) :
      cost(c), postprocessCost(0), label_values(label_values) {}

  // Assign
  spp_res_cont &operator=(const spp_res_cont &other) {
    if (this == &other) return *this;
    this->~spp_res_cont();
    new(this) spp_res_cont(other);
    return *this;
  }

  void dominate(const spp_res_cont &res);

  // Current cost
  double cost;
  double postprocessCost;

  // Current labels
  std::vector<int> label_values;

  // dominance lvl: is always strictly higher than any dominated label.
  // the parentDominanceLvl is fixed at the creation and is used for the
  // comparison in the queue.
  // The objective of this level is to ensure that labels of dominated paths
  // are all dominated:
  // the dominant label will be expanded until it reaches the end of the
  // dominated path (before the dominated path will be expanded).
  // Indeed, its parentDominanceLvl will always be higher than the one of the
  // dominated path, and thus will be expanded before.
  int dominanceLvl = 0, parentDominanceLvl = 0;

  int first_day = -1, day = -1;
#ifdef NS_DEBUG
  int pred_arc = -1;
  std::vector<int> shifts_;
#endif

  int label_value(int l) const {
    return label_values.at(l);
  }

  int size() const {
    return label_values.size();
  }

  void print(std::ostream &out = std::cout) const;
};

// Resources extension model (arc has cost + label consumptions)
class ref_spp {
 public:
  explicit ref_spp(const std::vector<LABEL> labels) : labels_(labels) {}
  bool operator()(const Graph &g,
                  spp_res_cont *new_cont,
                  const spp_res_cont &old_cont,
                  const edge &ed) const;

  // map the position in the resource container with the label
  // It's needed to fetch resource consumptions for an edge and
  // the bounds at a vertex
  const std::vector<LABEL> labels_;
};

// Dominance function model
class dominance_spp {
 public:
  dominance_spp(std::vector<LABEL> labels,
                const Penalties &penalties,
                double epsilon) :
      epsilon_(epsilon) {
    for (LABEL l : labels)
      penalties_.addLabel(l, penalties);
  }

  virtual bool operator()(spp_res_cont *res_cont_1,
                          const spp_res_cont *res_cont_2) const;

 protected:
  // Penalties that fit the order and the size of the labels considered
  Penalties penalties_;
  double epsilon_;

  double worstPenalty(const spp_res_cont &res_cont_1,
                      const spp_res_cont &res_cont_2) const;

  double worstPenalty(int i, int level1, int level2) const;
};

// default dominance class: just compare resource values
class default_dominance_spp : public dominance_spp {
 public:
  default_dominance_spp(std::vector<LABEL> labels,
                         const Penalties &penalties,
                         double epsilon) :
                         dominance_spp(labels, penalties, epsilon) {}

  bool operator()(spp_res_cont *res_cont_1,
                  const spp_res_cont *res_cont_2) const override;
};

// Comparator for the priority queue used by boost to process the labels.
// Note that the Compare parameter of a priority queue is defined such that
// it returns true if its first argument comes before its second argument in a
// weak ordering.
// But because the priority queue outputs largest elements first,
// the elements that "come before" are actually output last.
// That is, the front of the queue contains the "last" element according to
// the weak ordering imposed by Compare.
// r_c_shortest_paths_label struct
struct Label {
  Label(const int64_t n,
        const spp_res_cont& rc,
        const shared_ptr<Label> &pl,
        const typename boost::graph_traits< Graph >::edge_descriptor& ed,
        const typename boost::graph_traits< Graph >::vertex_descriptor& vd)
      : num(n)
      , cumulated_resource_consumption(rc)
      , p_pred_label(pl)
      , pred_edge(ed)
      , resident_vertex(vd)
      , b_is_dominated(false)
      , b_is_processed(false) {}

  Label& operator=(const Label& other) {
    if (this == &other)
      return *this;
    this->~Label();
    new (this) Label(other);
    return *this;
  }

  const int64_t num;
  spp_res_cont cumulated_resource_consumption;
  const shared_ptr<Label> p_pred_label;
  const typename
  boost::graph_traits< Graph >::edge_descriptor pred_edge;
  const typename
  boost::graph_traits< Graph >::vertex_descriptor resident_vertex;
  bool b_is_dominated;
  bool b_is_processed;
};  // r_c_shortest_paths_label

typedef shared_ptr<Label> Spplabel;

// Breath First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the smallest day (return true if day1 > day2)
// 2. break tie with cost (return true if cost1 > cost2)
struct SpplabelBreadthFirstComparator {
  bool operator()(const Spplabel &splabel1, const Spplabel &splabel2) const;
};

// Depth First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the biggest day (day1 < day2)
// 2. break tie with cost (return true if cost1 > cost2)
struct SpplabelDepthFirstComparator {
  bool operator()(const Spplabel &splabel1, const Spplabel &splabel2) const;
};

// Best First Comparator to order the processing of the nodes in boost rc spp:
// 1. starts with the ones with the smallest cost (return true if cost1 > cost2)
struct SpplabelBestFirstComparator {
  bool operator()(const Spplabel &splabel1, const Spplabel &splabel2) const;
};

// Dominant First Comparator to order the processing of the nodes in
// boost rc spp:
// 1. starts with the ones with the most parent dominant level (order the label
//    by dominance)
// 2. starts with the ones with the smallest cost (return true if cost1 > cost2)
struct SpplabelDominantFirstComparator {
  bool operator()(const Spplabel &splabel1, const Spplabel &splabel2) const;
};

// use default operator<
struct SpplabelDefaultComparator {
  bool operator()(const Spplabel &splabel1, const Spplabel &splabel2) const;
};

// rc_spp_visitor struct: derived from boost::default_r_c_shortest_paths_visitor
struct rc_spp_visitor {
  rc_spp_visitor(int nMax,
                 const std::vector<vertex> &sinks,
                 std::function<void(spp_res_cont *)> post_process_rc = nullptr,
                 double maxReducedCostBound = -EPSILON);

  void on_label_popped(Label &, const Graph &);

  void on_label_feasible(Label &, const Graph &);

  void on_label_not_feasible(Label &, const Graph &);

  void on_label_dominated(Label &, const Graph &);

  void on_label_not_dominated(Label &, const Graph &);

  template<typename Queue>
  bool on_enter_loop(const Queue &queue, const Graph &graph);

  const std::map<int, spp_res_cont> &paths() const {
    return paths_;
  }

  void printStats(std::ostream &out) const;

 private:
  int nMax_;
  std::vector<vertex> sinks_;
  std::map<int, spp_res_cont> paths_;
  std::function<void(spp_res_cont *)> postProcessRc_;
  double maxReducedCostBound_;
  int nPoppedLabels_ = 0, nFeasibleLabels_ = 0, nInfeasibleLabels_ = 0,
      nDominatedLabels_ = 0, nNotDominatedLabels_ = 0;
};

// Solver
class BoostRCSPPSolver {
 public:
  BoostRCSPPSolver(RCGraph *rcg,
                   double maxReducedCostBound,
                   int verbose,
                   double epsilon,
                   SPSearchStrategy strategy,
                   int nb_max_paths,
                   std::function<void(spp_res_cont *)> post_process_rc);

  std::vector<RCSolution> solve(std::vector<LABEL> labels,
                                const Penalties &penalties,
                                std::vector<vertex> sinks,
                                const PScenario &pScenario);

 protected:
  RCGraph *rcg_;
  int verbose_;
  double maxReducedCostBound_, epsilon_;
  SPSearchStrategy strategy_;
  int nb_max_paths_;
  std::function<void(spp_res_cont *)> post_process_rc_;
  Tools::Timer timer_;

  void printPath(std::ostream &out,
                 const std::vector<edge> &path,
                 const spp_res_cont &resource) const;

  bool processPath(
      std::vector<edge> *path,
      spp_res_cont *rc,
      const ref_spp &ref,
      std::function<void(spp_res_cont *)> post_process_rc = nullptr,
      bool printBadPath = true) const;

  RCSolution solution(const std::vector<edge> &path,
                      const spp_res_cont &resource,
                      const PScenario &pScenario) const;

  // modified boost::check_r_c_path function
  bool check_r_c_path(
      std::vector<edge> ed_vec_path,
      // if true, computed accumulated final resource levels must
      // be equal to desired_final_resource_levels
      // if false, computed accumulated final resource levels must
      // be less than or equal to desired_final_resource_levels
      const spp_res_cont &found_final_resource_levels,
      spp_res_cont *checked_final_resource_levels,
      const ref_spp &ref) const;

  // r_c_shortest_path class: derived from boost::r_c_shortest_paths_dispatch
  template<class Visitor>
  void r_c_shortest_paths_solve(const Graph &g,
                                vertex s,
                                const std::vector<vertex> &t,
                                vector2D<edge> *opt_solutions_spp,
                                std::vector<spp_res_cont> *pareto_opt_rcs_spp,
                                const spp_res_cont &rc,
                                const ref_spp &ref,
                                const dominance_spp &dominance,
                                Visitor *vis,
                                SPSearchStrategy strategy) {
    // solve with the right comparator
    switch (strategy) {
      case SP_BEST_FIRST:
        r_c_shortest_paths_dispatch(g, s, t,
                                    opt_solutions_spp, pareto_opt_rcs_spp,
                                    rc, ref, dominance, vis,
                                    SpplabelBestFirstComparator());
        break;
      case SP_BREADTH_FIRST:
        r_c_shortest_paths_dispatch(g, s, t,
                                    opt_solutions_spp, pareto_opt_rcs_spp,
                                    rc, ref, dominance, vis,
                                    SpplabelBreadthFirstComparator());
        break;
      case SP_DEPTH_FIRST:
        r_c_shortest_paths_dispatch(g, s, t,
                                    opt_solutions_spp, pareto_opt_rcs_spp,
                                    rc, ref, dominance, vis,
                                    SpplabelDepthFirstComparator());
        break;
      case SP_DOMINANT_FIRST:
        r_c_shortest_paths_dispatch(g, s, t,
                                    opt_solutions_spp, pareto_opt_rcs_spp,
                                    rc, ref, dominance, vis,
                                    SpplabelDominantFirstComparator());
        break;
      default:
        r_c_shortest_paths_dispatch(g, s, t,
                                    opt_solutions_spp, pareto_opt_rcs_spp,
                                    rc, ref, dominance, vis,
                                    SpplabelDefaultComparator());
    }
  }

// modified boost::r_c_shortest_paths_dispatch function (body/implementation)
  template<class Visitor,
      class Spplabel_Comparator>
  void r_c_shortest_paths_dispatch(
      const Graph &g,
      vertex s,
      const std::vector<vertex> &t,
      vector2D<edge> *opt_solutions_spp,
      std::vector<spp_res_cont> *pareto_opt_rcs_spp,
      const spp_res_cont &rc,
      const ref_spp &ref,
      const dominance_spp &dominance,
      Visitor *vis,
      Spplabel_Comparator) {
    // definition of the label and the priority queue that stores the labels
    // left to process
    std::priority_queue<Spplabel, std::vector<Spplabel>, Spplabel_Comparator>
        unprocessed_labels;

    // Allocator for the Splabel
    size_t i_label_num = 0;

    // build first label

    Spplabel spplabel_first_label = std::make_shared<Label>(
        i_label_num++, rc, nullptr, edge(), s);
    unprocessed_labels.push(spplabel_first_label);

    // build a vector to store a list of labels residing at a given vertex
    std::vector<std::list<Spplabel> > vec_vertex_labels(num_vertices(g));
    vec_vertex_labels[s].push_back(spplabel_first_label);

    // vector of iterator: store the last valid position in the list for
    // dominance:
    // every elements before the iterator cannot dominate each other,
    // the ones after haven't been tested for the moment
    typedef std::vector<std::list<Spplabel>::iterator>
        vec_last_pos_for_dominance_type;
    vec_last_pos_for_dominance_type vec_last_pos_for_dominance(num_vertices(g));
    // store in the map an iterator starting at the beginning of the vector of
    // labels for each vertex
    for (unsigned int v = 0; v < vec_vertex_labels.size(); ++v)
      vec_last_pos_for_dominance[v] = vec_vertex_labels[v].begin();

    /**
     * Dominance helper function: if spplabel1 dominates spplabel2,
     * erase it from the list, delete spplabel2 or mark it as dominated,
     * and move it2 to the next label.
     * return true if dominated.
     */
    auto dominance_func = [&dominance, &vis, &g](
        Spplabel &spplabel1, Spplabel &spplabel2,
        std::list<Spplabel>::iterator &it2,
        std::list<Spplabel> &list_labels_cur_vertex) {
      if (dominance(&spplabel1->cumulated_resource_consumption,
                    &spplabel2->cumulated_resource_consumption)) {
        it2 = list_labels_cur_vertex.erase(it2);
        if (spplabel2->b_is_processed) {
          vis->on_label_dominated(*spplabel2, g);
        } else {
          // mark as dominated: will be deleted when popped of the
          // unprocessed list
          // if current_label, it will be deleted latter on.
          spplabel2->b_is_dominated = true;
        }
        return true;
      }
      return false;
    };

    /**
      * MAIN LOOP THAT RUNS UNTIL ALL LABELS HAS BEEN PROCESSED OR
      * VIS.ON_ENTER_LOOP RETURNS FALSE
      */
    while (!unprocessed_labels.empty()
        && vis->on_enter_loop(unprocessed_labels, g)) {
      /**
        * 1 - Pop cur_label from the priority queue unprocessed_labels
        * (choose the one with the higher priority)
        */
      Spplabel cur_label = unprocessed_labels.top();
      unprocessed_labels.pop();
      vis->on_label_popped(*cur_label, g);

      // an Spplabel object in unprocessed_labels and the respective Spplabel
      // object in the respective list<Spplabel> of vec_vertex_labels share
      // their embedded r_c_shortest_paths_label object to avoid memory leaks.
      // dominated r_c_shortest_paths_label objects are marked and deleted
      // when popped from unprocessed_labels, as they can no longer be deleted
      // at the end of the function;
      // only the Spplabel object in unprocessed_labels still references
      // the r_c_shortest_paths_label object.
      // this is also for efficiency, because the else branch is executed only
      // if there is a chance that extending the
      // label leads to new undominated labels, which in turn is possible only
      // if the label to be extended is undominated

      /**
        * 2.a - If cur_label is non dominated,
        * check with dominance  all the label residing at
        * cur_label->resident_vertex
        */
      if (!cur_label->b_is_dominated) {
        // current resident vertex
        auto i_cur_resident_vertex = cur_label->resident_vertex;
        // list of labels at current vertex
        std::list<Spplabel>
            &list_labels_cur_vertex = vec_vertex_labels[i_cur_resident_vertex];
        auto last_pos_for_dominance =
            vec_last_pos_for_dominance[i_cur_resident_vertex];

        // if more than just current label at current resident vertex
        // and last_pos_for_dominance is not the last element of the list
        // (do not rechecked labels that were already at the vertex at the
        // previous check)
        // otherwise, go to 2.b
        if (list_labels_cur_vertex.size() > 1
            && last_pos_for_dominance != --list_labels_cur_vertex.end()) {
          // Two iterators will be used to check the dominance:
          // - it1 will move from the first label to the last one
          // - it2 will move from * to the last one:
          //    *: the next label after it1 if it1 is before last_position;
          //    *: just after the last label checked (will be the first one if
          //       first time).

          // start at beginning for it1
          auto it1 = list_labels_cur_vertex.begin();
          bool it1_at_or_beyond_last_pos_for_dominance = false;
          // if last_pos_for_dominance has not been used for the moment,
          // it is set to the end [a list returns the end() for begin() when
          // empty] -> set it to begin()
          if (last_pos_for_dominance == list_labels_cur_vertex.end())
            last_pos_for_dominance = list_labels_cur_vertex.begin();

          // while begin_it hasn't reached the end, compare for dominance with
          // the other labels (end_it)
          while (it1 != list_labels_cur_vertex.end()) {
            // fetch spplabel corresponding to the iterator
            Spplabel spplabel1 = *it1;

            // if it1 is exactly at the last position -> mark the iterator as
            // passed the last position
            if (!it1_at_or_beyond_last_pos_for_dominance
                && it1 == last_pos_for_dominance)
              it1_at_or_beyond_last_pos_for_dominance = true;

            // initialize second iterator. Will be used to compare labels to
            // begin_iter
            auto it2 = it1;

            // if not first time and it1 before last position -> set it2 to the
            // last position checked
            if (!it1_at_or_beyond_last_pos_for_dominance)
              it2 = last_pos_for_dominance;
            // just move it2 to the next position:
            // one position after last or it1
            ++it2;

            // check dominance for current position of it1 while non dominated
            bool it1_erased = false;
            while (it2 != list_labels_cur_vertex.end()) {
              Spplabel spplabel2 = *it2;

              // a - if spplabel1 dominates spplabel2
              // -> erase it, move to the next position, and continue current
              //    loop on it2
              if (dominance_func(spplabel1,
                                 spplabel2,
                                 it2,
                                 list_labels_cur_vertex))
                continue;
              else
                // otherwise move it2 to the next label
                ++it2;

              // b - if spplabel2 dominates spplabel1
              // -> erase it, move to the next position, and break current
              //    loop on it2
              if (dominance_func(spplabel2,
                                 spplabel1,
                                 it1,
                                 list_labels_cur_vertex)) {
                it1_erased = true;
                break;
              }
            }
            // if iterator not erased (i.e., already moved to the next label),
            // move to next label
            if (!it1_erased)
              ++it1;
          }  // end while loop for dominance

          // set the iterator at the last element of the list
          vec_last_pos_for_dominance[i_cur_resident_vertex] =
              --list_labels_cur_vertex.end();
        }
      }

      /**
        * 2.b - If cur_label is still non dominated, expand the label on all
        * successors
        */
      if (!cur_label->b_is_dominated) {
        // mark as processed and call vis
        cur_label->b_is_processed = true;
        vis->on_label_not_dominated(*cur_label, g);

        // iterate on all out going edges from the current resident vertex
        vertex cur_vertex = cur_label->resident_vertex;
        boost::graph_traits<Graph>::out_edge_iterator oei, oei_end;
        for (boost::tie(oei, oei_end) = out_edges(cur_vertex, g);
             oei != oei_end;
             ++oei) {
          // allocate a new label
          Spplabel new_sp_label = std::make_shared<Label>(
              i_label_num++,
              cur_label->cumulated_resource_consumption,
              cur_label,
              *oei,
              target(*oei, g));

          // populate the label and check its feasibility
          bool b_feasible =
              ref(g,
                  &new_sp_label->cumulated_resource_consumption,
                  new_sp_label->p_pred_label->cumulated_resource_consumption,
                  new_sp_label->pred_edge);

          // if infeasible, destroy it
          if (!b_feasible) {
            vis->on_label_not_feasible(*new_sp_label, g);
          } else {
            // otherwise, add it to the unprocessed priority queue
            vis->on_label_feasible(*new_sp_label, g);
            vec_vertex_labels[new_sp_label->resident_vertex].push_back(
                new_sp_label);
            unprocessed_labels.push(new_sp_label);
          }
        }
      } else {
        /**
        * 2.c - If cur_label is dominated, delete the label
        */
        vis->on_label_dominated(*cur_label, g);
      }
    }

    /**
      * 3 - The main loop has ended, retrieve all the paths that have reached
      *     one of the sink in t.
      *     The code has been modified as the boost version implied that all
      *     predecessors labels must be non-dominated, which is not the case it
      *     the process is stopped before the end.
      *     Instead, we reconstruct a non-dominated path, based on the
      *     non-dominated labels residing at each certain vertex.
      *     The resulting label at the sink will dominate the one found and could
      *     have a different path.
      */
    for (int sink : t) {
      const std::list<Spplabel> &dspplabels = vec_vertex_labels[sink];
      // if the sink could be reached from the source
      for (auto csi = dspplabels.begin(); csi != dspplabels.end(); ++csi)
        backtrack(g,
                  vec_vertex_labels,
                  *csi,
                  *csi,
                  ref,
                  dominance,
                  opt_solutions_spp,
                  pareto_opt_rcs_spp,
                  std::vector<edge>());
    }
    // --------------------------------------------------- END MODIFICATION
  }  // r_c_shortest_paths_dispatch

  void backtrack(const Graph &g,
                 const std::vector<std::list<Spplabel> > &vec_vertex_labels,
                 const Spplabel &p_original_label,
                 const Spplabel &p_cur_label,
                 const ref_spp &ref,
                 const dominance_spp &dominance,
                 vector2D<edge> *opt_solutions_spp,
                 std::vector<spp_res_cont> *pareto_opt_rcs_spp,
                 std::vector<edge> path) const;
};

}  // namespace boostRCSPP

#endif  // SRC_SOLVERS_MP_SP_RCSPP_BOOST_RCSPP_H_
