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

#ifndef SRC_SOLVERS_MP_MODELER_BCPBRANCHINGCANDIDATES_H_
#define SRC_SOLVERS_MP_MODELER_BCPBRANCHINGCANDIDATES_H_

#include <cstdio>
#include <list>
#include <numeric>
#include <utility>
#include <vector>

#include "solvers/mp/modeler/Modeler.h"

#include "CoinWarmStart.hpp"  // NOLINT (suppress cpplint error)
#include "CoinTime.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_math.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_enum.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_matrix.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_warmstart.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_result.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_node.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_user.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_functions.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_pool.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_branch.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp.hpp"  // NOLINT (suppress cpplint error)

class BcpModeler;

class BcpBranchingCandidates {
 public:
  explicit BcpBranchingCandidates(BcpModeler *pModel);

  BCP_lp_branching_object* selectCandidates(
      const std::vector<MyPBranchingCandidate> &candidates,
      BCP_lp_prob *p,
      // the variables in the current formulation.
      const BCP_vec<BCP_var *> &vars,
      // the cuts in the current formulation.
      const BCP_vec<BCP_cut *> &cuts,
      int lpIt);  // current lp iteration

  void updateTree(BCP_presolved_lp_brobj *best);

  void reset(double objValue);

  // build the candidate from my candidate
  BCP_lp_branching_object* buildCandidate(
      const MyBranchingCandidate &candidate,
                      const BCP_vec<BCP_var *> &vars,
                      const BCP_vec<BCP_cut *> &cuts,
                      int *nbNewVar,
                      int *nbNewCut);

  std::pair<int, int> BCP_add_branching_objects(
      BCP_lp_prob *p,
      OsiSolverInterface *lp,
      const BCP_vec<BCP_lp_branching_object*>& candidates);

  bool hasCandidate() const { return candidate_ != nullptr; }

  bool candidateIncludeVariable(MyVar *pVar) const;

  void
  BCP_mark_result_of_strong_branching(BCP_lp_prob *p,
                                      const BCP_lp_branching_object* can,
                                      const int added_col_num,
                                      const int added_row_num);

  MyPBranchingCandidate BCP_lp_perform_strong_branching(
      BCP_lp_prob *p,
      const std::vector<MyPBranchingCandidate> &candidates);

  BCP_branching_object_relation compare_branching_candidates(
      BCP_presolved_lp_brobj* new_solved,
      BCP_presolved_lp_brobj* old_solved);

 protected:
  BcpModeler *pModel_;
  // store best expected increase
  double min_expected_obj_;
  MyPBranchingCandidate candidate_;
  std::recursive_mutex mutex_;
};



#endif  // SRC_SOLVERS_MP_MODELER_BCPBRANCHINGCANDIDATES_H_
