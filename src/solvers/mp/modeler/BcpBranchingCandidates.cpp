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

#include "BcpBranchingCandidates.h"

#include <algorithm>
#include <list>
#include <map>

#include "solvers/mp/modeler/BcpModeler.h"

// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

//#############################################################################
BcpBranchingCandidates::BcpBranchingCandidates(BcpModeler *pModel):
    pModel_(pModel), min_expected_obj_(1e10), candidate_(nullptr) {}


BCP_lp_branching_object * BcpBranchingCandidates::selectCandidates(
    const std::list<MyPBranchingCandidate> &candidates,
    BCP_lp_prob *p,
    // the variables in the current formulation.
    const BCP_vec<BCP_var *> &vars,
    // the cuts in the current formulation.
    const BCP_vec<BCP_cut *> &cuts) {
  int nbNewVar = 0, nbNewCut = 0;  // counters for new columns and cuts
  BCP_vec<BCP_lp_branching_object*> cands;
  std::map<BCP_lp_branching_object*, MyPBranchingCandidate> candidatesMap;
  for (const auto &candidate : candidates) {
    BCP_lp_branching_object *can =
        buildCandidate(*candidate, vars, cuts, &nbNewVar, &nbNewCut);
    cands.push_back(can);
    candidatesMap[can] = candidate;
  }

  // perform strong branching if more than one candidate
  if (cands.size() <= 1) {
    candidate_ = candidatesMap.at(cands[0]);
    return cands[0];
  }

  // perform strong branching if more than one candidate
  BCP_lp_branching_object* best_candidate =
      BCP_lp_perform_strong_branching(p, cands);
  candidate_ = candidatesMap.at(best_candidate);
  return best_candidate;
}

void BcpBranchingCandidates::updateTree(BCP_presolved_lp_brobj *best) {
  for (int n=0; n < best->candidate()->child_num; ++n) {
    const MyPBranchingNode & pNode = candidate_->getChild(n);
    pModel_->addNode(pNode, best->lpres(n).objval());
  }
  candidate_ = nullptr;
}

void BcpBranchingCandidates::reset(double objValue) {
  min_expected_obj_ = objValue;
  candidate_ = nullptr;
}

BCP_lp_branching_object * BcpBranchingCandidates::buildCandidate(
    const MyBranchingCandidate &candidate,
    const BCP_vec<BCP_var *> &vars,
    const BCP_vec<BCP_cut *> &cuts,
    int *nbNewVar,
    int *nbNewCut) {
  BCP_vec<BCP_var *> new_vars;  // add a branching cut for the set of arcs
  BCP_vec<BCP_cut *> new_cuts;  // add a branching cut for the set of arcs
  BCP_vec<int> vpos;  // positions of the variables
  BCP_vec<double> vbd;  // new bounds for each variable and for each children
  BCP_vec<int> cpos;  // positions of the cuts
  BCP_vec<double> cbd;  // bounds of the cuts

  /*
   * Branching variables
   */

  BCP_vec<double> generalVarBounds;
  for (MyVar *var : candidate.getBranchingVars()) {
    // search the var in the vars
    int index = pModel_->getActiveIndex(var);
    if (index == -1) {  // new var
      index = vars.size() + *nbNewVar;
      ++*nbNewVar;
    }
    vpos.push_back(index);
  }

  for (MyVar *newVar : candidate.getNewBranchingVars())
    new_vars.push_back(dynamic_cast<BCP_var *>(newVar));

  // bounds
  for (const MyPBranchingNode &node : candidate.getChildren()) {
    auto lbIt = node->getLb().begin();
    for (auto ubIt = node->getUb().begin(); ubIt != node->getUb().end();
         ++ubIt, ++lbIt) {
      vbd.push_back(*lbIt);
      vbd.push_back(*ubIt);
    }
  }

  /*
   * Branching cuts
   */
  for (MyCons *cons : candidate.getBranchingCons()) {
    // search the var in the vars
    auto *cut = dynamic_cast<BcpBranchCons *>(cons);
    int index = cons->getIndex();
    if (cut) {
      index = cuts.size() + *nbNewCut;
      ++*nbNewCut;
    }
    cpos.push_back(index);
  }

  for (MyCons *newCut : candidate.getNewBranchingCons())
    new_cuts.push_back(dynamic_cast<BCP_cut *>(newCut));

  // bounds
  for (const MyPBranchingNode &node : candidate.getChildren()) {
    auto lhsIt = node->getLhs().begin();
    for (auto rhsIt = node->getRhs().begin(); rhsIt != node->getRhs().end();
         ++rhsIt) {
      cbd.push_back(*lhsIt);
      cbd.push_back(*rhsIt);
      ++lhsIt;
    }
  }

  return new BCP_lp_branching_object(
      candidate.getChildren().size(),
      new_vars.empty() ? nullptr : &new_vars,
      new_cuts.empty() ? nullptr : &new_cuts, /* vars/cuts_to_add */
      vpos.empty() ? nullptr : &vpos,
      cpos.empty() ? nullptr : &cpos,
      vbd.empty() ? nullptr : &vbd,
      cbd.empty() ? nullptr : &cbd, /* forced parts */
      nullptr,
      nullptr,
      nullptr,
      nullptr /* implied parts */);
}

// Decide which branching object is preferred for branching.
// Based on the member fields of the two presolved candidate branching
// objects decide which one should be preferred for really branching on it.
// Possible return values are: BCP_OldPresolvedIsBetter,
// BCP_NewPresolvedIsBetter and BCP_NewPresolvedIsBetter_BranchOnIt.
// This last value (besides specifying which candidate is preferred) also
// indicates that no further candidates should be examined, branching
// should be done on this candidate.
//
// Default: The behavior of this method is governed by the
// BranchingObjectComparison parameter in BCP_lp_par.
BCP_branching_object_relation
BcpBranchingCandidates::compare_branching_candidates(
    BCP_presolved_lp_brobj* new_solved,
    BCP_presolved_lp_brobj* old_solved) {
  // if a column node, do not need to do strong branching
  int num = new_solved->candidate()->child_num;
  if (num <= 1)
    return BCP_NewPresolvedIsBetter;

  double min_incr_new = DBL_MAX, avg = 0;
  if (pModel_->getVerbosity() > 1) std::cout << "Strong branching results:";
  // do not take into account the first node when 3 children as it's a
  // column node shared by all candidates
  int min_num = num < 3 ? 0 : 1;
  for (int n = min_num; n < num; n++) {
    double o = new_solved->lpres(n).objval();
    avg += o;
    if (pModel_->getVerbosity() > 1) std::cout << " " << o;
    if (o < min_incr_new - pModel_->epsilon())
      min_incr_new = o;
  }
  avg /= (num - min_num);
  min_incr_new += avg * pModel_->epsilon();  // add a small part to break tie
  if (pModel_->getVerbosity() > 1) std::cout << std::endl;

  std::lock_guard<std::recursive_mutex> lock(mutex_);
  // if better candidate
  if (min_incr_new > min_expected_obj_ - pModel_->epsilon()) {
    min_expected_obj_ = min_incr_new;
    return BCP_NewPresolvedIsBetter;
  }
  // otherwise keep previous one
  return BCP_OldPresolvedIsBetter;
}

std::pair<int, int> BcpBranchingCandidates::BCP_add_branching_objects(
    BCP_lp_prob *p,
    const BCP_vec<BCP_lp_branching_object*> &candidates) {
  // to make things short
  if (candidates.empty())
    return {0, 0};

  BCP_lp_branching_object* can;
  BCP_var_set& vars = p->node->vars;
  BCP_cut_set& cuts = p->node->cuts;

  // first count the number of cols/rows to add
  int newvar_num = 0;
  int newcut_num = 0;
  for (auto cani = candidates.begin(); cani != candidates.end(); ++cani) {
    can = *cani;
    can->init_pos_for_added(vars.size() + newvar_num,
                            cuts.size() + newcut_num);
    if (can->vars_to_add)
      newvar_num += can->vars_to_add->size();
    if (can->cuts_to_add)
      newcut_num += can->cuts_to_add->size();
  }

  const int orig_col_num = vars.size();
  const int orig_row_num = cuts.size();

  OsiSolverInterface* lp = p->lp_solver;

  // deal with the vars
  if (newvar_num > 0) {
    BCP_vec<BCP_var*> new_vars;
    new_vars.reserve(newvar_num);
    for (auto cani = candidates.begin(); cani != candidates.end(); ++cani) {
      can = *cani;
      if (can->vars_to_add)
        new_vars.append(*can->vars_to_add);
    }

    BCP_vec<BCP_col*> cols;
    cols.reserve(newvar_num);
    p->user->vars_to_cols(cuts, new_vars, cols,
                         *p->lp_result, BCP_Object_Branching, false);
    BCP_lp_add_cols_to_lp(cols, lp);
    purge_ptr_vector(cols);

    for (int i = 0; i < newvar_num; ++i) {
      new_vars[i]->set_bcpind(-BCP_lp_next_var_index(*p));
    }
    vars.append(new_vars);
  }

  // now add the rows
  if (newcut_num > 0) {
    BCP_vec<BCP_cut*> new_cuts;
    new_cuts.reserve(newcut_num);
    for (auto cani = candidates.begin(); cani != candidates.end(); ++cani) {
      can = *cani;
      if (can->cuts_to_add)
        new_cuts.append(*can->cuts_to_add);
    }

    BCP_vec<BCP_row*> rows;
    rows.reserve(newcut_num);
    BCP_fatal_error::abort_on_error = false;
    try {
      p->user->cuts_to_rows(vars, new_cuts, rows,
                           *p->lp_result, BCP_Object_Branching, false);
    } catch (...) {
    }
    BCP_fatal_error::abort_on_error = true;
    BCP_lp_add_rows_to_lp(rows, lp);
    purge_ptr_vector(rows);

    for (int i = 0; i < newcut_num; ++i) {
      new_cuts[i]->set_bcpind(-BCP_lp_next_cut_index(*p));
    }
    cuts.append(new_cuts);
    p->node->lb_at_cutgen.insert(p->node->lb_at_cutgen.end(), newcut_num,
                                p->lp_result->objval());
  }

  // mark the newly added vars as fixed to 0, and the newly added cuts as
  // free. (var_indices and cut_indices simply contain the indices of the
  // added vars/cuts.)

  if (newvar_num > 0) {
    for (int i = orig_col_num; i < orig_col_num + newvar_num; ++i)
      lp->setColBounds(i, 0.0, 0.0);
  }
  if (newcut_num > 0) {
    const double inf = lp->getInfinity();
    for (int i = orig_row_num; i < orig_row_num + newcut_num; ++i)
      lp->setRowBounds(i, -inf, inf);
  }

  return {newvar_num, newcut_num};
}

//#############################################################################

void BcpBranchingCandidates::BCP_mark_result_of_strong_branching(
    BCP_lp_prob* p,
    const BCP_lp_branching_object* can,
    const int added_col_num,
    const int added_row_num) {
  BCP_var_set& vars = p->node->vars;
  if (can->forced_var_pos) {
    BCP_vec<int>::const_iterator ii = can->forced_var_pos->begin();
    BCP_vec<int>::const_iterator lastii = can->forced_var_pos->end();
    for ( ; ii != lastii; ++ii)
      vars[*ii]->make_non_removable();
  }
  if (can->implied_var_pos) {
    // just to make sure that these vars are not deleted before the implied
    // bound change takes place...
    BCP_vec<int>::const_iterator ii = can->implied_var_pos->begin();
    BCP_vec<int>::const_iterator lastii = can->implied_var_pos->end();
    for ( ; ii != lastii; ++ii)
      vars[*ii]->make_active();
  }

  if (added_col_num) {
    BCP_var_set::iterator vari = vars.entry(vars.size() - added_col_num);
    BCP_var_set::const_iterator lastvari = vars.end();
    for ( ; vari != lastvari; ++vari) {
      if (!(*vari)->is_non_removable())
        (*vari)->make_to_be_removed();
    }
  }

  BCP_cut_set& cuts = p->node->cuts;
  if (can->forced_cut_pos) {
    BCP_vec<int>::const_iterator ii = can->forced_cut_pos->begin();
    BCP_vec<int>::const_iterator lastii = can->forced_cut_pos->end();
    for ( ; ii != lastii; ++ii)
      cuts[*ii]->make_non_removable();
  }
  if (can->implied_cut_pos) {
    // just to make sure that these cuts are not deleted before the implied
    // bound change takes place...
    BCP_vec<int>::const_iterator ii = can->implied_cut_pos->begin();
    BCP_vec<int>::const_iterator lastii = can->implied_cut_pos->end();
    for ( ; ii != lastii; ++ii)
      cuts[*ii]->make_active();
  }
  if (added_row_num > 0) {
    BCP_cut_set::iterator cuti = cuts.entry(cuts.size() - added_row_num);
    BCP_cut_set::const_iterator lastcuti = cuts.end();
    for ( ; cuti != lastcuti; ++cuti) {
      if (!(*cuti)->is_non_removable())
        (*cuti)->make_to_be_removed();
    }
  }
}

//#############################################################################

BCP_lp_branching_object*
BcpBranchingCandidates::BCP_lp_perform_strong_branching(
    BCP_lp_prob *p,
    const BCP_vec<BCP_lp_branching_object*>& candidates) {
  BCP_var_set& vars = p->node->vars;
  BCP_cut_set& cuts = p->node->cuts;

  const int orig_colnum = vars.size();

  const std::pair<int, int> added_object_num =
      BCP_add_branching_objects(p, candidates);

  const int added_colnum = added_object_num.first;
  const int added_rownum = added_object_num.second;

  const int colnum = vars.size();
  const int rownum = cuts.size();

  int i, j;  // loop variable

  const CoinWarmStart * ws = p->lp_solver->getWarmStart();

  // prepare for strong branching
  p->lp_solver->markHotStart();

  // save the lower/upper bounds of every var/cut
  BCP_vec<double> rowbounds(2 * rownum, 0.0);
  BCP_vec<double> colbounds(2 * colnum, 0.0);

  const int maxind = std::max<int>(rownum, colnum);
  BCP_vec<int> all_indices(maxind, 0);
  for (i = 0; i < maxind; ++i)
    all_indices[i] = i;

  const double * rlb_orig = p->lp_solver->getRowLower();
  const double * rub_orig = p->lp_solver->getRowUpper();
  for (j = -1, i = 0; i < rownum; ++i) {
    rowbounds[++j] = rlb_orig[i];
    rowbounds[++j] = rub_orig[i];
  }

  const double * clb_orig = p->lp_solver->getColLower();
  const double * cub_orig = p->lp_solver->getColUpper();
  for (j = -1, i = 0; i < colnum; ++i) {
    colbounds[++j] = clb_orig[i];
    colbounds[++j] = cub_orig[i];
  }

  // local thread pool (use all available threads)
  Tools::ThreadsPool pool;
  std::list<OsiSolverInterface*> lps;

  // Look at the candidates one-by-one and presolve them.
  p->user->print(p->param(BCP_lp_par::LpVerb_StrongBranchResult),
                 "\nLP: Starting strong branching:\n\n");
  const OsiBabSolver* babSolver = p->user->getOsiBabSolver();
  BCP_presolved_lp_brobj *best_presolved = nullptr;
  BCP_vec<BCP_lp_branching_object*>::iterator cani;
  for (auto cani = candidates.begin(); cani != candidates.end(); ++cani) {
   /**
   * Start of the job definition (part run in parallel)
   */
    Tools::Job job = [&, cani, this]() {
      // fetch an LP solver
      OsiSolverInterface* lp = nullptr;
      std::unique_lock<std::recursive_mutex> l(mutex_);
      if (!lps.empty()) {
        lp = lps.back();
        lps.pop_back();
      }
      l.unlock();
      // clone the lp if none is available
      if (lp == nullptr) {
        lp = p->lp_solver->clone();
        // prepare for strong branching
        lp->setWarmStart(ws);
        lp->markHotStart();
      }

      // Create a temporary branching object to hold the current results
      BCP_lp_branching_object *can = *cani;
      auto * tmp_presolved = new BCP_presolved_lp_brobj(can);
      for (int i = 0; i < can->child_num; ++i) {
        can->apply_child_bd(lp, i);
        // bound changes always imply that primal feasibility is lost.
        p->user->modify_lp_parameters(lp, 1, true);
        lp->solveFromHotStart();
        tmp_presolved->get_results(*lp, i);
        l.lock();
        BCP_lp_test_feasibility(*p, tmp_presolved->lpres(i));
        l.unlock();
        // reset the bounds of the affected vars/cuts
        if (can->cuts_affected() > 0)
          lp->setRowSetBounds(all_indices.begin(), all_indices.entry(rownum),
                              rowbounds.begin());
        if (can->vars_affected() > 0)
          lp->setColSetBounds(all_indices.begin(), all_indices.entry(colnum),
                              colbounds.begin());

//          if (p->param(BCP_lp_par::LpVerb_PresolveResult)) {
//            p->user->print(true, "LP:   Presolving:");
//            if (p->param(BCP_lp_par::LpVerb_PresolvePositions)) {
//              can->print_branching_info(orig_colnum,
//                                        p->lp_result->x(),
//                                        p->lp_solver->getObjCoefficients());
//            }
//            for (i = 0; i < can->child_num; ++i) {
//              const BCP_lp_result& res = tmp_presolved->lpres(i);
//              const double lb = res.objval();
//              p->user->print(true,
//                          (lb>BCP_DBL_MAX/10 ? " [%e,%i,%i]":" [%.4f,%i,%i]"),
//                            lb, res.termcode(), res.iternum());
//            }
//            p->user->print(true, "\n");
//          }
//      if (babSolver) {
//        p->user->generate_cuts_in_lp(tmp_presolved->lpres(i),
//                                    p->node->vars, p->node->cuts,
//                                    tmp_presolved->get_new_cuts()[i],
//                                    tmp_presolved->get_new_rows()[i]);
//      }
      }
//        for (int i = 0; i < can->child_num; ++i) {
//          const BCP_lp_result& res = tmp_presolved->lpres(i);
//          const double lb = res.objval();
//          p->user->print(true,
//                        (lb>BCP_DBL_MAX/10 ? " [%e,%i,%i]":" [%.4f,%i,%i]"),
//                        lb, res.termcode(), res.iternum());
//        }
//        p->user->print(true, "\n");

      // Compare the current one with the best so far
      l.lock();
      switch (compare_branching_candidates(tmp_presolved,
                                           best_presolved)) {
        case BCP_OldPresolvedIsBetter:
          break;
        case BCP_NewPresolvedIsBetter:
          std::swap(tmp_presolved, best_presolved);
          break;
        case BCP_NewPresolvedIsBetter_BranchOnIt:
          Tools::throwError(
              "BCP_NewPresolvedIsBetter_BranchOnIt not implemented.");
          // Free the remaining candidates if there are any. This also resets
          // candidates.end(), thus
//            purge_ptr_vector(candidates, cani + 1, candidates.end());
//            std::swap(tmp_presolved, best_presolved);
          break;
      }
      // push back the lp used
      lps.push_back(lp);
      l.unlock();
      delete tmp_presolved;
    };
    /**
    * End of the job definition (part run in parallel)
    */
    // The job will be run in parallel.
    pool.run(job);
  }
  // wait for the threads to be finished
  pool.wait();

  // delete lps
  for (auto *lp : lps) {
    lp->unmarkHotStart();
    delete lp;
  }
  lps.clear();

  // indicate to the lp solver that strong branching is done
  p->lp_solver->unmarkHotStart();
  p->lp_solver->setWarmStart(ws);
  delete ws;

  // store objval and termcode in candidate
  BCP_lp_branching_object* can = best_presolved->candidate();
  BCP_vec<double> objval(can->child_num);
  BCP_vec<int> termcode(can->child_num);
  for (i = 0; i < can->child_num; i++) {
    objval[i] = best_presolved->lpres(i).objval();
    termcode[i] = best_presolved->lpres(i).termcode();
  }
  can->set_presolve_result(objval, termcode);
  delete best_presolved;

  // delete all the candidates but the selected one (candidates will just
  // silently go out of scope and we'll be left with a pointer to the final
  // candidate in best_presolved).
  for (auto cani=candidates.begin(); cani != candidates.end(); ++cani) {
    if (*cani != can)
      delete *cani;
  }

  // Mark the cols/rows of the OTHER candidates as removable
  BCP_mark_result_of_strong_branching(p, can, added_colnum, added_rownum);

  // WARNING: cause some segfault as variables can be deleted but still
  // referenced in the variables vector used in BcpModeler
//  // Delete whatever cols/rows we want to delete. This function also updates
//  // var/cut_positions !!!
//  BCP_lp_delete_cols_and_rows(p, can, added_colnum, added_rownum,
//                              false /* not from fathom */,
//                              true /* to force deletion */);

  return can;
}
