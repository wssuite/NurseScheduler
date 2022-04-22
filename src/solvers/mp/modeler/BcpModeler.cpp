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

#include "solvers/mp/modeler/BcpModeler.h"

#include <string>

#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"
#include "solvers/mp/MasterProblem.h"
#include "solvers/mp/modeler/BcpBranchingCandidates.h"
#include "solvers/mp/modeler/BcpHeuristics.h"

#include "OsiClpSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "CoinTime.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_branch.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_lp_node.hpp"  // NOLINT (suppress cpplint error)

#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "cplex.h"  // NOLINT
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_CBC
#include "CbcModeler.h"
#endif

using std::string;
using std::vector;
using std::map;
using std::pair;

/*
 * BCP_lp_user methods
 */

BcpLpModel::BcpLpModel(BcpModeler *pModel) :
    pModel_(pModel),
    candidates_(pModel),
    currentNodelpIteration_(0),
    lpIteration_(0),
    last_node(-1),
    backtracked_(true),
    pHeuristics_(new BcpHeuristics(pModel, pModel->getSolverType())),
    nbNodesSinceLastHeuristic_(0),
    nbCurrentNodeGeneratedColumns_(0),
    nbGeneratedColumns_(0),
    nbCurrentNodeSPSolved_(0),
    // create the timer that records the life time of the solver and start it
    timerTotal_(true) {
  const auto &param = pModel_->getParameters();
  // Initialization of nb_dives_to_wait_before_branching_on_columns_
  if (param.nbDiveIfBranchOnColumns_) {
    for (int i = 1; i < 1000; ++i)
      nb_dives_to_wait_before_branching_on_columns_.push_back(pow(i, 2));
  }
  if (param.performMIPHeuristic_ && param.performHeuristicAfterXNode_ == 0 &&
      Tools::ThreadsPool::getMaxGlobalThreads() == 1) {
    Tools::throwError("The MIP heuristic is run sequentially "
                      "at each iteration of column generation "
                      "as multi-threading is not enable. "
                      "You must either enable multi-threading, "
                      "set performHeuristicAfterXNode_ "
                      "to a non-null value (-1 or >0), "
                      "or set performMIPHeuristic_ to false.");
  }
  if (param.performMIPHeuristic_ && param.performDiveHeuristic_)
    Tools::throwError("You cannot enable both performDiveHeuristic_ "
                      "and performMIPHeuristic_.");
  else if (param.performHeuristicAfterXNode_ >= 0 &&
      !param.performMIPHeuristic_ &&
      !param.performDiveHeuristic_)
    Tools::throwError("You need to enable one of the two heuristics "
                      "(performDiveHeuristic_ and performMIPHeuristic_) "
                      "in order to run a heuristic.");
}

BcpLpModel::~BcpLpModel() {
  // delete BcpHeuristics
  delete pHeuristics_;
  // kill the timer
  timerTotal_.stop();
  // record stats before being deleted
  pModel_->addTimeStats(getTimeStats());
  pModel_->addNbLpIterations(getNbLpIterations());
}

// Initialize the lp parameters and the OsiSolver
OsiSolverInterface *BcpLpModel::initialize_solver_interface() {
  for (const auto &entry : pModel_->getLpParameters())
    set_param(entry.first, entry.second);
  // set the optimal absolute gap for BCP
  set_param(BCP_lp_par::Granularity,
            pModel_->getParameters().optimalAbsoluteGap_ - pModel_->epsilon());
  if (pModel_->getParameters().nCandidates_ == 1)
    set_param(pModel_->strong_branching.first, -1);
  else if (pModel_->getParameters().nCandidates_ > 1)
    set_param(pModel_->strong_branching.first,
              pModel_->strong_branching.second);
  else
    Tools::throwError("The number of candidates needs to be positive "
                      "instead of %d.", pModel_->getParameters().nCandidates_);

  OsiSolverInterface *solver = getNewSolver(pModel_->getSolverType());
  if (solver == nullptr)
    Tools::throwError("%s is not linked to the executable or implemented.",
                      getNameForEnum(SolverTypesByName,
                                     pModel_->getSolverType()).c_str());

  int verbosity = std::max(0, pModel_->getVerbosity() - 1);
  solver->messageHandler()->setLogLevel(verbosity);
  pModel_->setInfinity(solver->getInfinity());

  return solver;
}

// Try to generate a heuristic solution (or return one generated during
// cut/variable generation).
// Return a pointer to the generated solution or return a NULL pointer.
BCP_solution *BcpLpModel::generate_heuristic_solution(
    const BCP_lp_result &lpres,
    const BCP_vec<BCP_var *> &vars,
    const BCP_vec<BCP_cut *> &cuts) {
  BCP_solution_generic *sol = nullptr;

  // if no integer solution is needed, don't run the heuristic
  if (pModel_->getParameters().performHeuristicAfterXNode_ == -1
      || pModel_->getParameters().stopAfterXSolution_ == 0)
    return sol;

  // if heuristic has already been run in these node or
  // it has not been long enough since the last run or
  // the objective of the sub-problem is too negative
  if (nbNodesSinceLastHeuristic_
      < pModel_->getParameters().performHeuristicAfterXNode_)
    return sol;

  // use a criterion on the fraction of the solution that is already integer or
  // on the level in the tree (the latter criterion should depend on the
  // instance though)
  if (pModel_->getParameters().heuristicMinIntegerPercent_ > 0) {
    // count the fraction of current solution that is integer
    MasterProblem *pMaster = pModel_->getMaster();
    double fractionInteger =
        pMaster->computeFractionOfIntegerInCurrentSolution();

    if (pModel_->getParameters().printFractionOfInteger_)
      std::cout << "FRACTION OF INTEGER ROSTERS= " << fractionInteger
                << std::endl;

    if (100.0 * fractionInteger
        < pModel_->getParameters().heuristicMinIntegerPercent_) {
      return sol;
    }
    if (pModel_->getParameters().printFractionOfInteger_)
      std::cout << "RUN THE HEURISTIC" << std::endl;
  }

  if (pModel_->getParameters().performDiveHeuristic_) {
    // diving heuristic. Only perform it if:
    // not a column node (already diving)
    // node has been solved close to optimality (no column generated)
    if (!pModel_->is_columns_node() && nbGeneratedColumns_ == 0)
      sol = pHeuristics_->rounding(getLpProblemPointer()->lp_solver, vars);
  } else if (pModel_->getParameters().performMIPHeuristic_) {
    // solve a mip with all active columns in parallel if multi-threading enable
    sol = pHeuristics_->mip_solve(getLpProblemPointer()->lp_solver, vars);
  }

  nbNodesSinceLastHeuristic_ = 0;

  return sol;
}

// Modify parameters of the LP solver before optimization.
// This method provides an opportunity for the user to change parameters of
// the LP solver before optimization in the LP solver starts.
// The second argument indicates whether the optimization is a "regular"
// optimization or it will take place in strong branching.
// Default: empty method.
void BcpLpModel::modify_lp_parameters(OsiSolverInterface *lp,
                                      const int changeType,
                                      bool in_strong_branching) {
  // DBG - write initial model
#ifdef  DBG
  if (last_node == -1)
    pModel_->writeProblem("master_model");
#endif

  if (current_index() != last_node) {
    last_node = current_index();
    currentNodeStartTime_ = CoinWallclockTime();

    // print a line as it is the first iteration of this node
    if (pModel_->getParameters().printBcpSummary_) {
      std::cout
          << "============================================================"
             "============================================================"
             "==========================="
          << std::endl;
      printSummaryLineHeaders();
    }

    if (pModel_->getParameters().printRelaxationLp_) {
      lp->writeLp("outfiles/test");
    }

#if DBG
    // writeLP("model_"+std::to_string(current_index()));
#endif
    // modify dual tolerance // DBG
    // double dualTol = std::min(0.1,
    // -pModel_->getParameters().sp_max_reduced_cost_bound_+pModel_->epsilon());
    // lp->setDblParam( OsiDualTolerance,dualTol);
    if (pModel_->getSolverType() == Cplex) {
#ifdef USE_CPLEX
      // do not let cplex use more than one thread, otherwise it will, and it
      // it will also keep resetting the parameter to the default value for some
      // unexplained reason
      if (auto * pSolverTmp = dynamic_cast<OsiCpxSolverInterface *>(lp)) {
        CPXsetintparam(pSolverTmp->getEnvironmentPtr(), CPX_PARAM_THREADS, 1);
      }
      // modifiy the value of cplex random seed otherwise the results will never
      // be reproduced
      if (auto * pSolverTmp = dynamic_cast<OsiCpxSolverInterface *>(lp)) {
        CPXsetintparam(pSolverTmp->getEnvironmentPtr(), CPXPARAM_RandomSeed, 1);
        std::cout << "Cplex random seed set to " << 1 << std::endl;
      }
      // use barrier optimization instead of dual simplex for reoptimization
      // after branching,because the solution parameter is just too degenerate
      if (changeType % 2 == 1) {
        if (auto * pSolverTmp = dynamic_cast<OsiCpxSolverInterface *>(lp)) {
          CPXLPptr cpxlp =
              pSolverTmp->getLpPtr(OsiCpxSolverInterface::FREECACHED_RESULTS);
          CPXbaropt(pSolverTmp->getEnvironmentPtr(), cpxlp);
        }
      }
#endif
    }

    // STAB
    // set the cost and upper bounds of every core variables back to their
    // values at the time of solution
    if (pModel_->getParameters().isStabilization_) {
      for (MyVar *pVar : pModel_->getCoreVars()) {
        int varind = pVar->getIndex();
        double cost = lp->getObjCoefficients()[varind];
        double ub = lp->getColUpper()[varind];
        pVar->setCost(cost);
        pVar->setUB(ub);
      }
    }
  } else {
    // DBG
    // double dualTol = std::min(0.1,
    // -pModel_->getParameters().sp_max_reduced_cost_bound_+pModel_->epsilon());
    // lp->setDblParam( OsiDualTolerance,dualTol);
  }
//  if (in_strong_branching)
//    lp->set

  // use default one
  BCP_lp_user::modify_lp_parameters(lp, changeType, in_strong_branching);
}

// print in cout a line summary headers and of the current solver state
void BcpLpModel::printSummaryLineHeaders() const {
  if (pModel_->getVerbosity() == 0) return;

  FILE *pFile =
      pModel_->logfile().empty() ?
      stdout :
      fopen(pModel_->logfile().c_str(), "a");

  fprintf(pFile,
          "BCP: %13s %5s | %10s  %10s  %10s | "
          "%8s %14s %14s %13s %10s | %14s %5s %5s \n",
          "Node",
          "Lvl",
          "BestUB",
          "RootLB",
          "BestLB",
          "#It",
          "Obj",
          "LagLB",
          "#Frac",
          "#Active",
          "ObjSP",
          "#SP",
          "#Col");
//  fprintf(pFile,
//          "BCP: %5d / %5d %5d | %10.0f  %10.2f  %10.2f | "
//          "%8s %14s %13s %10s | %14s %5s %5s \n",
//          current_index(),
//          pModel_->getTreeSize(),
//          current_level(),
//          pModel_->getObjective(),
//          pModel_->getRootLB(),
//          pModel_->getBestLB(),
//          "-",
//          "-",
//          "-",
//          "-",
//          "-",
//          "-",
//          "-");

  if (!pModel_->logfile().empty()) fclose(pFile);
}

// print in cout a line summary of the current solver state and of the tree
void BcpLpModel::printSummaryLine(bool printNode,
                                  const BCP_vec<BCP_var *> &vars) const {
  if (pModel_->getVerbosity() > 0 ||
      pModel_->getParameters().printBcpSummary_) {
    printSummaryLine(vars);
    if (printNode) printNodeSummaryLine();
  }
}

// print in cout a line summary of the current solver state
void BcpLpModel::printSummaryLine(const BCP_vec<BCP_var *> &vars) const {
  FILE *pFile =
      pModel_->logfile().empty() ?
      stdout :
      fopen(pModel_->logfile().c_str(), "a");

  double ub = pModel_->getObjective(), root = pModel_->getRootLB(),
  lb = pModel_->getBestLB(), lagLb = pModel_->getNodeLastLagLB();
  if (ub >= LARGE_SCORE - 1) ub = INFINITY;
  if (root >= LARGE_SCORE - 1) root = INFINITY;
  if (lb >= LARGE_SCORE - 1) lb = INFINITY;
  if (lagLb <= -LARGE_SCORE + 1) lagLb = -INFINITY;

  // if has some variables
  if (!vars.empty()) {
    std::pair<int, int> p = pModel_->getFractionalAndPositiveColumns();

    // if at least a subproblem has been solved
    if (pModel_->getLastNbSubProblemsSolved() > 0) {
      fprintf(pFile,
              "BCP: %5d / %5d %5d | %10.0f  %10.2f  %10.2f | "
              "%8d %14.2f %14.2f %5d / %5d %10ld | %14.2f %5d %5d \n",
              current_index(),
              pModel_->getTreeSize(),
              current_level(),
              ub,
              root,
              lb,
              lpIteration_,
              pModel_->getLastObj(),
              lagLb,
              p.first,
              p.second,
              vars.size() - pModel_->getCoreVars().size(),
              pModel_->getLastMinReducedCost(),
              pModel_->getLastNbSubProblemsSolved(),
              nbGeneratedColumns_);
    } else {
      fprintf(pFile,
              "BCP: %5d / %5d %5d | %10.0f  %10.2f  %10.2f | "
              "%8s %14.2f %14.2f %5d / %5d %10ld | %14s %5s %5s \n",
              current_index(),
              pModel_->getTreeSize(),
              current_level(),
              ub,
              root,
              lb,
              "-",
              pModel_->getLastObj(),
              lagLb,
              p.first,
              p.second,
              vars.size() - pModel_->getCoreVars().size(),
              "-",
              "-",
              "-");
    }
  } else {
    fprintf(pFile,
            "BCP: %5d / %5d %5d | %10.0f  %10.2f  %10.2f | "
            "%8d %14.2f %14.2f %5s / %5s %10s | %14s %5s %5s \n",
            current_index(),
            pModel_->getTreeSize(),
            current_level(),
            ub,
            root,
            lb,
            lpIteration_,
            pModel_->getLastObj(),
            lagLb,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-");
  }

  if (!pModel_->logfile().empty()) fclose(pFile);
}

void BcpLpModel::printNodeSummaryLine(int nbChildren) const {
  FILE *pFile =
      pModel_->logfile().empty() ?
      stdout :
      fopen(pModel_->logfile().c_str(), "a");

  fprintf(pFile,
          "BCP %6s: Node %5d processed; %5d nodes lefts; %5d nodes processed "
          "-- Node's totals: #It=%3d; #SP=%5d; #Col=%8d; Time:%8.2f s.\n",
          nbChildren ? "BRANCH" : "FATHOM",
          current_index(),
          pModel_->getTreeSize() - 1
              + nbChildren,  // node left -1 (current node) + nb children
          pModel_->getNbNodesProcessed() + 1,  // +1 for current node
          getNbCurrentNodeLpIterations(),
          nbCurrentNodeSPSolved_,
          nbCurrentNodeGeneratedColumns_,
          CoinWallclockTime() - currentNodeStartTime_);

  if (!pModel_->logfile().empty()) fclose(pFile);
}

// stop this node or BCP
bool BcpLpModel::doStop(const BCP_vec<BCP_var *> &vars) {
  if (pModel_->doStop(vars))
    return true;

  // fathom if the true lower bound greater than current upper bound
  return over_ub(getLpProblemPointer()->node->true_lower_bound);
}

//// store the LP results
void BcpLpModel::process_lp_result(const BCP_lp_result &lpres,
                                   const BCP_vec<BCP_var *> &vars,
                                   const BCP_vec<BCP_cut *> &cuts,
                                   const double old_lower_bound,
                                   double &true_lower_bound,
                                   BCP_solution *&sol,
                                   BCP_vec<BCP_cut *> &new_cuts,
                                   BCP_vec<BCP_row *> &new_rows,
                                   BCP_vec<BCP_var *> &new_vars,
                                   BCP_vec<BCP_col *> &new_cols) {
  // update LP solution
  pModel_->setLPSol(lpres, vars, cuts, lpIteration_);

  // print once at the beginning of each node
  if (currentNodelpIteration_ == 0)
    printSummaryLine(false, vars);

  // update the total number of LP solutions (from the beginning)
  ++currentNodelpIteration_;
  ++lpIteration_;

  // run the default method to set user_has_lp_result_processing to false
  // otherwise, BCP won't run the methods generate_vars_in_lp for example,
  // as it will suppose it has been run here.
  BCP_lp_user::process_lp_result(
      lpres, vars, cuts, old_lower_bound, true_lower_bound,
      sol, new_cuts, new_rows, new_vars, new_cols);
}

// This method provides an opportunity for the user to tighten the bounds of
// variables. The method is invoked after reduced cost fixing. The results
// are returned  in the last two parameters.
// Parameters:
// lpres    the result of the most recent LP optimization,
// vars  the variables in the current formulation,
// status   the stati of the variables as known to the system,
// var_bound_changes_since_logical_fixing    the number of variables whose
// bounds have changed (by reduced cost fixing) since the most recent
// invocation of this method that has actually forced changes returned
// something in the last two arguments,
// changed_pos   the positions of the variables whose bounds should be changed
// new_bd   the new bounds (lb/ub pairs) of these variables.
///////////////////////////////////////////////////////////////////////////////////////
// WARNING: JUST FIX BOUNDS THAT ARE COMPULSORY IN AN INTEGER POINT OF VIEW
// AS THE FINAL LP BOUND OF THE NODE IS GIVEN TO ITS CHILDREN,
// FIXING WRONG BOUNDS CAN LEAD TO A FALSE LOWER BOUND,
// THUS FATHOMING WRONG NODES.
///////////////////////////////////////////////////////////////////////////////////////
void BcpLpModel::logical_fixing(
    const BCP_lp_result &lpres,
    const BCP_vec<BCP_var *> &vars,
    const BCP_vec<BCP_cut *> &cuts,
    const BCP_vec<BCP_obj_status> &var_status,
    const BCP_vec<BCP_obj_status> &cut_status,
    const int var_bound_changes_since_logical_fixing,
    BCP_vec<int> &changed_pos,  // NOLINT
    BCP_vec<double> &new_bd) {  // NOLINT
  // DBG
  // pModel_->checkActiveColumns(vars);
}

////////////////////////////////////////////////////////////////////////////////
// WARNING: THIS METHOD NEEDS TO BE REIMPLEMENTED AS EMPTY, BECAUSE IT CAUSES
// CYCLING IN THE COLUMN GENERATION WHEN THE LOWER BOUND IS TOO CLOSE FROM THE
// UPPER BOUND
////////////////////////////////////////////////////////////////////////////////

// Restoring feasibility.
// This method is invoked before fathoming a search tree node that has been
// found infeasible and the variable pricing did not generate any new variables.
void BcpLpModel::restore_feasibility(const BCP_lp_result &lpres,
                                     const std::vector<double *> dual_rays,
                                     const BCP_vec<BCP_var *> &vars,
                                     const BCP_vec<BCP_cut *> &cuts,
                                     BCP_vec<BCP_var *> &vars_to_add,
                                     BCP_vec<BCP_col *> &cols_to_add) {
  writeMPS("infeasible");
  std::cout << "Infeasible LP: LP model written in infeasible.lp" << std::endl;
  // dive is finished
  // DBG dive_ = false;
  if (pModel_->getParameters().printBranchStats_) {
    std::cout << "DIVING STOPPED WITH INFEASIBLE SOLUTION" << std::endl;
    pModel_->printStats();
  }
}

// Convert a set of variables into corresponding columns for the current
// LP relaxation.
void BcpLpModel::vars_to_cols(
    const BCP_vec<BCP_cut *> &cuts,  // on what to extend
    BCP_vec<BCP_var *> &vars,  // what to extend   NOLINT
    BCP_vec<BCP_col *> &cols,  // the expanded cols  NOLINT
    // few things that the user can use for lifting vars if allowed
    const BCP_lp_result &lpres,
    BCP_object_origin origin,
    bool allow_multiple) {
  const int varnum = vars.size();
  if (varnum == 0)
    return;
  cols.reserve(varnum);

  for (int i = 0; i < varnum; ++i) {
    auto *var = dynamic_cast<CoinVar *>(vars[i]);
    if (!var)
      Tools::throwError("Bad variable casting.");

    // Copy the vectors var->getIndexRows() and var->getCoeffRows() in arrays
    const int size = var->getNbRows();

    // create a new array which will be deleted by ~BCP_col()
    int *indexRows = new int[size];
    const vector<int> &index = var->getIndexRows();
    copy(index.begin(), index.end(), indexRows);

    // create a new array which will be deleted by ~BCP_col()
    double *coeffRows = new double[size];
    const vector<double> &coeffs = var->getCoeffRows();
    copy(coeffs.begin(), coeffs.end(), coeffRows);

    cols.unchecked_push_back(
        new BCP_col(size,
                    indexRows,
                    coeffRows,
                    var->getCost(),
                    var->getLB(),
                    var->getUB()));
  }
}

// Generate variables within the LP process.
void BcpLpModel::generate_vars_in_lp(const BCP_lp_result &lpres,
                                     const BCP_vec<BCP_var *> &vars,
                                     const BCP_vec<BCP_cut *> &cuts,
                                     const bool before_fathom,
                                     BCP_vec<BCP_var *> &new_vars,  // NOLINT
                                     BCP_vec<BCP_col *> &new_cols) {  // NOLINT
  // if must stop, return immediately
  if (doStop(vars))
    return;

//  auto p = getLpProblemPointer();
//  auto f = std::setprecision(10);
//  auto print = [f](OsiSolverInterface *lp) {
//    std::stringstream ss;
//    ss << "Objective: " << f << lp->getObjValue()
//       << " " << lp->getIterationCount()
//       << " " << lp->isProvenOptimal()
//       << " " << lp->isProvenPrimalInfeasible()
//       << " " << lp->isProvenDualInfeasible()
//       << " " << lp->isAbandoned()
//       << " " << lp->isIterationLimitReached()
//       << " " << lp->isPrimalObjectiveLimitReached()
//       << " " << lp->isDualObjectiveLimitReached();
//    return ss;
//  };
//  auto ss = print(p->lp_solver);
//  p->lp_solver->resolve();
//  double diff = lpres.objval() - p->lp_solver->getObjValue();
//  if (fabs(diff) > 1e-2) {
//    int l, hl;
//    p->lp_solver->getIntParam(OsiMaxNumIteration, l);
//    p->lp_solver->getIntParam(OsiMaxNumIterationHotStart, hl);
//    double pb, db;
//    p->lp_solver->getDblParam(OsiPrimalObjectiveLimit, pb);
//    p->lp_solver->getDblParam(OsiDualObjectiveLimit, db);
//    std::cout << "Diff: " << f << diff << " " << f << lpres.primalTolerance()
//              << " " << f << pb << " " << f << db
//              << " it limit " << l << " hot it limit " << hl <<std::endl;
//    std::cout << ss.str() << std::endl;
//    std::cout << "Resolve LP and " << print(p->lp_solver).str() << std::endl;
//    auto lp = p->lp_solver->clone(true);
//    std::cout << "New LP and " << print(lp).str() << std::endl;
//    if (p->node->warmstart) {
//      auto ws = p->node->warmstart->convert_to_CoinWarmStart();
//      p->lp_solver->setWarmStart(ws);
//      p->lp_solver->resolve();
//      std::cout << "Warm LP and " << print(p->lp_solver).str() << std::endl;
//    }
//
//    std::cout << writeMPS("model_" + std::to_string(current_index()));
//    Tools::throwException("When resolving LP, got a different solution.");
//  }

  // Stop the algorithm if the objective value of the relaxation is
  // close enough to the current LB
  // WARNING: true if stabilization is inactive
  bool stabInactive = !pModel_->getParameters().isStabilization_ ||
      pModel_->stab().stabCheckStoppingCriterion();
  if (current_index() > 0 && stabInactive &&
      lpres.objval() < pModel_->getCurrentLB() + pModel_->epsilon())
    return;

  // call the rotation pricer to find columns that should be added to the LP
  bool after_fathom = (currentNodelpIteration_ == 1);
  // max reduced cost of a rotation that would be added to MP (a tolerance is
  // substracted in the SP)
  double maxReducedCost =
      pModel_->getParameters().spParam_.sp_max_reduced_cost_bound_;
  vector<MyVar *> generatedColumns = pModel_->pricing(maxReducedCost,
                                                      before_fathom,
                                                      after_fathom,
                                                      backtracked_);
  nbGeneratedColumns_ = generatedColumns.size();
  nbCurrentNodeGeneratedColumns_ += nbGeneratedColumns_;
  nbCurrentNodeSPSolved_ += pModel_->getLastNbSubProblemsSolved();

  // check if new columns add been added since the last time
  // if there are some, add all of them in new_vars
  new_vars.reserve(nbGeneratedColumns_);
  for (MyVar *var : generatedColumns) {
    auto *col = dynamic_cast<BcpColumn *>(var);
    // the BcpColumn which will be deleted by BCP (needs to be owned by BCP)
    new_vars.unchecked_push_back(col);
    // initialize the counter of active iteration for this new variable
    col->addActiveIteration(lpIteration_);
  }

  // debug
  // pModel_->checkActiveColumns(new_vars);

  // must be before fathoming. we need vars with red cost below the
  // negative of (lpobj-ub)/ks_num otherwise we can really fathom.
  //    const double rc_bound =
  //   (lpres.dualTolerance() + (lpres.objval() - upper_bound()))/kss.ks_num;
  //    generate_vars(lpres, vars, rc_bound, new_vars);
}

/*
 * BCP_DoNotBranch_Fathomed: The node should be fathomed without even trying to branch.
 * BCP_DoNotBranch: BCP should continue to work on this node.
 * BCP_DoBranch: branch on one of the candidates cands
 *
 */

BCP_branching_decision BcpLpModel::select_branching_candidates(
    // the result of the most recent LP optimization.
    const BCP_lp_result &lpres,
    // the variables in the current formulation.
    const BCP_vec<BCP_var *> &vars,
    // the cuts in the current formulation.
    const BCP_vec<BCP_cut *> &cuts,
    // the local pool that holds variables with negative reduced cost.
    const BCP_lp_var_pool &local_var_pool,
    // In case of continuing with the node the best so many variables
    // will be added to the formulation
    // (those with the most negative reduced cost).
    // the local pool that holds violated cuts.
    const BCP_lp_cut_pool &local_cut_pool,
    // In case of continuing with the node the best so many cuts will be
    // added to the formulation (the most violated ones).
    // the generated branching candidates.
    BCP_vec<BCP_lp_branching_object *> &cands,  // NOLINT
    // indicate whether to force branching regardless of the size of
    // the local cut/var pools
    bool force_branch) {
  // Print a line summary of the solver state
  pModel_->setCurrentTreeLevel(current_level());

  // select branching decision
  BCP_branching_decision decision = selectBranchingDecision(lpres,
                                                            vars,
                                                            cuts,
                                                            local_var_pool,
                                                            local_cut_pool,
                                                            cands);

  printSummaryLine(false, vars);

  // reset iteration sp data
  nbGeneratedColumns_ = 0;
  // WARNING: important to be reset, otherwise selectBranchingDecision
  // could throw an error on the following iteration
  pModel_->setLastMinReducedCost(0);
  pModel_->setLastNbSubProblemsSolved(0);

  // if fathoming or branching -> perform other operations
  int nbChildren = 0;
  switch (decision) {
    case BCP_DoNotBranch:return BCP_DoNotBranch;
    case BCP_DoNotBranch_Fathomed:backtracked_ = true;
      // if last node explored
      if (pModel_->getTreeSize() <= 1)
        pModel_->stop();
      break;
    case BCP_DoBranch:nbChildren = cands.front()->child_num;
      backtracked_ = false;
      if (pModel_->getVerbosity() >= 3) {
        std::cout << pModel_->getMaster()->costsConstrainstsToString()
                  << std::endl;
        std::cout << pModel_->getMaster()->currentSolToString()
                  << std::endl;
      }
      break;
    default:Tools::throwError("Decision not recognized.");
  }

  // node summary
  if (pModel_->getParameters().printBcpSummary_)
    printNodeSummaryLine(nbChildren);

  // print the current solution
  if (pModel_->getParameters().printRelaxationSol_) {
    std::cout << pModel_->getParameters().saveFunction_->currentSolToString()
              << std::endl;
  }

  // reset cg counters
  currentNodelpIteration_ = 0;
  nbCurrentNodeGeneratedColumns_ = 0;
  nbCurrentNodeSPSolved_ = 0;
  pModel_->setNbDegenerateIt(0);
  // update heuristic flags
  nbNodesSinceLastHeuristic_++;
  // increase the number of nodes
  pModel_->incrementNbNodes();
  // deactivate stabilization, feasibility, and dual UB
  feasible_ = false;
  approximatedDualUB_ = -LARGE_SCORE;
  if (pModel_->getParameters().isStabilization_)
    pModel_->stab().stabDeactivateBoundAndCost(
        getLpProblemPointer()->lp_solver);

  return decision;
}

BCP_branching_decision BcpLpModel::selectBranchingDecision(
    // the result of the most recent LP optimization.
    const BCP_lp_result &lpres,
    // the variables in the current formulation.
    const BCP_vec<BCP_var *> &vars,
    // the cuts in the current formulation.
    const BCP_vec<BCP_cut *> &cuts,
    // the local pool that holds variables with negative reduced cost.
    const BCP_lp_var_pool &local_var_pool,
    // In case of continuing with the node the best so many variables will be
    // added to the formulation (those with the most negative reduced cost).
    // the local pool that holds violated cuts.
    const BCP_lp_cut_pool &local_cut_pool,
    // In case of continuing with the node the best so many cuts will be
    // added to the formulation (the most violated ones).
    // the generated branching candidates.
    BCP_vec<BCP_lp_branching_object *> &cands) {  // NOLINT
  // if some variables have been generated, do not branch
  bool column_generated = !local_var_pool.empty();
  double lb = lpres.objval();

  // throw an error if no columns have been generated and the min reduced cost
  // is negative
  if (!column_generated
      && pModel_->getLastMinReducedCost() < -pModel_->epsilon())
    Tools::throwError("Column generation has finished with a negative reduced "
                      "cost (%.2f). There is a problem with the pricing.",
                      pModel_->getLastMinReducedCost());


  // STAB
  // activate and update stabilization variables if enable.
  // Also, check stopping dual stabilization criterion
  bool isStabActive = false;
  if (pModel_->getParameters().isStabilization_) {
    double dualUB =
        pModel_->getMaster()->computeApproximateDualUB(lpres.objval());

    isStabActive = !pModel_->stab().stabCheckStoppingCriterion();
    // if was infeasible at the previous iteration and not anymore:
    // activate stabilization
    if (becomeFeasible())
      pModel_->stab().stabInitializeBoundAndCost(
          getLpProblemPointer()->lp_solver);
      // Update the stabilization variables if:
      // feasible and the approximated dual UB has improved
    else if (feasible_)
      pModel_->stab().stabUpdate(
          getLpProblemPointer()->lp_solver,
          dualUB > approximatedDualUB_ || !column_generated);
    // update approximated dual UB.
    // Used as a criteria to decide when we obtain a better dual solution
    if (dualUB > approximatedDualUB_) approximatedDualUB_ = dualUB;
    // Do not branch if some stabilization variables are positive
    if (!column_generated && isStabActive)
      return BCP_DoNotBranch;
  }

  // STAB: compute the Lagrangian bound
  // It can also be used in general to fathom nodes when the the Lagrangian
  // bound is larger than the best UB
  // To be able to compute this bound, three conditions needs to be met:
  // 1. A lagrangian bound is available
  // 2. The subproblem were solved at optimality
  // 3. There are no stabilization variables in the solution
  if (pModel_->getMaster()->lagrangianBoundAvailable() &&
      pModel_->isLastPricingOptimal() &&
      !isStabActive) {
    double lagLb = pModel_->getMaster()->computeLagrangianBound(lpres.objval());
    pModel_->updateNodeLagLB(lagLb);
    double bestLagLb = pModel_->getNodeBestLagLB();
    if (column_generated) {
      // fathom only if column generation would continue
      // (otherwise would be fathom later in this function)
      if (pModel_->getParameters().isLagrangianFathom_ &&
          (pModel_->getParameters().isLagrangianFathomRootNode_
              || current_index() > 0)
          && over_ub(bestLagLb)) {
        // update lb with worst bound as fathoming
        pModel_->updateNodeLB(pModel_->getBestUB());
        if (pModel_->getParameters().printBcpSummary_)
          std::cout << "Forcibly fathom, because Lagrangian bound is exceeded."
                    << std::endl;
        return BCP_DoNotBranch_Fathomed;
      }
      // if lagrangian bound is really small -> stop column generation also
      // use lagrangian LB instead of current objective
      double diff = lb - bestLagLb;
      if (diff < pModel_->epsilon()) {
        column_generated = false;
        std::cout << "Forcibly stopped CG, because Lagrangian bound is very "
                     "closed to the current objective." << std::endl;
        lb = std::max(bestLagLb, pModel_->getCurrentLB());
      }
    }
  }

  // STAB:
  // Detect when the column generation is stalling
  // If stabilization is used this will determine when the stabilization costs
  // are updated
  // Otherwise, an option can be set on to stop column generation and branch
  // after a given number of degenerate iterations
  //
  bool isStalling = false;
  if (column_generated && getObjVariation() < pModel_->epsilon()) {
    isStalling = true;
    pModel_->incrementNbDegenerateIt();
    // stop column generation if not root node and too many iteration,
    // and current solution is fractional
    if (current_index() > 0
        && pModel_->getNbDegenerateIt()
            >= pModel_->getParameters().stopAfterXDegenerateIt_
        && pModel_->hasFractionalColumns()) {
      if (pModel_->getVerbosity() > 0)
        std::cout << "Branch with column generation stalling "
                     "(stop column generation)" << std::endl;
      column_generated = false;
      lb = std::max(pModel_->getNodeBestLagLB(), pModel_->getCurrentLB());
    }
  } else {
    // reset counter to 0 as soon as a significant step has been made
    pModel_->setNbDegenerateIt(0);
  }

  // check if continue column generation
  if (column_generated)
    return BCP_DoNotBranch;

  // update LB (as either branching or fathoming)
  // update node if stabilization not active.
  // If return false, the lb has decreased !
  // Bug needs to be found
  if (!isStabActive && !pModel_->updateNodeLB(lb)) {
    auto lp = getLpProblemPointer()->lp_solver;
    lp->resolve();
    std::cout << "Resolve LP and got objective: "
              << lp->getObjValue() << std::endl;
    std::cout << writeMPS("model_" + std::to_string(current_index()));
    std::cerr << "The best LB has decreased." << std::endl;
//    Tools::throwException("The best LB has decreased.");
  }
  // update true_lower_bound, as we reach the end of the column generation
  getLpProblemPointer()->node->true_lower_bound = lb;

#ifdef DBG
//   std::cout << pModel_->pMaster_->costsConstrainstsToString() << std::endl;
#endif

  // if root and a variable with the obj LARGE_SCORE is positive -> INFEASIBLE
  // otherwise, record the root solution for future use
  if (current_index() == 0) {
    if (pModel_->getTimeFirstRoot() < 0)
      pModel_->setTimeFirstRoot(timerTotal_.dSinceStart());
    // Check if infeasible
    if (!pModel_->isFeasible()) {
      pModel_->getMaster()->status(INFEASIBLE);
      std::cerr << "Feasibility core variable is still present "
                   "in the solution" << std::endl;
      pModel_->printInfeasibleVars();
      return BCP_DoNotBranch_Fathomed;
    }
  }

  // stop this process for BCP or the node
  if (doStop(vars))
    return BCP_DoNotBranch_Fathomed;

  // fathom if greater than current upper bound
  if (over_ub(lb))
    return BCP_DoNotBranch_Fathomed;

  // build a candidate
  bool generate = true;
  std::list<MyPBranchingCandidate> candidates;
  // if a column node, continue the dive
  if (pModel_->is_columns_node()) {
    generate = pModel_->branch_on_column(&candidates);
  } else {
    // branching decisions: rest on a day, branch on shifts ...
    generate = pModel_->branching_decisions(&candidates);

    // do we branch on columns ?
    if (pModel_->getNbDives()
        < pModel_->getParameters().nbDiveIfBranchOnColumns_) {
      pModel_->branch_on_column(&candidates);
    } else if (!nb_dives_to_wait_before_branching_on_columns_.empty() &&
        pModel_->getNbDives()
            >= nb_dives_to_wait_before_branching_on_columns_.front()) {
      // after a given number of nodes since last dive,
      // prepare to go for a new dive
      nb_dives_to_wait_before_branching_on_columns_.pop_front();
      pModel_->resetNbNodesSinceLastDive();
    }
  }

  if (!generate) {
    if (testPartialFeasibility())
      return BCP_DoNotBranch_Fathomed;
    // throw an error here. Should never happened
    find_infeasibility(lpres, vars);
    Tools::throwError("Solution should be fractional as no branching "
                      "candidate has been found.");
  }

  // create BCP candidates
  candidates_.reset(lpres.objval());
  BCP_lp_branching_object *can = candidates_.selectCandidates(
      candidates, getLpProblemPointer(), vars, cuts);
  cands.push_back(can);

  // check if BCP should fathom in case BCP has stopped
  // or a solution has been found during strong branching
  if (pModel_->isStopped() || over_ub(lb)) {
    for (auto *can : cands) delete can;
    cands.clear();
    return BCP_DoNotBranch_Fathomed;
  }

  // else return the potential candidate
  return BCP_DoBranch;
}

// Decide what to do with the children of the selected branching object.
// Fill out the _child_action field in best.
// This will specify for every child what to do with it.
// Possible values for each individual child are:
// BCP_PruneChild, BCP_ReturnChild and BCP_KeepChild.
// There can be at most child with this last action specified.
// It means that in case of diving this child will be processed
// by this LP process as the next search tree node.
// Default: Every action is BCP_ReturnChild.
// However, if BCP dives then one child will be mark with BCP_KeepChild.
// The decision which child to keep is based on the ChildPreference
// parameter in BCP_lp_par.
// Also, if a child has a presolved lower bound that is higher than
// the current upper bound then that child is mark as BCP_FathomChild.
void BcpLpModel::set_actions_for_children(BCP_presolved_lp_brobj *best) {
  if (best->candidate()->child_num == 0)
    Tools::throwError("No candidate has been generated.");

  // copy the nodes to the tree and clear candidates
  candidates_.updateTree(best);

  // set the lower bound of the children (very important when
  // performing strong branching)
  best->initialize_lower_bound(pModel_->getCurrentLB());

  // if only one child -> keep child as diving
  // if not looking for optimality -> keep diving
  // if the gap is less than the min relative gap and
  // the node's LP gap is less than integrality gap/2 or
  // maxRelativeLPGapToKeepChild_ -> keep child
  // The goal is to dive while the gap is more than
  // 2*maxRelativeLPGapToKeepChild_ and then explore the tree breadth
  double minGap = std::max(
      pModel_->getIntegralityGap()/2,
      pModel_->getParameters().maxRelativeLPGapToKeepChild_);
  if (best->action().size() == 1 ||
      !pModel_->getParameters().solveToOptimality_ ||
      pModel_->getCurrentNode()->getLPGap() < minGap) {
    best->action()[0] = BCP_KeepChild;
    // tell the tree that the first child if kept
    pModel_->keepFirstChild(best->action().size());
  }
}

// rerun the code use to test the integer feasibility of a solution and
// find why a solution is not feasible
void BcpLpModel::find_infeasibility(
    // the result of the most recent LP optimization.
    const BCP_lp_result &lpres,
    const BCP_vec<BCP_var *> &vars) {
  // check if the solution is feasible for BCP in case of exception
  auto p = getLpProblemPointer();
  // Do anything only if the termination code is sensible
  const int tc = lpres.termcode();
  if (!(tc & BCP_ProvenOptimal)) {
    std::cout << "Termination code: " << tc << " & " << BCP_ProvenOptimal
              << std::endl;
    return;
  }

  const double etol = p->param(BCP_lp_par::IntegerTolerance);
  const double *x = lpres.x();
  const int varnum = vars.size();
  const double etol1 = 1 - etol;
  int i;
  for (i = 0; i < varnum; ++i)
    switch (vars[i]->var_type()) {
      case BCP_BinaryVar: {
        const double val = x[i];
        if (val > etol && val < etol1) {
          std::cout << "Binary var " << i << ", bcp value=" << val << std::endl;
          auto var = dynamic_cast<MyVar *>(vars[i]);
          if (i < pModel_->getCoreVars().size())
            std::cout << "Core variable: ";
          else
            std::cout << "Column variable: ";
          std::cout << i << ", " << var->name_ << "model value="
                    << pModel_->getVarValue(var) << ", pattern :";
          for (int j : var->getPattern()) std::cout << " " << j;
          std::cout << std::endl;
          std::cout
              << "-----------------------------------------------------------"
              << std::endl;
        }
      }
        break;
      case BCP_IntegerVar: {
        const double val = x[i] - floor(x[i]);
        if (val > etol && val < etol1) {
          std::cout << "Integer var " << i << ": " << val << std::endl;
          auto *var = dynamic_cast<MyVar *>(vars[i]);
          if (i < pModel_->getCoreVars().size())
            std::cout << "Core variable: ";
          else
            std::cout << "Column variable: ";
          std::cout << i << ", " << var->name_ << " = "
                    << pModel_->getVarValue(var) << ", pattern :";
          for (int j : var->getPattern()) std::cout << " " << j;
          std::cout << std::endl;
          std::cout
              << "-----------------------------------------------------------"
              << std::endl;
        }
      }
        break;
      case BCP_ContinuousVar:break;
    }

  // print roster
  auto roster = pModel_->getMaster()->fractionalRoster();
  for (int n = 0; n < roster.size(); n++) {
    std::cout << "N" << n << ":";
    for (int k = 0; k < roster[n].size(); k++) {
      std::cout << " " << k << "(";
      double vLeft = 1;
      for (int s = 0; s < roster[n][k].size(); s++) {
        double v = roster[n][k][s];
        if (v < pModel_->epsilon()) continue;
        vLeft -= v;
        std::cout << s << ":" << v << ",";
      }
      if (vLeft > pModel_->epsilon())
        std::cout << "0:" << vLeft;
      std::cout << ")";
    }
    std::cout << std::endl;
  }
}

// check if the current solution is feasible according to the days that
// have been relaxed
bool BcpLpModel::testPartialFeasibility() {
  MasterProblem *pMaster = pModel_->getMaster();
  // check if the solution could be feasible without being integer
  if (!pMaster->isPartialRelaxed()) return false;

  // if partially relaxed, check that the unrelaxed days are integer
  vector3D<double> fracRosters = pMaster->fractionalRoster();
  for (const auto &vec2D : fracRosters)
    for (int k = 0; k < pMaster->nDays(); ++k) {
      if (pMaster->isRelaxDay(k))
        continue;
      // if should be integer, check integrality
      for (double v : vec2D[k])
        if (!pModel_->isInteger(v))
          return false;
    }

  // if partially integer, store solution and return true
  pModel_->addCurrentBcpSol(false);
  return true;
}

// Here, we generate a cut to branch on a set of variables
void BcpLpModel::cuts_to_rows(
    // the variables currently in the relaxation (IN)
    const BCP_vec<BCP_var *> &vars,
    // the cuts to be converted (IN/OUT)
    BCP_vec<BCP_cut *> &cuts,  // NOLINT
    // the rows into which the cuts are converted (OUT)
    BCP_vec<BCP_row *> &rows,  // NOLINT
    // solution to the current LP relaxation (IN)
    const BCP_lp_result &lpres,
    // where the cuts come from (IN)
    BCP_object_origin origin,
    // whether multiple expansion, i.e., lifting, is allowed (IN)
    bool allow_multiple) {
  rows.reserve(cuts.size());
  for (BCP_cut *cut : cuts) {
    auto *branchingCut = dynamic_cast<BcpBranchCons *>(cut);
    if (!cut)
      Tools::throwError("Should be a branching cut.");

    // create new arrays which will be deleted by ~BCP_row()
    const int size = branchingCut->getIndexCols().size();
    int *elementIndices = new int[size];
    copy(branchingCut->getIndexCols().begin(),
         branchingCut->getIndexCols().end(),
         elementIndices);
    double *elementValues = new double[size];
    copy(branchingCut->getCoeffCols().begin(),
         branchingCut->getCoeffCols().end(),
         elementValues);
    rows.unchecked_push_back(new BCP_row(size,
                                         elementIndices,
                                         elementValues,
                                         branchingCut->getLhs(),
                                         branchingCut->getRhs()));
  }
}

void BcpLpModel::select_vars_to_delete(const BCP_lp_result &lpres,
                                       const BCP_vec<BCP_var *> &vars,
                                       const BCP_vec<BCP_cut *> &cuts,
                                       const bool before_fathom,
                                       BCP_vec<int> &deletable) {
  if (before_fathom
      && getLpProblemPointer()->param(BCP_lp_par::NoCompressionAtFathom))
    return;

  // code from BCP_lp_user::select_vars_to_delete
  const int varnum = vars.size();
  deletable.reserve(varnum);
  for (int i = getLpProblemPointer()->core->varnum(); i < varnum; ++i) {
    BCP_var *var = vars[i];
    // Do not delete these variables as they are being branched on
    if (var->is_non_removable())
      continue;

    if (var->is_to_be_removed() || (var->lb() == 0 && var->ub() == 0)) {
      deletable.unchecked_push_back(i);
      continue;
    }

    auto *col = dynamic_cast<BcpColumn *>(vars[i]);
    if (col->getNbConsInactiveIteration(lpIteration_)
        > pModel_->getParameters().max_inactive_iteration_ &&
        col->getActivityRate(lpIteration_)
            < pModel_->getParameters().min_activity_rate_) {
      deletable.unchecked_push_back(i);
    }
  }

  if (pModel_->getParameters().printBranchStats_ && before_fathom) {
    std::cout << "ABOUT TO FATHOM CURRENT NODE" << std::endl;
    pModel_->printStats();
  }
}

/*
 * BcpBranchingTree
 */

BcpBranchingTree::BcpBranchingTree(BcpModeler *pModel) :
    pModel_(pModel), nbInitialColumnVars_(pModel->getActiveColumns().size()) {}

// setting the base
// Create the core of the problem by filling out the last three arguments.
void BcpBranchingTree::initialize_core(BCP_vec<BCP_var_core *> &vars,  // NOLINT
                                       BCP_vec<BCP_cut_core *> &cuts,  // NOLINT
                                       BCP_lp_relax *&matrix) {
  // initialize tm parameters
  set_param(BCP_tm_par::TmVerb_SingleLineInfoFrequency,
            pModel_->getFrequency());
  // always dive
  set_param(BCP_tm_par::UnconditionalDiveProbability, 1);
  // pModel_->getParameters().maxSolvingTimeSeconds_); ?
  set_param(BCP_tm_par::MaxRunTime, LARGE_TIME);
  for (const auto &entry : pModel_->getTmParameters())
    set_param(entry.first, entry.second);

  // build the problem
  BcpProblem pb(pModel_);

  // copy of the core variables
  vars.reserve(pb.colnum);
  for (int i = 0; i < pb.colnum; ++i)
    vars.push_back(new BcpCoreVar(*pb.vars[i]));

  // copy of the core cuts
  cuts.reserve(pb.rownum);
  for (int i = 0; i < pb.rownum; ++i)
    cuts.push_back(new BcpCoreCons(*pb.cons[i]));

  // copy matrix to the matrix of the solver
  matrix = new BCP_lp_relax;
  matrix->copyOf(*pb.matrix, &pb.obj[0], &pb.lb[0], &pb.ub[0],
                 &pb.lhs[0], &pb.rhs[0]);
}

// create the root node
// Create the set of extra variables and cuts that should be added
// to the formulation in the root node.
void BcpBranchingTree::create_root(BCP_vec<BCP_var *> &added_vars,
                                   BCP_vec<BCP_cut *> &added_cuts,
                                   BCP_user_data *&user_data) {
  added_vars.reserve(pModel_->getInitialColumns().size());
  for (MyVar *col : pModel_->getInitialColumns()) {
    auto *var = dynamic_cast<BcpColumn *>(col);
    if (!var)
      Tools::throwError("Bad variable casting.");
    // add BcpColumn which will be deleted by BCP
    added_vars.unchecked_push_back(var);
  }
  // clear vector as ownership as been transferred to BCP
  pModel_->clearInitialColumns();
}

BCP_solution *BcpBranchingTree::unpack_feasible_solution(BCP_buffer &buf) {  // NOLINT
  BCP_solution *sol = BCP_tm_user::unpack_feasible_solution(buf);

  if (pModel_->getVerbosity() >= 3 && sol->objective_value() <= LARGE_SCORE)
    std::cout << pModel_->getMaster()->costsConstrainstsToString() << std::endl;

  // store the solution
  pModel_->addBcpSol(sol);

  return sol;
}

// various initializations before a new phase (e.g., branching strategy)
void BcpBranchingTree::init_new_phase(int phase,
                                      BCP_column_generation &colgen,
                                      CoinSearchTreeBase *&candidates) {
  colgen = BCP_GenerateColumns;
  //   CoinSearchTreeCompareHighestIncrease::pModel = pModel_;
  switch (pModel_->getSearchStrategy()) {
    case BestFirstSearch:set_param(BCP_tm_par::TreeSearchStrategy, 0);
      pTree_ = new MyCompCoinSearchTree<CoinSearchTreeCompareBest>(pModel_);
      break;
    case BreadthFirstSearch:set_param(BCP_tm_par::TreeSearchStrategy, 1);
      pTree_ = new MyCompCoinSearchTree<CoinSearchTreeCompareBreadth>(pModel_);
      break;
    case DepthFirstSearch:set_param(BCP_tm_par::TreeSearchStrategy, 2);
      pTree_ = new MyCompCoinSearchTree<CoinSearchTreeCompareDepth>(pModel_);
      break;
    default:
      pTree_ = new MyCompCoinSearchTree<CoinSearchTreeCompareBest>(pModel_);
      break;
  }

  candidates = pTree_;
}


BcpProblem::BcpProblem(BcpModeler * pModel):
    rownum(pModel->getCoreCons().size()),
    colnum(pModel->getCoreVars().size()),
    lb(colnum), ub(colnum), obj(colnum),
    rhs(rownum), lhs(rownum), vars(colnum), cons(rownum) {
  // copy of the core variables
  for (int i = 0; i < colnum; ++i) {
    auto *var = dynamic_cast<BcpCoreVar *>(pModel->getCoreVars()[i]);
    if (!var)
      Tools::throwError("Bad variable casting.");
    // create a new BcpCoreVar which will be deleted by BCP
    vars[i] = var;
    lb[i] = var->getLB();
    ub[i] = var->getUB();
    obj[i] = var->getCost();
  }

  // copy of the core cuts
  for (int i = 0; i < rownum; ++i) {
    auto *cut = dynamic_cast<BcpCoreCons *>(pModel->getCoreCons()[i]);
    if (!cut)
      Tools::throwError("Bad constraint casting.");
    // create a new BcpCoreCons which will be deleted by BCP
    cons[i] = cut;
    lhs[i] = cut->getLhs();
    rhs[i] = cut->getRhs();
  }

  // build matrix
  matrix = pModel->buildCoinMatrix();

  // check columns size
  if (matrix->getNumCols() != colnum) {
    std::cerr << "Number of columns for CoinPackedMatrix="
              << matrix->getNumCols()
              << " and for the modeler=" << colnum << std::endl;
    Tools::throwError("CoinPackedMatrix does not have the same number "
                      "of columns than the modeler. ");
  }
  // check rows size
  if (matrix->getNumRows() != rownum) {
    std::cerr << "Number of rows for CoinPackedMatrix=" << matrix->getNumRows()
              << " and for the modeler=" << rownum << std::endl;
    Tools::throwError("CoinPackedMatrix does not have the same "
                      "number of rows than the modeler. ");
  }
}

BcpProblem::~BcpProblem() {
  delete matrix;
}

//-----------------------------------------------------------------------------
//
//  C l a s s   B c p M o d e l e r
//
// Specific class of modeler designed to use BCP.
// In particular, may methods defined in the BCP library need to be implemented
// here
//
//-----------------------------------------------------------------------------

BcpModeler::BcpModeler(MasterProblem *pMaster,
                       const char *name,
                       SolverType type) :
    CoinModeler(),
    pMaster_(pMaster),
    pBcp_(nullptr),
    primalValues_(0),
    dualValues_(0),
    reducedCosts_(0),
    lhsValues_(0),
    lastNbSubProblemsSolved_(0),
    lastMinReducedCost_(0),
    nbNodes_(0),
    solverType_(type) {
  pBcp_ = new BcpInitialize(this);
}

// destroy all the column in the solutions
BcpModeler::~BcpModeler() {
  clear();
  delete pBcp_;
}

void BcpModeler::clear() {
  deleteSolutions();
  CoinModeler::clear();
}

// copy the current active columns to keep a track of them after deleting BCP
void BcpModeler::copyActiveToInitialColumns() {
  initialColumnVars_.reserve(
      initialColumnVars_.size()+activeColumnVars_.size());
  for (MyVar *v : activeColumnVars_) {
    auto col = dynamic_cast<BcpColumn *>(v);
    col->resetBounds();
    initialColumnVars_.emplace_back(new BcpColumn(*col));
  }
  clearActiveColumns();
}

void BcpModeler::deleteSolutions() {
  bcpSolutions_.clear();
  for (MyVar *v : columnsInSolutions_)
    delete v;
  columnsInSolutions_.clear();
}

// solve the model
int BcpModeler::solve(bool relaxation) {
  // set stopAfterXSolution_ to 0 when solving relaxation
  if (relaxation) parameters_.stopAfterXSolution_ = 0;

  // create the root
  pTree_->createRootNode();

  // stopped to false
  stopped_ = false;

  // solve
  char **argv = NULL;
  int value = -1;
  try {
    value = bcp_main(0, argv, pBcp_);
    // set status to optimal if the status hasn't been set  for the moment
    // -> it implies that BCP exited normally after exploring the whole tree.
    if (getMaster()->status() == UNSOLVED)
      getMaster()->status(OPTIMAL);
  } catch (const Tools::NSException &e) {
    std::cerr << "Current LP solution:" << std::endl;
    std::cerr << pMaster_->currentSolToString() << std::endl;
    writeMPS("model");
    std::cerr << e.what() << std::endl;
    if (bcpSolutions_.empty()) getMaster()->status(INFEASIBLE);
    else
      getMaster()->status(FEASIBLE);
    throw;
  }

  // set to nullptr as deleted
  pBcp_->pLpModel_ = nullptr;

  // clear tree
  pTree_->clear();
  treeMapping_.clear();

  // clear columns variables as they have been deleted by BCP
  clearActiveColumns();

  return value;
}

// reinitialize all parameters and clear vectors
void BcpModeler::reset() {
  // reset parent model
  CoinModeler::reset();

  // reset all solutions related objects
  lastNbSubProblemsSolved_ = 0;
  lastMinReducedCost_ = 0;
  solHasChanged_ = false;

  obj_history_.clear();
  primalValues_.clear();
  dualValues_.clear();
  reducedCosts_.clear();
  lhsValues_.clear();

  // delete them
  deleteSolutions();
}

/*
 * Create core variable:
 *    var is a pointer to the pointer of the variable
 *    var_name is the name of the variable
 *    lb, ub are the lower and upper bound of the variable
 *    vartype is the type of the variable:
 *    VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
 */
int BcpModeler::createVar(MyVar **var,
                          const char *var_name,
                          int index,
                          double objCoeff,
                          double lb,
                          double ub,
                          VarType vartype,
                          const std::vector<double> &pattern,
                          double score) {
  *var = new BcpCoreVar(var_name, index, objCoeff, vartype, lb, ub, pattern);
  return 1;
}

int BcpModeler::createColumnVar(MyVar **var,
                                const char *var_name,
                                int index,
                                double objCoeff,
                                const std::vector<double> &pattern,
                                double dualObj,
                                double lb,
                                double ub,
                                VarType vartype,
                                double score) {
  *var = new BcpColumn(var_name,
                       index,
                       objCoeff,
                       pattern,
                       dualObj,
                       vartype,
                       lb,
                       ub);
  return 1;
}

/*
 * Create linear constraint:
 *   con is a pointer to the pointer of the constraint
 *   con_name is the name of the constraint
 *   lhs, rhs are the lower and upper bound of the constraint
 *   nonZeroVars is the number of non-zero coefficients to add to the constraint
 */
int BcpModeler::createCoinConsLinear(MyCons **con,
                                     const char *con_name,
                                     int index,
                                     double lhs,
                                     double rhs) {
  *con = new BcpCoreCons(con_name, index, lhs, rhs);
  return 1;
}

int BcpModeler::createCoinCutLinear(MyCons **con,
                                    const char *con_name,
                                    int index,
                                    double lhs,
                                    double rhs,
                                    const std::vector<int> &indexVars,
                                    const std::vector<double> &coeffs) {
  *con = new BcpBranchCons(con_name, index, lhs, rhs, indexVars, coeffs);
  return 1;
}

/*
 * Set the solution: this is where we set the active column variables
 */
void BcpModeler::setLPSol(const BCP_lp_result &lpres,
                          const BCP_vec<BCP_var *> &vars,
                          const BCP_vec<BCP_cut *> &cuts,
                          int lpIteration) {
  obj_history_.push_back(lpres.objval());
  solHasChanged_ = false;

  // copy the new arrays in the vectors for the core vars
  const int nbCoreVar = coreVars_.size();
  const int nbVar = vars.size();
  const int nbCons = cuts.size();

  // clear all
  primalValues_.clear();
  dualValues_.clear();
  reducedCosts_.clear();
  lhsValues_.clear();

  // assign value for core variables
  primalValues_.assign(lpres.x(), lpres.x() + nbVar);
  dualValues_.assign(lpres.pi(), lpres.pi() + nbCons);
  reducedCosts_.assign(lpres.dj(), lpres.dj() + nbVar);
  lhsValues_.assign(lpres.lhs(), lpres.lhs() + nbCons);

  clearActiveColumns();
  for (int i = nbCoreVar; i < nbVar; ++i) {
    auto *var = dynamic_cast<BcpColumn *>(vars[i]);
    if (primalValues_[i] > epsilon())
      var->addActiveIteration(lpIteration);  // update the different counters
    addActiveColumn(var, i);
  }
  // debug
  // checkActiveColumns(vars);
}

void BcpModeler::checkActiveColumns(const BCP_vec<BCP_var *> &vars) const {
  auto *shiftNode = dynamic_cast<ShiftNode*>(pTree_->getCurrentNode().get());
  if (shiftNode == nullptr) return;

  for (unsigned int i = coreVars_.size(); i < vars.size(); ++i) {
    auto *var = dynamic_cast<BcpColumn *>(vars[i]);
//    PPattern pat = pMaster_->getPattern(var);
    if (var->getUB() == 0 || var->ub() == 0
        || shiftNode->pNurse_->num_ != Pattern::nurseNum(var))
      continue;
    for (int k = Pattern::firstDay(var); k <= Pattern::lastDay(var); ++k)
      if (shiftNode->day_ == k &&
          find(shiftNode->forbiddenShifts_.begin(),
               shiftNode->forbiddenShifts_.end(),
               Pattern::shift(var, k))
              != shiftNode->forbiddenShifts_.end()) {
        PPattern pat = pMaster_->getPattern(var);
        std::cout << "problem: active column " << var->bcpind()
                  << " with forbidden shift " << pat->shift(k) << std::endl;
        std::cout << pat->toString() << std::endl;
        getchar();
      }
  }
}

void BcpModeler::addCurrentBcpSol(bool isInteger) {
  std::vector<MyVar *> vars;
  std::vector<double> values;

  for (MyVar *var : getCoreVars()) {
    double v = getVarValue(var);
    if (v > epsilon()) {
      vars.push_back(var);
      values.push_back(v);
    }
  }

  bool foundColumn = false;
  for (MyVar *var : getActiveColumns()) {
    double v = getVarValue(var);
    if (v > epsilon()) {
      vars.push_back(var);
      values.push_back(v);
      foundColumn = true;
    }
  }

  if (!foundColumn) {
    for (size_t i = 0; i < vars.size(); i++)
      std::cout << vars[i]->name_ << ": " << values[i] << std::endl;
    std::cout << pMaster_->coverageToString() << std::endl;
    std::cout << pMaster_->allocationToString() << std::endl;
    std::cout << pMaster_->currentSolToString() << std::endl;
    Tools::throwError("BCP found a solution without an active column !");
  }

  addBcpSol(getLastObj(), vars, values, isInteger);
}

void BcpModeler::addBcpSol(const BCP_solution *sol) {
  // if no integer solution is needed, don't store the solutions
  if (parameters_.stopAfterXSolution_ == 0)
    return;

  auto *sol2 = dynamic_cast<const BCP_solution_generic *>(sol);
  std::vector<BCP_var *> vars(sol2->_vars.begin(), sol2->_vars.end());
  std::vector<double> values(sol2->_values.begin(), sol2->_values.end());
  addBcpSol(sol2->objective_value(), vars, values);

  // check if should stop
  doStop();
}

template<typename T>
void BcpModeler::addBcpSol(double objValue,
                           const std::vector<T *> &vars,
                           const std::vector<double> &values,
                           bool isInteger) {
  MyBCPSolution mySol(isInteger);
  bool isArtificialSol = false;
  for (unsigned int i = 0; i < vars.size(); ++i) {
    auto *col = dynamic_cast<BcpColumn *>(vars[i]);
    if (col) {
      col = new BcpColumn(*col);
      // store pointers to be able to delete them
      columnsInSolutions_.push_back(col);
      mySol.add_entry(col, values[i]);
    } else {
      BcpCoreVar *var = getCoreVar(vars[i]);
      mySol.add_entry(var, values[i]);
      if (var->getCost() >= LARGE_SCORE - epsilon()) {
        if (values[i] > epsilon()) isArtificialSol = true;
      }
    }
  }

  if (!isArtificialSol) {
    bcpSolutions_.push_back(mySol);
    // if the solution improves the upper bound, record the new upper bound and
    // load the integer solution
    if (objValue + epsilon() < pTree_->getBestUB()) {
      pTree_->setBestUB(objValue);
      // print the solution in a text file
      if (parameters_.printEverySolution_ && isInteger) {
        solHasChanged_ = true;
        parameters_.saveFunction_->save(parameters_.weekIndices_,
                                        parameters_.outfile_);
      }
    }
  }

  // check if should stop
  doStop();
}

// Get the index of the best solution in the vector of solutions of BCP
int BcpModeler::getBestSolIndex(bool integer) const {
  int i = 0, index = -1;
  double bestObj = LARGE_SCORE;
  for (const MyBCPSolution &sol : bcpSolutions_) {
    if ((!integer || sol.isInteger_) && sol.objective_value() < bestObj) {
      bestObj = sol.objective_value();
      index = i;
    }
    ++i;
  }
  return index;
}

bool BcpModeler::loadBestSol(bool integer) {
  int index = getBestSolIndex(integer);
  if (index == -1)
    return false;

  loadBcpSol(index);
  return true;
}

// Clear the active column, set the active columns with those of the solution
// and set the primal values accordingly
void BcpModeler::loadBcpSol(int index) {
  indexBcpSol_ = index;
  const MyBCPSolution &sol = bcpSolutions_[index];
  // type of sol._vars[i] is either BcpColumn or BcpCoreVar
  // Get BcpColumn (after BcpCoreVar) and
  // map their index to their value in the primal vector
  clearActiveColumns();
  int coreSize = coreVars_.size(), colInd = coreSize;
  vector<double> primal(coreSize);  // fill primal vector with 0 for core vars
  for (int i = 0; i < sol._vars.size(); ++i) {
    // BcpCoreVar
    if (sol._vars[i]->bcpind() < coreSize) {
      primal[sol._vars[i]->bcpind()] = sol._values[i];
    } else {
      // BcpColumnVar
      primal.push_back(sol._values[i]);
      addActiveColumn(dynamic_cast<BcpColumn *>(sol._vars[i]), colInd++);
    }
  }
  setPrimal(primal);
}

/*
 * get/set the primal values
 */
double BcpModeler::getVarValue(MyVar *var) const {
  if (primalValues_.empty())
    Tools::throwError("Primal solution has not been initialized.");

  int index = getActiveIndex(var);
  if (index == -1) return 0;
  return primalValues_[index];
}

void BcpModeler::setVarValue(MyVar *var, double value) {
  if (primalValues_.empty())
    Tools::throwError("Primal solution has not been initialized.");

  int index = getActiveIndex(var);
  if (index == -1) return;
  primalValues_[index] = value;
}

/*
 * Get the dual variables
 */
double BcpModeler::getDual(MyCons *cons, bool transformed) const {
  if (dualValues_.empty())
    Tools::throwError("Dual solution has been initialized.");
  return dualValues_[cons->getIndex()];
}

/*
 * Get the reduced cost
 */
double BcpModeler::getReducedCost(MyVar *var) const {
  if (reducedCosts_.empty())
    Tools::throwError("Reduced cost solution has been initialized.");

  int index = getActiveIndex(var);
  if (index == -1) return 0;
  return reducedCosts_[index];
}

/**************
 * Parameters *
 *************/
int BcpModeler::setVerbosity(int v) {
  verbosity_ = v;

  if (v >= 1) {
    tm_parameters[BCP_tm_par::VerbosityShutUp] = 1;
    tm_parameters[BCP_tm_par::TmVerb_First] = 1;
    tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolutionValue] = 1;
    tm_parameters[BCP_tm_par::TmVerb_BestFeasibleSolution] = 1;
    tm_parameters[BCP_tm_par::TmVerb_NewPhaseStart] = 0;
    tm_parameters[BCP_tm_par::TmVerb_Last] = 1;
    tm_parameters[BCP_tm_par::DebugLpProcesses] = 1;

    // Just a marker for the last LpVerb
    lp_parameters[BCP_lp_par::LpVerb_Last] = 1;
    // Print information related to fathoming. (BCP_lp_main_loop,
    // BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
    lp_parameters[BCP_lp_par::LpVerb_FathomInfo] = 1;
  }

  if (v >= 2) {
    TmVerb_SingleLineInfoFrequency = 1;

    // Print the "Starting iteration x" line. (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_IterationCount] = 1;
    // Print the number of variables generated during this iteration.
    // (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_GeneratedVarCount] = 1;
    // Print information if receiving variables is timed out.
    // (BCP_lp_generate_vars)
    lp_parameters[BCP_lp_par::LpVerb_ReportVarGenTimeout] = 1;
    // Similar as above for variables. (BCP_lp_generate_vars)
    lp_parameters[BCP_lp_par::LpVerb_ReportLocalVarPoolSize] = 1;
    // Print the number of variables added from the local variable pool in the
    // current iteration. (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_AddedVarCount] = 1;
  }

  if (v >= 3) {
    tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolutionValue] = 1;
    tm_parameters[BCP_tm_par::TmVerb_PrunedNodeInfo] = 1;

    // After a branching object is selected print what happens to the presolved
    // children (e.g., fathomed). (BCP_print_brobj_stat)
    lp_parameters[BCP_lp_par::LpVerb_ChildrenInfo] = 1;
    // Print the number of variables generated before resolving the Lp
    // ir fathoming a node. (BCP_lp_fathom)
    lp_parameters[BCP_lp_par::LpVerb_ColumnGenerationInfo] = 1;
    // Print information related to fathoming.
    // (BCP_lp_main_loop, BCP_lp_perform_fathom, BCP_lp_branch) (BCP_lp_fathom)
    lp_parameters[BCP_lp_par::LpVerb_FathomInfo] = 1;
    // Print the size of the problem matrix and the LP solution value after
    // resolving the LP. (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_LpSolutionValue] = 1;
  }

  if (v >= 4) {
    tm_parameters[BCP_tm_par::TmVerb_BetterFeasibleSolution] = 1;
    tm_parameters[BCP_tm_par::TmVerb_AllFeasibleSolution] = 1;
    tm_parameters[BCP_tm_par::TmVerb_TimeOfImprovingSolution] = 1;
    tm_parameters[BCP_tm_par::TmVerb_TrimmedNum] = 1;
    tm_parameters[BCP_tm_par::ReportWhenDefaultIsExecuted] = 1;

    // Turn on the user hook "display_lp_solution" for the last LP relaxation
    // solved at a search tree node. (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_FinalRelaxedSolution] = 1;
    // Print out a message when the default version of an overridable method
    // is executed. Default: 1.
    lp_parameters[BCP_lp_par::ReportWhenDefaultIsExecuted] = 1;
    // Print the number of columns and rows that were deleted during matrix
    // compression. (BCP_lp_delete_cols_and_rows)
    lp_parameters[BCP_lp_par::LpVerb_MatrixCompression] = 1;
    // Print detailed information about all the branching candidates during
    // strong branching. LpVerb_PresolveResult must be set for this parameter
    // to have an effect. (BCP_lp_perform_strong_branching)
    lp_parameters[BCP_lp_par::LpVerb_PresolvePositions] = 1;
    // Print information on the presolved branching candidates during strong
    // branching. (BCP_lp_perform_strong_branching)
    lp_parameters[BCP_lp_par::LpVerb_PresolveResult] = 1;
    // Print the "Processing NODE x on LEVEL y" line. (BCP_lp-main_loop)
    lp_parameters[BCP_lp_par::LpVerb_ProcessedNodeIndex] = 1;
    // Print the number of variables whose bounds have been changed by reduced
    // cost fixing or logical fixing. (BCP_lp_fix_vars)
    lp_parameters[BCP_lp_par::LpVerb_VarTightening] = 1;
    // Print the number of ineffective rows in the current problem. The
    // definition of what rows are considered ineffective is determined by
    // the paramter IneffectiveConstraints. (BCP_lp_adjust_row_effectiveness)
    lp_parameters[BCP_lp_par::LpVerb_RowEffectivenessCount] = 1;
    // Print detailed information on the branching candidate selected by strong
    // branching. LpVerb_StrongBranchResult must be set fo this parameter to
    // have an effect. (BCP_print_brobj_stat)
    lp_parameters[BCP_lp_par::LpVerb_StrongBranchPositions] = 1;
    // Print information on the branching candidate selected by strong
    // branching. (BCP_print_brobj_stat)
    lp_parameters[BCP_lp_par::LpVerb_StrongBranchResult] = 1;
  }
  if (v >= 5) {
    // Turn on the user hook "display_lp_solution". (BCP_lp_main_loop)
    lp_parameters[BCP_lp_par::LpVerb_RelaxedSolution] = 1;
  }

  return 1;
}

// Check if BCP must stop
bool BcpModeler::doStop(const BCP_vec<BCP_var *> &vars) {
  // if already stooped
  if (isStopped())
    return true;
  // continue if doesn't have a lb
  if (pTree_->getBestLB() >= LARGE_SCORE - epsilon())
    return false;

  // check stopping criteria
  if (pTree_->getBestUB() - pTree_->getBestLB()
      < parameters_.optimalAbsoluteGap_ - epsilon()) {
//    // update with ub as optimal if solution is integer
//    if (getObjective() > std::max(getCurrentLB(), getBestUB()) - epsilon())
//      updateNodeLB(getBestUB());
    pMaster_->status(OPTIMAL);
    printf("BCP STOPPED: optimal absolute gap < %.2f.\n",
           parameters_.optimalAbsoluteGap_);
  } else if (getMaster()->timerTotal()->dSinceStart()
      > getParameters().maxSolvingTimeSeconds_) {
    pMaster_->status(TIME_LIMIT);
    printf("BCP STOPPED: Time has run out after %.2f s.\n",
           pMaster_->timerTotal()->dSinceStart());
  } else if (parameters_.solveToOptimality_) {
    // all the other criteria are not optimal -> stop here and return false
    return false;
  } else if (nbSolutions() >= parameters_.stopAfterXSolution_) {
    // check the number of solutions
    pMaster_->status(FEASIBLE);
    printf("BCP STOPPED: %d solutions have been founded.\n", nbSolutions());
  } else if (pTree_->getBestUB() - pTree_->getBestLB()
      < parameters_.maxAbsoluteGap_ - epsilon()) {
    pMaster_->status(FEASIBLE);
    printf("BCP STOPPED: absolute gap < %.2f.\n", parameters_.maxAbsoluteGap_);
  } else if (pTree_->getBestUB() - pTree_->getBestLB()
      < parameters_.minRelativeGap_ * pTree_->getBestLB() - epsilon()) {
    pMaster_->status(FEASIBLE);
    printf("BCP STOPPED: relative gap < %.2f.\n", parameters_.minRelativeGap_);
  } else if (pTree_->getBestUB() - pTree_->getBestLB()
      < parameters_.relativeGap_ * pTree_->getBestLB() - epsilon()) {
    // if the relative gap is small enough and if same incumbent
    // since the last dive, stop
    if (pTree_->getNbNodesSinceLastIncumbent()
        > parameters_.nbDiveIfMinGap_ * pTree_->getDiveLength()) {
      pMaster_->status(FEASIBLE);
      printf("BCP STOPPED: relative gap < %.2f and more than %d nodes without "
             "new incumbent.\n",
             parameters_.relativeGap_,
             parameters_.nbDiveIfMinGap_ * pTree_->getDiveLength());
    } else {
      return false;
    }
  } else if (nbSolutions() > 0) {
    // if(pTree_->getBestUB() - pTree_->getBestLB() <
    // 10.0 * pTree_->getBestLB()) {
    // if the relative gap is too big, wait 2 dives before stopping
    if (pTree_->getNbNodesSinceLastIncumbent()
        > parameters_.nbDiveIfRelGap_ * pTree_->getDiveLength()) {
      pMaster_->status(FEASIBLE);
      printf("BCP STOPPED: relative gap > %.2f and more than %d nodes without "
             "new incumbent.\n",
             parameters_.relativeGap_,
             parameters_.nbDiveIfRelGap_ * pTree_->getDiveLength());
    } else {
      return false;
    }
  } else {
    return false;
  }

  // stop
  stop();
  return true;
}

void BcpModeler::stop() {
  // copy active columns before being deleted
  copyActiveToInitialColumns();
  // remove all the nodes
  pBcp_->pTree_->clear();
  stopped_ = true;
}

/**************
 * Outputs *
 *************/

int BcpModeler::writeProblem(string fileName) const {
  return writeMPS(fileName) + writeLP(fileName);
}

int BcpModeler::writeLP(string fileName) const {
  return pBcp_->writeLP(fileName);
}

int BcpModeler::writeMPS(string fileName) const {
  return pBcp_->writeMPS(fileName);
}
