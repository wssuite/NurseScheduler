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

#include "BcpHeuristics.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "CoinBuild.hpp"  // NOLINT (suppress cpplint error)


#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "cplex.h"  // NOLINT
#endif

#ifdef USE_GUROBI
#include "OsiGrbSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_CBC
#include "OsiCbcSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

using std::string;
using std::vector;
using std::map;
using std::pair;

BcpHeuristics::BcpHeuristics(
    BcpModeler *pModel, SolverType type, int verbosity):
  pModel_(pModel), mip_(std::make_shared<HeuristicMIP>()) {
  mip_->pInitialSolver_ = getNewSolver(type);
  if (mip_->pInitialSolver_ == nullptr)
    Tools::throwError(
        "%s is not linked to the executable or implemented.",
        getNameForEnum(SolverTypesByName, type).c_str());
//  if (type == CLP)
//    std::cout << "WARNING: using CLP for branch&bound. "
//                 "You should use at least CBC." << std::endl;
  mip_->pInitialSolver_->messageHandler()->setLogLevel(verbosity);
  BcpProblem pb(pModel);
  mip_->pInitialSolver_->loadProblem(
      *pb.matrix, &pb.lb[0], &pb.ub[0], &pb.obj[0],
      &pb.lhs[0], &pb.rhs[0]);
}

BcpHeuristics::~BcpHeuristics() {
  std::lock_guard<std::recursive_mutex> l(mip_->mutex_);
  mip_->deleted = true;
}

static inline bool compareCol(const pair<int, double> &p1,
                              const pair<int, double> &p2) {
  return (p1.second < p2.second);
}

BCP_solution_generic * BcpHeuristics::rounding(
    OsiSolverInterface *solver, const BCP_vec<BCP_var*> &vars) {
  // store the basis
  const CoinWarmStart *ws = solver->getWarmStart();

  // define different size
  const size_t size = vars.size(),
        coreSize = pModel_->getCoreVars().size();

  // store lower bounds
  map<int, double> indexColLbChanged;

  // REMARK
  // another possible, but more costly heuristic is simply to solve the MIP
  // involving all the columns currently in the problem, for this simply run
  // solver->branchAndBound()

  // while the solution is feasible
//  solver->resolve();
  BCP_solution_generic *sol = nullptr;
  while (solver->isProvenOptimal()) {
    // find the best not integer columns
    vector<pair<int, double>> candidates;
    for (size_t i = coreSize; i < size; ++i) {
      double value = solver->getColSolution()[i];
      if (value < pModel_->epsilon()
          || solver->getColLower()[i] == 1)  // value > 1 - pModel_->epsilon())
        continue;
      candidates.emplace_back(i, 1 - value);
    }
    stable_sort(candidates.begin(), candidates.end(), compareCol);

    // if we have found a column
    if (!candidates.empty()) {
      double valueLeft = .99;
      for (const pair<int, double> &p : candidates) {
        if (p.second > valueLeft)
          break;
        if (p.second > .2)
          valueLeft -= p.second;
        indexColLbChanged[p.first] = solver->getColLower()[p.first];
        solver->setColLower(p.first, 1);
      }
      solver->resolve();
    } else {
      std::vector<BcpColumn*> columns; columns.reserve(size - coreSize);
      for (int i = coreSize; i < size; ++i)
          columns.push_back(dynamic_cast<BcpColumn *>(vars[i]));
      sol = buildSolution(solver, columns);
      break;
    }
  }

  // restore bounds
  for (const auto &p : indexColLbChanged)
    solver->setColLower(p.first, p.second);

  //   // indicate to the lp solver that the heuristic branching is done
  //   solver->unmarkHotStart();
  solver->setWarmStart(ws);

  delete ws;

  return sol;
}

BCP_solution_generic * BcpHeuristics::mip_solve(
    OsiSolverInterface *solver,
    const BCP_vec<BCP_var*> &vars,
    SolverType type) {
  // lock and retrieve the current solution
  std::lock_guard<std::recursive_mutex> l(mip_->mutex_);
  BCP_solution_generic *sol = mip_->sol_;
  mip_->sol_ = nullptr;

  // check if solution improved Bcp
  if (sol && pModel_->getObjective() > sol->objective_value())
    std::cout << "The MIP heuristic has generated a solution of value "
              << sol->objective_value() << std::endl;

  // if running, return current solution
  stop();
  if (mip_->pSolver_) return sol;

  // else, start a new mip
  // Copy the active columns
  int colSize = pModel_->getActiveColumns().size();
  std::vector<BcpColumn*> columns; columns.reserve(colSize);
  for (auto * pVar : pModel_->getActiveColumns())
    columns.push_back(new BcpColumn(*dynamic_cast<BcpColumn *>(pVar)));

  double maxObj =
      pModel_->getObjective() - pModel_->getParameters().optimalAbsoluteGap_;
  auto mip = mip_;
  Tools::Job job = [columns, mip, maxObj, this]() {
    // get a clone of the solver
    mip_->clone();

    // stop as soon as a better solution has been found
    mip->pSolver_->setDblParam(OsiPrimalObjectiveLimit, maxObj);
    mip->pSolver_->setIntParam(OsiMaxNumIteration, 1e6);

    // add the columns to the new solver
    CoinBuild helper(1);
    std::vector<int> colIndices;
    colIndices.reserve(columns.size());
    int i = mip->pSolver_->getNumCols();
    for (auto *pCol : columns) {
      pCol->resetBounds();
      colIndices.push_back(i++);
      // add it to the solver
      helper.addCol(pCol->getNbRows(),
                    &pCol->getIndexRows().front(),
                    &pCol->getCoeffRows().front(),
                    pCol->lb(), pCol->ub(),
                    pCol->obj());
    }
    mip->pSolver_->addCols(helper);
    mip->pSolver_->setInteger(&*colIndices.begin(), columns.size());

#ifdef USE_CBC
    auto cbcSolver = dynamic_cast<OsiCbcSolverInterface *>(mip->pSolver_);
    if (cbcSolver) {
      auto *cbcModel = cbcSolver->getModelPtr();
      cbcModel->setLogLevel(mip->pSolver_->messageHandler()->logLevel());
    }
#endif

    // solve the MIP
    mip->pSolver_->branchAndBound();

    // build the solution if feasible
    BCP_solution_generic *sol = nullptr;
    std::lock_guard<std::recursive_mutex> l(mip->mutex_);
    if (!mip->deleted && !mip->pSolver_->isProvenPrimalInfeasible() &&
        !mip->pSolver_->isProvenDualInfeasible())
      sol = buildSolution(mip->pSolver_, columns);

    for (auto *pCol : columns) delete pCol;

    // if BcpHeuristics has been destroyed, delete everything
    delete mip->pSolver_;
    mip->pSolver_ = nullptr;
    mip->sol_ = sol;
  };

  // launch the job
  Tools::ThreadsPool::runOneJob(job);

  return sol;
}

void BcpHeuristics::stop() {
  std::lock_guard<std::recursive_mutex> l(mip_->mutex_);
  if (mip_->pSolver_) mip_->pSolver_->setIntParam(OsiMaxNumIteration, 0);
}


BCP_solution_generic * BcpHeuristics::buildSolution(
    OsiSolverInterface *solver, const std::vector<BcpColumn*> &columns) const {
  // define different size
  const size_t colSize = columns.size(),
        coreSize = pModel_->getCoreVars().size();
  // create a BCP_solution_generic to return
  auto *sol = new BCP_solution_generic();
  for (size_t i = 0; i < coreSize; ++i)
    if (solver->getColSolution()[i] > pModel_->epsilon()) {
        // create new var that will be deleted by the solution sol
        auto *var0 = dynamic_cast<BcpCoreVar *>(pModel_->getCoreVars()[i]);
        sol->add_entry(new BcpCoreVar(*var0), solver->getColSolution()[i]);
      }

  for (int i = 0; i < colSize; ++i)
    if (solver->getColSolution()[coreSize+i] > pModel_->epsilon()) {
      // create new var that will be deleted by the solution sol
        double v = solver->getColSolution()[coreSize+i];
      sol->add_entry(new BcpColumn(*columns[i]),
                     solver->getColSolution()[coreSize+i]);
    }

  return sol;
}

HeuristicMIP::~HeuristicMIP() {
  delete pSolver_;
  delete pInitialSolver_;
  delete sol_;
}

