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
#include <utility>
#include <vector>
#include <chrono>  // NOLINT (suppress cpplint error)

#include "solvers/mp/RotationMP.h"
#include "solvers/mp/RosterMP.h"

#include "CoinBuild.hpp"  // NOLINT (suppress cpplint error)

#ifdef USE_CBC
#include "OsiCbcSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "CbcEventHandler.hpp"  // NOLINT (suppress cpplint error)
#endif

#ifdef USE_CPLEX
#include "OsiCpxSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "cplex.h"  // NOLINT
SolverType BcpHeuristics::defaultSolverType = Cplex
#endif

#ifdef USE_GUROBI
#include "gurobi_c++.h"  // NOLINT
#include "OsiGrbSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#endif

using std::vector;
using std::map;
using std::pair;


void HeuristicThread::safe_solve(const std::vector<BcpColumn *> &columns) {
  if (!start()) return;
  // solve the roster MP
  try {
    solve(columns);
  } catch (const std::exception &e) {
    std::cout << "HeuristicThread::safe_solve() caught an exception=: "
              << e.what() << std::endl;
  }
  stop();
  // delete
  for (auto *pCol : columns) delete pCol;
}

BcpHeuristics::BcpHeuristics(
    BcpModeler *pModel, bool useRotationModel, SolverType type, int verbosity):
    pModel_(pModel),
    mip_() {
  if (pModel->getParameters().performMIPHeuristic_) {
    if (useRotationModel && pModel->getParameters().spType_ == ROSTER)
      mip_ = std::make_shared<HeuristicRotation>(
              pModel->pMaster_, type, verbosity);
    else
      mip_ = std::make_shared<HeuristicMIP>(
              pModel->pMaster_, type, verbosity);
  }
  if (pModel->getParameters().performLNSHeuristic_ &&
      !pModel_->pMaster_->pHeuristics().empty())
    lns_ = std::make_shared<HeuristicSolver>(
            pModel->pMaster_, pModel_->pMaster_->pHeuristics().front().get(),
            "LNS", type, verbosity);
}

BcpHeuristics::~BcpHeuristics() {
  stop();
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
      sol = buildSolution(solver, pModel_, columns);
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

BCP_solution_generic * BcpHeuristics::mip_solve() {
  std::unique_lock<std::mutex> l(mutex_);
  BCP_solution_generic *sol;

  // retrieve the current solution
  sol = mip_->getEraseSolution();

  // if running, return current solution
  if (mip_->running()) return sol;

  // else, start a new mip
  // Copy the active columns
  std::vector<BcpColumn *> columns;
  size_t colSize = pModel_->getActiveColumns().size();
  columns.reserve(colSize);
  for (auto *pVar : pModel_->getActiveColumns())
    columns.push_back(new BcpColumn(*dynamic_cast<BcpColumn *>(pVar)));

  auto mip = mip_;
  // take copies so they continue to live outside the scope
  Tools::Job job([columns, mip](Tools::Job job) {
    mip->safe_solve(columns);
  }, pModel_->getParameters().MIPHeuristicNThreads_);

  // launch the job
  Tools::ThreadsPool::runOneJob(job);

  // wait for the resolution to start
  cRunning_.wait_for(l, std::chrono::milliseconds(1), [mip]() {
    return mip->running();
  });

  return sol;
}

BCP_solution_generic * BcpHeuristics::lns_solve() {
  // one static field causing an issue is _cols_are_valid.
  std::cerr << "2 instances of BCP cannot run simultaneously because of many "
               "static fields. Therefore, BcpHeuristics::lns_solve() "
               "is ignored." << std::endl;
  return nullptr;
//  std::unique_lock<std::mutex> l(gMutex_);
//  // lock and retrieve the current solution
//  BCP_solution_generic *sol = lns_->getEraseSolution();
//
//  // if running, return current solution if any, otherwise continue
//  if (lns_->running()) {
//    // check if current solution is better than the one of the solver
//    if (lns_->objValue() >= pModel_->getObjective() + 1e-1)
//      lns_->askStop();  // stop current resolution to start a new one
//    return sol;
//  }
//
//  // else, start the solver
//  // store the current solution if any
//  if (pModel_->nbSolutions() > 0) {
//    // store current primal solution
//    auto primal = pModel_->getPrimal();
//    auto activeColumns = pModel_->getActiveColumns();
//    // store best solutions in the nurse rosters
//    pModel_->loadBestSol(true);
//    pModel_->pMaster_->storeSolution();
//    // reload current primal solution
//    pModel_->setPrimal(primal);
//    int colInd = pModel_->getCoreVars().size();
//    for (auto *pCol : activeColumns)
//      pModel_->addActiveColumn(pCol, colInd++);
//
//    // take copies so they continue to live outside the scope
//    auto lns = lns_;
//    Tools::Job job([lns]() { lns->safe_solve({}); });
//    lns_->attachJob(job);
//
//    // launch the job
//    Tools::ThreadsPool::runOneJob(job);
//
//    // wait for the resolution to start
//    cRunning_.wait_for(l, std::chrono::milliseconds(1), [lns]() {
//      return lns->running();
//    });
//  }
//
//  return sol;
}

void HeuristicMIP::buildBcpSolution(const std::vector<Stretch> &solution) {
  lock();
  // check if pMaster still alive
  if (shouldStop()) return;
  // fetch columns
  std::vector<BcpColumn*> columns;
  for (int n=0; n < solution.size(); ++n) {
    std::vector<Column> cols;
    if (pMaster()->pModel()->getParameters().spType_ == ROSTER) {
      RosterColumn rosterCol(RCSolution(solution[n]), n);
      pMaster()->computeColumnCost(&rosterCol);
      columns.push_back(dynamic_cast<BcpColumn *>(
                            pMaster()->createColumn(rosterCol, "roster")));
    } else {
      auto rotations = computeWorkedRotations(solution[n], n);
      for (auto &rot : rotations) {
        pMaster()->computeColumnCost(&rot);
        columns.push_back(dynamic_cast<BcpColumn *>(
                              pMaster()->createColumn(rot, "rotation")));
      }
    }
  }
  unlock();
  // solve the MIP with these columns
  solveMP(columns);
  // delete
  for (auto *pCol : columns) delete pCol;
}


void HeuristicMIP::createInitialSolver(
    MasterProblem *pMaster,
    SolverType type,
    int verbosity,
    OsiSolverInterface **pSolver,
    std::vector<double> *obj) {
  *pSolver = getNewSolver(type);
  if (*pSolver == nullptr)
    Tools::throwError(
        "%s is not linked to the executable or implemented.",
        namesBySolverType.at(type).c_str());
  (*pSolver)->messageHandler()->setLogLevel(verbosity);
  BcpProblem pb(dynamic_cast<BcpModeler *>(pMaster->pModel()));
// add objective as a cut
  *obj = pb.obj;
  (*pSolver)->loadProblem(
      *pb.matrix, &pb.lb[0], &pb.ub[0], &pb.obj[0],
      &pb.lhs[0], &pb.rhs[0]);
}

void HeuristicMIP::solve(const std::vector<BcpColumn *> &columns) {
  solveMP(columns);
}

void HeuristicMIP::solveMP(const std::vector<BcpColumn*> &columns) {
  // select rosters
  std::vector<BcpColumn *> selectedCols = selectColumns(columns);
  // get a clone of the solver
  OsiSolverInterface *pSolver = clone(pInitialSolver_);
  // solve the MIP with the columns
  bool feasible = solveMIP(pSolver, obj_, selectedCols);
  // build the solution if feasible
  lock();
  if (!shouldStop() && feasible) {
    BCP_solution_generic * sol = BcpHeuristics::buildSolution(
        pSolver, pMaster()->pModel(), selectedCols);
    std::cout << "The MIP heuristic has generated a solution of value "
              << sol->objective_value() << std::endl;
    setSolution(sol);
  }
  unlock();
  delete pSolver;
}

std::vector<BcpColumn*> HeuristicMIP::selectColumns(
    const std::vector<BcpColumn*> &columns) {
  // WARNING: do not use master without locking the mutex
  // select a subset of columns for each nurse
  vector2D<BcpColumn*> colsPerNurse(nNurses_);
  vector<BcpColumn*> selectedCols;
  selectedCols.reserve(nNurses_ * nColsToSelectPerNurses_);
  for (auto pCol : columns)
    colsPerNurse.at(Column::nurseNum(pCol)).push_back(pCol);
  for (auto & cols : colsPerNurse) {
    if (cols.size() > nColsToSelectPerNurses_) {
      std::stable_sort(
          cols.begin(), cols.end(),
          [](BcpColumn *pCol1, BcpColumn *pCol2) {
            // keep the most recent one used
            if (pCol1->getLastActive() > pCol2->getLastActive()) return true;
            if (pCol2->getLastActive() > pCol1->getLastActive()) return false;
            // then the cheapest
            if (pCol1->getCost() <= pCol2->getCost() - 1e-5)
              return true;
            if (pCol2->getCost() <= pCol1->getCost() - 1e-5)
              return false;
            // keep the latest one created
            return pCol1->getIndex() > pCol2->getIndex();
          });
      cols.resize(nColsToSelectPerNurses_);
    }
    selectedCols = Tools::appendVectors(selectedCols, cols);
  }
  return selectedCols;
}

#ifdef USE_CBC
// create a callback to stop if the solution is below the UB
class SolCallback : public CbcEventHandler {
 public:
  explicit SolCallback(HeuristicMIP *pSolver, CbcModel *model = nullptr):
  CbcEventHandler(model), pSolver_(pSolver) {}

  CbcAction event(CbcEvent whichEvent) override {
    double ub = pSolver_->safeComputeObjUB();
    if (!model_->parentModel()) {
      if (whichEvent == solution || whichEvent == heuristicSolution) {
        if (model_->getObjValue() <= ub)
          return stop;  // say finished
      }
    }
    double lb = model_->getBestPossibleObjValue();
    if (lb >= ub)
      return stop;  // say finished
    return noAction;  // carry on
  }

  CbcEventHandler *clone() const override {
    return new SolCallback(pSolver_, model_);
  }

 private:
  HeuristicMIP *pSolver_;
};
#endif

bool HeuristicMIP::solveMIP(OsiSolverInterface *pSolver,
                            std::vector<double> obj,
                            const std::vector<BcpColumn*> &columns) {
  try {
    // add the columns to the new solver
    CoinBuild helper(1);
    std::vector<int> colIndices;
    colIndices.reserve(columns.size());
    // add objective as a cut
    obj.reserve(obj.size() + columns.size());
    int i = pSolver->getNumCols();
    for (auto pCol : columns) {
      pCol->resetBounds();
      colIndices.push_back(i++);
      // add it to the solver
      helper.addCol(pCol->getNbRows(),
                    &pCol->getIndexRows().front(),
                    &pCol->getCoeffRows().front(),
                    pCol->lb(), pCol->ub(),
                    pCol->obj());
      obj.push_back(pCol->obj());
    }
    pSolver->addCols(helper);
    pSolver->setInteger(&*colIndices.begin(), columns.size());

    // set parameters and objective limit if necessary
    lock();
    if (shouldStop()) return false;
    auto *pModel = dynamic_cast<BcpModeler*>(pMaster()->pModel());
    double ub = computeObjUB();
    double maxUb = pModel->getLastObj() *
        (1 + pModel->getParameters().MIPHeuristicGapLimit_);
    if (ub >= maxUb) ub = maxUb;

    // add a constraint to bound the objective
    std::vector<int> elements;
    elements.reserve(obj.size());
    std::vector<double> obj2;
    obj2.reserve(obj.size());
    for (int j = 0; j < obj.size(); ++j)
      if (abs(obj[j]) > pModel->epsilon()) {
        elements.push_back(j);
        obj2.push_back(obj[j]);
      }
    double lhs = -pSolver->getInfinity();
    pSolver->addRow(obj2.size(), &elements[0], &obj2[0], lhs, ub);

    // stop as soon as a better solution has been found
    if (pModel->getParameters().MIPHeuristicObjLimit_)
      pSolver->setDblParam(OsiPrimalObjectiveLimit, ub);
    // maxNIterations_
    int maxIt = maxNIterations_ >= 0 ? maxNIterations_ :
        10 * nNurses_ * nColsToSelectPerNurses_;
    pSolver->setIntParam(OsiMaxNumIteration, maxIt);

#ifdef USE_CBC
    auto cbcSolver = dynamic_cast<OsiCbcSolverInterface *>(pSolver);
    if (cbcSolver) {
      auto *cbcModel = cbcSolver->getModelPtr();
      cbcModel->setLogLevel(pSolver->messageHandler()->logLevel());
      // use a callback to check the LB/UB of CBC vs the current UB
      // it's not necessary, but allows to stop CBC as soon as necessary
      std::shared_ptr<SolCallback> pCall = std::make_shared<SolCallback>(this);
      cbcModel->passInEventHandler(pCall.get());
    }
#endif

#ifdef USE_GUROBI
    auto grbSolver = dynamic_cast<OsiGrbSolverInterface *>(pSolver);
    if (grbSolver) {
      GRBmodel *grbMod =
          grbSolver->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL);
      GRBenv *grbEnv = GRBgetenv(grbMod);
      GRBsetintparam(grbEnv, GRB_INT_PAR_THREADS,
                     pModel->getParameters().MIPHeuristicNThreads_);
      if (pModel->getParameters().MIPHeuristicObjLimit_)
        GRBsetdblparam(grbEnv, GRB_DBL_PAR_BESTOBJSTOP, ub);
      GRBsetdblparam(grbEnv, GRB_DBL_PAR_MIPGAPABS,
                     pModel->getParameters().absoluteGap_);
      // use a callback to check the LB of gurobi vs the current UB
      // it's not necessary, but allows to stop gurobi as soon as necessary
      auto cb = [](GRBmodel *model, void *cbdata, int where, void *usrdata) {
        auto *pH = reinterpret_cast<HeuristicRotation *>(usrdata);
        if (pH->shouldStop()) return GRB_ERROR_CALLBACK;
        if (where == GRB_CB_MIPNODE) {
          double lb;
          if (GRBcbget(cbdata, where, GRB_CB_MIPNODE_OBJBND,
                       reinterpret_cast<void *>(&lb)))
            return 0;  // do nothing as an error has been encountered
          double ub = pH->safeComputeObjUB();
          if (lb >= ub) return GRB_ERROR_CALLBACK;
          return 0;
        }
        return 0;
      };
      GRBsetcallbackfunc(grbMod, cb, this);
    }
#else
    // ensure termination of the solver by setting a maximum number of
    // iterations
    if (pModel->getParameters().MIPHeuristicMaxIteration_ > 10e9)
        pSolver->setIntParam(OsiMaxNumIteration, 10e6);
#endif
    pModel = nullptr;
    unlock();

    // solve the MIP
    pSolver->branchAndBound();
    bool feasible =  !pSolver->isProvenPrimalInfeasible() &&
        !pSolver->isProvenDualInfeasible() && (pSolver->getObjValue() <= ub);
    if (feasible) {
      nConsInfeasible_ = 0;
      return true;
    }
  } catch (...) {}

  // the solution is infeasible if reaches here
  // increase the number of rosters per nurses if necessary
  if (nNurses_ * nColsToSelectPerNurses_ <= columns.size()) {
    nColsToSelectPerNurses_ *= 2;
  } else if (pSolver->isIterationLimitReached()) {
    maxNIterations_ *= 2;
  } else if (nConsInfeasible_++ > maxConsInfeasible_) {
    // update nConsInfeasible_
    nConsInfeasible_ = 0;
    nSolveToSkip_ = maxSolveToSkip_;
    maxSolveToSkip_ *= 2;
  }
  return false;
}

std::vector<RotationColumn> HeuristicMIP::computeWorkedRotations(
    const Stretch &st, int nurseNum) const {
  std::vector<RotationColumn> rotations;
  Stretch rot;
  int k = st.firstDayId();
  for (const PShift &pS : st.pShifts()) {
    // push back shift
    if (pS->isWork()) {
      // init if first one
      if (rot.nDays() == 0)
        rot = Stretch(k, pS);
      else
        rot.pushBack(pS);
    } else if (rot.nDays() > 0) {
      // save rotation
      rotations.emplace_back(RCSolution(rot), nurseNum);
      // clear rot
      rot = Stretch();
    }
    ++k;
  }

#ifdef NS_DEBUG
  if (k != st.nDays())
    Tools::throwError("Not all days have been visited as k=%d "
                      "for this roster: %s", k, st.toString().c_str());
#endif

  // save last rotation
  if (rot.nDays() > 0)
    rotations.emplace_back(RCSolution(rot), nurseNum);

  return rotations;
}

double HeuristicMIP::safeComputeObjUB() {
  lock();
  double ub = computeObjUB();
  unlock();
  return ub;
}

double HeuristicMIP::computeObjUB() {
  if (pMaster() == nullptr)
    return INFEAS_COST;
  auto *pModel = pMaster()->pModel();
  double ub = pModel->getObjective()
      - pModel->getParameters().absoluteGap_
      + 10 * pModel->epsilon();
  return ub;
}

HeuristicRotation::HeuristicRotation(
    MasterProblem *pMaster, SolverType type, int verbosity):
    HeuristicMIP(pMaster, type, verbosity) {
  // create RotationMP
  if (type == FirstAvailable) type = getFirstSolverTypeAvailable();
  auto param = pMaster->pModel()->getParameters();
  param.spType_ = ALL_ROTATION;
  param.spParam_.rcsppBidirectional_ = false;
  pRotMP_ = new RotationMP(pMaster->pScenario(), type, param);
  createInitialSolver(pRotMP_, type, verbosity,
                      &pInitialRotSolver_, &objRot_);
}

HeuristicRotation::~HeuristicRotation() {
  delete pInitialRotSolver_;
  delete pRotMP_;
}

void HeuristicRotation::solve(const std::vector<BcpColumn *> &columns) {
  // select the roster columns
  std::vector<BcpColumn*> allRosters = selectColumns(columns);
  // solve the rotation MP
  std::vector<Stretch> solution = solveRotationMP(allRosters);

  if (!solution.empty()) buildBcpSolution(solution);
}

std::vector<Stretch> HeuristicRotation::solveRotationMP(
    const std::vector<BcpColumn*> &columns) {
  // build rotation columns
  auto *pModel = pRotMP_->pModel();
  auto param = pModel->getParameters();
  std::vector<BcpColumn*> rotColumns = buildRotationColumns(columns);

  // solve the roster MP problem
  OsiSolverInterface *pSolver = clone(pInitialRotSolver_);
  bool feasible = solveMIP(pSolver, objRot_, rotColumns);

  // check if infeasible or should stop
  if (!feasible || shouldStop()) {
    for (auto pCol : rotColumns) delete pCol;
    delete pSolver;
    return {};
  }

  // fetch the rotation in the solution for each nurse
  std::vector<vector<PColumn>> nurseColumns(pRotMP_->nNurses());
  size_t coreSize = pModel->getCoreVars().size();
  // create a BCP_solution_generic to return
  int nCol = pSolver->getNumCols();
  bool infeasible = false;
  for (size_t i = coreSize; i < nCol; ++i) {
    double vCol = pSolver->getColSolution()[i];
    MyVar *pCol = rotColumns.at(i - coreSize);
    if (vCol > pModel->epsilon()) {
      if (abs(1 - vCol) > pModel->epsilon()) {
#ifdef NS_DEBUG
        Tools::throwError("HeuristicRotation::solveRotationMP: Value of the "
                          "rotation column should be equal to 1 and not %.2f",
                          vCol);
#endif
        infeasible = true;
      }
      auto pRot = pRotMP_->getPColumn(pCol);
      nurseColumns.at(pRot->nurseNum()).push_back(pRot);
    }
  }
  for (auto pCol : rotColumns) delete pCol;
  delete pSolver;

  // check if need to stop
  if (infeasible || shouldStop()) return {};

  // transform the rotations back to a roster
  std::vector<Stretch> solution;
  solution.reserve(nurseColumns.size());
  for (auto &rotations : nurseColumns)
    solution.push_back(computeRoster(rotations));

  return solution;
}

std::vector<BcpColumn*> HeuristicRotation::buildRotationColumns(
    const std::vector<BcpColumn*> &columns) {
  // separate rosters into rotations and keep unique rotations for each nurse
  std::vector<BcpColumn*> rotColumns;
  vector4D<RotationColumn> rotationCols;  // per nurse, first day, last day
  Tools::initVector3D(
      &rotationCols, pRotMP_->nNurses(), pRotMP_->nDays(), pRotMP_->nDays());
  for (BcpColumn *pCol : columns) {
    PColumn roster =
        std::make_shared<RosterColumn>(pCol, pRotMP_->pScenario());
    std::vector<RotationColumn> rotations =
        computeWorkedRotations(*roster, roster->nurseNum());
    for (auto &rot : rotations) {
      // check if rotation already present
      bool present = false;
      auto &rotCols = rotationCols.at(rot.nurseNum()).at(
          rot.firstDayId()).at(rot.nDays() - 1);
      for (const RotationColumn &rot0 : rotCols) {
        present = true;
        auto it = rot.pShifts().begin();
        for (const auto &pS : rot0.pShifts()) {
          if (!pS->equals(**it)) {
            present = false;
            break;
          }
          ++it;
        }
        if (present)
          break;
      }
      if (!present) {
        pRotMP_->computeColumnCost(&rot);
        MyVar *pRotCol = pRotMP_->createColumn(rot, "rotation");
        rotColumns.push_back(dynamic_cast<BcpColumn *>(pRotCol));
        rotCols.push_back(rot);
      }
    }
  }

  return rotColumns;
}

Stretch HeuristicRotation::computeRoster(std::vector<PColumn> rotations) const {
  std::stable_sort(
      rotations.begin(), rotations.end(),
      [](const PColumn &rot1, const PColumn &rot2) {
#ifdef NS_DEBUG
        if (rot1->nurseNum() != rot2->nurseNum())
          Tools::throwError("Two rotations for different nurses are used to "
                            "generate a roster: %s and %s",
                            rot1->toString().c_str(), rot2->toString().c_str());
        if (rot1->firstDayId() <= rot2->firstDayId() &&
            rot1->lastDayId() >= rot2->firstDayId())
          Tools::throwError("Two overlapping rotations are used to generate a "
                            "roster: %s and %s",
                            rot1->toString().c_str(), rot2->toString().c_str());
#endif
        return rot1->firstDayId() <= rot2->firstDayId();
      });

  int nDays = pRotMP_->nDays();
  const PShift &pRest = pRotMP_->pScenario()->pRestShift();
  int nextRotationFirstDay = rotations.empty() ?
                             nDays : rotations.front()->firstDayId();
  // stretch for the future roster containing the first shift
  std::vector<PShift> shifts(nextRotationFirstDay, pRest);
  shifts.reserve(nDays);
  for (auto itRot = rotations.begin(); itRot != rotations.end();) {
    // add work days shifts
#ifdef NS_DEBUG
    if (shifts.size() != (*itRot)->firstDayId())
      Tools::throwError("Rotation should start on %d", shifts.size());
#endif
    shifts = Tools::appendVectors(shifts, (*itRot)->pShifts());

    // insert rest days shifts between rotations or until the end
    itRot++;
    nextRotationFirstDay =
        itRot == rotations.end() ? nDays : (*itRot)->firstDayId();
#ifdef NS_DEBUG
    if (shifts.size() != nDays && shifts.size() >= nextRotationFirstDay)
      Tools::throwError("Rotation should start on at least one day after the "
                        "previous one (%d) instead of before (%d)",
                        shifts.size(), nextRotationFirstDay);
#endif
    std::vector<PShift> restShifts(nextRotationFirstDay - shifts.size(), pRest);
    shifts = Tools::appendVectors(shifts, restShifts);
  }

#ifdef NS_DEBUG
  if (shifts.size() != nDays)
      Tools::throwError("Rotation should end on %d", nDays - 1);
#endif

  return {0, shifts};
}

void HeuristicSolver::solve(const std::vector<BcpColumn *> &columns) {
  // build an initial solution
  std::vector<Roster> rosters;
  lock();
  if (shouldStop()) {
    for (const PLiveNurse &pN : pMaster()->pLiveNurses()) {
      // if no valid solution, give an empty solution
      if (pN->roster_.nDays() != pMaster()->nDays() ||
          pN->roster_.duration() == 0) {
        rosters.clear();
        break;
      }
      rosters.push_back(pN->roster_);
    }
    unlock();
    // solve
    double obj = pSolver_->solve(rosters);
  } else {
    unlock();
  }
}

BCP_solution_generic * HeuristicSolver::getEraseSolution() {
  lock();
  // check if has a better solution
  if (pMaster() && pSolver_->objValue() >=
      pMaster()->pModel()->getObjective() - 1e-1) {
    unlock();
    return nullptr;
  }
  std::cout << "HeuristicMIP " << name_ << " found a better solution of cost: "
            << pSolver_->objValue() << std::endl;
  // retrieve solution
  std::vector<Roster> solution = pSolver_->solution();
  if (solution.empty()) return nullptr;
  // get BCP solution
  std::vector<Stretch> stretches;
  stretches.resize(solution.size());
  for (const Roster &r : solution) stretches.push_back(r);
  buildBcpSolution(stretches);
  unlock();
  // run parent method
  return HeuristicMIP::getEraseSolution();
}

BCP_solution_generic * BcpHeuristics::buildSolution(
    OsiSolverInterface *solver,
    Modeler* pModel,
    const std::vector<BcpColumn*> &columns) {
  // check feasibility
  for (MyVar *var : pModel->getFeasibilityCoreVars())
    if (solver->getColSolution()[var->getIndex()] > pModel->epsilon()) {
      std::cerr << "HeuristicMIP found a infeasible solution of value: "
                << solver->getObjValue() << std::endl;
      return nullptr;
    }

  // define different size
  const size_t coreSize = pModel->getCoreVars().size();
  // create a BCP_solution_generic to return
  int nCol = solver->getNumCols();
  auto *sol = new BCP_solution_generic();
  for (size_t i = 0; i < nCol; ++i) {
    double vCol = solver->getColSolution()[i];
    if (vCol > pModel->epsilon()) {
      // create new var that will be deleted by the solution sol
      BCP_var* pV;
      if (i < coreSize)
        pV = new BcpCoreVar(*dynamic_cast<BcpCoreVar*>(
            pModel->getCoreVars().at(i)));
      else
        pV = new BcpColumn(*columns.at(i - coreSize));
#ifdef NS_DEBUG
      if (abs(vCol - round(vCol)) > pModel->epsilon())
        Tools::throwError("BcpHeuristics: Column value should be integer "
                          "and not equals to %.2f", vCol);
#endif
      sol->add_entry(pV, vCol);
    }
  }

  return sol;
}
