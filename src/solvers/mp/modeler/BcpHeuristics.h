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

#ifndef SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_
#define SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_

#include <memory>
#include <string>
#include <vector>
#include <condition_variable>  // NOLINT (suppress cpplint error)
#include <mutex>  // NOLINT (suppress cpplint error)

#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RosterMP.h"
#include "solvers/mp/RotationMP.h"
#include "solvers/Solver.h"

#include "OsiSolverInterface.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_solution.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_vector.hpp"  // NOLINT (suppress cpplint error)
#include "BCP_var.hpp"  // NOLINT (suppress cpplint error)

class HeuristicThread {
 public:
  explicit HeuristicThread(MasterProblem *pMaster) :
  pMaster_(pMaster), run_lock_(run_mutex_) {
    run_lock_.unlock();
  }

  virtual ~HeuristicThread() {
    std::lock_guard<std::recursive_mutex> l(sol_mutex_);
    delete sol_;
  }

  virtual void safe_solve(const std::vector<BcpColumn *> &columns);

  bool running() {
    std::lock_guard<std::recursive_mutex> l(run_mutex_);
    return running_;
  }

  bool lock() {
    if (run_lock_.owns_lock())
      return false;
    run_lock_.lock();
    return true;
  }

  void unlock() { run_lock_.unlock(); }

  MasterProblem* pMaster() {
    lock();
    // cannot access master when not running or should be stopping
    if (!running_ || needStop_) return nullptr;
    return pMaster_;
  }

  virtual BCP_solution_generic *getEraseSolution() {
    std::lock_guard<std::recursive_mutex> l(sol_mutex_);
    BCP_solution_generic *sol = sol_;
    sol_ = nullptr;
    return sol;
  }

  virtual void setSolution(BCP_solution_generic *sol) {
    std::lock_guard<std::recursive_mutex> l(sol_mutex_);
    delete sol_;
    sol_ = sol;
  }

  virtual void askStop() {
    std::lock_guard<std::recursive_mutex> l(run_mutex_);
    needStop_ = true;
  }

  bool shouldStop() {
    std::lock_guard<std::recursive_mutex> l(run_mutex_);
    return needStop_;
  }

 private:
  std::recursive_mutex sol_mutex_, run_mutex_;
  std::unique_lock<std::recursive_mutex> run_lock_;
  MasterProblem *pMaster_;
  BCP_solution_generic *sol_ = nullptr;
  bool running_ = false;
  bool needStop_ = false;

 protected:
  virtual void solve(const std::vector<BcpColumn *> &columns) = 0;

  virtual bool start() {
    std::lock_guard<std::recursive_mutex> l(run_mutex_);
    if (running_) {
      std::cerr << "Heuristic can't be started, as it's already running."
                << std::endl;
      return false;
    }
    // if need to stop, do not start
    if (needStop_) return false;
    // start
    running_ = true;
    return true;
  }

  void stop() {
    std::lock_guard<std::recursive_mutex> l(run_mutex_);
    if (!running_)
      Tools::throwError("Heuristic can't be stopped, as it's not running or "
                        "is already stopped.");
    running_ = false;
    needStop_ = false;
  }
};

class HeuristicMIP : public HeuristicThread {
 public:
  explicit HeuristicMIP(
      MasterProblem *pMaster, SolverType type, int verbosity) :
      HeuristicThread(pMaster), nNurses_(pMaster->nNurses()),
      nColsToSelectPerNurses_(
      pMaster->pModel()->getParameters().MIPHeuristicNColumnsPerNurse_),
      maxNIterations_(
          pMaster->pModel()->getParameters().MIPHeuristicMaxIteration_) {
    // create initial solver
    createInitialSolver(pMaster, type, verbosity, &pInitialSolver_, &obj_);
  }

  virtual ~HeuristicMIP() {
    delete pInitialSolver_;
  }

  void buildBcpSolution(const std::vector<Stretch> &solution);

  std::vector<RotationColumn> computeWorkedRotations(
      const Stretch &st, int nurseNum) const;

  static void createInitialSolver(MasterProblem *pMaster,
                                  SolverType type,
                                  int verbosity,
                                  OsiSolverInterface **pSolver,
                                  std::vector<double> *obj);

  double safeComputeObjUB();  // unlock thread at the end

 protected:
  OsiSolverInterface *pInitialSolver_;
  std::vector<double> obj_;
  int nNurses_;
  int nColsToSelectPerNurses_;
  int maxNIterations_;
  // count the number of consecutive infeasible solutions
  // if too many (reaches maxConsInfeasible_),
  // skip some solve (maxSolveToSkip_) and retry
  // everytime maxSolveToSkip_ is reached,
  // restart counters and multiply by 2 the max
  int nConsInfeasible_ = 0;
  int maxConsInfeasible_ = 5;
  int nSolveToSkip_ = 0;
  int maxSolveToSkip_ = 1;

  void solve(const std::vector<BcpColumn *> &columns) override;

  bool start() override {
    if (!HeuristicThread::start()) return false;
    // if should skip some solving round, do not start
    if (nSolveToSkip_ > 0) {
      --nSolveToSkip_;
      stop();
      return false;
    }
    return true;
  }

  OsiSolverInterface *clone(OsiSolverInterface *pInitSolver) {
    // get a copy of the solver
    OsiSolverInterface *pSolver = pInitSolver->clone(true);
    pSolver->messageHandler()->logLevel(
        pInitSolver->messageHandler()->logLevel());
    return pSolver;
  }

  std::vector<BcpColumn*> selectColumns(const std::vector<BcpColumn*> &columns);

  void solveMP(const std::vector<BcpColumn*> &columns);

  bool solveMIP(OsiSolverInterface *pSolver,
                std::vector<double> obj,
                const std::vector<BcpColumn *> &columns);

  double computeObjUB();  // do not unlock thread at the end
};

class HeuristicRotation : public HeuristicMIP {
 public:
  HeuristicRotation(MasterProblem *pMaster, SolverType type, int verbosity);
  ~HeuristicRotation();

  void solve(const std::vector<BcpColumn *> &columns) override;

 protected:
  MasterProblem *pRotMP_ = nullptr;
  OsiSolverInterface *pInitialRotSolver_ = nullptr;
  std::vector<double> objRot_;

  std::vector<Stretch> solveRotationMP(
      const std::vector<BcpColumn*> &columns);

  std::vector<BcpColumn*> buildRotationColumns(
      const std::vector<BcpColumn*> &columns);

  Stretch computeRoster(std::vector<PColumn> rotations) const;

//  void addRotations(const std::vector<PColumn> &rotations);
};

class HeuristicSolver : public HeuristicMIP {
 public:
  HeuristicSolver(MasterProblem *pMaster,
                  Solver *pSolver,
                  const std::string &name,
                  SolverType type, int verbosity) :
      HeuristicMIP(pMaster, type, verbosity), pSolver_(pSolver), name_(name) {}

  void solve(const std::vector<BcpColumn *> &columns) override;

  BCP_solution_generic * getEraseSolution() override;

  void attachJob(const Tools::Job &job) {
    pSolver_->attachJob(job);
  }

  void askStop() override {
    HeuristicMIP::askStop();
    pSolver_->getJob().askStop();
  }

  double objValue() const { return pSolver_->objValue(); }

 protected:
  Solver *pSolver_;
  const std::string name_;
};

class BcpHeuristics {
 public:
  explicit BcpHeuristics(
      BcpModeler *pModel, bool useRotationModel = true,
      SolverType type = FirstAvailable, int verbosity = 0);

  ~BcpHeuristics();

  BCP_solution_generic * rounding(OsiSolverInterface *solver,
                                  const BCP_vec<BCP_var*> &vars);

  BCP_solution_generic * mip_solve();

  BCP_solution_generic * lns_solve();

  void stop() {
    if (mip_) mip_->askStop();
    if (lns_) lns_->askStop();
  }

  static BCP_solution_generic * buildSolution(
      OsiSolverInterface *solver,
      Modeler* pModel,
      const std::vector<BcpColumn*> &columns);

 protected:
  std::mutex mutex_;
  std::condition_variable cRunning_;
  BcpModeler *pModel_;
  shared_ptr<HeuristicMIP> mip_;
  shared_ptr<HeuristicSolver> lns_;
};

#endif  // SRC_SOLVERS_MP_MODELER_BCPHEURISTICS_H_
