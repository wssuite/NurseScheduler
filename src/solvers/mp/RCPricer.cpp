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

#include "solvers/mp/RCPricer.h"

#include <algorithm>
#include <memory>
#include <string>

#include "solvers/mp/sp/CyclicRosterSP.h"
#include "solvers/mp/sp/RosterSP.h"
#include "solvers/mp/sp/RotationSP.h"
#include "solvers/mp/modeler/BcpModeler.h"
#ifdef BOOST
#include "solvers/mp/sp/rcspp/boost/LongRotationSP.h"
#include "solvers/mp/sp/rcspp/boost/RosterSP.h"
#endif


/* namespace usage */
using std::pair;
using std::set;
using std::vector;
using std::lock_guard;
using std::unique_lock;
using std::recursive_mutex;


//////////////////////////////////////////////////////////////
//
// R C   P R I C E R
//
//////////////////////////////////////////////////////////////

/* Constructs the pricer object. */
RCPricer::RCPricer(MasterProblem *pMaster,
                   const char *name,
                   const SolverParam &param) :
    MyPricer(name),
    pMaster_(pMaster),
    pScenario_(pMaster->pScenario()),
    nbDays_(pMaster->nDays()),
    pModel_(pMaster->pModel()),
    nursesToSolve_(pMaster->sortedLiveNurses()),
    minReducedCost_(0),
    rdm_(Tools::getANewRandomGenerator()) {
  /* sort the nurses */
  std::shuffle(nursesToSolve_.begin(), nursesToSolve_.end(), rdm_);
}

/* Destructs the pricer object. */
RCPricer::~RCPricer() {}

/******************************************************
 * Perform pricing
 ******************************************************/
// at root node, after_fathom and backtracked are true
vector<MyVar *> RCPricer::pricing(double bound,
                                  bool before_fathom,
                                  bool after_fathom,
                                  bool backtracked) {
  // reset the current strategies if needed
  // modify spMaxSolvingTimeRatioForRelaxation_ ->
  // sp parameters would be updated at the beginning of a node at minimum
  int index = pModel_->pTree()->getCurrentNode()->getIndex();
  SubProblemParam spParam = pModel_->getParameters().spParam_;
  double timeRatio =
      index == 0 ? spParam.spMaxSolvingTimeRatioForRelaxationAtRoot_ :
      spParam.spMaxSolvingTimeRatioForRelaxation_;
  spParam.spMaxSolvingTimeRatioForRelaxation_ = timeRatio;  // used as default
  for (const auto &p : subProblems_) {
    p.second->setSpMaxSolvingTimeRatioForRelaxation(timeRatio);
    if (after_fathom)  // first pricing for a new node
      p.second->updateParameters(pModel_->isInfeasible());
  }

  // get duals
  auto pDualCosts = std::make_shared<DualCosts>(pMaster_);
  pDualCosts->updateDuals();

  // Reset all columns, counters, etc.
  resetSolutions();

  // reshuffle from time to time based on dual solutions
  if (pModel_->getParameters().spParam_.spSortNursesBasedOnDuals_ &&
      rdm_() % 10 == 0) {
    vector<double> duals = pDualCosts->getMaxDualValues();
    std::sort(nursesToSolve_.begin(), nursesToSolve_.end(),
              [&duals](const PLiveNurse &pN1, const PLiveNurse &pN2) {
                return duals[pN1->num_] > duals[pN2->num_];
              });
  }

  // count and store the nurses whose subproblems produced rotations.
  vector<PLiveNurse> nursesSolved,
      nursesNoSolution, nursesComputeLB, nursesComputingLB,
      nursesSolvedOrder = nursesToSolve_;
  // flag if any shift has been  forbidden because of column disjoint feature
  bool disjointForbidden = false;
  // local thread pool (use all available thread)
  Tools::PThreadsPool pPool = Tools::ThreadsPool::newThreadsPool();

  int nSPBeingSolved = 0;
  vector2D<RCSolution> solutionsPerNurses(pMaster_->nNurses());
  vector<Tools::Job> allJobs;
  bool stop = false;
  while (true) {
    unique_lock<recursive_mutex> mLock(m_subproblem_);

    // check if any nurse's subproblem is waiting to be solved
    if (nursesToSolve_.empty()) {
      mLock.unlock();
      // wait for 1 thread to finish
      bool threadWasRunning = pPool->wait(1);
      // if still empty, solve the LB if no columns have been generated
      // otherwise, solve one nurse sub-problem
      mLock.lock();
      if (nursesToSolve_.empty()) {
        if (nursesComputeLB.empty() || nSPSolvedWithSuccess_ > 0) {
          // it there was still a thread running -> continue
          if (threadWasRunning)
            continue;
          // break if there is no thread running anymore
          else
            break;
        }
        // compute LB as no subproblem provided any solutions and
        // some were not solved until optimality
        nursesToSolve_ = nursesComputeLB;
        nursesComputingLB =
            Tools::appendVectors(nursesComputingLB, nursesComputeLB);
        nursesComputeLB.clear();
      }
    }

    // RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
    PLiveNurse pNurse = nursesToSolve_.front();

    // try next nurse if forbidden
    if (isNurseForbidden(pNurse->num_)) {
      nursesToSolve_.erase(nursesToSolve_.begin());
      nursesNoSolution.push_back(pNurse);
      continue;
    }

    // check if there is all sub problem pricing columns have been solved
    // and some columns have been found
    // -> break (no need to compute LBs for the moment)
    if (nbSPTried_ == pMaster_->nNurses() && nSPSolvedWithSuccess_ > 0) {
      stop = true;
      break;
    }

    // increase the number of nurses being solved
    ++nSPBeingSolved;

    // if there are already enough subproblems being solved, wait
    // Decrease a little efficiency when close to the max, but allows to be
    // determistic as the same number of subproblems is always solved whatever
    // the number of threads used
    mLock.unlock();
    while (nSPBeingSolved > pModel_->getParameters().nSubProblemsToSolve_) {
      // wait for 1 thread to finish before checking again
      if (!pPool->wait(1))
        // if there is no thread running anymore, break
        break;
    }

    // if the maximum number of subproblem solved is reached, break.
    if (nSPSolvedWithSuccess_ >= pModel_->getParameters().nSubProblemsToSolve_)
      break;

    // remove the nurse as being solved
    mLock.lock();
    nursesToSolve_.erase(nursesToSolve_.begin());
    mLock.unlock();

    // define function to be run by the pool
    // The Job that is defined below will be run in parallel.
    // We create a job for each nurse and then run it in the pool
    /**
     * Start of the job definition (part run in parallel)
     */
    Tools::Job job([pNurse, pDualCosts, spParam, &nSPBeingSolved, &nursesSolved,
        &nursesNoSolution, &nursesComputeLB, &nursesComputingLB,
        &solutionsPerNurses, &disjointForbidden, &bound, this](Tools::Job job) {
      // lock the pricer
      unique_lock<recursive_mutex> lock(m_subproblem_);

      // UPDATE NURSE FORBIDDEN SHIFTS
      bool locDisjointForbidden = disjointForbidden;
      set<pair<int, int> > nurseForbiddenShifts(
          nursesForbiddenShifts_[pNurse->num_]);
      nurseForbiddenShifts.insert(forbiddenShifts_.begin(),
                                  forbiddenShifts_.end());
      pModel_->addForbiddenShifts(pNurse, &nurseForbiddenShifts);

      // check if trying to compute an LB
      bool computingLB = nursesComputingLB.end() !=
          std::find(nursesComputingLB.begin(), nursesComputingLB.end(), pNurse);

      // increment tried only if not computing an LB
      if (!computingLB)
        ++nbSPTried_;
      lock.unlock();  // unlock before solving

      // BUILD OR RE-USE THE SUBPROBLEM
      SubProblem *subProblem = retrieveSubproblem(pNurse, spParam);

      // attach job
      subProblem->attachJob(job);

      // SOLVE THE PROBLEM
      subProblem->solve(pDualCosts, nurseForbiddenShifts, bound, computingLB);

      // RETRIEVE THE GENERATED ROTATIONS
      std::vector<RCSolution> solutions = subProblem->getSolutions();

#ifdef NS_DEBUG
//      for (RCSolution &sol : solutions)
//        std::cout << sol.toString() << std::endl;
#endif

#ifdef CTR
      for (RCSolution &sol : solutions) {
        // check if any shifts are forbidden
        for (const auto &p : nurseForbiddenShifts)
          if (p.first >= sol.firstDayId() && p.first <= sol.lastDayId()
              && sol.shift(p.first) == p.second) {
            std::cerr << "Generated pattern does not respect forbidden "
                         "(day, shift): " << p.first << "," << p.second
                      << " for nurse " << pNurse->name_ << std::endl;
            std::cerr << sol.toString() << std::endl;
          }
        subProblem->computeCost(pMaster_, &sol);
      }
      if (pModel_->getParameters().rcsppType_ == LABEL_SETTING &&
          subProblem->isLastRunOptimal()) {
        SubProblemParam par2(pMaster_->pModel()->getParameters());
        SubProblem *sub2 = nullptr;
#ifdef BOOST
        par2.strategyLevel_ =
            boostRCSPP::SubProblem::maxSubproblemStrategyLevel_;
        if (pModel_->getParameters().spType_ == ALL_ROTATION)
          sub2 =
              new boostRCSPP::RotationSP(pScenario_, nbDays_, pNurse, par2);
        else if (pModel_->getParameters().spType_ == ROSTER)
          sub2 =
              new boostRCSPP::RosterSP(pScenario_, nbDays_, pNurse, par2);
#endif
//        par2.rcsppBidirectional_ = !par2.rcsppBidirectional_;
//        sub2 = new RosterSP(pScenario_, 0, nbDays_, pNurse,
//                            pMaster_->getSPResources(pNurse), par2);

        if (sub2 != nullptr) {
          sub2->build();
          sub2->solve(pDualCosts,
                      nurseForbiddenShifts,
                      bound);
          std::vector<RCSolution> sols = sub2->getSolutions();
          delete sub2;

          RCSolution::sort(&solutions);
          RCSolution::sort(&sols);
          if (!sols.empty() || !solutions.empty()) {
            subProblem->computeCost(pMaster_, &sols.front());
            if (sols.empty() ^ solutions.empty()) {
              std::cout << "For nurse " << pNurse->name_ <<std::endl;
              if (!solutions.empty())
                std::cout << "Test: " << solutions.front().toString()
                          << std::endl;
              else
                std::cout << "CTR: " << sols.front().toString() << std::endl;
              Tools::throwError("One of the subproblem has found "
                                "a solution and the other not.");
            }
            double diff = sols.front().reducedCost() -
                solutions.front().reducedCost();
            if (diff > pMaster_->pModel()->epsilon()
                || diff < -pMaster_->pModel()->epsilon()) {
              std::cout << "All solutions found for nurse "
                        << pNurse->name_ << ":" << std::endl;
              for (const RCSolution &sol : solutions)
                std::cout << sol.toString() << std::endl;
              std::cout << std::endl << "Both best solutions for nurse "
                                     << pNurse->name_ << ":" << std::endl;
              std::cout << solutions.front().toString() << std::endl;
              std::cout << "CTR: " << sols.front().toString() << std::endl;
              std::cerr << "The subproblems haven't found the same best "
                           "reduced costs." << std::endl;
            }
          }
        }
      }
#endif

      // Lock the pricer
      lock.lock();

      // CHECK IF THE SUBPROBLEM GENERATED NEW RC SOLUTIONS
      // If yes, store the nurses
      if (!solutions.empty()) {
        ++nSPSolvedWithSuccess_;

        // update the reduced costs
        updateRedCost(subProblem, solutions, locDisjointForbidden);

        // store the solutions
        solutionsPerNurses[pNurse->num_] = solutions;

        // UPDATE FORBIDDEN SHIFTS
        if (pModel_->getParameters().isColumnDisjoint_)
          // set disjointForbidden to true,
          // if addForbiddenShifts returns true once
          if (addForbiddenShifts(solutions))
            disjointForbidden = true;

        // add the nurse to the nurses solved
        nursesSolved.push_back(pNurse);
      } else {
        // update reduced costs LB
        reducedCostLBs_[pNurse->num_] = subProblem->minRedCostLB();
        // decrease nSPBeingSolved as no solutions were produced
        --nSPBeingSolved;
        // if solved to optimality, put the nurse on the side
        if (subProblem->isLastRunOptimal() || computingLB) {
          // update the reduced costs
          updateRedCost(subProblem, solutions, locDisjointForbidden);
          nursesNoSolution.push_back(pNurse);
          // if computing LB, remove nurse from the vector
          if (computingLB)
            nursesComputingLB.erase(std::find(
                nursesComputingLB.begin(), nursesComputingLB.end(), pNurse));
        } else if (!subProblem->hasANextExecutionAvailable()) {
          // else if solver does not have any other execution available
          // if necessary to compute LB and is not a column node
          if (pModel_->getParameters().spParam_.spComputeLB_ &&
              !pModel_->isColumnsNode()) {
            nursesComputeLB.push_back(pNurse);
          } else {
            // otherwise, put the nurse on the side
            nursesNoSolution.push_back(pNurse);
            // update the reduced costs
            updateRedCost(subProblem, solutions, locDisjointForbidden);
          }
        } else {
          // otherwise, add it to being solved again if necessary
          // if not computing LBs -> add it to the back
          if (nursesComputingLB.empty())
            nursesToSolve_.push_back(pNurse);
          else  // add it to the front
            nursesToSolve_ = Tools::appendVectors({pNurse}, nursesToSolve_);
          --nbSPTried_;
        }
      }
    });  // END JOB
    /**
     * End of the job definition (part run in parallel)
     */

    // store job
    allJobs.push_back(job);

    // The job will be run in parallel.
    pPool->run(job);

    // if the maximum number of subproblems solved is reached, break.
    if (nSPSolvedWithSuccess_ >= pModel_->getParameters().nSubProblemsToSolve_)
      break;
  }

  // ask all jobs to stop
  if (stop)
    for (Tools::Job &job : allJobs) job.askStop();

  // wait for the threads to be finished
  pPool->wait();

  // ADD THE SOLUTIONS TO THE MASTER PROBLEM
  for (const auto &pN : nursesSolvedOrder)
    addColumnsToMaster(pN->num_, &solutionsPerNurses[pN->num_]);

  /* Add the nurses back into the nursesToSolve_ list
   * Reverse the vector before putting it back in
   * */

  // At least one subproblem unsolved
  if (!nursesToSolve_.empty()) optimal_ = false;

  // check if a lower bound has been computed for every sub problem
  lowerBounded_ = true;
  for (double v : reducedCostLBs_)
    if (v == -DBL_MAX) {
      lowerBounded_ = false;
      break;
    }

  // nurses order
  std::map<PLiveNurse, int> nursesOrder;
  int n = 0;
  for (const auto &pN : nursesSolvedOrder)
    nursesOrder[pN] = n++;

  nursesNoSolution = Tools::appendVectors(nursesNoSolution, nursesComputeLB);
  if (disjointForbidden) {
    // 1- Add the nurses in nursesNoSolution:
    // we hope that they will generate columns on next run as
    // some shifts were forbidden because of the column disjoint feature.
    reverseOrderPushBackNurses(&nursesNoSolution, nursesOrder);

    // 2- Add the nurses in nursesSolved
    reverseOrderPushBackNurses(&nursesSolved, nursesOrder);
  } else {
    // 1- Add the nurses in nursesSolved
    // (we hope that they will keep generating columns)
    reverseOrderPushBackNurses(&nursesSolved, nursesOrder);

    // 2- Add the nurses in nursesSolved
    // (we hope that they will keep not generating columns)
    reverseOrderPushBackNurses(&nursesNoSolution, nursesOrder);
  }

  // set statistics
  auto *pModel = dynamic_cast<BcpModeler *>(pModel_);
  pModel->setLastNbSubProblemsSolved(nbSPTried_);
  pModel->setLastMinReducedCost(minReducedCost_);

  return allNewColumns_;
}

void RCPricer::updateParameters(bool useMoreTime) {
  for (const auto &p : subProblems_) {
    p.second->setPreferences(pModel_->getParameters().spParam_);
    p.second->updateParameters(useMoreTime);
  }
}


void RCPricer::updateRedCost(
    SubProblem *pSP,
    const std::vector<RCSolution> &solutions,
    bool disjointForbidden) {
  lock_guard<recursive_mutex> lock(m_subproblem_);

  // if not on last level and disjointForbidden = false -> not optimal
  if (disjointForbidden || !pSP->isLastRunOptimal())
    optimal_ = false;

  // update reduced cost and min
  double redCost = solutions.empty() ? 0 : solutions.front().reducedCost();
  if (redCost < minReducedCost_)
    minReducedCost_ = redCost;

  // if not an optimal run, use LB (if any)
  if (!pSP->isLastRunOptimal())
    redCost = pSP->minRedCostLB();
  reducedCostLBs_[pSP->pLiveNurse()->num_] = redCost;
}

/******************************************************
 * add some forbidden shifts
 ******************************************************/
void RCPricer::initNursesAvailabilities() {
  Tools::initVector(&nursesForbiddenShifts_,
                    pMaster_->nNurses(),
                    std::set<std::pair<int, int>>());
  for (int n = 0; n < pMaster_->nNurses(); ++n)
    for (int k = 0; k < pMaster_->nDays(); ++k)
      for (int s = 0; s < pMaster_->nShifts(); ++s)
        if (!pMaster_->isNurseAvailableOnDayShift(n, k, s))
          nursesForbiddenShifts_[n].insert({k, s});
}

bool RCPricer::addForbiddenShifts(const std::vector<RCSolution> &solutions) {
  // search best rotation
  RCSolution bestSolution;
  double bestRedcost = DBL_MAX;
  for (const RCSolution &sol : solutions)
    if (sol.reducedCost() < bestRedcost) {
      bestRedcost = sol.reducedCost();
      bestSolution = sol;
    }

  // forbid working shifts of the best column if small enough
  if (bestRedcost > pModel_->getParameters().minReducedCostDisjoint_)
    return false;

  lock_guard<recursive_mutex> lock(m_subproblem_);
  bool shiftAdded = false;
  int k = bestSolution.firstDayId();
  for (const PShift &pS : bestSolution.pShifts()) {
    forbiddenShifts_.insert(pair<int, int>(k, pS->id));
    shiftAdded = true;
    k++;
  }
  return shiftAdded;
}

// Returns a pointer to the right subproblem
SubProblem *RCPricer::retrieveSubproblem(const PLiveNurse &pNurse,
                                         const SubProblemParam &spParam) {
  // lock the pricer
  unique_lock<recursive_mutex> lock(m_subproblem_);
  // Each nurse has a subproblem. If null, create a new one.
  SubProblem *subProblem = subProblems_[pNurse].get();
  if (subProblem == nullptr) {
    subProblem = buildSubproblem(pNurse, spParam);
    subProblems_[pNurse] = std::unique_ptr<SubProblem>(subProblem);
  }
  return subProblem;
}

SubProblem *RCPricer::buildSubproblem(const PLiveNurse &pNurse,
                                      const SubProblemParam &spParam) const {
  SubProblem *subProblem;
  switch (pModel_->getParameters().spType_) {
    case LONG_ROTATION:
      if (pModel_->getParameters().rcsppType_ == BOOST_LABEL_SETTING) {
#ifdef BOOST
        subProblem = new boostRCSPP::LongRotationSP(
                pScenario_, nbDays_, pNurse, spParam);
#else
        Tools::throwError("Boost subproblems are not compiled, "
                          "thus rcspp BOOST_LABEL_SETTING cannot be used.");
#endif
      } else {
        subProblem = new RotationSP(pScenario_, 0, nbDays_, pNurse,
                                    pMaster_->getSPResources(pNurse),
                                    spParam);
      }
      break;
    case ALL_ROTATION: {
      if (pModel_->getParameters().rcsppType_ == BOOST_LABEL_SETTING) {
#ifdef BOOST
        subProblem = new boostRCSPP::RotationSP(
                pScenario_, nbDays_, pNurse, spParam);
#else
        Tools::throwError("Boost subproblems are not compiled, "
                          "thus rcspp BOOST_LABEL_SETTING cannot be used.");
#endif
      } else {
        subProblem = new RotationSP(pScenario_, 0, nbDays_, pNurse,
                                    pMaster_->getSPResources(pNurse),
                                    spParam);
      }
      break;
    }
    case ROSTER: {
      if (pModel_->getParameters().rcsppType_ == BOOST_LABEL_SETTING) {
#ifdef BOOST
        subProblem = new boostRCSPP::RosterSP(
                pScenario_, nbDays_, pNurse, spParam);
#else
        Tools::throwError("Boost subproblems are not compiled, "
                          "thus rcspp BOOST_LABEL_SETTING cannot be used.");
#endif
      } else if (pScenario_->isCyclic()) {
        subProblem = new CyclicRosterSP(pScenario_, 0, nbDays_, pNurse,
                                        pMaster_->getSPResources(pNurse),
                                        spParam);
      } else {
        subProblem = new RosterSP(pScenario_, 0, nbDays_, pNurse,
                                  pMaster_->getSPResources(pNurse),
                                  spParam);
      }
      break;
    }
    default:
      Tools::throwError(
          "There is no subproblem defined associated to this type");
  }
  // then build the rcspp
  subProblem->build();
  return subProblem;
}

void RCPricer::computeCost(Column *col) {
  // get the live nurse of the column
  const PLiveNurse &pNurse = pMaster_->pLiveNurses()[col->nurseNum()];

  // get the subproblem associated to the live nurse
  SubProblem *pSP =
      retrieveSubproblem(pNurse, pModel_->getParameters().spParam_);

  // compute the cost of the column from scratch
  pSP->computeCost(pMaster_, col);
}

// Add the rotations to the master problem
int RCPricer::addColumnsToMaster(int nurseNum,
                                 std::vector<RCSolution> *solutions) {
  // SORT THE SOLUTIONS
  sortGeneratedSolutions(solutions);

  // SECOND, ADD THE ROTATIONS TO THE MASTER PROBLEM
  // (in the previously computed order)
  lock_guard<recursive_mutex> lock(m_subproblem_);
  int nbcolumnsAdded = 0;
  for (const RCSolution &sol : *solutions) {
    allNewColumns_.push_back(pMaster_->addColumn(nurseNum, sol));
#ifdef NS_DEBUG
//      std::cout << sol.toString() << std::endl;
    if (dynamic_cast<BcpColumn *>(allNewColumns_.back()) == nullptr) {
      std::cerr << "Should be a BcpColumn." << std::endl;
    }
#endif
    ++nbcolumnsAdded;
    if (nbcolumnsAdded >= pModel_->getParameters().spParam_.nbMaxColumnsToAdd_)
      break;
  }

#ifdef NS_DEBUG
  //  std::cerr << "--------------------------------"
//               "--------------------------------" << std::endl;
#endif

  return nbcolumnsAdded;
}

// Sort the rotations that just were generated for a nurse.
// Default option is sort by increasing reduced cost, but we
// could try something else (involving disjoint columns or shuffle for ex.)
void RCPricer::sortGeneratedSolutions(std::vector<RCSolution> *solutions) {
  RCSolution::sort(solutions);
  // shuffle the solutions with the best ones to ensure diversity of columns
  if (!solutions->empty()) {
    auto it = solutions->begin();
    double maxDualCost = solutions->front().reducedCost() * .9;
    while (it != solutions->end() && it->reducedCost() <= maxDualCost)
      it++;
    std::shuffle(solutions->begin()+1, it, rdm_);
  }
}

// ------------------------------------------
//
// PRINT functions
//
// ------------------------------------------
void RCPricer::printStatSPSolutions() {
  double tMeanSubproblems = timeInExSubproblems_ / nbExSubproblems_,
      tMeanS = timeForS_ / nbS_,
      tMeanNL = timeForNL_ / nbNL_,
      tMeanN = timeForN_ / nbN_;
  std::string sepLine = "+-----------------+------------------------"
                        "+-----------+\n";
  printf("\n");
  printf("%s", sepLine.c_str());
  printf("| %-15s |%10s %12s |%10s |\n", "type", "time", "number", "mean time");
  printf("%s", sepLine.c_str());
  printf("| %-15s |%10.2f %12d |%10.4f |\n",
         "Ex. Subproblems",
         timeInExSubproblems_,
         nbExSubproblems_,
         tMeanSubproblems);
  printf("| %-15s |%10.2f %12d |%10.4f |\n",
         "Short rotations",
         timeForS_,
         nbS_,
         tMeanS);
  printf("| %-15s |%10.2f %12d |%10.4f |\n",
         "NL rotations",
         timeForNL_,
         nbNL_,
         tMeanNL);
  printf("%s", sepLine.c_str());
  printf("| %-15s |%10.2f %12d |%10.4f |\n",
         "N rotations",
         timeForN_,
         nbN_,
         tMeanN);
  printf("%s", sepLine.c_str());
  printf("\n");
}
