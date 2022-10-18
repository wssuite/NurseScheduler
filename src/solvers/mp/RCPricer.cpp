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
#include "solvers/mp/sp/rcspp/boost/LongRotationSP.h"
#include "solvers/mp/sp/rcspp/boost/RosterSP.h"


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
    nb_int_solutions_(0) {
  /* sort the nurses */
  std::shuffle(nursesToSolve_.begin(), nursesToSolve_.end(),
               Tools::getANewRandomGenerator());
}

/* Destructs the pricer object. */
RCPricer::~RCPricer() {
  for (auto &p : subProblems_)
    delete p.second;
}

/******************************************************
 * Perform pricing
 ******************************************************/
// at root node, after_fathom and backtracked are true
vector<MyVar *> RCPricer::pricing(double bound,
                                  bool before_fathom,
                                  bool after_fathom,
                                  bool backtracked) {
  // reset the current strategies at the beginning of a node
  if (after_fathom)  // first pricing for a new node
    for (const auto &p : subProblems_)
      p.second->updateParameters(pModel_->isFeasible());

  // get duals
  auto pDualCosts = std::make_shared<DualCosts>(pMaster_);
  pDualCosts->updateDuals();

  // Reset all rotations, columns, counters, etc.
  resetSolutions();

  // count and store the nurses whose subproblems produced rotations.
  vector<PLiveNurse> nursesSolved,
      nursesNoSolution,
      nursesSolvedOrder = nursesToSolve_;
  // flag if any shift has been  forbidden because of column disjoint feature
  bool disjointForbidden = false;
  // local thread pool (use all available thread)
  Tools::PThreadsPool pPool = Tools::ThreadsPool::newThreadsPool();

  int nSPBeingSolved = 0;
  vector2D<RCSolution> solutionsPerNurses(pMaster_->nNurses());
  while (!nursesToSolve_.empty()) {
    unique_lock<recursive_mutex> mLock(m_subproblem_);
    // RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
    PLiveNurse pNurse = nursesToSolve_.front();

    // try next nurse if forbidden
    if (isNurseForbidden(pNurse->num_)) {
      nursesToSolve_.erase(nursesToSolve_.begin());
      nursesNoSolution.push_back(pNurse);
      continue;
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
        break;  // break if there is no thread running anymore
    }

    // if the maximum number of subproblem solved is reached, break.
    if (nSPSolvedWithSuccess_ >= pModel_->getParameters().nSubProblemsToSolve_)
      break;

    // remove the nurse as being solved
    nursesToSolve_.erase(nursesToSolve_.begin());

    // define function to be run by the pool
    // The Job that is defined below will be run in parallel.
    // We create a job for each nurse and then run it in the pool
    /**
     * Start of the job definition (part run in parallel)
     */
    Tools::Job job([pNurse, pDualCosts, &nSPBeingSolved, &nursesSolved,
        &nursesNoSolution, &solutionsPerNurses,
        &disjointForbidden, &bound, this]() {
      // Add this part again, if we'd like to solve as many sub problem
      // as possible on the available threads
//      // before starting, verify if really need to be processed.
//      // if the maximum number of subproblem solved is reached, return.
//      if (nSPSolvedWithSuccess_
//          >= pModel_->getParameters().nSubProblemsToSolve_) {
//        lock_guard<recursive_mutex> lock(m_subproblem_);
//        nursesToSolve_.insert(nursesToSolve_.begin(),
//                              pNurse);  // add back the nurse
//        return;
//      }

      // lock the pricer
      unique_lock<recursive_mutex> lock(m_subproblem_);

      // UPDATE NURSE FORBIDDEN SHIFTS
      bool locDisjointForbidden = disjointForbidden;
      set<pair<int, int> > nurseForbiddenShifts(
          nursesForbiddenShifts_[pNurse->num_]);
      nurseForbiddenShifts.insert(forbiddenShifts_.begin(),
                                  forbiddenShifts_.end());
      pModel_->addForbiddenShifts(pNurse, &nurseForbiddenShifts);

      ++nbSPTried_;
      lock.unlock();  // unlock before solving

      // BUILD OR RE-USE THE SUBPROBLEM
      SubProblem *subProblem =
          retrieveSubproblem(pNurse, pModel_->getParameters().spParam_);

      // SOLVE THE PROBLEM
      subProblem->solve(pDualCosts,
                        nurseForbiddenShifts,
                        bound);
      /*if (!subProblem->solve(pDualCosts,
                        nurseForbiddenShifts,
                        bound)) {
        subProblem->setDefaultDomination(true);
        subProblem->updateParameters(true);
        bool isColumns = subProblem->solve(pDualCosts,
                          nurseForbiddenShifts,
                          bound);
        subProblem->setDefaultDomination(false);
        subProblem->updateParameters(true);
        std::cout << "with default domination, has generated columns ? " <<
        isColumns << std::endl;
      }*/


      // RETRIEVE THE GENERATED ROTATIONS
      std::vector<RCSolution> solutions = subProblem->getSolutions();

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
      if (pModel_->getParameters().rcspp_type_ == LABEL_SETTING &&
          subProblem->isLastRunOptimal()) {
        SubProblemParam par2(pMaster_->pModel()->getParameters());
        par2.strategyLevel_ =
            boostRCSPP::SubProblem::maxSubproblemStrategyLevel_;
        SubProblem *sub2 = nullptr;
        if (pModel_->getParameters().sp_type_ == ALL_ROTATION)
          sub2 =
              new boostRCSPP::RotationSP(pScenario_, nbDays_, pNurse, par2);
        else if (pModel_->getParameters().sp_type_ == ROSTER)
          sub2 =
              new boostRCSPP::RosterSP(pScenario_, nbDays_, pNurse, par2);

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

      // update the reduced costs
      updateRedCost(subProblem, solutions, locDisjointForbidden);

      // CHECK IF THE SUBPROBLEM GENERATED NEW ROTATIONS
      // If yes, store the nurses
      if (!solutions.empty()) {
        ++nSPSolvedWithSuccess_;
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
        // decrease nSPBeingSolved as no solutions were produced
        --nSPBeingSolved;
        nursesNoSolution.push_back(pNurse);
      }
    });  // END JOB
    /**
     * End of the job definition (part run in parallel)
     */

    // The job will be run in parallel.
    pPool->run(job);

    // if the maximum number of subproblems solved is reached, break.
    if (nSPSolvedWithSuccess_ >= pModel_->getParameters().nSubProblemsToSolve_)
      break;
  }

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

  // nurses order
  std::map<PLiveNurse, int> nursesOrder;
  int n = 0;
  for (const auto &pN : nursesSolvedOrder)
    nursesOrder[pN] = n++;

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
  auto *model = dynamic_cast<BcpModeler *>(pModel_);
  model->setLastNbSubProblemsSolved(nbSPTried_);
  model->setLastMinReducedCost(minReducedCost_);

  return allNewColumns_;
}

void RCPricer::updateRedCost(
    SubProblem *pSP,
    const std::vector<RCSolution> &solutions,
    bool disjointForbidden) {
  lock_guard<recursive_mutex> lock(m_subproblem_);

  // if not on on last level and disjointForbidden = false -> not optimal
  if (disjointForbidden || !pSP->isLastRunOptimal())
    optimal_ = false;

  // update reduced cost and min
  double redCost = solutions.empty() ? 0 : solutions.front().reducedCost();
  minReducedCosts_[pSP->pLiveNurse()->num_] = redCost;
  if (redCost < minReducedCost_)
    minReducedCost_ = redCost;
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
  SubProblem *&subProblem = subProblems_[pNurse];
  if (subProblem == nullptr)
    subProblem = buildSubproblem(pNurse, spParam);
  return subProblem;
}

SubProblem *RCPricer::buildSubproblem(const PLiveNurse &pNurse,
                                      const SubProblemParam &spParam) const {
  SubProblem *subProblem;
  switch (pModel_->getParameters().sp_type_) {
    case LONG_ROTATION:
      if (pModel_->getParameters().rcspp_type_ == BOOST_LABEL_SETTING)
        subProblem = new boostRCSPP::LongRotationSP(
            pScenario_, nbDays_, pNurse, spParam);
      else
        subProblem = new RotationSP(pScenario_, 0, nbDays_, pNurse,
                                    pMaster_->getSPResources(pNurse),
                                    spParam);
      break;
    case ALL_ROTATION: {
      if (pModel_->getParameters().rcspp_type_ == BOOST_LABEL_SETTING)
        subProblem = new boostRCSPP::RotationSP(
            pScenario_, nbDays_, pNurse, spParam);
      else
        subProblem = new RotationSP(pScenario_, 0, nbDays_, pNurse,
                                    pMaster_->getSPResources(pNurse),
                                    spParam);
      break;
    }
    case ROSTER: {
      if (pModel_->getParameters().rcspp_type_ == BOOST_LABEL_SETTING)
        subProblem = new boostRCSPP::RosterSP(
            pScenario_, nbDays_, pNurse, spParam);
      else if (pScenario_->isCyclic())
        subProblem = new CyclicRosterSP(pScenario_, 0, nbDays_, pNurse,
                                        pMaster_->getSPResources(pNurse),
                                        spParam);
      else
        subProblem = new RosterSP(pScenario_, 0, nbDays_, pNurse,
                                  pMaster_->getSPResources(pNurse),
                                  spParam);
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
#ifdef DBG
//      std::cout << sol.toStringINRC2() << std::endl;
if (dynamic_cast<BcpColumn*>(allNewColumns_.back()) == nullptr) {
  std::cerr << "Should be a BcpColumn." << std::endl;
}
#endif
    ++nbcolumnsAdded;
    if (nbcolumnsAdded >= pModel_->getParameters().spParam_.nbMaxColumnsToAdd_)
      break;
  }

#ifdef DBG
  //  std::cerr << "--------------------------------"
//               "--------------------------------" << std::endl;
#endif

  return nbcolumnsAdded;
}

// Sort the rotations that just were generated for a nurse.
// Default option is sort by increasing reduced cost but we
// could try something else (involving disjoint columns for ex.)
void RCPricer::sortGeneratedSolutions(
    std::vector<RCSolution> *solutions) const {
  RCSolution::sort(solutions);
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
