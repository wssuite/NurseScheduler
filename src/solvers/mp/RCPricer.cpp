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

#include <string>

#include "solvers/mp/sp/CyclicRosterSP.h"
#include "solvers/mp/sp/RosterSP.h"
#include "solvers/mp/sp/RotationSP.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/sp/rcspp/boost/LongRotationSP.h"
#include "solvers/mp/sp/rcspp/boost/RosterSP.h"
#include "solvers/mp/RosterMP.h"


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
RCPricer::RCPricer(MasterProblem *master,
                   const char *name,
                   const SolverParam &param) :
    MyPricer(name),
    pMaster_(master),
    pScenario_(master->pScenario()),
    nbDays_(master->nDays()),
    pModel_(master->getModel()),
    nursesToSolve_(master->sortedLiveNurses()),
    minReducedCost_(0),
    nb_int_solutions_(0),
    rand_(Tools::getANewRandomGenerator()) {
  /* sort the nurses */
  std::shuffle(nursesToSolve_.begin(), nursesToSolve_.end(), rand_);
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
  if (after_fathom) {  // first pricing for a new node
    // if backtracking -> restart at 0
    int diff = SubproblemParam::maxSubproblemStrategyLevel_
        - pModel_->getParameters().sp_default_strategy_;
    // was using (backtracked || diff<=0), but seems more efficient
    if (!pModel_->isFeasible() || diff <= 0)
      Tools::initVector(&currentSubproblemStrategy_, pMaster_->nNurses(),
                        pModel_->getParameters().sp_default_strategy_);
    else if (diff > 1)
      // otherwise, just try previous level if more than 2 levels
      Tools::initVector(&currentSubproblemStrategy_, pMaster_->nNurses(),
                        SubproblemParam::maxSubproblemStrategyLevel_ - 1);
  }

  // Reset all rotations, columns, counters, etc.
  resetSolutions();

  // count and store the nurses whose subproblems produced rotations.
  vector<PLiveNurse> nursesSolved,
      nursesIncreasedStrategyAndNoSolution,
      nursesNoSolution;
  // flag if any shift has been  forbidden because of column disjoint feature
  bool disjointForbidden = false;
  // local thread pool (use all available thread)
  Tools::ThreadsPool pool;

  while (!nursesToSolve_.empty()) {
    // RETRIEVE THE NURSE AND CHECK THAT HE/SHE IS NOT FORBIDDEN
    PLiveNurse pNurse = nursesToSolve_.front();
    nursesToSolve_.erase(nursesToSolve_.begin());

    // try next nurse if forbidden
    if (isNurseForbidden(pNurse->num_)) {
      nursesNoSolution.push_back(pNurse);
      continue;
    }

    // define function to be run by the pool
    // The Job that is defined below will be run in parallel.
    // We create a job for each nurse and then run it in the pool
    /**
     * Start of the job definition (part run in parallel)
     */
    Tools::Job job =
        [pNurse, &nursesSolved, &nursesIncreasedStrategyAndNoSolution,
            &nursesNoSolution, &disjointForbidden, &bound, this]() {
          // before starting, verify if really need to be processed.
          // if the maximum number of subproblem solved is reached, return.
          if (nbSPSolvedWithSuccess_
              >= pModel_->getParameters().nbSubProblemsToSolve_) {
            lock_guard<recursive_mutex> lock(m_subproblem_);
            nursesToSolve_.insert(nursesToSolve_.begin(),
                                  pNurse);  // add back the nurse
            return;
          }

          // RETRIEVE DUAL VALUES
          PDualCosts pDualCosts = pMaster_->buildDualCosts(pNurse);

          // lock the pricer
          unique_lock<recursive_mutex> lock(m_subproblem_);

          // UPDATE NURSE FORBIDDEN SHIFTS
          bool locDisjointForbidden = disjointForbidden;
          set<pair<int, int> > nurseForbiddenShifts(
              nursesForbiddenShifts_[pNurse->num_]);
          nurseForbiddenShifts.insert(forbiddenShifts_.begin(),
                                      forbiddenShifts_.end());
          pModel_->addForbiddenShifts(pNurse, &nurseForbiddenShifts);

          // SET SOLVING OPTIONS
          SubproblemParam sp_param(currentSubproblemStrategy_[pNurse->num_],
                                   pNurse,
                                   pMaster_->getModel()->getParameters());

          // BUILD OR RE-USE THE SUBPROBLEM
          SubProblem *subProblem = retrieveSubproblem(pNurse, sp_param);

          // SOLVE THE PROBLEM
          ++nbSPTried_;
          lock.unlock();  // unlock before solving
          subProblem->solve(pNurse,
                            pDualCosts,
                            sp_param,
                            nurseForbiddenShifts,
                            bound);

          // RETRIEVE THE GENERATED ROTATIONS
          std::vector<RCSolution> solutions = subProblem->getSolutions();

#ifdef DBG
//          if (pModel_->getParameters().sp_type_ == ALL_ROTATION &&
//              pModel_->getParameters().rcspp_type_ == LABEL_SETTING) {
//            SubProblem *sub2 =
//                new boostRCSPP::RotationSP(pScenario_,
//                                           nbDays_,
//                                           pNurse->pContract_,
//                                           pMaster_->pInitialStates());
//            SubproblemParam par2(SubproblemParam::maxSubproblemStrategyLevel_,
//                                 pNurse,
//                                 pMaster_->getModel()->getParameters());
//            sub2->build();
//            sub2->solve(pNurse,
//                        pDualCosts,
//                        par2,
//                        nurseForbiddenShifts,
//                        bound);
//            std::vector<RCSolution> sols = sub2->getSolutions();
//            delete sub2;
//
//            sortGeneratedSolutions(&solutions);
//            sortGeneratedSolutions(&sols);
//            if (!sols.empty() || !solutions.empty()) {
//              if (sols.empty() ^ solutions.empty()) {
//                if (!solutions.empty())
//                  std::cout << solutions.front().toString() << std::endl;
//                else
//                  std::cout << sols.front().toString() << std::endl;
//                Tools::throwError("One of the subproblem has found "
//                                  "a solution and the other not.");
//              }
//              double diff = sols.front().cost() - solutions.front().cost();
//              if (diff > pMaster_->getModel()->epsilon()
//                  || diff < -pMaster_->getModel()->epsilon()) {
//                std::cout << solutions.front().toString() << std::endl;
//                std::cout << sols.front().toString() << std::endl;
//                Tools::throwError(
//                    "The subproblems haven't find solutions of same cost.");
//              }
//            }
//          }
#endif

          // Lock the pricer
          lock.lock();

          // ADD THE ROTATIONS TO THE MASTER PROBLEM
          addColumnsToMaster(pNurse->num_, &solutions);

          // update the strategy and the reduced costs
          bool updated = updateCurrentStrategyAndRedCost(pNurse,
                                                         solutions,
                                                         locDisjointForbidden);

          // CHECK IF THE SUBPROBLEM GENERATED NEW ROTATIONS
          // If yes, store the nurses
          if (!solutions.empty()) {
            ++nbSPSolvedWithSuccess_;

            // UPDATE FORBIDDEN SHIFTS
            if (pModel_->getParameters().isColumnDisjoint_)
              // set disjointForbidden to true,
              // if addForbiddenShifts returns true once
              if (addForbiddenShifts(solutions, pDualCosts))
                disjointForbidden = true;

            // add the nurse to the nurses solved
            nursesSolved.push_back(pNurse);
          } else if (updated) {
            // Otherwise (no solution generated),
            // if strategy has been updated
            nursesIncreasedStrategyAndNoSolution.push_back(pNurse);
          } else {
            nursesNoSolution.push_back(pNurse);
          }
        };
    /**
     * End of the job definition (part run in parallel)
     */

    // begin parallel. The job will be run in parallel.
    // However, if no thread are available, the function run waits (i.e. blocks)
    // until a thread becomes available.
    pool.run(job);

    // if the maximum number of subproblem solved is reached, break.
    if (nbSPSolvedWithSuccess_
        >= pModel_->getParameters().nbSubProblemsToSolve_)
      break;

    // all nurses have been processed AND  we have increased strategy level
    // -> try to solve again for these nurses
    if (nursesToSolve_.empty()) {
      // wait for either:
      // a. have other nurses to try
      // b. all the threads to be finished
      // retest every time a thread has finished
      while (nursesIncreasedStrategyAndNoSolution.empty() && pool.wait(1))
        continue;
      // lock before copying and clearing
      lock_guard<recursive_mutex> lock(m_subproblem_);
      // then update the nurse to solve to continue the loop
      nursesToSolve_ = nursesIncreasedStrategyAndNoSolution;
      // remove all the nurses that have just been added back
      nursesIncreasedStrategyAndNoSolution.clear();
    }
  }

  // wait for the threads to be finished
  pool.wait();

  /* Add the nurses back into the nursesToSolve_ list
   * Reverse the vector before putting it back in
   * */

  // At least one subproblem unsolved
  if (!nursesToSolve_.empty()) optimal_ = false;

  // 1- Add the nurses in nursesIncreasedStrategyAndNoSolution:
  // first to be retried on next run as the subproblem will be more effective
  // at generating new columns
  reversePushBackNurses(&nursesIncreasedStrategyAndNoSolution);

  if (disjointForbidden) {
    // 2- Add the nurses in nursesNoSolution:
    // we hope that they will generate columns on next run as
    // some shifts were forbidden because of the column disjoint feature.
    reversePushBackNurses(&nursesNoSolution);

    // 3- Add the nurses in nursesSolved
    reversePushBackNurses(&nursesSolved);
  } else {
    // 2- Add the nurses in nursesSolved
    // (we hope that they will keep generating columns)
    reversePushBackNurses(&nursesSolved);

    // 3- Add the nurses in nursesSolved
    // (we hope that they will keep not generating columns)
    reversePushBackNurses(&nursesNoSolution);
  }

  // set statistics
  BcpModeler *model = static_cast<BcpModeler *>(pModel_);
  model->setLastNbSubProblemsSolved(nbSPTried_);
  model->setLastMinReducedCost(minReducedCost_);

  return allNewColumns_;
}

bool RCPricer::updateCurrentStrategyAndRedCost(
    PLiveNurse pNurse,
    const std::vector<RCSolution> &solutions,
    bool disjointForbidden) {
  lock_guard<recursive_mutex> lock(m_subproblem_);

  // if not on on last level and disjointForbidden = false -> not optimal
  if (disjointForbidden || currentSubproblemStrategy_[pNurse->num_]
      < SubproblemParam::maxSubproblemStrategyLevel_)
    optimal_ = false;

  // check if enough columns have been generated and no shifts forbidden
  // by disjoint column strategy,
  // otherwise, increase strategy
  bool increased = false;
  if (!disjointForbidden &&
      solutions.size() < pModel_->getParameters().nbMaxColumnsToAdd_ *
          pModel_->getParameters().sp_min_columns_ratio_for_increase_ &&
      currentSubproblemStrategy_[pNurse->num_]
          < SubproblemParam::maxSubproblemStrategyLevel_) {
    currentSubproblemStrategy_[pNurse->num_]++;
#ifdef DBG
    //    std::cout << "nurse " << pNurse->num_ << ": lvl "
    //              << currentSubproblemStrategy_[pNurse->num_] << std::endl;
#endif
    increased = true;
  }

  // update reduced cost and min
  double redCost = solutions.empty() ? 0 : solutions.front().cost();
  minReducedCosts_[pNurse->num_] = redCost;
  if (redCost < minReducedCost_)
    minReducedCost_ = redCost;

  return increased;
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

bool RCPricer::addForbiddenShifts(const std::vector<RCSolution> &solutions,
                                  const PDualCosts &pDuals) {
  // search best rotation
  RCSolution bestSolution;
  double bestRedcost = DBL_MAX;
  for (const RCSolution &sol : solutions)
    if (sol.cost() < bestRedcost) {
      bestRedcost = sol.cost();
      bestSolution = sol;
    }

  // forbid working shifts of the best column if small enough
  if (bestRedcost > pModel_->getParameters().minReducedCostDisjoint_)
    return false;

  lock_guard<recursive_mutex> lock(m_subproblem_);
  bool shiftAdded = false;
  int k = bestSolution.firstDay();
  for (const PShift &pS : bestSolution.pShifts()) {
    // work shift with a dual cost high enough
    // (have enough diversity in the solutions generated)
    if (pDuals->workedDayShiftCost(k, pS->id) > pModel_->epsilon()) {
      forbiddenShifts_.insert(pair<int, int>(k, pS->id));
      shiftAdded = true;
    }
    k++;
  }
  return shiftAdded;
}

// Returns a pointer to the right subproblem
SubProblem *RCPricer::retrieveSubproblem(const PLiveNurse &pNurse,
                                         const SubproblemParam &spParam) {
  // lock the pricer
  unique_lock<recursive_mutex> lock(m_subproblem_);
  // Each nurse has a subproblem. If null, create a new one.
  SubProblem *&subProblem = subProblems_[pNurse];
  if (subProblem == nullptr) {
    switch (pModel_->getParameters().sp_type_) {
      case LONG_ROTATION:
        subProblem = new boostRCSPP::LongRotationSP(pScenario_,
                                                    nbDays_,
                                                    pNurse->pContract_,
                                                    pMaster_->pInitialStates());
        break;
      case ALL_ROTATION: {
        if (pModel_->getParameters().rcspp_type_ == BOOST_LABEL_SETTING)
          subProblem = new boostRCSPP::RotationSP(pScenario_,
                                                  nbDays_,
                                                  pNurse->pContract_,
                                                  pMaster_->pInitialStates());
        else
          subProblem = new RotationSP(pScenario_, nbDays_, pNurse,
                                      pMaster_->createPResources(pNurse),
                                      spParam);
        break;
      }
      case ROSTER: {
        if (pModel_->getParameters().rcspp_type_ == BOOST_LABEL_SETTING)
          subProblem = new boostRCSPP::RosterSP(pScenario_,
                                                nbDays_,
                                                pNurse->pContract_,
                                                pMaster_->pInitialStates());
        else if (pScenario_->isCyclic())
          subProblem = new CyclicRosterSP(pScenario_, nbDays_, pNurse,
                                          pMaster_->createPResources(pNurse),
                                          spParam);
        else
          subProblem = new RosterSP(pScenario_, nbDays_, pNurse,
                                    pMaster_->createPResources(pNurse),
                                    spParam);
        break;
      }
      default:
        Tools::throwError(
            "There is no subproblem defined associated to this type");
    }
    // then build the rcspp
    subProblem->build();
  }
  return subProblem;
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
#ifdef DBG
//  std::cout << sol.toString(pScenario_->shiftIDToShiftTypeID_) << std::endl;
#endif
    allNewColumns_.emplace_back(pMaster_->addColumn(nurseNum, sol));
    ++nbcolumnsAdded;
    if (nbcolumnsAdded >= pModel_->getParameters().nbMaxColumnsToAdd_)
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
  std::stable_sort(solutions->begin(), solutions->end(),
                   [](const RCSolution &sol1, const RCSolution &sol2) {
                     return sol1.cost() < sol2.cost();
                   });
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
