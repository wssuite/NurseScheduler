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

#include "RCSPPSubProblem.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <set>
#include <utility>

#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

RCSPPSubProblem::RCSPPSubProblem(
    PScenario scenario,
    int nbDays,
    PLiveNurse nurse,
    std::vector<PResource> pResources,
    SubProblemParam param) :
    SubProblem(std::move(scenario),
               nbDays,
               std::move(nurse),
               std::move(param)),
    pRCGraph_(std::make_shared<RCGraph>(nbDays, pScenario_->nShifts())),
    pEnumGraph_(nullptr),
    pResources_(std::move(pResources)),
    pRcsppSolver_(nullptr),
    timerEnumerationOfSubPath_(),
    timerComputeMinCostFromSink_() {
  Tools::initVector2D<PRCNode>(&pNodesPerDay_, nDays_,
                               pScenario_->nShifts(), nullptr);
  fixParameters();
}

void RCSPPSubProblem::build(const PRCGraph &pRCGraph) {
  pRCGraph->reset();
  // initialize the graph structure and base costs
  createNodes(pRCGraph);
  createArcs(pRCGraph);
}

void RCSPPSubProblem::build() {
  SubProblem::initStructuresForSolve();

  // set the status of all arcs to authorized
  Tools::initVector2D(&dayShiftStatus_, nDays_, pScenario_->nShifts(), true);
  nPathsMin_ = 0;

  //  Initialization of the vectors that will contain all information
  //  corresponding to bounds and to costs associated with the violations of
  //  soft bounds of each shift
  Tools::initVector(&consShiftsUbs_, pScenario_->nShiftTypes(), 0);
  Tools::initVector(&consShiftsLbs_, pScenario_->nShiftTypes(), 0);
  Tools::initVector(&consShiftsUbCosts_, pScenario_->nShiftTypes(), .0);
  Tools::initVector(&consShiftsLbCosts_, pScenario_->nShiftTypes(), .0);

  for (int st = 0; st < pScenario_->nShiftTypes(); st++) {
    shared_ptr<AbstractShift> absShift = std::make_shared<AnyOfTypeShift>(st);
    if (absShift->isWork()) {
      consShiftsLbs_.at(st) = pScenario_->minConsShiftsOfType(st);
      consShiftsUbs_.at(st) = pScenario_->maxConsShiftsOfType(st);
      consShiftsLbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_SHIFTS;
      consShiftsUbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_SHIFTS;
    } else if (absShift->isRest()) {
      consShiftsLbs_.at(st) = pLiveNurse_->minConsDaysOff();
      consShiftsUbs_.at(st) = pLiveNurse_->maxConsDaysOff();
      consShiftsLbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
      consShiftsUbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
    }
  }

  // build the graph
  this->build(pRCGraph_);

  // preprocess graph
  preprocessRCGraph();

  // create solver
  createRCSPPSolver();
}

void RCSPPSubProblem::enumerateSubPaths(const PRCGraph &pRCGraph) {
  // We enumerate all of the sub paths starting from any node of the network
  for (const PRCNode &pOrigin : pRCGraph->pNodes())
      enumerateConsShiftType(pRCGraph, pOrigin);

  // Finally, the costs of the arcs which existed before the enumeration
  // process are modified
  updateOfExistingArcsCost(pRCGraph);

  // set the resources that has not been enumerate
  vector<PResource> pNonEnumResources;
  for (const PResource &pR : pResources_) {
    auto pSR = std::dynamic_pointer_cast<SoftConsShiftResource>(pR);
    // if a soft const shift resource and the abstract is indeed a shift,
    // the resource has been enumerated (i.e., do not count AnyWorkShift
    // for example).
    if (pSR && dynamic_cast<AnyOfTypeShift*>(pSR->pAShift().get())) continue;
    pNonEnumResources.push_back(pR);
  }
  pRCGraph->pNonEnumResources(pNonEnumResources);
}

void RCSPPSubProblem::enumerateConsShiftType(
    const PRCGraph &pRCGraph, const PRCNode &pOrigin) {
  // Recovery of the last shift performed by the current nurse.
  // This is the 'initial Shift'.
  PShift pShiftIni = std::dynamic_pointer_cast<Shift>(pOrigin->pAShift);
  if (!pShiftIni) Tools::throwError("Enumerate works only with graph where "
                                    "the nodes represent a shift.");

  // Recovery of the number of times the last shift was completed by the
  // current nurse (It's either a worked shift or a rest shift). This value
  // is stored in the variable 'initialConsumption'.
  int initialConsumption(0);
  double initialBaseCost(0);
  if (pOrigin->day == -1) {
    if (pShiftIni->isRest())
      initialConsumption = pLiveNurse_->pStateIni_->consDaysOff_;
    else if (pShiftIni->isWork())
      initialConsumption = pLiveNurse_->pStateIni_->consShifts_;
    // Arc cost due to penalty of soft lower bound of the initial shift
    // will be added only if working on a different type of shift
    initialBaseCost = consShiftsLbCosts_.at(pShiftIni->id) *
        std::max(0, (consShiftsLbs_.at(pShiftIni->id) - initialConsumption));
  }

  // All following arcs added will have pOrigin node as their origin
  // Iterate through all the possible successor shifts of the initial shift
  for (auto indSuccessorShift : pShiftIni->successors) {
    // The successor shift is either of the same type as the initial shift
    // if starting from day -1 or of a different type.
    // 1- Compute the initial consumption of this shift type:
    // either 0 or initialConsumption if from source.
    int initialC(0);
    double initialBCost(0);
    if (indSuccessorShift == pShiftIni->id) {
      // Can work on the same shift only if starting from day -1.
      // Otherwise, the stretch must change of shift type
      if (pOrigin->day != -1) continue;
      // set the initial consumption
      initialC = initialConsumption;
    } else {
      // Arc cost due to penalty of soft lower bound of the initial
      // shift if on a different shift and if was starting from day -1
      initialBCost = initialBaseCost;
    }

    // 2- compute all possible stretch length of size <= UB or
    // UB - initialC if starting from day -1 and on the same initial shift
    PRCNode pTarget = pOrigin;
    vector<PShift> vecShift;
    vecShift.reserve(consShiftsUbs_.at(pShiftIni->id));
    for (int d{1}; d <= consShiftsUbs_.at(indSuccessorShift) - initialC &&
        d < pRCGraph_->nDays() - pOrigin->day; d++) {
      // if sink, stop
      if (pTarget->type == SINK_NODE) break;
      // Target node of the arc that will be added
      PShift pShift;
      bool succFound = false;
      for (const PRCArc &pArc : pTarget->outArcs) {
        pShift = std::dynamic_pointer_cast<Shift>(pArc->target->pAShift);
#ifdef DBG
        if (!pShift)
          Tools::throwError("Enumerate works only with graph where"
                            " the nodes represent a shift.");
#endif
        if (pShift->id == indSuccessorShift) {
          pTarget = pRCGraph->pNode(pArc->target->id);
          succFound = true;
          break;
        }
      }

      // if successor not found, break
      if (!succFound) break;

      // Creation of the stretch of the arc that will be added
      vecShift.push_back(pShift);  // add one more same shift to the stretch
      Stretch stretch(pOrigin->day+1, vecShift);

      // we do not add the first sub path as it already exists in the graph
      if (d == 1) continue;

      // Recovery of the base cost of the arc stretch (taking account the
      // previous shift, i.e. the initial shift)
      double bCost = baseCost(stretch, pOrigin->pAShift);

      // Arc cost due to penalty of soft lower bound of the corresponding
      // shift (added only if the target node is not a sink node).
      // Also, add initialC in case of starting from day -1, generally it's 0.
      if (pTarget->type != SINK_NODE)
        bCost += consShiftsLbCosts_.at(indSuccessorShift)
            * std::max(0, consShiftsLbs_.at(indSuccessorShift)
                - stretch.nDays() - initialC);

      // Arc cost due to penalty of soft lower bound of the initial shift :
      // generally 0
      bCost += initialBCost;

      // The new arc is added to the rcGraph with its corresponding costs
      PRCArc pArc = pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
    }
  }
}


void RCSPPSubProblem::updateOfExistingArcsCost(const PRCGraph &pRCGraph) {
  // Recovery of the last shift performed by the current nurse (It can be a
  // worked shift or a rest shift). This is the 'initial Shift'.
  State* pInitialState = this->pLiveNurse_->pStateIni_;
  PShift pShiftIni = pInitialState->pShift_;

  // Recovery of the number of times the last shift was completed by the
  // current nurse (It's either a worked shift or a rest shift). This value
  // is stored in the variable 'initialConsumption'.
  int initialConsumption(0);
  if (pShiftIni->isWork())
    initialConsumption = pInitialState->consShifts_;
  else if (pShiftIni->isRest())
    initialConsumption = pInitialState->consDaysOff_;

  // Iterate through all the arcs of the rcGraph which weren't added during
  // the enumeration process (only the 'non-enumerated' arcs).
  for (const auto &pArc : pRCGraph->pArcs()) {
    if (pArc->isEnumeratedArc()) continue;

    PShift pOShift = std::dynamic_pointer_cast<Shift>(pArc->origin->pAShift),
        pTShift = std::dynamic_pointer_cast<Shift>(pArc->target->pAShift);
    // Arcs coming from the source node are treated separately from the
    // others.
    if (pArc->origin->day != -1) {
      if (pOShift->id == pTShift->id) {
        // If the arc links two nodes with the same shift type, the cost due
        // to the violation of soft upper bound of the corresponding shifts
        // is added to the base cost of the arc.
        pArc->addBaseCost(consShiftsUbCosts_.at(pOShift->id));
      } else {
        // Otherwise, the cost due to the violation of the soft lower bound
        // of the target node's shifts is added to the base cost of the arc
        // (Only if the target node is not a sink node)
        if (pArc->target->type != SINK_NODE &&
            consShiftsLbs_.at(pTShift->id) >= 2)
          pArc->addBaseCost(consShiftsLbCosts_.at(pTShift->id) *
              (consShiftsLbs_.at(pTShift->id) - 1));
      }
      continue;
    }

    // Otherwise, starting from day -1
    if (pTShift->id == pShiftIni->id) {
      // If the arc links two nodes with the same shift type,
      // the costs due to the violation of soft bounds of
      // the initial shift (taking account the initial consumption) are
      // added to the base cost of the arc.
      if (initialConsumption + 1 < consShiftsLbs_.at(pShiftIni->id) )
        pArc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
            (consShiftsLbs_.at(pShiftIni->id) - initialConsumption - 1));
      if (initialConsumption + 1 > consShiftsUbs_.at(pShiftIni->id))
        pArc->addBaseCost(consShiftsUbCosts_.at(pShiftIni->id) *
            (initialConsumption+1 - consShiftsUbs_.at(pShiftIni->id)));
    } else {
      // Otherwise, the cost due to the violation of the soft lower bound
      // of the target node's shifts and the cost due to the violation of
      // the soft lower bound of the initial shift (taking account the
      // initial consumption) is added to the base cost of the arc
      if (consShiftsLbs_.at(pTShift->id) >= 2)
        pArc->addBaseCost(consShiftsLbCosts_.at(pTShift->id) *
            (consShiftsLbs_.at(pTShift->id)- 1));
      if (initialConsumption < consShiftsLbs_.at(pShiftIni->id))
        pArc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
            (consShiftsLbs_.at(pShiftIni->id)-initialConsumption));
    }
  }
}

PRCArc RCSPPSubProblem::addSingleArc(const PRCGraph &pRCGraph,
                                     const PRCNode &pOrigin,
                                     const PRCNode &pTarget,
                                     const PShift &pS,
                                     int day) {
  Stretch stretch(day, pS);
  double bCost = baseCost(stretch, pOrigin->pAShift);
  return pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
}

double RCSPPSubProblem::baseCost(
    const Stretch &stretch, PAbstractShift pPrevShift) {
  double cost = 0;
  int curDay = stretch.firstDay();
  // iterate through the shift to update the cost
  for (const auto& pS : stretch.pShifts()) {
    cost += preferencesCosts_[curDay][pS->id];
    if (pS->isWork()) {
      // add complete weekend cost if first worked shift of a rotation
      if (pPrevShift->isRest())
        cost += startWeekendCosts_[curDay];
    } else if (pS->isRest() && curDay > 0 && pPrevShift->isWork()) {
      cost += endWeekendCosts_[curDay - 1];
    }
    curDay++;
    pPrevShift = pS;
  }
  return cost;
}

double RCSPPSubProblem::dualCost(const PRCArc &pArc) {
  return pCosts_->getCost(
      pLiveNurse_->num_, pArc->stretch, pArc->origin->pAShift);
}


void RCSPPSubProblem::createResources(const PRCGraph &pRCGraph) {
  pRCGraph->clearResources();
  for (const PResource& pR : pResources_)
    pRCGraph->addResource(pR);
}

void RCSPPSubProblem::createRCSPPSolver() {
  pRcsppSolver_ = std::make_shared<RCSPPSolver>(pRCGraph_, param_);
}

void RCSPPSubProblem::preprocessRCGraph() {
  // resources
  createResources(pRCGraph_);

  if (param_.rcsppEnumSubpaths_) {
    timerEnumerationOfSubPath_.start();
    // Enumeration of sub paths in the rcGraph
    enumerateSubPaths(pRCGraph_);
    timerEnumerationOfSubPath_.stop();
//    std::cout << " - Total time spent in enumeration of sub-paths:  " <<
//              timerEnumerationOfSubPath_.dSinceStart() << std::endl;
  }

  if (param_.rcsppEnumSubpathsForMinCostToSinks_) {
    pEnumGraph_ = std::make_shared<RCGraph>(
        pRCGraph_->nDays(), pRCGraph_->nShifts());
    build(pEnumGraph_);
    createResources(pEnumGraph_);
    timerEnumerationOfSubPath_.start();
    // Enumeration of sub paths in the rcGraph
    enumerateSubPaths(pEnumGraph_);
    timerEnumerationOfSubPath_.stop();
//      std::cout << " - Total time spent in enumeration of sub-paths:  " <<
//                timerEnumerationOfSubPath_.dSinceStart() << std::endl;
  }

  // Initialization of each expander of each arc corresponding to each
  // resource
  pRCGraph_->initializeExpanders();

  // initialize the label
  initializeResources();
}

void RCSPPSubProblem::initializeResources() {
  pRCGraph_->initializeDominance();
}

void RCSPPSubProblem::updateParameters(bool masterFeasible) {
  // re-create a solver with the initial parameters
  createRCSPPSolver();
}

void RCSPPSubProblem::fixParameters() {
  if (param_.rcsppEnumSubpathsForMinCostToSinks_ &&
      (param_.rcsppEnumSubpaths_  || !param_.rcsppMinCostToSinks_)) {
    std::cout << "Warning: rcsppEnumSubpathsForMinCostToSinks has no point"
                 " if rcsppEnumSubpaths is on or rcsppMinCostToSinks is "
                 "off: it has been turned off" << std::endl;
    param_.rcsppEnumSubpathsForMinCostToSinks_ = false;
  }
  if (param_.rcsppNbToExpand_ > 0 && param_.rcsppDssr_) {
    std::cout << "Warning: rcsppNbToExpand_ and rcsppDssr_ do not work really "
                 "well together." << std::endl;
  }
}

bool RCSPPSubProblem::preprocess() {
  // update dual costs
  updateArcDualCosts();

  if (param_.rcsppEnumSubpathsForMinCostToSinks_) {
      if (pEnumGraph_ == nullptr)
        Tools::throwError("Enum graph has not been built.");
      for (const auto &pA : pEnumGraph_->pArcs())
        updateArcDualCost(pA);
      timerComputeMinCostFromSink_.start();
      // Computation of the costs of the shortest paths from each sink node
      // to all the other nodes in the rcGraph
      pRcsppSolver_->setMinimumCostToSinks(
          minCostPathToSinksAcyclic(pEnumGraph_));
      timerComputeMinCostFromSink_.stop();
//     std::cout << " - Total time spent in computing minimum paths costs from "
//                  "sinks: " << timerComputeMinCostFromSink_.dSinceStart() <<
//                std::endl;
  } else if (param_.rcsppMinCostToSinks_) {
    timerComputeMinCostFromSink_.start();
    // Computation of the costs of the shortest paths from each sink node
    // to all the other nodes in the rcGraph
    pRcsppSolver_->setMinimumCostToSinks(minCostPathToSinksAcyclic(pRCGraph_));
    timerComputeMinCostFromSink_.stop();
//    std::cout << " - Total time spent in computing minimum paths costs from "
//                 "sinks: " << timerComputeMinCostFromSink_.dSinceStart() <<
//              std::endl;
  }
  if (param_.verbose_ >= 4) {
    // print graph with updated costs
    pRCGraph_->printSummaryOfGraph();
  }

  // resetLabels solver in case it is not the first time it is called
  pRcsppSolver_->reset();
  pRcsppSolver_->initializeLabels();

  return true;
}

// TODO(JO): verify the computation with enumerated subpath once the
//  enumeration function is functional again

vector<double> RCSPPSubProblem::minCostPathToSinksAcyclic(
    const PRCGraph &pRCGraph) {
  // Initialization of all the costs of the shortest paths at infinity
  double inf = std::numeric_limits<double>::infinity();
  vector<double> shortestPathToSinks(pRCGraph->nNodes(), inf);
  // Only the cost of the shortest path to the sink node is set to 0
  for (const auto &pN : pRCGraph->pSinks())
    shortestPathToSinks.at(pN->id) = 0.0;

  // Iteration through all the nodes of the acyclic graph in the reverse
  // order of the topological order
  vector<PRCNode> sortedNodes = pRCGraph->sortNodes();
  auto itN = sortedNodes.rbegin();
  for (; itN != sortedNodes.rend(); itN++) {
    double sp = shortestPathToSinks[(*itN)->id];
    for (const auto &pArc : (*itN)->inArcs) {
      if (pArc->forbidden) continue;  // do not use forbidden arcs
      double cost = pArc->cost;
      int day = pArc->stretch.firstDay();
      auto it = pArc->stretch.pShifts().begin();
      vector<shared_ptr<Shift>> pShifts = pArc->stretch.pShifts();
      PAbstractShift pPrevShift = pArc->origin->pAShift;
      int predecessor = pArc->origin->id;
      for (; it != pArc->stretch.pShifts().end(); ++it, ++day) {
        // update shortest path if needed
        if (shortestPathToSinks[predecessor] > sp + cost) {
          shortestPathToSinks[predecessor] = sp + cost;
        }
        // when enumerating subpath only for computing minimum cost to sinks,
        // we need to compute what is the cost to each node of each stretch
        // of the enumerated arcs, otherwise, this is not necessary
        if (!param_.rcsppEnumSubpathsForMinCostToSinks_) break;
        if (day == pArc->stretch.lastDay()) break;

        // next, we will focus on the following node of the arc's stretch, so
        // update the cost to remove the base and dual cost associated with
        // current node
        // TODO(JO): what is done below is not general, I am actually
        //  re-computing the cost of an arc of the original graph. The
        //  simplest way forward would be to keep pointers to the original
        //  arcs in the enumerated arcs, but it raises several other issues.
        //  Rather create a method that computes a transition from day-shift1
        //  to day+1-shift2.
        cost -= preferencesCosts_[day][(*it)->id];
        if ((*it)->isWork()) {
          // add complete weekend cost if first worked shift of a rotation
          if (pPrevShift->isRest())
            cost -= startWeekendCosts_[day];
        } else if ((*it)->isRest() && day > 0 && pPrevShift->isWork()) {
          cost -= endWeekendCosts_[day - 1];
        }
        cost -= pCosts_->getCost(pLiveNurse_->num_,
                                 Stretch(day, *it),
                                 pPrevShift);

        pPrevShift = (*it);
        predecessor = pNodesPerDay_[day][(*it)->id]->id;
      }
    }
  }
  return shortestPathToSinks;
}

//// TODO(JO): verify the computation with enumerated subpath once the
////  enumeration function is functional again
// vector<double> RCSPPSubProblem::minCostPathFromGivenDayAcyclic(
//    const PRCGraph &pRCGraph, int sourceDay) {
//  // Initialization of all the costs of the shortest paths at infinity
//  double inf = std::numeric_limits<double>::infinity();
//  vector<double> minimumCost(pRCGraph->nNodes(), inf);
//  // get a topological order of the nodes (i.e. by increasing day)
//  vector<PRCNode> sortedNodes = pRCGraph->sortNodes();
//  // Only the cost of the shortest path to the sink node is set to 0
//  for (const auto &pN : sortedNodes) {
//    if (pN->day > sourceDay) break;
//    if (pN->day == sourceDay) minimumCost.at(pN->id) = 0.0;
//  }
//
//  // Iteration through all the nodes of the acyclic graph in topological
//  // order, starting from the first node of sourceDay
//  auto itN = sortedNodes.begin();
//  while ((*itN)->day < sourceDay) itN++;
//  for (; itN != sortedNodes.end(); itN++) {
//    double sp = minimumCost[(*itN)->id];
//    for (const auto &pArc : (*itN)->outArcs) {
//      if (pArc->forbidden) continue;  // do not use forbidden arcs
//      double cost = pArc->cost;
//      int day = pArc->stretch.lastDay();
//      auto it = pArc->stretch.pShifts().rbegin();
//      vector<shared_ptr<Shift>> pShifts = pArc->stretch.pShifts();
//      PAbstractShift pNextShift = pArc->target->pAShift;
//      int successor = pArc->target->id;
//      for (; it != pArc->stretch.pShifts().rend(); ++it, --day) {
//        // update shortest path if needed
//        if (minimumCost[successor] > sp + cost) {
//          minimumCost[successor] = sp + cost;
//        }
//        // when enumerating subpath only for computing minimum cost from a
//        // source day, we need to compute what is the cost to each node of
//      // each stretch of the enumerated arcs, otherwise, this is not necessary
//        if (!param_.rcsppEnumSubpathsForMinCostToSinks_) break;
//        if (day == pArc->stretch.firstDay()) break;
//
//        // next, we will focus on the previous node of the arc's stretch, so
//        // update the cost to remove the base and dual cost associated with
//        // current node
//        // TODO(JO): what is done below is not general, I am actually
//        //  re-computing the cost of an arc of the original graph. The
//        //  simplest way forward would be to keep pointers to the original
//        //  arcs in the enumerated arcs, but it raises several other issues.
//        //  Rather create a method that computes a transition from day-shift1
//        //  to day+1-shift2.
//        // TODO(AL): I don't understand what is happening here ...
//        //  I have the impression that pNextShift, should be previous
//        cost -= preferencesCosts_[day][(*it)->id];
//        if ((*it)->isWork()) {
//          // add complete weekend cost if first worked shift of a rotation
//          if (pNextShift->isRest())
//            cost -= startWeekendCosts_[day];
//        } else if ((*it)->isRest() && day > 0 && pNextShift->isWork()) {
//          cost -= endWeekendCosts_[day - 1];
//        }
//        cost -= pCosts_->getDualCost(pLiveNurse_->nurseNum_,
//                                     Stretch(day, *it),
//                                     pPrevShift);
//
//        pNextShift = (*it);
//        successor = pNodesPerDay_[day][(*it)->id]->id;
//      }
//    }
//  }
//  return minimumCost;
//}

void RCSPPSubProblem::updateArcDualCosts() {
  for (const auto& pA : pRCGraph_->pArcs())
    updateArcDualCost(pA);
}

void RCSPPSubProblem::updateArcDualCost(const PRCArc &pA) {
  pA->resetDualCost();
  pA->addDualCost(dualCost(pA));
}

bool RCSPPSubProblem::solveRCGraph() {
  // fix parameters
  fixParameters();

  // create the first label on the source node
  createInitialLabels();

  // Solution of the RCSPP obtained with the RCSPP Solver
  theSolutions_ = pRcsppSolver_->solve(maxReducedCostBound_);

  // Extract the best reduced cost
  if (!theSolutions_.empty())
    bestReducedCost_ = theSolutions_.front().reducedCost();

#ifdef DBG
  /*
  for (const RCSolution &sol : solutions) {
    theSolutions_.push_back(sol);
    nPaths_++;
    nFound_++;
    // we need to recompute the cost of each solution here
    RosterPattern pat(
        sol.shifts, pLiveNurse_->originalNurseId_, DBL_MAX, sol.cost);
    pat.computeCost(pScenario_, pLiveNurse_);
    double cost = pat.cost_;
    for (int day = 0; day < nDays_; ++day) {
      if (sol.shifts[day] > 0)
        cost -= pCosts_->workedDayShiftCost(day, sol.shifts[day]);
    }
    if (bestReducedCost_ > cost)
      bestReducedCost_ = cost;
  }*/
#endif

  return !theSolutions_.empty();
}

// Forbids a day-shift couple
void RCSPPSubProblem::forbidDayShift(int k, int s) {
  SubProblem::forbidDayShift(k, s);
  pRCGraph_->forbidDayShift(k, s);
}

// (re)Authorizes the day-shift couple
void RCSPPSubProblem::authorizeDayShift(int k, int s) {
  SubProblem::authorizeDayShift(k, s);
  pRCGraph_->authorizeDayShift(k, s);
}

// Reset all authorizations to true
void RCSPPSubProblem::resetAuthorizations() {
  SubProblem::resetAuthorizations();
  pRCGraph_->resetAuthorizationsArcs();
}

void RCSPPSubProblem::computeResourcesCosts(
    const State &initialState,
    MasterProblem *pMaster,
    RCSolution *rcSol) const {
  // reset costs
  rcSol->resetCosts();
  // create origin, destination and arc
  PRCNode pSource = std::make_shared<RCNode>(
      0, SOURCE_NODE, rcSol->firstDay() - 1, initialState.pShift_),
      pSink = std::make_shared<RCNode>(
      1, SINK_NODE, rcSol->lastDay(), rcSol->pShifts().back());
  PRCArc pArc = std::make_shared<RCArc>(
      0, pSource, pSink, *rcSol, rcSol->firstDay(), TO_SINK);
  // create expander for stretch
  double c = pArc->cost;
  for (const auto &pR : pResources_) {
    pR->initialize(*pSource->pAShift, *rcSol, pArc);
    rcSol->addCost(pArc->cost - c, pMaster->resourceCostType(pR, pLiveNurse_));
    c = pArc->cost;
  }
  // expand
  PRCLabel pL = std::make_shared<RCLabel>(pResources_, initialState);
  c = pL->cost();
  for (const auto &pE : pArc->expanders) {
    ResourceValues &v = pL->getResourceValues(pE->resourceId);
    pE->expand(pL, &v);
    rcSol->addCost(pL->cost() - c, pMaster->resourceCostType(
        pResources_[pE->resourceId], pLiveNurse_));
    c = pL->cost();
  }
}
