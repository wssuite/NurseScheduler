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

#include "solvers/mp/sp/RosterSP.h"

#include <algorithm>
#include <memory>
#include <set>
#include <utility>
#include <limits>

#include "solvers/mp/RosterMP.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

RosterSP::RosterSP(PScenario scenario, int nbDays, PLiveNurse nurse,
    std::vector<PResource> pResources,
    const SubproblemParam &param) :
    SubProblem(std::move(scenario),
               nbDays,
               std::move(nurse),
               param),
    rcGraph_(nbDays, pScenario_->nShifts()),
    pResources_(std::move(pResources)),
    pRcsppSolver_(nullptr),
    timerEnumerationOfSubPath_(),
    timerComputeMinCostFromSink_() {
  Tools::initVector2D<PRCNode>(&pNodesPerDayShift_, nDays_,
                               pScenario_->nShifts(), nullptr);
}

void RosterSP::build() {
  SubProblem::initStructuresForSolve();
  this->build(&rcGraph_);
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
      consShiftsLbs_.at(st) = pScenario_->minConsShiftsOf(st);
      consShiftsUbs_.at(st) = pScenario_->maxConsShiftsOf(st);
      consShiftsLbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_SHIFTS;
      consShiftsUbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_SHIFTS;
    } else if (absShift->isRest()) {
      consShiftsLbs_.at(st) = pLiveNurse_->minConsDaysOff();
      consShiftsUbs_.at(st) = pLiveNurse_->maxConsDaysOff();
      consShiftsLbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
      consShiftsUbCosts_.at(st) = pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
    }
  }
  // resources
  createResources(&rcGraph_);

  if (param_.rcsppEnumSubpaths_) {
    timerEnumerationOfSubPath_.start();
    // Enumeration of sub paths in the rcGraph
    enumerateSubPaths(&rcGraph_);
    timerEnumerationOfSubPath_.stop();
    std::cout << " - Total time spent in enumeration of sub-paths:  " <<
              timerEnumerationOfSubPath_.dSinceStart() << std::endl;
  }

  // Initialization of each expander of each arc corresponding to each
  // resource
  rcGraph_.initializeExpanders();

  // create solver
  createRCSPPSolver();

  // preprocess graph
  preprocessRCGraph();
}

void RosterSP::build(RCGraph *pRCGraph) {
  pRCGraph->reset();
  // initialize the graph structure and base costs
  createNodes(pRCGraph);
  createArcs(pRCGraph);
}


void RosterSP::createNodes(RCGraph *pRCGraph) {
  pRCGraph->addSingleNode(SOURCE_NODE, -1, pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d < nDays_ - 1; ++d)
    for (const auto& pShift : pScenario_->pShifts())
      pNodesPerDayShift_[d][pShift->id] =
          pRCGraph->addSingleNode(PRINCIPAL_NETWORK, d, pShift);

  // every shift on the last day is a sink
  for (const auto& pShift : pScenario_->pShifts())
    pNodesPerDayShift_[nDays_ - 1][pShift->id] =
        pRCGraph->addSingleNode(SINK_NODE, nDays_ - 1, pShift);
}


void RosterSP::createArcs(RCGraph* pRCGraph) {
  // arcs from source to first day
  PShift pShiftIni = pScenario_->pShift(pLiveNurse_->pStateIni_->shift_);
  for (auto shiftId : pShiftIni->successors) {
    PRCNode pN = pNodesPerDayShift_[0][shiftId];
    addSingleArc(
        pRCGraph, pRCGraph->pSource(), pN, pScenario_->pShift(shiftId), 0);
  }

  // arcs from the previous day to current day;
  // those from nDays-2 to nDays-1 are the arcs to the sinks
  for (int d = 1; d < nDays_; ++d) {
    for (const PShift &pS : pScenario_->pShifts()) {
      PRCNode pOrigin = pNodesPerDayShift_[d-1][pS->id];
      for (int succId : pS->successors) {
        PRCNode pTarget = pNodesPerDayShift_[d][succId];
        addSingleArc(pRCGraph, pOrigin, pTarget,
                     pScenario_->pShift(succId), d);
      }
    }
  }
}

void RosterSP::enumerateSubPaths(RCGraph *pRcGraph) {
  // First, we enumerate sub paths from the source node taking account the
  // history of the current nurse
  enumerateConsShiftTypeFromSource(pRcGraph);

  // Then, we enumerate sub paths from a PRINCIPAL_NETWORK node in the rcGraph
  // corresponding to a given day and to a given shift
  for (int d = 0; d < nDays_ - 1; ++d)
    for (const auto& pS : pScenario_->pShifts())
      enumerateConsShiftType(pRcGraph, pS, d);

  // Finally, the costs of the arcs which existed before the enumeration
  // process are modified
  updateOfExistingArcsCost(pRcGraph);
}


void RosterSP::enumerateConsShiftTypeFromSource(RCGraph *pRCGraph) {
/*  // Recovery of the last shift performed by the current nurse (It can be a
  // worked shift or a rest shift). This is the 'initial Shift'.
  State* initialState = this->pLiveNurse_->pStateIni_;
  PShift pShiftIni = pScenario_->pShifts_[initialState->shift_];

  // Recovery of the number of times the last shift was completed by the
  // current nurse (It's either a worked shift or a rest shift). This value
  // is stored in the variable 'initialConsumption'.
  int initialConsumption(0);
  if (pShiftIni->isRest())
    initialConsumption = initialState->consDaysOff_;
  else if (pShiftIni->isWork())
    initialConsumption = initialState->consShifts_;

  // All following arcs added will have the source node as their origin
  PRCNode pOrigin = pRCGraph->pSource();
  // Iterate through all the possible successor shifts of the initial shift
  for (auto indSuccessorShift : pShiftIni->successors) {
    // The successor shift is either of the same type as the initial shift or
    // of a different type. The two cases are treated separately in order to
    // take account the initial consumption of the last shift performed.
    if (indSuccessorShift == pShiftIni->id) {
      for (int d{1}; d <= consShiftsUbs_.at(pShiftIni->id) -
          initialConsumption - 1 && d < pScenario_->nbDays() ; d++) {
        // Target node of the arc that will be added
        PRCNode pTarget = pNodesPerDayShift_[d][indSuccessorShift];

        // Creation of the stretch of the arc that will be added
        vector<PShift> vecShift(d+1, nullptr);
        for (size_t i = 0; i < d+1; i++)
          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pAShift);
        Stretch stretch(vecShift, 0);

        // Recovery of the base cost of the arc stretch (taking account the
        // previous shift)
        double bCost = baseCost(stretch, pOrigin->pAShift);

        // Arc cost due to penalty of soft lower bound of the corresponding
        // shift (this initial shift)
        bCost += consShiftsLbCosts_.at(pShiftIni->id) *
            std::max(0, consShiftsLbs_.at(pShiftIni->id) -
                initialConsumption - 1 - d);

        // The new arc is added to the rcGraph with its corresponding costs
        PRCArc pArc =
            pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
      }
    } else {
      for (int d{1}; d <= consShiftsUbs_.at(indSuccessorShift)-1 && d <
          pScenario_->nbDays(); d++) {
        // Target node of the arc that will be added
        PRCNode pTarget = pNodesPerDayShift_[d][indSuccessorShift];

        // Creation of the stretch of the arc that will be added
        vector<PShift> vecShift(d+1, nullptr);
        for (size_t i = 0; i < d+1; i++)
          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pAShift);
        Stretch stretch(vecShift, 0);

        // Recovery of the base cost of the arc stretch (taking account the
        // previous shift)
       double bCost = baseCost(stretch, pOrigin->pAShift);

        // Arc cost due to penalty of soft lower bound of the corresponding
        // shift
        bCost += std::max(0, consShiftsLbs_.at(indSuccessorShift) - d - 1)
            * consShiftsLbCosts_.at(indSuccessorShift);

        // Arc cost due to penalty of soft lower bound of the initial
        // shift
        bCost += consShiftsLbCosts_.at(pShiftIni->id) *
           std::max(0, (consShiftsLbs_.at(pShiftIni->id)-initialConsumption));

        // The new arc is added to the rcGraph with its corresponding costs
        PRCArc pArc =
            pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
      }
    }
  }*/
}

void RosterSP::enumerateConsShiftType(RCGraph *pRCGraph,
                                      const PShift& pS,
                                      int day) {
//  // All arcs that will be added will have as origin the node corresponding to
//  // the day and to the shift given in parameter
//  PRCNode pOrigin = pNodesPerDayShift_[day][pS->id];
//
//  // Iterate through all the possible successor shifts of the origin node's
//  // shift
//  for (auto indSuccessorShift : pS->successors) {
//    // New arcs will be added only toward nodes whose shift is different from
//    // the origin node's one
//    if (pS->id != indSuccessorShift) {
//      for (int d = day+2; d <= day + consShiftsUbs_.at(indSuccessorShift)
//          && d < pScenario_->nbDays(); d++) {
//        // Target node of the arc that will be added
//        PRCNode pTarget = pNodesPerDayShift_[d][indSuccessorShift];
//
//        // Creation of the stretch of the arc that will be added
//        vector<PShift> vecShift(d-day, nullptr);
//        for (size_t i = 0; i < d-day; i++)
//          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pShift);
//        Stretch stretch(vecShift, day+1, d-day);
//
//        // Recovery of the base cost of the arc stretch (taking  into account
//        // the previous shift)
//        double bCost = baseCost(stretch, pOrigin->pShift);
//
//        // The arc cost due to penalty of soft lower bound of the
//        // corresponding shift is added only if the target node is not
//        // a sink node
//        if (pTarget->type != SINK_NODE) {
//          int missingDays = consShiftsLbs_.at(indSuccessorShift) - (d - day);
//          if (missingDays > 0)
//            bCost += consShiftsLbCosts_.at(indSuccessorShift) * missingDays;
//        }
//
//        // The new arc is added to the rcGraph with its corresponding costs
//        PRCArc pArc =
//            pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
//      }
//    }
//  }
}


void RosterSP::updateOfExistingArcsCost(RCGraph *pRCGraph) {
//  // Recovery of the last shift performed by the current nurse (It can be a
//  // worked shift or a rest shift). This is the 'initial Shift'.
//  State* initialState = this->pLiveNurse_->pStateIni_;
//  PShift pShiftIni = pScenario_->pShifts_[initialState->shift_];
//
//  // Recovery of the number of times the last shift was completed by the
//  // current nurse (It's either a worked shift or a rest shift). This value
//  // is stored in the variable 'initialConsumption'.
//  int initialConsumption(0);
//  if (pShiftIni->isWork())
//    initialConsumption = initialState->consShifts_;
//  else if (pShiftIni->isRest())
//    initialConsumption = initialState->consDaysOff_;
//
//  // Iterate through all the arcs of the rcGraph which weren't added during
//  // the enumeration process (only the 'non-enumerated' arcs).
//  for (const auto &arc : pRCGraph->pArcs()) {
//    if (!arc->isEnumeratedArc()) {
//      // Arcs coming from the source node are treated separately from the
//      // others.
//      if (arc->origin->id != 0) {
//        if (arc->origin->pShift->id == arc->target->pShift->id) {
//          // If the arc links two nodes with the same shift type, the cost due
//          // to the violation of soft upper bound of the corresponding shifts
//          // is added to the base cost of the arc.
//          arc->addBaseCost(consShiftsUbCosts_.at(arc->origin->pShift->id));
//        } else {
//          // Otherwise, the cost due to the violation of the soft lower bound
//          // of the target node's shifts is added to the base cost of the arc
//          // (Only if the target node is not a sink node)
//          if (arc->target->type != SINK_NODE &&
//              consShiftsLbs_.at(arc->target->pShift->id) >= 2)
//            arc->addBaseCost(consShiftsLbCosts_.at(arc->target->pShift->id) *
//                (consShiftsLbs_.at(arc->target->pShift->id)-1));
//        }
//      } else {
//        if (arc->target->pShift->id == pShiftIni->id) {
//          // If the arc links two nodes with the same shift type,
//          // the costs due to the violation of soft bounds of
//          // the initial shift (taking account the initial consumption) are
//          // added to the base cost of the arc.
//          if (initialConsumption + 1 < consShiftsLbs_.at(pShiftIni->id) )
//            arc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
//                (consShiftsLbs_.at(pShiftIni->id) - initialConsumption - 1));
//          if (initialConsumption + 1 > consShiftsUbs_.at(pShiftIni->id))
//            arc->addBaseCost(consShiftsUbCosts_.at(pShiftIni->id) *
//                (initialConsumption+1 - consShiftsUbs_.at(pShiftIni->id)));
//        } else {
//          // Otherwise, the cost due to the violation of the soft lower bound
//          // of the target node's shifts and the cost due to the violation of
//          // the soft lower bound of the initial shift (taking account the
//          // initial consumption) is added to the base cost of the arc
//          if (consShiftsLbs_.at(arc->target->pShift->id) >= 2)
//            arc->addBaseCost(consShiftsLbCosts_.at(arc->target->pShift->id) *
//                (consShiftsLbs_.at(arc->target->pShift->id)- 1));
//          if (initialConsumption < consShiftsLbs_.at(pShiftIni->id))
//            arc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
//                (consShiftsLbs_.at(pShiftIni->id)-initialConsumption));
//        }
//      }
//    }
//  }
}

PRCArc RosterSP::addSingleArc(RCGraph* pRCGraph,
                              const PRCNode &pOrigin,
                              const PRCNode &pTarget,
                              const PShift &pS,
                              int day) {
  Stretch stretch(pS, day);
  double bCost = baseCost(stretch, pOrigin->pAShift);
  return pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
}

double RosterSP::dualCost(const Stretch &stretch) {
  double dualCost = 0;
  int curDay = stretch.firstDay();
  // if start, add constant
  if (curDay == 0)
    dualCost -= pCosts_->constant();
  // iterate through the shift to update the cost
  for (const auto& pS : stretch.pShifts()) {
    if (pS->isWork())
      dualCost -= pCosts_->workedDayShiftCost(curDay, pS->id);
    curDay++;
  }
  return dualCost;
}

double RosterSP::baseCost(
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

void RosterSP::createResources(RCGraph *pRCGraph) {
  pRCGraph->clearResources();
  for (const PResource& pR : pResources_)
    pRCGraph->addResource(pR);
}

void RosterSP::createRCSPPSolver() {
  pRcsppSolver_ = std::make_shared<MyRCSPPSolver>(&rcGraph_, param_);
}

void RosterSP::preprocessRCGraph() {
  initializeResources();
}

void RosterSP::initializeResources() {
  rcGraph_.initializeDominance();
}

bool RosterSP::preprocess() {
  // update dual costs
  updateArcDualCosts();

  // TODO(JO): create a method that checks the validity and compatibility of
  //  options
  if (param_.rcsppEnumSubpathsForMinCostToSinks_) {
    if (!param_.rcsppEnumSubpaths_ && param_.rcsppMinCostToSinks_) {
      RCGraph enumGraph(rcGraph_.nDays(), rcGraph_.nShifts());
      build(&enumGraph);
      createResources(&enumGraph);
      for (const auto &pA : enumGraph.pArcs())
        updateArcDualCost(pA);
      timerEnumerationOfSubPath_.start();
      // Enumeration of sub paths in the rcGraph
      enumerateSubPaths(&enumGraph);
      timerEnumerationOfSubPath_.stop();
      std::cout << " - Total time spent in enumeration of sub-paths:  " <<
                timerEnumerationOfSubPath_.dSinceStart() << std::endl;

      timerComputeMinCostFromSink_.start();
      // Computation of the costs of the shortest paths from each sink node
      // to all the other nodes in the rcGraph
      pRcsppSolver_->setMinimumCostToSinks(
          shortestPathToSinksAcyclic(&enumGraph));
      timerComputeMinCostFromSink_.stop();
      std::cout << " - Total time spent in computing minimum paths costs from "
                   "sinks: " << timerComputeMinCostFromSink_.dSinceStart() <<
                std::endl;
    } else {
      std::cout << "Warning: rcsppEnumSubpathsForMinCostToSinks has no point"
                   " if rcsppEnumSubpaths is on or rcsppMinCostToSinks is "
                   "off: it has been turned off" << std::endl;
      param_.rcsppEnumSubpathsForMinCostToSinks_ = false;
    }
  } else if (param_.rcsppMinCostToSinks_) {
    timerComputeMinCostFromSink_.start();
    // Computation of the costs of the shortest paths from each sink node
    // to all the other nodes in the rcGraph
    pRcsppSolver_->setMinimumCostToSinks(shortestPathToSinksAcyclic(&rcGraph_));
    timerComputeMinCostFromSink_.stop();
    std::cout << " - Total time spent in computing minimum paths costs from "
                 "sinks: " << timerComputeMinCostFromSink_.dSinceStart() <<
              std::endl;
  }
  if (param_.verbose_ >= 4) {
    // print graph with updated costs
    rcGraph_.printSummaryOfGraph();
  }

  // resetLabels solver in case it is not the first time it is called
  pRcsppSolver_->reset(param_);
  pRcsppSolver_->initializeLabels();

  return true;
}

// TODO(JO): verify the computation with enumerated subpath once the
//  enumeration function is functional again
vector<double> RosterSP::shortestPathToSinksAcyclic(
    const RCGraph *pRCGraph) {
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
        cost -= preferencesCosts_[day][(*it)->id];
        if ((*it)->isWork()) {
          // add complete weekend cost if first worked shift of a rotation
          if (pPrevShift->isRest())
            cost -= startWeekendCosts_[day];
        } else if ((*it)->isRest() && day > 0 && pPrevShift->isWork()) {
          cost -= endWeekendCosts_[day - 1];
        }
        if (day == 0) cost -= pCosts_->constant();
        if ((*it)->isWork())
            cost += pCosts_->workedDayShiftCost(day, (*it)->id);

        pPrevShift = (*it);
        predecessor = pNodesPerDayShift_[day][(*it)->id]->id;
      }
    }
  }
  return shortestPathToSinks;
}

void RosterSP::updateArcDualCosts() {
  for (const auto& pA : rcGraph_.pArcs())
    updateArcDualCost(pA);
}

void RosterSP::updateArcDualCost(const PRCArc &pA) {
  pA->resetDualCost();
  pA->addDualCost(dualCost(pA->stretch));
}

bool RosterSP::solveRCGraph() {
  // create the first label on the source node
  createInitialLabel();

  // Solution of the RCSPP obtained with the RCSPP Solver
  theSolutions_ = pRcsppSolver_->solve(maxReducedCostBound_,
                                       param_.rcsppMinNegativeLabels_);

  // Extract the best reduced cost
  if (!theSolutions_.empty())
    bestReducedCost_ = theSolutions_[0].cost;

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

void RosterSP::createInitialLabel() {
  auto pL = std::make_shared<RCLabel>(rcGraph_.pResources(),
                                      *pLiveNurse_->pStateIni_);
  pL->setNode(rcGraph_.pSource());
  pRcsppSolver_->setSourceLabel(pL);
}

// Forbids a day-shift couple
void RosterSP::forbidDayShift(int k, int s) {
  SubProblem::forbidDayShift(k, s);
  rcGraph_.forbidDayShift(k, s);
}

// (re)Authorizes the day-shift couple
void RosterSP::authorizeDayShift(int k, int s) {
  SubProblem::authorizeDayShift(k, s);
  rcGraph_.authorizeDayShift(k, s);
}

// Reset all authorizations to true
void RosterSP::resetAuthorizations() {
  SubProblem::resetAuthorizations();
  rcGraph_.resetAuthorizationsArcs();
}
