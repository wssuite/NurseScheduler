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

#include <limits>
#include <memory>
#include <set>
#include <utility>
#include <map>

#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"

RCSPPSubProblem::RCSPPSubProblem(
    PScenario scenario,
    int firstDayId,
    int nbDays,
    PLiveNurse nurse,
    std::vector<PResource> pResources,
    SubProblemParam param) :
    SubProblem(std::move(scenario),
               firstDayId,
               nbDays,
               std::move(nurse),
               param),
    pRCGraph_(std::make_shared<RCGraph>(
        firstDayId, nbDays, pScenario_->nShifts(),
        pScenario_->shiftsFactory().pAnyTypeShifts())),
    pEnumGraph_(nullptr),
    pResources_(std::move(pResources)),
    pRcsppSolver_(nullptr),
    timerEnumerationOfSubPath_("Sub path enumeration"),
    timerComputeMinCostFromSink_("Min cost from sink") {
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

  // build the graph
  build(pRCGraph_);

  // presolve the graph that will be used to solve the RCSPP
  preprocessRCGraph(pRCGraph_, param_.rcsppEnumSubpaths_);

  // if we do not enumerate the consecutive shift resources in the RCSPP, but
  // wish to enumerate them to compute better bounds on the minimum cost from
  // each node to sink, we build an enumerated graph
  if (!param_.rcsppEnumSubpaths_ &&
      param_.rcsppEnumSubpathsForMinCostToSinks_) {
    pEnumGraph_ = std::make_shared<RCGraph>(pRCGraph_->firstDayId(),
                                            pRCGraph_->nDays(),
                                            pRCGraph_->nShifts(),
                                            pRCGraph_->pAShifts());
    build(pEnumGraph_);
    preprocessRCGraph(pEnumGraph_, true);
  }

  // create solver
  createRCSPPSolver();
}

// void RCSPPSubProblem::enumerateSubPaths(const PRCGraph &pRCGraph) {
//  // We enumerate all of the sub paths starting from any node of the network
//  for (const PRCNode &pOrigin : pRCGraph->pNodes())
//      enumerateConsShiftType(pRCGraph, pOrigin);
//
//  // Finally, the costs of the arcs which existed before the enumeration
//  // process are modified
//  updateOfExistingArcsCost(pRCGraph);
//
//  // set the resources that has not been enumerate
//  vector<PResource> pNonEnumResources;
//  for (const PResource &pR : pResources_) {
//    auto pSR = std::dynamic_pointer_cast<SoftConsShiftResource>(pR);
//    // if a soft const shift resource and the abstract is indeed a shift,
//    // the resource has been enumerated (i.e., do not count AnyWorkShift
//    // for example).
//    if (pSR && dynamic_cast<AnyOfTypeShift*>(pSR->pAShift().get())) continue;
//    pNonEnumResources.push_back(pR);
//  }
//  pRCGraph->pNonEnumResources(pNonEnumResources);
// }
//
// void RCSPPSubProblem::enumerateConsShiftType(
//    const PRCGraph &pRCGraph, const PRCNode &pOrigin) {
//  // Recovery of the last shift performed by the current nurse.
//  // This is the 'initial Shift'.
//  PShift pShiftIni = std::dynamic_pointer_cast<Shift>(pOrigin->pAShift);
//  if (!pShiftIni) Tools::throwError("Enumerate works only with graph where "
//                                    "the nodes represent a shift.");
//
//  // Recovery of the number of times the last shift was completed by the
//  // current nurse (It's either a worked shift or a rest shift). This value
//  // is stored in the variable 'initialConsumption'.
//  int initialConsumption(0);
//  double initialBaseCost(0);
//  if (pOrigin->day == -1) {
//    if (pShiftIni->isRest())
//      initialConsumption = pLiveNurse_->pStateIni_->consDaysOff_;
//    else if (pShiftIni->isWork())
//      initialConsumption = pLiveNurse_->pStateIni_->consShifts_;
//    // Arc cost due to penalty of soft lower bound of the initial shift
//    // will be added only if working on a different type of shift
//    initialBaseCost = consShiftsLbCosts_.at(pShiftIni->id) *
//        std::max(0, (consShiftsLbs_.at(pShiftIni->id) - initialConsumption));
//  }
//
//  // All following arcs added will have pOrigin node as their origin
//  // Iterate through all the possible successor shifts of the initial shift
//  for (auto indSuccessorShift : pShiftIni->successors) {
//    // The successor shift is either of the same type as the initial shift
//    // if starting from day -1 or of a different type.
//    // 1- Compute the initial consumption of this shift type:
//    // either 0 or initialConsumption if from source.
//    int initialC(0);
//    double initialBCost(0);
//    if (indSuccessorShift == pShiftIni->id) {
//      // Can work on the same shift only if starting from day -1.
//      // Otherwise, the stretch must change of shift type
//      if (pOrigin->day != -1) continue;
//      // set the initial consumption
//      initialC = initialConsumption;
//    } else {
//      // Arc cost due to penalty of soft lower bound of the initial
//      // shift if on a different shift and if was starting from day -1
//      initialBCost = initialBaseCost;
//    }
//
//    // 2- compute all possible stretch length of size <= UB or
//    // UB - initialC if starting from day -1 and on the same initial shift
//    PRCNode pTarget = pOrigin;
//    vector<PShift> vecShift;
//    vecShift.reserve(consShiftsUbs_.at(pShiftIni->id));
//    for (int d{1}; d <= consShiftsUbs_.at(indSuccessorShift) - initialC &&
//        d < pRCGraph_->nDays() - pOrigin->day; d++) {
//      // if sink, stop
//      if (pTarget->type == SINK_NODE) break;
//      // Target node of the arc that will be added
//      PShift pShift;
//      bool succFound = false;
//      for (const PRCArc &pArc : pTarget->outArcs) {
//        pShift = std::dynamic_pointer_cast<Shift>(pArc->target->pAShift);
// #ifdef NS_DEBUG
//        if (!pShift)
//          Tools::throwError("Enumerate works only with graph where"
//                            " the nodes represent a shift.");
// #endif
//        if (pShift->id == indSuccessorShift) {
//          pTarget = pRCGraph->pNode(pArc->target->id);
//          succFound = true;
//          break;
//        }
//      }
//
//      // if successor not found, break
//      if (!succFound) break;
//
//      // Creation of the stretch of the arc that will be added
//      vecShift.push_back(pShift);  // add one more same shift to the stretch
//      Stretch stretch(pOrigin->day+1, vecShift);
//
//      // we do not add the first sub path as it already exists in the graph
//      if (d == 1) continue;
//
//      // Recovery of the base cost of the arc stretch (taking account the
//      // previous shift, i.e. the initial shift)
//      double bCost = baseCost(stretch, pOrigin->pAShift);
//
//      // Arc cost due to penalty of soft lower bound of the corresponding
//      // shift (added only if the target node is not a sink node).
//      // Also, add initialC in case of starting from day -1, generally it's 0.
//      if (pTarget->type != SINK_NODE)
//        bCost += consShiftsLbCosts_.at(indSuccessorShift)
//            * std::max(0, consShiftsLbs_.at(indSuccessorShift)
//                - stretch.nDays() - initialC);
//
//      // Arc cost due to penalty of soft lower bound of the initial shift :
//      // generally 0
//      bCost += initialBCost;
//
//      // The new arc is added to the rcGraph with its corresponding costs
//      PRCArc pArc = pRCGraph->addSingleArc(pOrigin, pTarget, stretch, bCost);
//    }
//  }
//}
//
//
// void RCSPPSubProblem::updateOfExistingArcsCost(const PRCGraph &pRCGraph) {
//  // Recovery of the last shift performed by the current nurse (It can be a
//  // worked shift or a rest shift). This is the 'initial Shift'.
//  State* pInitialState = this->pLiveNurse_->pStateIni_;
//  PShift pShiftIni = pInitialState->pAShift_;
//
//  // Recovery of the number of times the last shift was completed by the
//  // current nurse (It's either a worked shift or a rest shift). This value
//  // is stored in the variable 'initialConsumption'.
//  int initialConsumption(0);
//  if (pShiftIni->isWork())
//    initialConsumption = pInitialState->consShifts_;
//  else if (pShiftIni->isRest())
//    initialConsumption = pInitialState->consDaysOff_;
//
//  // Iterate through all the arcs of the rcGraph which weren't added during
//  // the enumeration process (only the 'non-enumerated' arcs).
//  for (const auto &pArc : pRCGraph->pArcs()) {
//    if (pArc->isEnumeratedArc()) continue;
//
//    PShift pOShift = std::dynamic_pointer_cast<Shift>(pArc->origin->pAShift),
//        pTShift = std::dynamic_pointer_cast<Shift>(pArc->target->pAShift);
//    // Arcs coming from the source node are treated separately from the
//    // others.
//    if (pArc->origin->day != -1) {
//      if (pOShift->id == pTShift->id) {
//        // If the arc links two nodes with the same shift type, the cost due
//        // to the violation of soft upper bound of the corresponding shifts
//        // is added to the base cost of the arc.
//        pArc->addBaseCost(consShiftsUbCosts_.at(pOShift->id));
//      } else {
//        // Otherwise, the cost due to the violation of the soft lower bound
//        // of the target node's shifts is added to the base cost of the arc
//        // (Only if the target node is not a sink node)
//        if (pArc->target->type != SINK_NODE &&
//            consShiftsLbs_.at(pTShift->id) >= 2)
//          pArc->addBaseCost(consShiftsLbCosts_.at(pTShift->id) *
//              (consShiftsLbs_.at(pTShift->id) - 1));
//      }
//      continue;
//    }
//
//    // Otherwise, starting from day -1
//    if (pTShift->id == pShiftIni->id) {
//      // If the arc links two nodes with the same shift type,
//      // the costs due to the violation of soft bounds of
//      // the initial shift (taking account the initial consumption) are
//      // added to the base cost of the arc.
//      if (initialConsumption + 1 < consShiftsLbs_.at(pShiftIni->id) )
//        pArc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
//            (consShiftsLbs_.at(pShiftIni->id) - initialConsumption - 1));
//      if (initialConsumption + 1 > consShiftsUbs_.at(pShiftIni->id))
//        pArc->addBaseCost(consShiftsUbCosts_.at(pShiftIni->id) *
//            (initialConsumption+1 - consShiftsUbs_.at(pShiftIni->id)));
//    } else {
//      // Otherwise, the cost due to the violation of the soft lower bound
//      // of the target node's shifts and the cost due to the violation of
//      // the soft lower bound of the initial shift (taking account the
//      // initial consumption) is added to the base cost of the arc
//      if (consShiftsLbs_.at(pTShift->id) >= 2)
//        pArc->addBaseCost(consShiftsLbCosts_.at(pTShift->id) *
//            (consShiftsLbs_.at(pTShift->id)- 1));
//      if (initialConsumption < consShiftsLbs_.at(pShiftIni->id))
//        pArc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
//            (consShiftsLbs_.at(pShiftIni->id)-initialConsumption));
//    }
//  }
// }

PRCArc RCSPPSubProblem::addSingleArc(const PRCGraph &pRCGraph,
                                     const PRCNode &pOrigin,
                                     const PRCNode &pTarget,
                                     const PShift &pS,
                                     int day) {
  Stretch stretch(day, pS);
  // the base cost of the arc  is not computed here, it is computed at the
  // resource level when preprocessing the graph
  return pRCGraph->addSingleArc(pOrigin, pTarget, stretch, 0);
}

double RCSPPSubProblem::dualCost(const PRCArc &pArc) {
  return pCosts_->getCost(
      pLiveNurse_->num_, pArc->stretch, pArc->origin->pAShift);
}

void RCSPPSubProblem::createResources(const PRCGraph &pRCGraph) {
  pRCGraph->clearResources();
  for (const PResource& pR : pResources_) {
    if (!pR->isPreprocessed())
      pRCGraph->addResource(pR);
  }
}

void RCSPPSubProblem::createRCSPPSolver(int seed) {
  if (seed < 0) seed = rdm_();
  pRcsppSolver_ = std::make_shared<RCSPPSolver>(
      pRCGraph_, param_, seed, pRcsppSolver_ ? pRcsppSolver_.get() : nullptr);
}

void RCSPPSubProblem::preprocessRCGraph(const PRCGraph &pRCGraph,
                                        bool forceEnum) {
  // the preprocessing has three steps:
  // 1. add arcs to the RCGraph to deal with all the resources that are
  // enumerated delete them from the RCSPP
  // 2. modify the base cost of all the arcs to deal with the resources that
  // can be dealt with by modifying those costs instead of appearing in the
  // RCSPP
  // 3. Initialize the expanders of all the resources that will be considered
  // in the RCSPP

  // 1. enumeration of the resources
  timerEnumerationOfSubPath_.start();
  for (const PResource& pR : pResources_) {
    pR->enumerate(pRCGraph, forceEnum);
  }
  timerEnumerationOfSubPath_.stop();

  // 2. update the base costs of the arcs
  for (const PResource& pR : pResources_) {
    pR->preprocess(pRCGraph);
  }

  // 3. initialization of expanders and update of arcs base costs for the
  // resources of the RCSPP
  createResources(pRCGraph);
  pRCGraph->initializeExpanders();
}

void RCSPPSubProblem::updateParameters(bool useMoreTime) {
  // re-create a solver with the initial parameters
  createRCSPPSolver();
  if (useMoreTime)
    pRcsppSolver_->resetSearchLevel();
}

void RCSPPSubProblem::fixParameters() {
  if (param_.rcsppEnumSubpathsForMinCostToSinks_ &&
      (param_.rcsppEnumSubpaths_  || param_.rcsppMinCostToSinks_)) {
    std::cout << "Warning: rcsppEnumSubpathsForMinCostToSinks has no point"
                 " if either rcsppEnumSubpaths or rcsppMinCostToSinks is on: "
                 "it has been turned off" << std::endl;
    param_.rcsppEnumSubpathsForMinCostToSinks_ = false;
  }
  if (param_.rcsppToOptimality_ &&
      (param_.rcsppNbToExpand_ > 0 || param_.rcsppDssr_)) {
    std::cout << "Warning: rcsppNbToExpand_ or/and rcsppDssr_ are activated as "
                 "well as rcsppToOptimality_. Optimality will thus be "
                 "deactivated as will be reach if needed."  << std::endl;
    param_.rcsppToOptimality_ = false;
    if (param_.rcsppMinNegativeLabels_ <= 0) {
      std::cout << "Warning: rcsppMinNegativeLabels_ <= 0. It will be set to 1 "
                   "in order to reach optimality." << std::endl;
    }
  }
  if (!param_.rcsppToOptimality_ &&
      (param_.rcsppNbToExpand_ == 0 && !param_.rcsppDssr_)) {
    // optimality should be true as it's the real case
    param_.rcsppToOptimality_ = true;
  }
}

bool RCSPPSubProblem::presolve() {
  // update dual costs with the super method
  SubProblem::presolve();

  if (param_.rcsppResetParamAtEachIteration_) {
    // fix parameters
    fixParameters();

    // update parameters (create a new solver)
    updateParameters();
  }

  if (param_.rcsppEnumSubpathsForMinCostToSinks_) {
    if (pEnumGraph_ == nullptr)
      Tools::throwError("Enum graph has not been built.");
    for (const auto &pA : pEnumGraph_->pArcs())
      updateArcDualCost(pA);
    timerComputeMinCostFromSink_.start();
    // Computation of the costs of the shortest paths from each sink node
    // to all the other nodes in the rcGraph
    auto minEnumCosts = minCostPathToSinks(pEnumGraph_, pRCGraph_);
    pRcsppSolver_->setMinimumCostToSinks(minEnumCosts);
    timerComputeMinCostFromSink_.stop();
#ifdef NS_DEBUG
    auto minCosts = minCostPathToSinks(pRCGraph_);
    int equal = 0, bigger = 0;
    for (int n=0; n < pRCGraph_->nNodes(); ++n) {
      if (minCosts[n] >= minEnumCosts[n] + 1e-3)
        std::cerr << "Enum min costs (" << minEnumCosts[n]
                  << ") should never be lower than normal min cost ("
                  << minCosts[n] << ")" << std::endl;
      if (minCosts[n] <= minEnumCosts[n] - 1e-3)
        bigger++;
      else
        equal++;
    }
    std::cout << "Min enum costs is better " << bigger << " times and equal "
              << equal << " times." <<std::endl;
//     std::cout << " - Total time spent in computing minimum paths costs from "
//                  "sinks: " << timerComputeMinCostFromSink_.dSinceStart() <<
//                std::endl;
#endif
  } else if (param_.rcsppMinCostToSinks_) {
    timerComputeMinCostFromSink_.start();
    // Computation of the costs of the shortest paths from each sink node
    // to all the other nodes in the rcGraph
    pRcsppSolver_->setMinimumCostToSinks(minCostPathToSinks(pRCGraph_));
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

vector<double> RCSPPSubProblem::minCostPathToSinks(
    const PRCGraph &pRCGraph, const PRCGraph &pRCGraphOriginal) {
  // Initialization of all the costs of the shortest paths at infinity
  double inf = std::numeric_limits<double>::infinity();
  vector<double> minCostThrough(pRCGraph->nNodes(), inf);
  vector<double> minCostFrom(pRCGraph->nNodes(), inf);

  // Only the cost of the shortest path to the sink node is set to 0
  for (const auto &pN : pRCGraph->pSinks()) {
    minCostThrough.at(pN->id) = 0.0;
    minCostFrom.at(pN->id) = 0.0;
  }

  // Iteration through all the nodes of the acyclic graph in topological
  // order, starting from sourceDay
  for (int targetDay = this->nDays(); targetDay >= 0; targetDay--) {
    for (const PRCNode &pN : pRCGraph->pNodesPerDayId(targetDay)) {
      int targetId = pN->id;
      for (const auto &pArc : pN->inArcs) {
        if (pArc->forbidden)
          continue;  // do not use forbidden arcs

        // update shortest path if needed
        RCNode *pOrigin = pArc->origin;
        double cost = pArc->cost, t_cost = minCostFrom[targetId],
               o_cost = t_cost + cost;
        if (minCostFrom[pOrigin->id] > o_cost)
          minCostFrom[pOrigin->id] = o_cost;

        if (!pRCGraphOriginal) continue;

        // update the minimum cost to target when taking the full path
        if (minCostThrough[pOrigin->id] > o_cost)
          minCostThrough[pOrigin->id] = o_cost;

        // When disaggregating, we wish to be able to take a portion of any
        // long arc that goes through a given day/shift,
        // but when using the min cost in the RCSPP, we will complete it
        // with single day arcs whose costs should not be counted twice.
        // This needs to be done only for the last enumerated arc leading
        // to each given node.
        // All the intermediate arcs must be taken completely.
        // Therefore, we will go through each node of the original network
        // corresponding to the stretch of the last enumerated arc pArc.
        double inter_cost = t_cost;  // intermediary cost
        // find same target node in the original graph
        RCNode *pEndNode = pRCGraphOriginal->pNode(targetId).get();
        // iterate through the shift of the arc's stretch in reverse
        vector<PAbstractShift> originShifts = {pOrigin->pAShift};
        for (auto itS = pArc->stretch.pShifts().begin();
             itS != --pArc->stretch.pShifts().end(); ++itS)
          originShifts.push_back(*itS);
        auto itS = originShifts.end();
        int originDayId = pEndNode->dayId;
        while (itS != originShifts.begin()) {
          // update cost for all node except the origin one
          // (must take the full path for it)
          if (minCostThrough[pEndNode->id] > inter_cost)
            minCostThrough[pEndNode->id] = inter_cost;
          // move to previous day and find new origin
          --itS; --originDayId;
          // We iterate through the stretch in reverse,
          // find the previous node, and update the path cost
          for (const auto &pA : pEndNode->inArcs) {
            if (pA->origin->dayId == originDayId &&
                pA->origin->pAShift->equals(**itS)) {
              inter_cost += pA->cost;
              pEndNode = pA->origin;
              break;
            }
          }
        }
#ifdef NS_DEBUG
      if (inter_cost >= o_cost + 1e-3)
        std::cerr << "Through path cost (" << inter_cost
                  << ") should never be greater than enum min cost ("
                  << o_cost << ")" << std::endl;
#endif
      }
    }
  }
  return pRCGraphOriginal ? minCostThrough : minCostFrom;
}

void RCSPPSubProblem::updateArcDualCosts() {
  for (const auto& pA : pRCGraph_->pArcs())
    updateArcDualCost(pA);
}

void RCSPPSubProblem::updateArcDualCost(const PRCArc &pA) {
  pA->resetDualCost();
  pA->addDualCost(dualCost(pA));
}

bool RCSPPSubProblem::solveRCGraph(bool initialSolve, bool relaxation) {
  if (initialSolve) {
    // set the ids of the resources to make sure that the access to resource
    // values is correct
    int id = 0;
    for (const PResource &pR : pRCGraph_->pResources()) {
      pR->setId(id);
      id++;
    }

    // create the first label on the source node
    createInitialLabels();
  }

  // attach job to the solver
  pRcsppSolver_->attachJob(job_);

  // Solution of the RCSPP obtained with the RCSPP Solver
  theSolutions_ = pRcsppSolver_->solve(maxReducedCostBound_, relaxation);

  // Extract the best reduced cost
  if (!theSolutions_.empty())
    bestReducedCost_ = theSolutions_.front().reducedCost();

  return !theSolutions_.empty();
}


// Forbids a day-shift couple
void RCSPPSubProblem::forbidDayShift(int k, int s) {
  SubProblem::forbidDayShift(k, s);
  pRCGraph_->forbidDayShift(k, s);
  if (pEnumGraph_) pEnumGraph_->forbidDayShift(k, s);
}

// (re)Authorizes the day-shift couple
void RCSPPSubProblem::authorizeDayShift(int k, int s) {
  SubProblem::authorizeDayShift(k, s);
  pRCGraph_->authorizeDayShift(k, s);
  if (pEnumGraph_) pEnumGraph_->authorizeDayShift(k, s);
}

// Reset all authorizations to true
void RCSPPSubProblem::resetAuthorizations() {
  SubProblem::resetAuthorizations();
  pRCGraph_->resetAuthorizationsArcs();
  if (pEnumGraph_) pEnumGraph_->resetAuthorizationsArcs();
}

vector<PResource> RCSPPSubProblem::computeResourcesCosts(
    const State &initialState,
    const Stretch &stretch,
    std::map<CostType, double> *costsPerType) {
  RCSolution sol(stretch);
  auto pCosts = computeResourcesCosts(initialState, &sol);
  *costsPerType = sol.costs();
  return pCosts.first;
}

std::pair<vector<PResource>, PRCGraph> RCSPPSubProblem::computeResourcesCosts(
    const State &initialState, RCSolution *rcSol) const {
  // reset costs
  rcSol->resetCosts();
  // create origin, destination and arc
  PRCGraph pRCGraph = std::make_shared<RCGraph>(
          rcSol->firstDayId() - 1,
          rcSol->lastDayId() + 1,
          pScenario_->nShifts(),
          pScenario_->shiftsFactory().pAnyTypeShifts());
  PDay pFirstD = std::make_shared<Day>(rcSol->firstDayId() - 1),
       pLastD = std::make_shared<Day>(rcSol->lastDayId());
  PRCNode pSource =
          pRCGraph->addSingleNode(SOURCE_NODE, pFirstD, initialState.pShift_);
  PRCNode pSink =
          pRCGraph->addSingleNode(SINK_NODE, pLastD, rcSol->pShifts().back());
  PRCArc pArc = pRCGraph->addSingleArc(pSource, pSink, *rcSol, 0);

  // preprocess resources for arc
  for (const auto &pR : pResources_) {
    double c = 0;
    pR->preprocess(pArc, &c);
    rcSol->addCost(c, pR->costType());
  }

  // create expander for stretch
  double cost = pArc->cost;
  vector<PResource> pActiveResources;
  for (const auto &pR : pResources_) {
    if (pR->initialize(*pSource->pAShift, *rcSol, pArc,
                       static_cast<int>(pActiveResources.size())))
      pActiveResources.push_back(pR);
    rcSol->addCost(pArc->cost - cost, pR->costType());
    cost = pArc->cost;
  }

  // expand
  PRCLabel pLS = std::make_shared<RCLabel>(pActiveResources, initialState),
          pL = std::make_shared<RCLabel>(pActiveResources);
  pL->setAsNext(pLS, pArc);
  cost = pL->cost();
  for (const auto &pE : pArc->expanders) {
    ResourceValues &v = pL->getResourceValues(pE->indResource);
    pE->expand(pL, &v);
    rcSol->addCost(pL->cost() - cost, pE->costType);
    cost = pL->cost();
  }
#ifdef NS_DEBUG
  rcSol->pLabel_ = pL;
#endif

  return {pActiveResources, pRCGraph};
}
