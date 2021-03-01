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

#include "solvers/mp/sp/MyRosterSP.h"

#include <algorithm>
#include <memory>
#include <set>
#include <utility>

#include "solvers/mp/RosterMP.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

MyRosterSP::MyRosterSP(PScenario pScenario,
                       int nDays,
                       PLiveNurse pNurse,
                       const SubproblemParam &param,
                       bool enumerateSubPath) :
    SubProblem(std::move(pScenario),
               nDays,
               std::move(pNurse),
               param), rcGraph_(nDays), pRcsppSolver_(nullptr),
    mrsp_enumeratedSubPathOption_(enumerateSubPath) {
  mrsp_sortLabelsOption_ = param.sp_sortLabelsOption;
  mrsp_minimumCostFromSinksOption_ = param.sp_minimumCostFromSinksOption;
  mrsp_worstCaseCostOption_ = param.sp_worstCaseCostOption;
  mrsp_enumeratedSubPathOption_ = param.sp_enumeratedSubPathOption;

  Tools::initVector2D<PRCNode>(&pNodesPerDayShift_, nDays_,
                               pScenario_->nbShifts(), nullptr);
}

void MyRosterSP::build() {
  rcGraph_.reset();
  // initialize the graph structure and base costs
  createNodes();
  SubProblem::initStructuresForSolve();
  createArcs();

  // set the status of all arcs to authorized
  Tools::initVector2D(&dayShiftStatus_, nDays_, pScenario_->nbShifts_, true);
  nPathsMin_ = 0;

  //  Initialization of the vectors that will contain all information
  //  corresponding to bounds and to costs associated with the violations of
  //  soft bounds of each shift
  Tools::initVector(&consShiftsUbs_, pScenario_->nbShiftsType_, 0);
  Tools::initVector(&consShiftsLbs_, pScenario_->nbShiftsType_, 0);
  Tools::initVector(&consShiftsUbCosts_, pScenario_->nbShiftsType_, .0);
  Tools::initVector(&consShiftsLbCosts_, pScenario_->nbShiftsType_, .0);

  for (int st = 0; st < pScenario_->nbShiftsType_; st++) {
    shared_ptr<AbstractShift> absShift = std::make_shared<AnyOfTypeShift>(st);
    if (absShift->isWork()) {
      consShiftsLbs_.at(st) = pScenario_->minConsShiftsOf(st);
      consShiftsUbs_.at(st) = pScenario_->maxConsShiftsOf(st);
      consShiftsLbCosts_.at(st) = pScenario_->pWeights_->WEIGHT_CONS_SHIFTS;
      consShiftsUbCosts_.at(st) = pScenario_->pWeights_->WEIGHT_CONS_SHIFTS;
    } else if (absShift->isRest()) {
      consShiftsLbs_.at(st) = pLiveNurse_->minConsDaysOff();
      consShiftsUbs_.at(st) = pLiveNurse_->maxConsDaysOff();
      consShiftsLbCosts_.at(st) = pScenario_->pWeights_->WEIGHT_CONS_DAYS_OFF;
      consShiftsUbCosts_.at(st) = pScenario_->pWeights_->WEIGHT_CONS_DAYS_OFF;
    }
  }

  // resources
  createResources();

  // numerate subpaths for rsources if needed
  if (mrsp_enumeratedSubPathOption_) {
    timerEnumerationOfSubPath_->start();
    // Enumeration of sub paths in the rcGraph
    enumerationSubPaths();
    timerEnumerationOfSubPath_->stop();
    if (param_.verbose_ >= 3)
      std::cout << " - Total time spent in enumeration of sub-paths: "
                << timerEnumerationOfSubPath_->dSinceStart() << std::endl;
  }

  // Initialization of each expander of each arc corresponding to each
  // resource
  rcGraph_.initializeExpanders();

  // Initialization of the timers
  timerEnumerationOfSubPath_ = new Tools::Timer();
  timerComputeMinCostFromSink_ = new Tools::Timer();

  // create solver
  this->createRCSPPSolver();
}


void MyRosterSP::createNodes() {
  rcGraph_.addSingleNode(SOURCE_NODE, -1, pLiveNurse_->pStateIni_->pShift_);
  // principal network is from day 0 to day nDays_-2
  for (int d = 0; d < nDays_ - 1; ++d)
    for (auto pShift : pScenario_->pShifts_)
      pNodesPerDayShift_[d][pShift->id] =
          rcGraph_.addSingleNode(PRINCIPAL_NETWORK, d, pShift);

  // every shift on the last day is a sink
  for (auto pShift : pScenario_->pShifts_)
    pNodesPerDayShift_[nDays_ - 1][pShift->id] =
        rcGraph_.addSingleNode(SINK_NODE, nDays_ - 1, pShift);
}


void MyRosterSP::createArcs() {
  // arcs from source to first day
  PShift pShiftIni = pScenario_->pShifts_[pLiveNurse_->pStateIni_->shift_];
  for (auto shiftId : pShiftIni->successors) {
    PRCNode pN = pNodesPerDayShift_[0][shiftId];
    PRCArc pArc = rcGraph_.addSingleArc(
        rcGraph_.pSource(), pN, Stretch(pN->pShift, 0, 1), 0);
  }

  // arcs from each day to next day; those from nDays-2 to nDays-1 are the
  // arcs to the sinks
  for (int d = 0; d < nDays_ - 1; ++d) {
    for (int shiftId = 0; shiftId < pScenario_->nbShifts(); ++shiftId) {
      PRCNode pOrigin = pNodesPerDayShift_[d][shiftId];
      for (int succId : pOrigin->pShift->successors) {
        PRCNode pTarget = pNodesPerDayShift_[d + 1][succId];
        Stretch stretch(pTarget->pShift, d + 1, 1);
        double cost = baseCost(stretch, pOrigin->pShift);
        rcGraph_.addSingleArc(pOrigin, pTarget, stretch, cost);
      }
    }
  }
}


void MyRosterSP::enumerationSubPaths() {
  // First, we enumerate sub paths from the source node taking account the
  // history of the current nurse
  enumerateConsShiftTypeFromSource();

  // Then, we enumerate sub paths from a PRINCIPAL_NETWORK node in the rcGraph
  // corresponding to a given day and to a given shift
  for (int d = 0; d < nDays_-1; ++d)
    for (auto pS : pScenario_->pShifts_)
      enumerateConsShiftType(pS, d);

  // Finally, the costs of the arcs which existed before the enumeration
  // process are modified
  updateOfExistingArcsCost();
}


void MyRosterSP::enumerateConsShiftTypeFromSource() {
  // Recovery of the last shift performed by the current nurse (It can be a
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
  PRCNode pOrigin = rcGraph_.pNode(0);
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
          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pShift);
        Stretch stretch(vecShift, 0, d + 1);

        // Recovery of the base cost of the arc stretch (taking account the
        // previous shift)
        double bCost = baseCost(stretch, pOrigin->pShift);

        // Recovery of the dual cost of the arc stretch
        double dCost = dualCost(stretch);

        // Arc cost due to penalty of soft lower bound of the corresponding
        // shift (this initial shift)
        double cost = consShiftsLbCosts_.at
            (pShiftIni->id) * std::max(0,
                                       consShiftsLbs_.at(pShiftIni->id) -
                                           initialConsumption - 1 - d);

        // The new arc is added to the rcGraph with its corresponding costs
        PRCArc pArc = rcGraph_.addSingleArc(pOrigin, pTarget, stretch, cost);
        pArc->addBaseCost(bCost);
        pArc->addDualCost(dCost);
      }
    } else {
      for (int d{1}; d <= consShiftsUbs_.at(indSuccessorShift)-1 && d <
          pScenario_->nbDays(); d++) {
        // Target node of the arc that will be added
        PRCNode pTarget = pNodesPerDayShift_[d][indSuccessorShift];

        // Creation of the stretch of the arc that will be added
        vector<PShift> vecShift(d+1, nullptr);
        for (size_t i = 0; i < d+1; i++)
          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pShift);
        Stretch stretch(vecShift, 0, d+1);

        // Recovery of the base cost of the arc stretch (taking account the
        // previous shift)
        double bCost = baseCost(stretch, pOrigin->pShift);

        // Recovery of the dual cost of the arc stretch
        double dCost = dualCost(stretch);

        // Arc cost due to penalty of soft lower bound of the corresponding
        // shift
        double cost = std::max(0, consShiftsLbs_.at(indSuccessorShift) - d - 1)
            * consShiftsLbCosts_.at(indSuccessorShift);

        // Arc cost due to penalty of soft lower bound of the initial
        // shift
        cost +=
            std::max(0, (consShiftsLbs_.at(pShiftIni->id) - initialConsumption))
                * consShiftsLbCosts_.at(pShiftIni->id);

        // The new arc is added to the rcGraph with its corresponding costs
        PRCArc pArc = rcGraph_.addSingleArc(pOrigin, pTarget, stretch, cost);
        pArc->addBaseCost(bCost);
        pArc->addDualCost(dCost);
      }
    }
  }
}

void MyRosterSP::enumerateConsShiftType(PShift pS, int day) {
  // All arcs that will be added will have as origin the node corresponding to
  // the day and to the shift given in parameter
  PRCNode pOrigin = pNodesPerDayShift_[day][pS->id];

  // Iterate through all the possible successor shifts of the origin node's
  // shift
  for (auto indSuccessorShift : pS->successors) {
    // New arcs will be added only toward nodes whose shift is different from
    // the origin node's one
    if (pS->id != indSuccessorShift) {
      for (int d = day+2; d <= day + consShiftsUbs_.at(indSuccessorShift)
          && d < pScenario_->nbDays(); d++) {
        // Target node of the arc that will be added
        PRCNode pTarget = pNodesPerDayShift_[d][indSuccessorShift];

        // Creation of the stretch of the arc that will be added
        vector<PShift> vecShift(d-day, nullptr);
        for (size_t i = 0; i < d-day; i++)
          vecShift.at(i) = std::make_shared<Shift>(*pTarget->pShift);
        Stretch stretch(vecShift, day+1, d-day);

        // Recovery of the base cost of the arc stretch (taking account the
        // previous shift)
        double bCost = baseCost(stretch, pOrigin->pShift);

        // Recovery of the dual cost of the arc stretch
        double dCost = dualCost(stretch);

        // The arc cost due to penalty of soft lower bound of the corresponding
        // shift is added only if the target node is not a sink node
        double cost = 0;
        if (pTarget->type != SINK_NODE) {
          int missingDays = consShiftsLbs_.at(indSuccessorShift) - (d - day);
          if (missingDays > 0)
            cost = consShiftsLbCosts_.at(indSuccessorShift) * missingDays;
        }

        // The new arc is added to the rcGraph with its corresponding costs
        PRCArc pArc = rcGraph_.addSingleArc(pOrigin, pTarget, stretch, cost);
        pArc->addBaseCost(bCost);
        pArc->addDualCost(dCost);
      }
    }
  }
}


void MyRosterSP::updateOfExistingArcsCost() {
  // Recovery of the last shift performed by the current nurse (It can be a
  // worked shift or a rest shift). This is the 'initial Shift'.
  State* initialState = this->pLiveNurse_->pStateIni_;
  PShift pShiftIni = pScenario_->pShifts_[initialState->shift_];

  // Recovery of the number of times the last shift was completed by the
  // current nurse (It's either a worked shift or a rest shift). This value
  // is stored in the variable 'initialConsumption'.
  int initialConsumption(0);
  if (pShiftIni->isWork())
    initialConsumption = initialState->consShifts_;
  else if (pShiftIni->isRest())
    initialConsumption = initialState->consDaysOff_;

  // Iterate through all the arcs of the rcGraph which weren't added during
  // the enumeration process (only the 'non-enumerated' arcs).
  for (const auto &arc : rcGraph_.pArcs()) {
    if (!arc->isEnumeratedArc()) {
      // Arcs coming from the source node are treated separately from the
      // others.
      if (arc->origin->id != 0) {
        if (arc->origin->pShift->id == arc->target->pShift->id) {
          // If the arc links two nodes with the same shift type, the cost due
          // to the violation of soft upper bound of the corresponding shifts
          // is added to the base cost of the arc.
          arc->addBaseCost(consShiftsUbCosts_.at(arc->origin->pShift->id));
        } else {
          // Otherwise, the cost due to the violation of the soft lower bound
          // of the target node's shifts is added to the base cost of the arc
          // (Only if the target node is not a sink node)
          if (arc->target->type != SINK_NODE) {
            if (consShiftsLbs_.at(arc->target->pShift->id) >= 2)
              arc->addBaseCost(consShiftsLbCosts_.at(arc->target->pShift->id) *
                  (consShiftsLbs_.at(arc->target->pShift->id)-1));
          }
        }
      } else {
        if (arc->target->pShift->id == pShiftIni->id) {
          // If the arc links two nodes with the same shift type, the costs due
          // to the violation of soft lower bound and of soft upper bound of
          // the initial shift (taking account the initial consumption) are
          // added to the base cost of the arc.
          if (initialConsumption + 1 < consShiftsLbs_.at(pShiftIni->id) )
            arc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
                (consShiftsLbs_.at(pShiftIni->id) - initialConsumption - 1));
          if (initialConsumption + 1 > consShiftsUbs_.at(pShiftIni->id))
            arc->addBaseCost(consShiftsUbCosts_.at(pShiftIni->id) *
                (initialConsumption+1 - consShiftsUbs_.at(pShiftIni->id)));
        } else {
          // Otherwise, the cost due to the violation of the soft lower bound
          // of the target node's shifts and the cost due to the violation of
          // the soft lower bound of the initial shift (taking account the
          // initial consumption) is added to the base cost of the arc
          if (consShiftsLbs_.at(arc->target->pShift->id) >= 2)
            arc->addBaseCost(consShiftsLbCosts_.at(arc->target->pShift->id) *
                (consShiftsLbs_.at(arc->target->pShift->id)- 1));
          if (initialConsumption < consShiftsLbs_.at(pShiftIni->id))
            arc->addBaseCost(consShiftsLbCosts_.at(pShiftIni->id) *
                (consShiftsLbs_.at(pShiftIni->id)-initialConsumption));
        }
      }
    }
  }
}


double MyRosterSP::dualCost(const Stretch &stretch) {
  double dualCost = 0;
  int curDay = stretch.firstDay();
  // if start, add constant
  if (curDay == 0)
    dualCost -= pCosts_->constant();
  // iterate through the shift to update the cost
  for (auto pS : stretch.pShifts()) {
    if (pS->isWork())
      dualCost -= pCosts_->workedDayShiftCost(curDay, pS->id);
    curDay++;
  }
  return dualCost;
}

double MyRosterSP::baseCost(const Stretch &stretch, PShift pPrevShift) {
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

void MyRosterSP::createResources() {
  rcGraph_.clearResources();

  PNurse pN = this->pLiveNurse_;
  PScenario pS = this->pScenario_;
  PWeights pW = pS->pWeights_;

  // initialize resource on the total number of working days
  rcGraph_.addResource(std::make_shared<SoftTotalShiftResource>(
      pN->minTotalShifts(),
      pN->maxTotalShifts(),
      pW->WEIGHT_TOTAL_SHIFTS,
      pW->WEIGHT_TOTAL_SHIFTS,
      std::make_shared<AnyWorkShift>(),
      rcGraph_.nDays()));

  // initialize resource on the total number of week-ends
  rcGraph_.addResource(std::make_shared<SoftTotalWeekendsResource>(
      pN->maxTotalWeekends(),
      pW->WEIGHT_TOTAL_WEEKENDS,
      rcGraph_.nDays()));

  // initialize resource on the number of consecutive worked days
  rcGraph_.addResource(std::make_shared<SoftConsShiftResource>(
      pN->minConsDaysWork(),
      pN->maxConsDaysWork(),
      pW->WEIGHT_CONS_DAYS_WORK,
      pW->WEIGHT_CONS_DAYS_WORK,
      std::make_shared<AnyWorkShift>(),
      rcGraph_.nDays()));

  if (!mrsp_enumeratedSubPathOption_) {
    // initialize resources on the number of consecutive shifts of each type
    for (int st = 0; st < pS->nbShiftsType_; st++) {
      shared_ptr<AbstractShift> absShift = std::make_shared<AnyOfTypeShift>(st);
      if (absShift->isWork())
        rcGraph_.addResource(std::make_shared<SoftConsShiftResource>(
            consShiftsLbs_.at(st),
            consShiftsUbs_.at(st),
            consShiftsLbCosts_.at(st),
            consShiftsUbCosts_.at(st),
            absShift,
            rcGraph_.nDays()));
      else if (absShift->isRest())
        rcGraph_.addResource(std::make_shared<SoftConsShiftResource>(
            consShiftsLbs_.at(st),
            consShiftsUbs_.at(st),
            consShiftsLbCosts_.at(st),
            consShiftsUbCosts_.at(st),
            absShift,
            rcGraph_.nDays()));
    }
  }
}

void MyRosterSP::createRCSPPSolver() {
  pRcsppSolver_ = std::make_shared<MyRCSPPSolver>(&rcGraph_);
}

void MyRosterSP::preprocessRCGraph() {
  initializeResources();
}

void MyRosterSP::initializeResources() {
  rcGraph_.initializeDominance();
  // rcGraph_.initializeExpanders();
}

bool MyRosterSP::preprocess() {
  // add dual costs to the arcs costs
  updateArcDualCosts();

  // print graph with updated costs
  // rcGraph_.printSummaryOfGraph();

  // reset solver in case it is not the first time it is called
  pRcsppSolver_->reset();

  if (mrsp_minimumCostFromSinksOption_) {
    timerComputeMinCostFromSink_->start();
    // Computation of the costs of the shortest paths from each sink node
    // to all the other nodes in the rcGraph
    pRcsppSolver_->computeMinimumCostToSinks();
    timerComputeMinCostFromSink_->stop();
    if (param_.verbose_ >= 3)
      std::cout << " - Total time spent in computing minimum paths costs "
                << "from sinks: " << timerComputeMinCostFromSink_->dSinceStart()
                << std::endl;
  }

  return true;
}

void MyRosterSP::updateArcDualCosts() {
  for (const auto& pA : rcGraph_.pArcs())
    updateArcDualCost(pA);
}

void MyRosterSP::updateArcDualCost(PRCArc pA) {
  pA->resetDualCost();
  pA->addDualCost(dualCost(pA->stretch));
}

bool MyRosterSP::solveRCGraph() {
  // create the first label on the source node
  createInitialLabel();

  // Solution of the RCSPP obtained with the RCSPP Solver
  std::vector<RCSolution> solutions =
      pRcsppSolver_->solve(maxReducedCostBound_,
                           param_.verbose_,
                           param_.epsilon_,
                           param_.search_strategy_,
                           param_.nb_max_paths_,
                           pScenario_,
                           consShiftsUbs_,
                           consShiftsLbs_,
                           mrsp_enumeratedSubPathOption_,
                           mrsp_sortLabelsOption_,
                           mrsp_minimumCostFromSinksOption_,
                           mrsp_worstCaseCostOption_);

  // Extract the best reduced cost
  bestReducedCost_ = maxReducedCostBound_;
  for (const RCSolution &sol : solutions) {
    // only add solutions with a cost < maxReducedCostBound_
    if (sol.cost > maxReducedCostBound_)
      continue;
    theSolutions_.push_back(sol);
    nPaths_++;
    nFound_++;
    if (bestReducedCost_ > sol.cost)
      bestReducedCost_ = sol.cost;
  }

  return !solutions.empty();
}

void MyRosterSP::createInitialLabel() {
  auto pL = std::make_shared<RCLabel>(rcGraph_.pResources(),
                                      *pLiveNurse_->pStateIni_);
  pL->setNode(rcGraph_.pSource());
  pRcsppSolver_->addLabel(pL);
}


