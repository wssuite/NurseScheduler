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

#include "RotationSP.h"

#include <algorithm>

namespace boostRCSPP {

RotationSP::RotationSP(PScenario scenario,
                       int nDays,
                       PLiveNurse pNurse,
                       SubProblemParam param) :
    SubProblem(scenario, nDays, pNurse, param) {
  labels_ = {CONS_DAYS};
}

RotationSP::~RotationSP() {}


//--------------------------------------------
//
// Functions for the NODES of the rcspp
//
//--------------------------------------------

// Function that creates the nodes of the network
void RotationSP::createNodes() {
  // INITIALIZATION
  principalGraphs_.clear();

  // 1. SOURCE NODE
  int v = addSingleNode(SOURCE_NODE);
  g_.setSource(v);

  // 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
  // just add a dummy rcspp to have the right indices
  principalGraphs_.emplace_back(PrincipalGraph(0, nullptr));
  // For each possible worked shift
  for (int sh = 1; sh < pScenario_->nShiftTypes(); sh++)
    principalGraphs_.emplace_back(PrincipalGraph(sh, this));

  // 3. DAILY SINKS AND GLOBAL SINK
  for (int k = 0; k < nDays_; k++) {
    v = addSingleNode(SINK_NODE);
    g_.addSink(v);
  }
  v = addSingleNode(SINK_NODE);
  g_.addSink(v);
}

//--------------------------------------------
//
// Functions for the ARCS of the rcspp
//
//--------------------------------------------

// Create all arcs whose origin is the source nodes (all go to short
// rotations nodes)
void RotationSP::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (int sh=1; sh < pScenario_->nShiftTypes(); ++sh)
    for (int k = minConsDays_ - 1; k < nDays_; k++)
      for (int dest : principalGraphs_[sh].getDayNodes(k)) {
        std::vector<int> vec;
        for (int s : pScenario_->shiftTypeIDToShiftID(sh))
          vec.emplace_back(addSingleArc(
              origin, dest, 0, startConsumption(k, {pScenario_->pShift(s)}),
              SOURCE_TO_PRINCIPAL, k, pScenario_->pShift(s)));
        arcsFromSource_[sh][k].push_back(vec);
      }
}

// Create all arcs from one principal subgraph to another one
void RotationSP::createArcsPrincipalToPrincipal() {
  int nShiftsType = pScenario_->nShiftTypes();
  for (int sh = 1; sh < nShiftsType; sh++) {
    for (int newSh = 1; newSh < nShiftsType; newSh++) {
      // check if succession is allowed
      if (newSh != sh
          && !pScenario_->isForbiddenSuccessorShiftType_ShiftType(newSh, sh)) {
        // do not create arcs for the last day, as it won't be possible to work
        for (int k = 0; k < nDays_ - 1; k++) {
          // last level for day k
          int o = principalGraphs_[sh].exit(k);
          // entrance level for day k
          int d = principalGraphs_[newSh].entrance(k);
          principalToPrincipal_[sh][newSh][k] =
              g_.addSingleArc(o, d, 0, {0, 0, 0}, SHIFT_TO_NEWSHIFT, k);
        }
      }
    }
  }
}

void RotationSP::createArcsPrincipalToSink() {
  Tools::initVector2D(&arcsPrincipalToSink,
                      pScenario_->nShiftTypes(),
                      nDays_,
                      -1);

  // Create pricing arcs for the label CONS_DAYS
  Penalties penalties = initPenalties();
  std::vector<int> cons = {-10*maxRotationLength_, 0, 0};
  // Arcs to daily sinks
  for (int k = 0; k < nDays_; k++) {
    // last day: price just max
    if (k == nDays_ - 1) penalties.minLevel(CONS_DAYS, 0);
    int s = g_.sink(k);
    for (int sh = 1; sh < pScenario_->nShiftTypes(); sh++) {
      // incoming  arc
      int o = principalGraphs_[sh].exit(k);
      arcsPrincipalToSink[sh][k] =
          g_.addPricingArc(o, s, 0, cons, TO_SINK, k, {CONS_DAYS}, penalties);
    }

    // Arcs from daily sinks to global one
    int o = s;
    int d = g_.lastSink();
    arcsTosink_.push_back(
        g_.addSingleArc(o, d, 0, {0, 0, 0}, TO_SINK, nDays_ - 1));
  }
}

// Updates the costs depending on the reduced costs given for the nurse
void RotationSP::updateArcDualCosts() {
  // A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
  // B. ARCS : PRINCIPAL GRAPH
  SubProblem::updateArcDualCosts();

  // C. ARCS : PRINCIPAL_TO_SINK
  // starts at 1, as on 0 we rest (so no end work) and also not defined
  for (int s = 1; s < pScenario_->nShiftTypes(); s++)
    for (int a : arcsPrincipalToSink[s])
      g_.updateCost(a, endWorkCost(a));
}

double RotationSP::startWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // retrieve the work cost
  int start = arc_prop.day;
  // test is true when the arc does not match an assignment;
  // it may be then end of a rotation for instance
  if (arc_prop.pShifts.empty())
    ++start;  // start work the next day as no shift today
  double cost = shiftCost(a, nullptr);
  cost += startWeekendCosts_[start];

  return cost;
}

double RotationSP::shiftCost(int a, const PAbstractShift &prevS) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // initial cost of the arc is the part of the cost that is common to all
  // nurses with the same contract
  double cost = arc_prop.initialCost;

  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for (const PShift &pS : arc_prop.pShifts)
    cost += preferencesCosts_[k++][pS->id];
  Stretch st(arc_prop.day, arc_prop.pShifts);
  cost -= pCosts_->getCost(pLiveNurse_->num_, st, prevS);
  return cost;
}

double RotationSP::endWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  PShift pW = pScenario_->pAnyWorkShift();
  double cost = shiftCost(a, pW);
  int length = arc_prop.pShifts.size(), end = arc_prop.day;
  // compute the end of the sequence of shifts
  // length could be equal to 0, but still represents the end of the rotation
  if (length > 1) end += length - 1;
  if (arc_prop.pShifts.empty() || arc_prop.pShifts.back())
    cost += endWeekendCosts_[end];
  // create a stretch of 1 rest day at the end to get the end cost
  // if not last day
  if (++end <  nDays()) {
    Stretch st(end, pScenario_->pRestShift());
    cost -= pCosts_->getCost(pLiveNurse_->num_, st, pW);
  }
  return cost;
}

void RotationSP::computeCost(MasterProblem *, RCSolution *rcSol) const {
  /************************************************
   * Compute all the costs of a rotation:
   ************************************************/
#ifdef DBG
  double cost = rcSol->cost();
#endif
  rcSol->resetCosts();

  // if first day of the planning, check on the past, otherwise 0 (rest)
  int lastShiftType = (rcSol->firstDay() == 0) ?
      pLiveNurse_->pStateIni_->pShift_->type : -1;
  // nbConsShift = number of consecutive shift
  // if first day of the planning, check on the past, otherwise 0
  int nbConsShifts = (rcSol->firstDay() == 0) ?
      pLiveNurse_->pStateIni_->consShifts_ : 0;

  // nbConsWorked = number of consecutive worked days
  // if first day of the planning, check on the past , otherwise 0
  int nbConsDaysWorked =
      (rcSol->firstDay() == 0) ? pLiveNurse_->pStateIni_->consDaysWorked_ : 0;

  /*
   * Compute consShiftCost
   */

  // if the initial shift has already exceeded the max, substract now the cost
  // that will be readd later
  if ((rcSol->firstDay() == 0) && (lastShiftType > 0) &&
      (nbConsShifts > pScenario_->maxConsShiftsOfType(lastShiftType))) {
    rcSol->addCost(
        -(nbConsShifts - pScenario_->maxConsShiftsOfType(lastShiftType))
            * pScenario_->weights().WEIGHT_CONS_SHIFTS,
        CONS_SHIFTS_COST);
  }

  for (int k = rcSol->firstDay(); k <= rcSol->lastDay(); ++k) {
    const PShift &pS = rcSol->pShift(k);
    if (lastShiftType == pS->type) {
      nbConsShifts++;
      continue;
    }
    if (lastShiftType > 0) {
      int diff = std::max(
          pScenario_->minConsShiftsOfType(lastShiftType) - nbConsShifts,
          nbConsShifts - pScenario_->maxConsShiftsOfType(lastShiftType));
      if (diff > 0)
        rcSol->addCost(diff * pScenario_->weights().WEIGHT_CONS_SHIFTS,
                       CONS_SHIFTS_COST);
    }
    // initialize nbConsShifts and lastShift
    nbConsShifts = 1;
    lastShiftType = pS->type;
  }

  // compute consShiftsCost for the last shift
  if (lastShiftType > 0) {
    int diff =
        std::max((rcSol->lastDay()+1 == nDays()) ? 0 :
                 pScenario_->minConsShiftsOfType(lastShiftType) - nbConsShifts,
                 nbConsShifts - pScenario_->maxConsShiftsOfType(lastShiftType));
    if (diff > 0)
      rcSol->addCost(diff * pScenario_->weights().WEIGHT_CONS_SHIFTS,
                     CONS_SHIFTS_COST);
  }


  /*
   * Compute consDaysWorkedCost
   */

  // if already worked too much
  double diffDays =
      nbConsDaysWorked - pLiveNurse_->pContract_->maxConsDaysWork_;
  rcSol->addCost(diffDays > 0 ?
                 -diffDays * pScenario_->weights().WEIGHT_CONS_DAYS_WORK : 0,
                 CONS_WORK_COST);

  nbConsDaysWorked += rcSol->nDays();
  // check if nbConsDaysWorked < min, if finishes on last day, does not count
  // nbConsDaysWorked > 0, normally should always be the case,
  // but could for the initial variables of the rotation graph
  if (nbConsDaysWorked > 0 && nbConsDaysWorked < pLiveNurse_->minConsDaysWork()
      && rcSol->lastDay() + 1 < nDays())
    rcSol->addCost((pLiveNurse_->minConsDaysWork() - nbConsDaysWorked)
                       * pScenario_->weights().WEIGHT_CONS_DAYS_WORK,
                   CONS_WORK_COST);
  else if (nbConsDaysWorked > pLiveNurse_->maxConsDaysWork())
    // check if nbConsDaysWorked > max
    rcSol->addCost((nbConsDaysWorked - pLiveNurse_->maxConsDaysWork())
                       * pScenario_->weights().WEIGHT_CONS_DAYS_WORK,
                   CONS_WORK_COST);

  /*
   * Compute completeWeekendCost
   */
  if (pLiveNurse_->needCompleteWeekends()) {
    // if first day is a Sunday, the saturday is not worked
    if (Tools::isSunday(rcSol->firstDay()))
      rcSol->addCost(pScenario_->weights().WEIGHT_COMPLETE_WEEKEND,
          COMPLETE_WEEKEND_COST);
    // if last day + 1 is a Sunday, the sunday is not worked
    if (Tools::isSunday(rcSol->lastDay()+1))
      rcSol->addCost(pScenario_->weights().WEIGHT_COMPLETE_WEEKEND,
                     COMPLETE_WEEKEND_COST);
  }

  computePreferencesCost(rcSol);

  /*
   * Compute initial resting cost if rotation not empty
   * (used to price initial rest arcs)
   */
  if (rcSol->firstDay() == 0 && rcSol->nDays() > 0 &&
      pLiveNurse_->pStateIni_->pShift_->isRest()) {
    int diff =
        pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
    rcSol->addCost(
        (diff > 0) ? diff * pScenario_->weights().WEIGHT_CONS_DAYS_OFF : 0,
        CONS_REST_COST);
  }

  /*
   * Compute the sum of the cost and stores it in cost_
   */
#ifdef DBG
  if (cost < DBL_MAX - 1 && std::abs(cost - rcSol->cost()) > EPSILON) {
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: " << rcSol->cost() << " != " << cost
              << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    Tools::throwError("boostRCSPP::computeCost does not get the same cost.");
  }
#endif
}

}  // namespace boostRCSPP
