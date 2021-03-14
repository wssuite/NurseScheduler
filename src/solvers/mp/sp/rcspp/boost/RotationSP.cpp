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

namespace boostRCSPP {

RotationSP::RotationSP(PScenario scenario,
                                 int nbDays,
                                 PConstContract contract,
                                 std::vector<State> *pInitState) :
    SubProblem(scenario, nbDays, contract, pInitState) {
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
  double cost = shiftCost(a);
  // add weekend cost if first day is sunday
  if (Tools::isSunday(start))
    cost -= pCosts_->workedWeekendCost();
  cost -= pCosts_->startWorkCost(start);
  cost += startWeekendCosts_[start];

  return cost;
}

double RotationSP::shiftCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // initial cost of the arc is the part of the cost that is common to all
  // nurses with the same contract
  double cost = arc_prop.initialCost;

  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for (const PShift &pS : arc_prop.pShifts) {
    cost += preferencesCosts_[k][pS->id];
    if (pS->isWork()) {
      cost -= pCosts_->workedDayShiftCost(k, pS->id);
      // add weekend cost if working on saturday
      // sunday exception is managed by startWorkCost
      if (Tools::isSaturday(k))
        cost -= pCosts_->workedWeekendCost();
    }
    ++k;
  }
  return cost;
}

double RotationSP::endWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  double cost = shiftCost(a);
  int length = arc_prop.pShifts.size(), end = arc_prop.day;
  // compute the end of the sequence of shifts
  // length could be equal to 0, but still represents the end of the rotation
  if (length > 1) end += length - 1;
  if (arc_prop.pShifts.empty() || arc_prop.pShifts.back())
    cost += endWeekendCosts_[end];
  cost -= pCosts_->endWorkCost(end);
  return cost;
}

}  // namespace boostRCSPP
