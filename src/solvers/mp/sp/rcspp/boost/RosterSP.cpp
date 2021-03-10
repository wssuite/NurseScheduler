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

#include "RosterSP.h"

#include <string>
#include <algorithm>

#include "RCSPP.h"

using std::string;
using std::vector;
using std::map;

namespace boostRCSPP {

// Constructors and destructor
RosterSP::RosterSP(PScenario scenario,
                   int nbDays,
                   PConstContract contract,
                   vector<State> *pInitState) :
    SubProblem(scenario, nbDays, contract, pInitState) {
  labels_ = {CONS_DAYS, DAYS, WEEKEND};
  build();
}

RosterSP::~RosterSP() {}

BoostRCSPPSolver *RosterSP::initRCSSPSolver() {
  double constant = pCosts_->constant();
  Penalties penalties = initPenalties();
  // lambda expression to post process the solutions found by the RCSPP solver
  auto postProcess = [constant, penalties](spp_res_cont *res_cont) {
    res_cont->postprocessCost =
        penalties.penalty(DAYS, res_cont->label_value(DAYS)) +
            penalties.penalty(WEEKEND, res_cont->label_value(WEEKEND)) -
            constant;
    res_cont->cost += res_cont->postprocessCost;
  };
  return new BoostRCSPPSolver(&g_,
                              maxReducedCostBound_,
                              param_.verbose_,
                              param_.epsilon_,
                              param_.search_strategy_,
                              param_.rcsppMinNegativeLabels_,
                              postProcess);
}

// Function that creates the nodes of the network
void RosterSP::createNodes() {
  // INITIALIZATION
  principalGraphs_.clear();

  // 1. SOURCE NODE
  int v = addSingleNode(SOURCE_NODE);
  g_.setSource(v);

  // 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
  // For each possible worked shift
  for (int sh = 0; sh < pScenario_->nShiftTypes(); sh++)
    principalGraphs_.emplace_back(PrincipalGraph(sh, this));

  // 3. SINK
  v = addSingleNode(SINK_NODE);
  g_.addSink(v);
}

void RosterSP::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (PrincipalGraph &pg : principalGraphs_) {
    vector2D<int> vec2;
    for (int dest : pg.getDayNodes(0)) {
      std::vector<int> vec;
      for (int s : pScenario_->shiftTypeIDToShiftID(pg.shiftType()))
        vec.emplace_back(addSingleArc(origin,
                                      dest,
                                      0,
                                      {},
                                      SOURCE_TO_PRINCIPAL,
                                      0,
                                      s));
      vec2.push_back(vec);
    }
    arcsFromSource_[pg.shiftType()] = {vec2};
  }
}

// Create all arcs from one principal subgraph to another one
// add a pricing arc when going from work to rest
void RosterSP::createArcsPrincipalToPrincipal() {
  int nShiftsType = pScenario_->nShiftTypes();
  Penalties penalties = initPenalties();
  // consumption when pricing to reset CONS_DAYS
  std::vector<int> pcons = {-10*maxRotationLength_, 0, 0};
  for (int sh = 0; sh < nShiftsType; sh++) {
    for (int newSh = 0; newSh < nShiftsType; newSh++) {
      // check if succession is allowed
      if (newSh != sh
          && !pScenario_->isForbiddenSuccessorShiftType_ShiftType(newSh, sh)) {
        for (int k = 0; k < nDays_ - 1; k++) {
          // last level for day k
          int o = principalGraphs_[sh].exit(k);
          // entrance level for day k
          int d = principalGraphs_[newSh].entrance(k);
          // if ends a rotation, add a pricing arc
          int a;
          if (newSh == 0) {
            a = g_.addPricingArc(o, d, 0, pcons, SHIFT_TO_NEWSHIFT, k,
                                 {CONS_DAYS}, penalties);
          } else {
            // otherwise, just add a normal arc
            // if start work on sunday (rest to start shift)
            bool w = (sh == 0 && Tools::isSunday(k + 1));
            a = g_.addSingleArc(o, d, 0, {0, 0, w}, SHIFT_TO_NEWSHIFT, k);
          }
          principalToPrincipal_[sh][newSh][k] = a;
        }
      }
    }
  }
}

void RosterSP::createArcsPrincipalToSink() {
  Penalties penalties = initPenalties();
  // last day: price just max
  penalties.minLevel(CONS_DAYS, 0);
  // consumption to reset CONS_DAYS
  std::vector<int> cons = {-maxRotationLength_, 0, 0};
  int d = g_.lastSink();
  for (int sh = 0; sh < pScenario_->nShiftTypes(); sh++) {
    // incoming  arc
    int o = principalGraphs_[sh].exit(nDays_ - 1);
    arcsTosink_.push_back(
        g_.addPricingArc(o, d, 0, cons, TO_SINK, nDays_ - 1,
                         {CONS_DAYS}, penalties));
  }
}

void RosterSP::updateArcDualCosts() {
  // A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
  // B. ARCS : PRINCIPAL GRAPH
  SubProblem::updateArcDualCosts();

  // C. ARCS : PRINCIPAL GRAPH TO PRINCIPAL GRAPH
  // for rest principal graph
  for (int s = 1; s < pScenario_->nShiftTypes(); s++) {
    for (int a : principalToPrincipal_[0][s])
      if (a != -1) g_.updateCost(a, startWorkCost(a));
    for (int a : principalToPrincipal_[s][0])
      if (a != -1) g_.updateCost(a, endWorkCost(a));
  }
}

double RosterSP::startWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // retrieve the work cost
  double cost = shiftCost(a);
  int start = arc_prop.day;
  // test is true when the arc does not match an assignment;
  // it may be then end of a rotation for instance
  if (arc_prop.shifts.empty())
    ++start;  // start work the next day as no shift today
  cost += startWeekendCosts_[start];

  return cost;
}

double RosterSP::shiftCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // initial cost of the arc is the part of the cost that is common to all
  // nurses with the same contract
  double cost = arc_prop.initialCost;

  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for (int s : arc_prop.shifts) {
    cost += preferencesCosts_[k][s];
    if (pScenario_->isWorkShift(s))
      cost -= pCosts_->workedDayShiftCost(k, s);
    ++k;
  }
  return cost;
}

double RosterSP::endWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  double cost = shiftCost(a);
  int length = arc_prop.shifts.size(), end = arc_prop.day;
  // compute the end of the sequence of shifts
  // length could be equal to 0, but still represents the end of the rotation
  if (length > 1) end += length - 1;
  if (arc_prop.shifts.empty() || arc_prop.shifts.back())
    cost += endWeekendCosts_[end];
  return cost;
}

}  // namespace boostRCSPP
