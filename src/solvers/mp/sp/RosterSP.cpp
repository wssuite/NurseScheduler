/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 *  full license detail.
 */

#include "RosterSP.h"

#include <string>
#include <algorithm>

#include "solvers/mp/sp/rcspp/BoostRCSPP.h"

using std::string;
using std::vector;
using std::map;

// Constructors and destructor
RosterSP::RosterSP(PScenario scenario,
                   int nbDays,
                   PConstContract contract,
                   vector<State> *pInitState) :
    SubProblem(scenario, nbDays, contract, pInitState) {
  labels_ = {CONS_DAYS, DAYS, WEEKEND};
}

RosterSP::~RosterSP() {}

RCSPPSolver *RosterSP::initRCSSPSolver() {
  double constant = pCosts_->constant();
  Penalties penalties = initPenalties();
  // lambda expression to post process the solutions found by the RCSPP solver
  auto postProcess = [constant, penalties](spp_res_cont *res_cont) {
    res_cont->cost -= constant;
    res_cont->cost += penalties.penalty(DAYS, res_cont->label_value(DAYS));
    res_cont->cost +=
        penalties.penalty(WEEKEND, res_cont->label_value(WEEKEND));
  };
  return new BoostRCSPPSolver(&g_,
                              maxReducedCostBound_,
                              param_.verbose_,
                              param_.epsilon_,
                              param_.search_strategy_,
                              param_.nb_max_paths_,
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
  for (int sh = 0; sh < pScenario_->nbShiftsType_; sh++)
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
      for (int s : pScenario_->shiftTypeIDToShiftID_[pg.shiftType()])
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
  int nShiftsType = pScenario_->nbShiftsType_;
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
  for (int sh = 0; sh < pScenario_->nbShiftsType_; sh++) {
    // incoming  arc
    int o = principalGraphs_[sh].exit(nDays_ - 1);
    arcsTosink_.push_back(
        g_.addPricingArc(o, d, 0, cons, TO_SINK, nDays_ - 1,
                         {CONS_DAYS}, penalties));
  }
}

void RosterSP::updateArcCosts() {
  // A-B: call parent method first
  SubProblem::updateArcCosts();

  // C. ARCS : PRINCIPAL GRAPH TO PRINCIPAL GRAPH
  // for rest principal graph
  for (int s = 1; s < pScenario_->nbShiftsType_; s++) {
    for (int a : principalToPrincipal_[0][s])
      if (a != -1) g_.updateCost(a, startWorkCost(a));
    for (int a : principalToPrincipal_[s][0])
      if (a != -1) g_.updateCost(a, endWorkCost(a));
  }
}
