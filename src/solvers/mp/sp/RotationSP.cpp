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

#include "RotationSP.h"

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
  for (int sh = 1; sh < pScenario_->nbShiftsType_; sh++)
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
  for (PrincipalGraph &pg : principalGraphs_)
    for (int k = minConsDays_ - 1; k < nDays_; k++)
      for (int dest : pg.getDayNodes(k)) {
        std::vector<int> vec;
        for (int s : pScenario_->shiftTypeIDToShiftID_[pg.shiftType()])
          vec.emplace_back(addSingleArc(origin,
                                        dest,
                                        0,
                                        startConsumption(k, {s}),
                                        SOURCE_TO_PRINCIPAL,
                                        k,
                                        s));
        arcsFromSource_[pg.shiftType()][k].push_back(vec);
      }
}

// Create all arcs from one principal subgraph to another one
void RotationSP::createArcsPrincipalToPrincipal() {
  int nShiftsType = pScenario_->nbShiftsType_;
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
                      pScenario_->nbShiftsType_,
                      nDays_,
                      -1);

  // Create pricing arcs for the label CONS_DAYS
  Penalties penalties = initPenalties();
  std::vector<int> cons = {-maxRotationLength_, 0, 0};
  // Arcs to daily sinks
  for (int k = 0; k < nDays_; k++) {
    // last day: price just max
    if (k == nDays_ - 1) penalties.minLevel(CONS_DAYS, 0);
    int s = g_.sink(k);
    for (int sh = 1; sh < pScenario_->nbShiftsType_; sh++) {
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

void RotationSP::updateArcCosts() {
  // A-B: call parent method first
  SubProblem::updateArcCosts();

  // C. ARCS : PRINCIPAL_TO_SINK
  // starts at 1, as on 0 we rest (so no end work) and also not defined
  for (int s = 1; s < pScenario_->nbShiftsType_; s++)
    for (int a : arcsPrincipalToSink[s])
      g_.updateCost(a, endWorkCost(a));
}
