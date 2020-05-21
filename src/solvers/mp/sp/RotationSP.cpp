//
// Created by antoine legrain on 2020-05-21.
//

#include "RotationSP.h"

RotationSP::RotationSP(PScenario scenario, int nbDays, PConstContract contract, std::vector<State>* pInitState):
  SubProblem(scenario, nbDays, contract,  pInitState) {}

RotationSP::~RotationSP() {}


//--------------------------------------------
//
// Functions for the NODES of the rcspp
//
//--------------------------------------------

// Function that creates the nodes of the network
void RotationSP::createNodes(){
  // INITIALIZATION
  principalGraphs_.clear();
  priceLabelsGraphs_.clear();

  // 1. SOURCE NODE
  //
  int v = addSingleNode(SOURCE_NODE);
  g_.setSource(v);

  // 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
  //
  principalGraphs_.emplace_back(PrincipalGraph(0,  nullptr)); // just add a dummy rcspp to have the right indices
  for(int sh=1; sh<pScenario_->nbShiftsType_; sh++)		// For each possible worked shift
    principalGraphs_.emplace_back(PrincipalGraph(sh, this));

  // 3. ROTATION LENGTH CHECK
  //
  // For each of the days, do a rotation-length-checker
  // Deactivate the min cost for the last day
  for(int k=0; k<nDays_-1; k++) {
    priceLabelsGraphs_.emplace_back(std::vector<PriceLabelGraph>(
        {PriceLabelGraph(k, pContract_->maxConsDaysWork_, maxRotationLength_, MAX_CONS_DAYS, this),
         PriceLabelGraph(k, 0, CDMin_, MIN_CONS_DAYS, this) }));

    // link the sub graphs
    priceLabelsGraphs_.back().front().linkOutSubGraph(priceLabelsGraphs_.back().back());

    // Daily sink node
    g_.addSink(priceLabelsGraphs_.back().back().exit());
  }
  // last day: price just max
  priceLabelsGraphs_.emplace_back(std::vector<PriceLabelGraph>(
      {PriceLabelGraph(nDays_-1, pContract_->maxConsDaysWork_, maxRotationLength_, MAX_CONS_DAYS, this)}));
  // Daily sink node
  g_.addSink(priceLabelsGraphs_.back().back().exit());

  // 4. SINK NODE
  //
  v = addSingleNode(SINK_NODE);
  g_.addSink(v);
}

//--------------------------------------------
//
// Functions for the ARCS of the rcspp
//
//--------------------------------------------

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void RotationSP::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (PrincipalGraph &pg: principalGraphs_)
    for (int k = minConsDays_ - 1; k < nDays_; k++)
      for (int dest : pg.getDayNodes(k)) {
        std::vector<int> vec;
        for (int s: pScenario_->shiftTypeIDToShiftID_[pg.shiftType()])
          vec.emplace_back(addSingleArc(origin, dest, 0, startConsumption(k, {s}),
                                        SOURCE_TO_PRINCIPAL, k, s));
        arcsFromSource_[pg.shiftType()][k].push_back(vec);
      }
}

// Create all arcs that involve the rotation size checking subnetwork (incoming, internal, and exiting that subnetwork)
void RotationSP::createArcsAllPriceLabels(){
  for(int k=0; k<nDays_; k++){				// For all days
    for(int sh=1; sh<pScenario_->nbShiftsType_; sh++){		// For all shifts
      // incoming  arc
      int origin = principalGraphs_[sh].exit(k);
      int destin = priceLabelsGraphs_[k].front().entrance();
      arcsPrincipalToPriceLabelsIn_[sh][k] =
          g_.addSingleArc(origin, destin, 0, {0,0,0,0,0}, PRINCIPAL_TO_PRICE_LABEL, k);	// Allow to stop rotation that day
    }

    // outgoing  arcs
    int origin = priceLabelsGraphs_[k].back().exit();
    int destin = g_.lastSink();
    arcsTosink_.push_back(g_.addSingleArc(origin, destin, 0, {0,0,0,0,0}, PRICE_LABEL_OUT_TO_SINK, k));
  }
}