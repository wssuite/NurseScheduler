//
// Created by antoine legrain on 2020-04-19.
//

#include "RosterSP.h"

using std::string;
using std::vector;
using std::map;

// Constructors and destructor
RosterSP::RosterSP(): SubProblem() {}

RosterSP::RosterSP(Scenario* scenario, int nbDays, const Contract* contract, vector<State>* pInitState):
    SubProblem(scenario, nbDays, contract,  pInitState) {

  nLabels_ = 5;
}

RosterSP::~RosterSP(){}

std::function<void (spp_res_cont&)> RosterSP::postProcessResCont() const {
  double constant = pCosts_->constant();
  int max_days = pLiveNurse_->maxTotalShifts(),
      max_weekends = pLiveNurse_->maxTotalWeekends();
  return [constant, max_days, max_weekends](spp_res_cont& res_cont) {
    res_cont.cost -= constant;
    res_cont.cost += res_cont.label_value(MIN_DAYS) * WEIGHT_TOTAL_SHIFTS;
    res_cont.cost += std::max(0, res_cont.label_value(MAX_DAYS) - max_days) * WEIGHT_TOTAL_SHIFTS;
    res_cont.cost += std::max(0, res_cont.label_value(MAX_WEEKEND) - max_weekends) * WEIGHT_TOTAL_WEEKENDS;
  };
}

// Function that creates the nodes of the network
void RosterSP::createNodes() {
  // INITIALIZATION
  principalGraphs_.clear();
  priceLabelsGraphs_.clear();

  // 1. SOURCE NODE
  //
  int v = addSingleNode(SOURCE_NODE);
  g_.setSource(v);

  // 2. PRINCIPAL NETWORK(S) [ONE PER SHIFT TYPE]
  //
  for(int sh=0; sh<pScenario_->nbShiftsType_; sh++)		// For each possible worked shift
    principalGraphs_.emplace_back(PrincipalGraph(sh, this));

  // 3. PRICE LABELS
  //
  // For each of the days, price labels for min and max consecutive days and reset these labels
  for(int k=0; k<nDays_-1; k++) {
    priceLabelsGraphs_.emplace_back(vector<PriceLabelGraph>(
        {PriceLabelGraph(k, pContract_->maxConsDaysWork_, maxRotationLength_, MAX_CONS_DAYS, this, true),
         PriceLabelGraph(k, 0, CDMin_, MIN_CONS_DAYS, this, true) }));

    // link the sub graphs
    std::vector<PriceLabelGraph>& plgs = priceLabelsGraphs_.back();
    principalGraphs_[0].linkInSubGraph(plgs.back(), k); // MIN -> rest principal
    plgs.back().linkInSubGraph(plgs.front(), k); // MAX -> MIN
  }

  // last day: price max cons days, min/max days and max weekend
  priceLabelsGraphs_.emplace_back(vector<PriceLabelGraph>(
      {PriceLabelGraph(nDays_-1, pContract_->maxConsDaysWork_, maxRotationLength_, MAX_CONS_DAYS, this),
       // use a postprocess function instead
       //       PriceLabelGraph(nDays_-1, pContract_->maxTotalShifts_, maxTotalDuration_, MAX_DAYS, this),
//       PriceLabelGraph(nDays_-1, 0, pContract_->minTotalShifts_, MIN_DAYS, this),
//       PriceLabelGraph(nDays_-1, pContract_->maxTotalWeekends_, pScenario_->nbWeeks(), MAX_WEEKEND, this)
       }));
  // link subgraphs
  SubGraph* previousGraph = nullptr;
  for(PriceLabelGraph& plg: priceLabelsGraphs_.back()) {
    if(previousGraph) previousGraph->linkOutSubGraph(plg, nDays_-1);
    previousGraph = &plg;
  }

  // 4. SINK NODE
  //
  v = addSingleNode(SINK_NODE);
  g_.addSink(v);
}

void RosterSP::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (PrincipalGraph &pg: principalGraphs_) {
    vector2D<int> vec2;
    for (int dest : pg.getDayNodes(0)) {
      std::vector<int> vec;
      for (int s: pScenario_->shiftTypeIDToShiftID_[pg.shiftType()])
        vec.emplace_back(addSingleArc(origin, dest, 0, {}, SOURCE_TO_PRINCIPAL, 0, s));
      vec2.push_back(vec);
    }
    arcsFromSource_[pg.shiftType()] = {vec2};
  }
}

void RosterSP::createArcsAllPriceLabels() {
    for(int sh=0; sh<pScenario_->nbShiftsType_; sh++){		// For all shifts
      // incoming  arc
      int origin = principalGraphs_[sh].exit(nDays_-1);
      // do not take the first pricing graph (max cons) if rest
      int destin = priceLabelsGraphs_[nDays_-1].front().entrance();
      arcsPrincipalToPriceLabelsIn_[sh] =
          {g_.addSingleArc(origin, destin, 0, {0,0,0,0,0}, PRINCIPAL_TO_PRICE_LABEL, nDays_-1)};
    }

    // outgoing  arcs
    int origin = priceLabelsGraphs_[nDays_-1].back().exit();
    int destin = g_.lastSink();
    g_.addSingleArc(origin, destin, 0, {0,0,0,0,0}, PRICE_LABEL_OUT_TO_SINK, nDays_-1);
}