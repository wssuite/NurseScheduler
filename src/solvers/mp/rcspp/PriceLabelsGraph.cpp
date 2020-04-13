//
// Created by antoine legrain on 2020-04-11.
//

#include "PriceLabelsGraph.h"
#include "SubProblem.h"


PriceLabelsGraph::PriceLabelsGraph(int min_ub, int max_ub, SubProblem* sp):
    pSP_(sp), min_ub_(min_ub), max_ub_(max_ub) {
  if(sp) build();
}

PriceLabelsGraph::~PriceLabelsGraph() {}

void PriceLabelsGraph::build() {
  // CREATE THE NODES
  // Entrance
  entrance_ = pSP_->addSingleNode(ROTATION_LENGTH_ENTRANCE, {0,0}, {max_ub_, min_ub_});
  // Check nodes
  for(int l=0; l<=max_ub_; l++)	// Check nodes: from CD_max (longest free) to maximum rotation length, for each day
    checkNodes_.emplace_back(pSP_->addSingleNode(ROTATION_LENGTH, {0, 0}, {l, std::max(0, min_ub_-l)}));
  // exit day
  exit_ = pSP_->addSingleNode(ROTATION_LENGTH_EXIT, {0,0}, {max_ub_, min_ub_});

  // CREATE THE ARCS
  int l =  0;
  for(int v: checkNodes_) {
    int c = pSP_->contract()->consDaysCost(l++);
    in_arcs_.emplace_back(pSP_->addSingleArc(entrance_, v, c, {0,0}, ROTSIZEIN_TO_ROTSIZE));
    out_arcs_.emplace_back(pSP_->addSingleArc(v, exit_, 0, {0,0}, ROTSIZE_TO_ROTSIZEOUT));
  }
}

void PriceLabelsGraph::updateArcCosts() {

}