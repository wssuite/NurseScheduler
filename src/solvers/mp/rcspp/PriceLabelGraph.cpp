//
// Created by antoine legrain on 2020-04-11.
//

#include "PriceLabelGraph.h"
#include "SubProblem.h"


PriceLabelGraph::PriceLabelGraph(int lb, int ub, LABEL label, SubProblem* sp,
    bool reset_labels_after_pricing):
    SubGraph(), pSP_(sp), lb_(lb), ub_(ub), label_(label),
    reset_labels_after_pricing_(reset_labels_after_pricing),
    inSubGraph_(nullptr), outSubGraph_(nullptr) {
  if(sp) build();
}

PriceLabelGraph::~PriceLabelGraph() = default;

void PriceLabelGraph::build() {
  // CREATE THE NODES
  // dafault bounds
  std::vector<int> lbs = pSP_->defaultLBs(),
                   ubs = pSP_->defaultUBs();
  // Entrance
  entrance_ = pSP_->addSingleNode(PRICE_LABEL_ENTRANCE, lbs, ubs);
  // Check nodes
  for(int l=lb_; l<=ub_; l++)	{
    std::vector<int> ubs2 = ubs;
    ubs2[label_] = l;
    checkNodes_.emplace_back(pSP_->addSingleNode(PRICE_LABEL, lbs, ubs2));
  }
  // exit day
  exit_ = pSP_->addSingleNode(PRICE_LABEL_EXIT, lbs, ubs);

  // CREATE THE ARCS
  int l = lb_;
  for(int v: checkNodes_) {
    int c = minLabel() ? getLabelCost(ub_-l) : getLabelCost(l);
    std::vector<int> consumptions = {0,0,0,0,0};
    in_arcs_.emplace_back(pSP_->addSingleArc(entrance_, v, c, consumptions, PRICE_LABEL_IN_TO_PRICE_LABEL));
    if(reset_labels_after_pricing_) {
      std::vector<int> ubs = pSP_->defaultUBs();
      if (minLabel()) consumptions[label_] = ubs[label_] - l;
      else consumptions[label_] = -ubs[label_];
    }
    out_arcs_.emplace_back(pSP_->addSingleArc(v, exit_, 0, consumptions , PRICE_LABEL_TO_PRICE_LABEL_OUT));
    ++l;
  }
}

void PriceLabelGraph::updateArcCosts() {

}

int PriceLabelGraph::getLabelCost(int l) const {
  switch (label_) {
    case MAX_CONS_DAYS:
    case MIN_CONS_DAYS:
      return pSP_->contract()->consDaysCost(l);
    case MAX_DAYS:
    case MIN_DAYS:
      return pSP_->contract()->totalShiftCost(l);
    case MAX_WEEKEND:
      return pSP_->contract()->totalWeekendCost(l);
  }
}

void PriceLabelGraph::linkInSubGraph(SubGraph& inSubGraph, int day) {
  inSubGraph_ = &inSubGraph;
  pSP_->addSingleArc(inSubGraph.exit(day), entrance_, 0, {0,0,0,0,0}, NONE_ARC);
}

void PriceLabelGraph::linkOutSubGraph(SubGraph& outSubGraph, int day) {
  outSubGraph_ = &outSubGraph;
  pSP_->addSingleArc(exit_, outSubGraph.entrance(day), 0, {0, 0, 0, 0, 0}, NONE_ARC);
}