//
// Created by antoine legrain on 2020-04-11.
//

#include "PriceLabelGraph.h"
#include "SubProblem.h"


PriceLabelGraph::PriceLabelGraph(int day, int lb, int ub, LABEL label, SubProblem* sp,
    bool reset_labels_after_pricing):
    SubGraph(), pSP_(sp), day_(day), lb_(lb), ub_(ub), label_(label),
    reset_labels_after_pricing_(reset_labels_after_pricing),
    inSubGraph_(nullptr), outSubGraph_(nullptr) {
  if(lb > ub) {
    lb_ = ub;
    std::cout << "WARNING: price label lb(" << lb << ") is greater than ub("
              << ub << ") for label "  << labelName[label]<< std::endl;
  }
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
  for (int l = lb_; l <= ub_; l++) {
    // WARNING: it's important to set the LB = UB to ensure a valid dominance operator
    // Otherwise, some path could be wrongly dominated leading to segfault in boost.
    std::vector<int> lbs2 = lbs, ubs2 = ubs;
    lbs2[label_] = l;
    ubs2[label_] = l;
    checkNodes_.emplace_back(pSP_->addSingleNode(PRICE_LABEL, lbs2, ubs2));
  }
  // exit day
  if (reset_labels_after_pricing_ && minLabel()) // lb = ub, reach the UB as the label will decrease with consumption
    lbs[label_] = ubs[label_];
  exit_ = pSP_->addSingleNode(PRICE_LABEL_EXIT, lbs, ubs);

  // CREATE THE ARCS
  std::vector<int> in_consumptions = {0, 0, 0, 0, 0}, out_consumptions = {0, 0, 0, 0, 0};
  if (reset_labels_after_pricing_ && !minLabel())
    out_consumptions[label_] = -ubs[label_]; // decrease the label to a negative value -> will be set to 0 by the lb

  int l = lb_;
  for (int v: checkNodes_) {
    int c = minLabel() ? getLabelCost(ub_ - l) : getLabelCost(
        l); // min cost is decreasing with l, and max cost is increasing with l
    in_arcs_.emplace_back(pSP_->addSingleArc(entrance_, v, c, in_consumptions, PRICE_LABEL_IN_TO_PRICE_LABEL, day_));
    out_arcs_.emplace_back(pSP_->addSingleArc(v, exit_, 0, out_consumptions, PRICE_LABEL_TO_PRICE_LABEL_OUT, day_));
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
  pSP_->addSingleArc(inSubGraph.exit(day), entrance_, 0, {0,0,0,0,0}, NONE_ARC, day_);
}

void PriceLabelGraph::linkOutSubGraph(SubGraph& outSubGraph, int day) {
  outSubGraph_ = &outSubGraph;
  pSP_->addSingleArc(exit_, outSubGraph.entrance(day), 0, {0, 0, 0, 0, 0}, NONE_ARC, day_);
}