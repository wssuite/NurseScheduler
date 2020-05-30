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

#include "PrincipalGraph.h"

#include <string>

#include "solvers/mp/sp/SubProblem.h"

using std::string;
using std::vector;

PrincipalGraph::PrincipalGraph(int shift_type, SubProblem *sp) :
    SubGraph(), pSP_(sp), shift_type_(shift_type), max_cons_(-1) {
  if (sp) {
    max_cons_ = sp->maxCons(shift_type);
    int i = 0;
    for (int s : sp->scenario()->shiftTypeIDToShiftID_[shift_type_])
      shifts_to_indices_[s] = i++;
    Tools::initVector<SubGraph *>(&inSubGraphs_, pSP_->nDays(), nullptr);
    Tools::initVector<SubGraph *>(&outSubGraphs_, pSP_->nDays(), nullptr);

    build();
  }
}

PrincipalGraph::~PrincipalGraph() {}

void PrincipalGraph::build() {
  PScenario pScenario = pSP_->scenario();
  int nDays = pSP_->nDays();

  // CREATE THE NODES
  for (int k = 0; k < nDays; k++) {  // For each date
    vector<int> vec;
    // For each level of network
    for (int cons = 0; cons <= max_cons_; cons++)
      vec.emplace_back(pSP_->addSingleNode(PRINCIPAL_NETWORK));
    principalNetworkNodes_.push_back(vec);
  }


  // CREATE THE ARCS
  unsigned int nShifts = pScenario->shiftTypeIDToShiftID_[shift_type_].size();

  // initialization
  Tools::initVector3D(&arcsShiftToSameShift_,
                      nDays - 1,
                      max_cons_,
                      nShifts,
                      -1);
  Tools::initVector2D(&arcsRepeatShift_, nDays - 1, nShifts, -1);
  Tools::initVector2D(&arcsShiftToEndsequence_, nDays, max_cons_, -1);

  // FOR EACH OF THE DAYS
  int origin, destin;
  double cost = consCost(max_cons_ + 1);
  for (int k = 0; k < nDays - 1; k++) {
    // 1. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS NOT REACHED YET
    for (int nCons = 0; nCons < max_cons_; nCons++) {
      origin = principalNetworkNodes_[k][nCons];
      destin = principalNetworkNodes_[k + 1][nCons + 1];
      for (unsigned int s = 0; s < nShifts; s++) {
        int shiftID = pScenario->shiftTypeIDToShiftID_[shift_type_][s];
        int a = pSP_->addSingleArc(origin,
                                   destin,
                                   0,
                                   getConsumption(k + 1, shiftID),
                                   SHIFT_TO_SAMESHIFT,
                                   k + 1,
                                   shiftID);
        arcsShiftToSameShift_[k][nCons][s] = a;
      }
    }

    // 2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
    origin = principalNetworkNodes_[k][max_cons_];
    destin = principalNetworkNodes_[k + 1][max_cons_];
    for (unsigned int s = 0; s < nShifts; s++) {
      int shiftID = pScenario->shiftTypeIDToShiftID_[shift_type_][s];
      int a = pSP_->addSingleArc(origin,
                                 destin,
                                 cost,
                                 getConsumption(k + 1, shiftID),
                                 REPEATSHIFT,
                                 k + 1,
                                 shiftID);
      arcsRepeatShift_[k][s] = a;
    }

    // 3. END THE CONSECUTIVE SHIFT SEQUENCE
    // cannot end for the initial level: must perform at least one shift
    for (int nCons = 1; nCons < max_cons_; nCons++) {
      origin = principalNetworkNodes_[k][nCons];
      destin = principalNetworkNodes_[k][max_cons_];
      arcsShiftToEndsequence_[k][nCons] =
          pSP_->addSingleArc(origin, destin, consCost(nCons),
                             getConsumption(-1, -1), SHIFT_TO_ENDSEQUENCE, k);
    }
  }

  // SPECIAL CASE: LAST DAY
  // do not pay any penalty for the minimum
  for (int nCons = 1; nCons < max_cons_; nCons++) {
    origin = principalNetworkNodes_[nDays - 1][nCons];
    destin = principalNetworkNodes_[nDays - 1][max_cons_];
    arcsShiftToEndsequence_[nDays - 1][nCons] =
        pSP_->addSingleArc(origin, destin, 0, getConsumption(-1, -1),
                           SHIFT_TO_ENDSEQUENCE, nDays - 1);
  }
}

// return the right vector of consumption based on the day
// (if < 0, not performing any shift)
std::vector<int> PrincipalGraph::getConsumption(int day, int shift) const {
  // not performing any shift or the shift type is rest
  if (day < 0 || shift_type_ == 0) return {0, 0, 0, 0, 0};

  // otherwise work => consume one resource of each if needed
  int t = pSP_->scenario()->timeDurationToWork_[shift];
  return {1, -1, t, -t, Tools::isSaturday(day)};
}

// check if feasible to link this arc at this level
bool PrincipalGraph::checkFeasibilityEntranceArc(
    const Arc_Properties &arc_prop,
    int level) const {
  if (!checkIfShiftBelongsHere(arc_prop.shifts.back(), true))
    return false;

  // find which level should be reached
  int sh = -1, n = 0;
  if (arc_prop.day == 0) {
    sh = pSP_->liveNurse()->pStateIni_->shiftType_;
    n = pSP_->liveNurse()->pStateIni_->consShifts_;
    if (arc_prop.shifts.empty()) {
      std::cerr << "Arc must contain at least a shift when starting the first "
                   "to take into account the historical state." << std::endl;
      return false;
    }
  }

  for (int s : arc_prop.shifts) {
    int new_sh = pSP_->scenario()->shiftIDToShiftTypeID_[s];
    if (new_sh == sh) ++n;
    else
      n = 1;
    sh = new_sh;
  }
  if (n > max_cons_)
    n = max_cons_;

  return n == level;
}

void PrincipalGraph::updateArcCosts() {
  if (!pSP_ || shift_type_ == 0) return;

  for (int k = 0; k < pSP_->nDays() - 1; k++) {
    for (int n = 0; n < max_cons_; n++)
      for (int a : arcsShiftToSameShift_[k][n])
        pSP_->g().updateCost(a, pSP_->shiftCost(a));

    for (int a : arcsRepeatShift_[k])
      pSP_->g().updateCost(a, pSP_->shiftCost(a));
  }
}

double PrincipalGraph::consCost(int n) const {
  if (shift_type_ == 0)
    return pSP_->contract()->consDaysOffCost(n);
  return pSP_->scenario()->consShiftTypeCost(shift_type_, n);
}

void PrincipalGraph::linkInSubGraph(SubGraph *inSubGraph, int day) {
  inSubGraphs_[day] = inSubGraph;
  pSP_->addSingleArc(inSubGraph->exit(day),
                     getNode(day, 0),
                     0,
                     {0, 0, 0, 0, 0},
                     NONE_ARC);
}

void PrincipalGraph::linkOutSubGraph(SubGraph *outSubGraph, int day) {
  outSubGraphs_[day] = outSubGraph;
  pSP_->addSingleArc(getNode(day, max_cons_),
                     outSubGraph->entrance(day),
                     0,
                     {0, 0, 0, 0, 0},
                     NONE_ARC);
}

void PrincipalGraph::forbidDayShift(int k, int s) {
  if (!checkIfShiftBelongsHere(s, true))
    return;

  int i = shifts_to_indices_.at(s);
  if (k-- > 0) {  // if k>0, continue and do --k
    // forbid any ingoing arcs using shift s
    for (int n = 0; n < max_cons_; n++) {
//    pSP_->g().forbidNode(principalNetworkNodes_[k][n]);
      pSP_->g().forbidArc(arcsShiftToSameShift_[k][n][i]);
    }
    pSP_->g().forbidArc(arcsRepeatShift_[k][i]);
  }
}

void PrincipalGraph::authorizeDayShift(int k, int s) {
  if (!checkIfShiftBelongsHere(s, true))
    return;

  int i = shifts_to_indices_.at(s);
  if (k > 0) {  // authorize any ingoing arcs using shift s
    for (int n = 0; n < max_cons_; n++) {
//    pSP_->g().forbidNode(principalNetworkNodes_[k][n]);
      pSP_->g().authorizeArc(arcsShiftToSameShift_[k - 1][n][i]);
    }
    pSP_->g().authorizeArc(arcsRepeatShift_[k - 1][i]);
  }
}

// forbid any arc that authorizes the violation of a consecutive constraint
void PrincipalGraph::forbidViolationConsecutiveConstraints() {
  if (!pSP_)
    return;

  // forbid the arcs that allow to violate the min consecutive constraint
  int minCons = pSP_->minCons(shift_type_);
  for (auto &endArcs : arcsShiftToEndsequence_)
    for (int n = 1; n < minCons; n++)
      pSP_->g().forbidArc(endArcs[n]);

  // forbid the arcs that allow to violate the max consecutive constraint if
  // not unlimited
  if (!pSP_->isUnlimited(shift_type_))
    for (auto &repeatArc : arcsRepeatShift_)
      for (int a : repeatArc)
        pSP_->g().forbidArc(a);
}

bool PrincipalGraph::checkIfShiftBelongsHere(int s, bool print_err) const {
  if (pSP_ == nullptr) return false;
  int sh = pSP_->scenario()->shiftIDToShiftTypeID_[s];
  if (shift_type_ != sh) {
    if (print_err)
      std::cerr << "Shift " << s << " of type " << sh <<
                "does not belong to the principal graph associated to type "
                << shift_type_ << std::endl;
    return false;
  }
  return true;
}
