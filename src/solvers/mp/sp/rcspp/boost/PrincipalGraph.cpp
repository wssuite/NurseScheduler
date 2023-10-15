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

#include <memory>
#include <string>

#include "SubProblem.h"

using std::string;
using std::vector;

namespace boostRCSPP {

PrincipalGraph::PrincipalGraph(int shift_type, SubProblem *sp) :
    SubGraph(), pSP_(sp), shift_type_(shift_type), max_cons_(-1),
    pAS_(sp == nullptr ? nullptr :
         sp->pScenario()->pAnyTypeShift(shift_type_)) {
  if (sp) {
    max_cons_ = sp->maxCons(shift_type);
    int i = 0;
    for (const PShift &pS : sp->pScenario()->shiftsFactory().
         pAnyTypeShift(shift_type_)->pIncludedShifts())
      shifts_to_indices_[pS->id] = i++;
    Tools::initVector<SubGraph *>(&inSubGraphs_, pSP_->nDays(), nullptr);
    Tools::initVector<SubGraph *>(&outSubGraphs_, pSP_->nDays(), nullptr);

    build();
  }
}

PrincipalGraph::~PrincipalGraph() {}

void PrincipalGraph::build() {
  PScenario pScenario = pSP_->pScenario();
  int nDays = pSP_->nDays();

  // CREATE THE NODES
  for (int k = 0; k < nDays; k++) {  // For each date
    vector<int> vec;
    // For each level of network
    for (int cons = 0; cons <= max_cons_; cons++)
      vec.push_back(pSP_->addSingleNode(PRINCIPAL_NETWORK));
    principalNetworkNodes_.push_back(vec);
  }


  // CREATE THE ARCS
  unsigned int nShifts = pScenario->pShiftsOfType(shift_type_).size();

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
      int s = 0;
      for (const PShift &pS : pScenario->pShiftsOfType(shift_type_)) {
        int a = pSP_->addSingleArc(origin,
                                   destin,
                                   0,
                                   getConsumption(k + 1, pS->id),
                                   SHIFT_TO_SAMESHIFT,
                                   k + 1,
                                   pS);
        arcsShiftToSameShift_[k][nCons][s++] = a;
      }
    }

    // 2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
    origin = principalNetworkNodes_[k][max_cons_];
    destin = principalNetworkNodes_[k + 1][max_cons_];
    int s = 0;
    for (const PShift &pS : pScenario->pShiftsOfType(shift_type_)) {
      int a = pSP_->addSingleArc(origin,
                                 destin,
                                 cost,
                                 getConsumption(k + 1, pS->id),
                                 REPEATSHIFT,
                                 k + 1,
                                 pS);
      arcsRepeatShift_[k][s++] = a;
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
  if (day < 0 || shift_type_ == 0) return {0, 0, 0};

  // otherwise work => consume one resource of each if needed
  int t = pSP_->pScenario()->duration(shift);
  return {1, t, pSP_->pLiveNurse()->pContract()->isFirstWeekendDay(day)};
}

// check if feasible to link this arc at this level
bool PrincipalGraph::checkFeasibilityEntranceArc(
    const Arc_Properties &arc_prop,
    int level) const {
  if (!checkIfShiftBelongsHere(arc_prop.pShifts.back()->id, true))
    return false;

  // find which level should be reached
  int sh = -1, n = 0;
  if (arc_prop.day == 0) {
    sh = pSP_->pLiveNurse()->pStateIni_->pShift_->type;
    n = pSP_->pLiveNurse()->pStateIni_->consShifts_;
    if (arc_prop.pShifts.empty()) {
      std::cerr << "Arc must contain at least a shift when starting the first "
                   "to take into account the historical state." << std::endl;
      return false;
    }
  }

  for (const PShift &pS : arc_prop.pShifts) {
    if (pS->type == sh) ++n;
    else
      n = 1;
    sh = pS->type;
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
        pSP_->g().updateCost(a, pSP_->shiftCost(a, pAS_));

    for (int a : arcsRepeatShift_[k])
      pSP_->g().updateCost(a, pSP_->shiftCost(a, pAS_));
  }
}

double PrincipalGraph::consCost(int n) const {
  if (shift_type_ == 0)
    return pSP_->pContract()->consDaysOffCost(n);
  return pSP_->pScenario()->consShiftTypeCost(shift_type_, n);
}

void PrincipalGraph::linkInSubGraph(SubGraph *inSubGraph, int day) {
  checkInitialization();
  inSubGraphs_[day] = inSubGraph;
  pSP_->addSingleArc(inSubGraph->exit(day),
                     getNode(day, 0),
                     0,
                     {0, 0, 0},
                     NONE_ARC);
}

void PrincipalGraph::linkOutSubGraph(SubGraph *outSubGraph, int day) {
  checkInitialization();
  outSubGraphs_[day] = outSubGraph;
  pSP_->addSingleArc(getNode(day, max_cons_),
                     outSubGraph->entrance(day),
                     0,
                     {0, 0, 0},
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
vector<int> PrincipalGraph::forbidViolationConsecutiveConstraints() {
  if (pSP_ == nullptr) return {};

  // forbid the arcs that allow to violate the min consecutive constraint
  vector<int> forbiddenArcs;
  int minCons = pSP_->minCons(shift_type_);
  for (auto &endArcs : arcsShiftToEndsequence_)
    for (int n = 1; n < minCons; n++)
      if (pSP_->g().forbidArc(endArcs[n]))
        forbiddenArcs.push_back(endArcs[n]);

  // forbid the arcs that allow to violate the max consecutive constraint if
  // not unlimited
  if (!pSP_->isUnlimited(shift_type_))
    for (auto &repeatArc : arcsRepeatShift_)
      for (int a : repeatArc)
        if (pSP_->g().forbidArc(a))
          forbiddenArcs.push_back(a);
  return forbiddenArcs;
}

bool PrincipalGraph::checkIfShiftBelongsHere(int s, bool print_err) const {
  if (pSP_ == nullptr) return false;
  int sh = pSP_->pScenario()->shiftIDToShiftTypeID(s);
  if (shift_type_ != sh) {
    if (print_err)
      std::cerr << "Shift " << s << " of type " << sh <<
                "does not belong to the principal graph associated to type "
                << shift_type_ << std::endl;
    return false;
  }
  return true;
}

}  // namespace boostRCSPP