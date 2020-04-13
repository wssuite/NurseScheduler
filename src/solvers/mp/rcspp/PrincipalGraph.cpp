//
// Created by antoine legrain on 2020-04-11.
//

#include "PrincipalGraph.h"
#include "SubProblem.h"

using std::string;
using std::vector;

PrincipalGraph::PrincipalGraph(int shift_type,  SubProblem* sp):
  pSP_(sp), shift_type_(shift_type), max_cons_(0) {
  if(sp) {
    max_cons_ = sp->maxCons(shift_type);
    int i = 0;
    for (int s: sp->scenario()->shiftTypeIDToShiftID_[shift_type_])
      shifts_to_indices_[s] = i++;
    build();
  }
}

PrincipalGraph::~PrincipalGraph() {}

void PrincipalGraph::build() {
  Scenario* pScenario = pSP_->scenario();
  int nDays = pSP_->nDays();

  // CREATE THE NODES
  for(int k=0; k<nDays; k++){		// For each date
    vector<int> vec;
    for(int cons=0; cons<=max_cons_; cons++)		// For each level of network
      vec.emplace_back(pSP_->addSingleNode(PRINCIPAL_NETWORK));
    principalNetworkNodes_.push_back(vec);
  }


  // CREATE THE ARCS
  unsigned int nShifts = pScenario->shiftTypeIDToShiftID_[shift_type_].size();

  // initialization
  Tools::initVector3D(arcsShiftToSameShift_, nDays-1, max_cons_, nShifts, -1);
  Tools::initVector2D(arcsRepeatShift_, nDays-1, nShifts, -1);
  Tools::initVector2D(arcsShiftToEndsequence_, nDays, max_cons_, -1);

  // FOR EACH OF THE DAYS
  //
  int origin, destin;
  for(int k=0; k<nDays-1; k++){

    // 1. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS NOT REACHED YET
    //
    for(int nCons=0; nCons<max_cons_; nCons ++){
      origin = principalNetworkNodes_[k][nCons];
      int dest_day = k+(nCons>0);
      destin = principalNetworkNodes_[dest_day][nCons+1];  // stay on same day for the first level

      for (unsigned int s = 0; s < nShifts; s++) {
        int  shiftID = pScenario-> shiftTypeIDToShiftID_[shift_type_][s];
        int a = pSP_->addSingleArc(origin, destin, 0, {1,-1}, SHIFT_TO_SAMESHIFT, dest_day, {shiftID});
        arcsShiftToSameShift_[k][nCons][s] = a;
      }
    }

    // 2. WORK ONE MORE DAY ON THE SAME SHIFT WHEN MAXIMUM IS ALREADY REACHED
    //
    origin = principalNetworkNodes_[k][max_cons_];
    destin = principalNetworkNodes_[k+1][max_cons_];
    double cost = pSP_->isUnlimited(shift_type_) ? 0 : WEIGHT_CONS_SHIFTS;

    for (unsigned int s = 0; s < nShifts; s++) {
      int  shiftID = pScenario-> shiftTypeIDToShiftID_[shift_type_][s];
      int a = pSP_->addSingleArc(origin, destin, cost, {1,-1}, REPEATSHIFT, k+1, {shiftID});
      arcsRepeatShift_[k][s] = a;
    }

    // 3. END THE CONSECUTIVE SHIFT SEQUENCE
    // cannot end for the initial level: must perform at least one shift
    for(int nCons=1; nCons<max_cons_; nCons++){
      origin = principalNetworkNodes_[k][nCons];
      destin = principalNetworkNodes_[k][max_cons_];
      arcsShiftToEndsequence_[k][nCons] =
          pSP_->addSingleArc(origin, destin, pScenario->consShiftTypeCost(shift_type_, nCons),
              {0,0}, SHIFT_TO_ENDSEQUENCE, k);
    }
  }

  // SPECIAL CASE: LAST DAY
  // be able to work one shift to reach level 1 from level 0
  origin = principalNetworkNodes_[nDays-1][0];
  destin = principalNetworkNodes_[nDays-1][1];
  std::vector<int> last_arcs;
  for (unsigned int s = 0; s < nShifts; s++) {
    int  shiftID = pScenario-> shiftTypeIDToShiftID_[shift_type_][s];
    last_arcs.emplace_back(pSP_->addSingleArc(origin, destin, 0, {1,-1}, SHIFT_TO_SAMESHIFT, nDays-1, {shiftID}));
  }
  arcsShiftToSameShift_.push_back({last_arcs});
  // do not pay any penalty for the minimum
  for(int nCons=1; nCons<max_cons_; nCons++){
    origin = principalNetworkNodes_[nDays-1][nCons];
    destin = principalNetworkNodes_[nDays-1][max_cons_];
    arcsShiftToEndsequence_[nDays-1][nCons] =
        pSP_->addSingleArc(origin, destin, 0, {0,0}, SHIFT_TO_ENDSEQUENCE, nDays-1);
  }
}


// check if feasible to link this arc at this level
bool PrincipalGraph::checkFeasibilityEntranceArc(const Arc_Properties& arc_prop, int level) const {
  if( !checkIfShiftBelongsHere(arc_prop.shifts.back(), true) )
    return false;

  // find which level should be reached
  int sh = -1, n = 0;
  if(arc_prop.day == 0) {
    sh = pSP_->liveNurse()->pStateIni_->shiftType_;
    n = pSP_->liveNurse()->pStateIni_->consShifts_;
    if(arc_prop.shifts.empty()) {
      std::cerr << "Arc must contain at least a shift "
                << "when starting the first to take into account the historical state." << std::endl;
      return false;
    }
  }
  for(int s: arc_prop.shifts) {
    int new_sh = pSP_->scenario()->shiftIDToShiftTypeID_[s];
    if(new_sh == sh) ++n;
    else n = 1;
    sh = new_sh;
  }
  if(n > max_cons_)
    n = max_cons_;

  return n == level;
}


void PrincipalGraph::updateArcCosts() {
  if(!pSP_) return;

  for (int k = 0; k < pSP_->nDays() - 1; k++) {
    for (int n = 0; n < max_cons_; n++)
      for (int a: arcsShiftToSameShift_[k][n])
        pSP_->g().updateCost(a, pSP_->workCost(a));

    for(int a: arcsRepeatShift_[k])
      pSP_->g().updateCost(a, pSP_->workCost(a));
  }

  for (int a: arcsShiftToSameShift_[pSP_->nDays()-1][0])
    pSP_->g().updateCost(a, pSP_->workCost(a));
}

void PrincipalGraph::forbidDayShift(int k, int s) {
  if( !checkIfShiftBelongsHere(s, true) )
    return;

  int i = shifts_to_indices_.at(s);
  pSP_->g().forbidArc(arcsShiftToSameShift_[k][0][i]); // forbid first cons shift
  if(k>0) { // forbid any ingoing arcs using shift s
    for (int n = 1; n < max_cons_; n++) {
//    pSP_->g().forbidNode(principalNetworkNodes_[k][n]);
      pSP_->g().forbidArc(arcsShiftToSameShift_[k-1][n][i]);
    }
    pSP_->g().forbidArc(arcsRepeatShift_[k-1][i]);
  }
}

void PrincipalGraph::authorizeDayShift(int k, int s){
  if( !checkIfShiftBelongsHere(s, true) )
    return;

  int i = shifts_to_indices_.at(s);
  pSP_->g().authorizeArc(arcsShiftToSameShift_[k][0][i]); // authorize first cons shift
  if(k>0) { // authorize any ingoing arcs using shift s
    for (int n = 1; n < max_cons_; n++) {
//    pSP_->g().forbidNode(principalNetworkNodes_[k][n]);
      pSP_->g().authorizeArc(arcsShiftToSameShift_[k-1][n][i]);
    }
    pSP_->g().authorizeArc(arcsRepeatShift_[k-1][i]);
  }
}

bool PrincipalGraph::checkIfShiftBelongsHere(int s, bool print_err) const {
  if(pSP_ == nullptr) return false;
  int sh = pSP_->scenario()->shiftIDToShiftTypeID_[s];
  if(shift_type_ != sh) {
    if(print_err) std::cerr << "Shift " << s << " of type " << sh <<
                      "does not belong to the principal graph associated to type " << shift_type_ << std::endl;
    return false;
  }
  return true;
}