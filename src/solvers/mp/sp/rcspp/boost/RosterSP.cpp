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

#include "RosterSP.h"

#include <string>
#include <algorithm>

#include "RCSPP.h"

using std::string;
using std::vector;
using std::map;

namespace boostRCSPP {

// Constructors and destructor
RosterSP::RosterSP(PScenario scenario,
                   int nDays,
                   PLiveNurse pNurse,
                   SubProblemParam param) :
    SubProblem(scenario, nDays, pNurse, param) {
  labels_ = {CONS_DAYS, DAYS, WEEKEND};
  build();
}

RosterSP::~RosterSP() {}

BoostRCSPPSolver *RosterSP::initRCSSPSolver() {
  Penalties penalties = initPenalties();
  // lambda expression to post process the solutions found by the RCSPP solver
  auto postProcess = [penalties](spp_res_cont *res_cont) {
    res_cont->postprocessCost =
        penalties.penalty(DAYS, res_cont->label_value(DAYS)) +
            penalties.penalty(WEEKEND, res_cont->label_value(WEEKEND));
    res_cont->cost += res_cont->postprocessCost;
  };
  return new BoostRCSPPSolver(&g_,
                              maxReducedCostBound_,
                              param_.verbose_,
                              param_.epsilon_,
                              param_.search_strategy_,
                              param_.rcsppMaxNegativeLabels_,
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
  for (int sh = 0; sh < pScenario_->nShiftTypes(); sh++)
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
      for (int s : pScenario_->shiftTypeIDToShiftID(pg.shiftType()))
        vec.emplace_back(addSingleArc(origin,
                                      dest,
                                      0,
                                      {},
                                      SOURCE_TO_PRINCIPAL,
                                      0,
                                      pScenario_->pShift(s)));
      vec2.push_back(vec);
    }
    arcsFromSource_[pg.shiftType()] = {vec2};
  }
}

// Create all arcs from one principal subgraph to another one
// add a pricing arc when going from work to rest
void RosterSP::createArcsPrincipalToPrincipal() {
  int nShiftsType = pScenario_->nShiftTypes();
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
  for (int sh = 0; sh < pScenario_->nShiftTypes(); sh++) {
    // incoming  arc
    int o = principalGraphs_[sh].exit(nDays_ - 1);
    arcsTosink_.push_back(
        g_.addPricingArc(o, d, 0, cons, TO_SINK, nDays_ - 1,
                         {CONS_DAYS}, penalties));
  }
}

void RosterSP::updateArcDualCosts() {
  // A. ARCS : SOURCE_TO_PRINCIPAL [baseCost = 0]
  // B. ARCS : PRINCIPAL GRAPH
  SubProblem::updateArcDualCosts();

  // C. ARCS : PRINCIPAL GRAPH TO PRINCIPAL GRAPH
  // for rest principal graph
  for (int s = 1; s < pScenario_->nShiftTypes(); s++) {
    for (int a : principalToPrincipal_[0][s])
      if (a != -1) g_.updateCost(a, startWorkCost(a));
    for (int a : principalToPrincipal_[s][0])
      if (a != -1) g_.updateCost(a, endWorkCost(a));
  }
}

double RosterSP::startWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // retrieve the work cost
  // previous does not make any difference for roster duals
  double cost = shiftCost(a, nullptr);
  int start = arc_prop.day;
  // test is true when the arc does not match an assignment;
  // it may be then end of a rotation for instance
  if (arc_prop.pShifts.empty())
    ++start;  // start work the next day as no shift today
  cost += startWeekendCosts_[start];

  return cost;
}

double RosterSP::shiftCost(int a, const PAbstractShift &prevS) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // initial cost of the arc is the part of the cost that is common to all
  // nurses with the same contract
  double cost = arc_prop.initialCost;

  int k = arc_prop.day;
  // iterate through the shift to update the cost
  for (const PShift &pS : arc_prop.pShifts)
    cost += preferencesCosts_[k++][pS->id];
  Stretch st(arc_prop.day, arc_prop.pShifts);
  cost -= pCosts_->getCost(pLiveNurse_->num_, st, prevS);
  return cost;
}

double RosterSP::endWorkCost(int a) const {
  const Arc_Properties &arc_prop  = g_.arc(a);
  // previous does not make any difference for roster duals
  double cost = shiftCost(a, nullptr);
  int length = arc_prop.pShifts.size(), end = arc_prop.day;
  // compute the end of the sequence of shifts
  // length could be equal to 0, but still represents the end of the rotation
  if (length > 1) end += length - 1;
  if (arc_prop.pShifts.empty() || arc_prop.pShifts.back())
    cost += endWeekendCosts_[end];
  return cost;
}

void RosterSP::computeCost(MasterProblem *, RCSolution *rcSol) const {
  /************************************************
   * Compute all the costs of a roster:
   ************************************************/
#ifdef DBG
  double cost = rcSol->cost();
#endif
  /*
   * Compute initial cost
   */
  rcSol->resetCosts();
  // initial state of the nurse
  PShift lastPShift = pLiveNurse_->pStateIni_->pShift_;
  int nbConsShifts = pLiveNurse_->pStateIni_->consShifts_;
  int nbConsDaysWorked = pLiveNurse_->pStateIni_->consDaysWorked_;
  int nbConsDaysOff = pLiveNurse_->pStateIni_->consDaysOff_;

  // 1. if the initial shift has already exceeded the max,
  // fix these values to the UB to not pay it twice
  if (lastPShift->isWork()) {  // was working
    if (nbConsShifts > pScenario_->maxConsShiftsOfType(lastPShift->type))
      nbConsShifts = pScenario_->maxConsShiftsOfType(lastPShift->type);
    if (nbConsDaysWorked > pLiveNurse_->maxConsDaysWork())
      nbConsDaysWorked = pLiveNurse_->maxConsDaysWork();
  } else if (nbConsShifts > pLiveNurse_->maxConsDaysOff()) {
    nbConsDaysOff = pLiveNurse_->maxConsDaysOff();
  }

  // 2. compute the cost
  for (const PShift &pS : rcSol->pShifts()) {
    // a. same shift type -> increment the counters
    if (lastPShift->type == pS->type) {
      if (pS->isWork()) {
        nbConsShifts++;
        nbConsDaysWorked++;
      } else {
        nbConsDaysOff++;
      }
      continue;
    }
    // b. different shift type -> add the corresponding cost
    // i) goes to rest
    if (pS->isRest()) {
      // compute cost
      rcSol->addCost(
          pScenario_->consShiftTypeCost(lastPShift->type, nbConsShifts),
          CONS_SHIFTS_COST);
      rcSol->addCost(
          pLiveNurse_->consDaysCost(nbConsDaysWorked), CONS_WORK_COST);
      // update counters
      nbConsShifts = 0;
      nbConsDaysWorked = 0;
      nbConsDaysOff = 1;
    } else if (lastPShift->isRest()) {
      // ii) goes to work
      // compute cost
      rcSol->addCost(
          pLiveNurse_->consDaysOffCost(nbConsDaysOff), CONS_REST_COST);
      // update counters
      nbConsDaysOff = 0;
      nbConsDaysWorked = 1;
      nbConsShifts = 1;
    } else {
      // iii) continue to work on a different shift
      // compute cost
      rcSol->addCost(
          pScenario_->consShiftTypeCost(lastPShift->type, nbConsShifts),
          CONS_SHIFTS_COST);
      // update counters
      nbConsShifts = 1;
      nbConsDaysWorked++;
    }
    // update
    lastPShift = pS;
  }

  // pay the max for last day
  if (nbConsDaysOff > pLiveNurse_->maxConsDaysOff())
    rcSol->addCost(
        pLiveNurse_->consDaysOffCost(nbConsDaysOff), CONS_REST_COST);
  if (lastPShift->isWork()
      && nbConsShifts > pScenario_->maxConsShiftsOfType(lastPShift->type))
    rcSol->addCost(
        pScenario_->consShiftTypeCost(lastPShift->type, nbConsShifts),
        CONS_SHIFTS_COST);
  if (nbConsDaysWorked > pLiveNurse_->maxConsDaysWork())
    rcSol->addCost(
        pLiveNurse_->consDaysCost(nbConsDaysWorked), CONS_WORK_COST);

  computePreferencesCost(rcSol);

  /*
   * Compute time duration and complete weekend
   */
  int nWeekends = 0;
  int k = 0;
  bool rest = false;
  for (const PShift &pS : rcSol->pShifts()) {
    if (pS->isWork()) {
      if (Tools::isSaturday(k)) {
        nWeekends++;
      } else if (rest && Tools::isSunday(k)) {
        nWeekends++;
        if (pLiveNurse_->needCompleteWeekends())
          rcSol->addCost(pScenario_->weights().WEIGHT_COMPLETE_WEEKEND,
                         COMPLETE_WEEKEND_COST);
      }
      rest = false;
    } else {
      if (pLiveNurse_->needCompleteWeekends() && !rest && Tools::isSunday(k))
        rcSol->addCost(pScenario_->weights().WEIGHT_COMPLETE_WEEKEND,
                       COMPLETE_WEEKEND_COST);
      rest = true;
    }
    k++;
  }

  if (pLiveNurse_->minTotalShifts() - rcSol->duration() > 0)
    rcSol->addCost(pScenario_->weights().WEIGHT_TOTAL_SHIFTS
        * (pLiveNurse_->minTotalShifts() - rcSol->duration()),
        TOTAL_WORK_COST);
  if (rcSol->duration() - pLiveNurse_->maxTotalShifts() > 0)
    rcSol->addCost(pScenario_->weights().WEIGHT_TOTAL_SHIFTS
        * (rcSol->duration() - pLiveNurse_->maxTotalShifts()),
        TOTAL_WORK_COST);
  rcSol->addCost(pLiveNurse_->totalWeekendCost(nWeekends), TOTAL_WEEKEND_COST);

#ifdef DBG
  if (cost < DBL_MAX - 1 && std::abs(cost - rcSol->cost()) > EPSILON) {
    std::cerr << "# " << std::endl;
    std::cerr << "Bad cost: " << rcSol->cost() << " != " << cost
              << std::endl;
    std::cerr << rcSol->costsToString();
    std::cerr << rcSol->toString();
    std::cerr << "# " << std::endl;
    Tools::throwError("boostRCSPP::computeCost does not get the same cost.");
  }
#endif
}

}  // namespace boostRCSPP
