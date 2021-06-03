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

#include "AssignmentConstraints.h"

#include <memory>
#include <utility>
#include <vector>

#include "solvers/mp/RotationMP.h"


RosterAssignmentConstraint::RosterAssignmentConstraint(MasterProblem *pMaster):
    ConstraintMP(pMaster, true) {
  build();
}

// update the dual values of the constraints based on the current solution
void RosterAssignmentConstraint::updateDuals() {
  dualValues_ = pMaster_->pModel()->getDuals(assignmentCons_);
}

// update the dual values of the constraints randomly
void RosterAssignmentConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  dualValues_ = Tools::randomDoubleVector(
      pMaster_->nNurses(),
      -10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND,
      10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
}

// return the dual cost of a stretch based on its consumption of
// the constraints
double RosterAssignmentConstraint::getDualCost(
    int nurseNum,
    const Stretch &st,
    const PAbstractShift &prevS) const {
  // if stretch empty, does not count
  if (st.nDays() == 0) return 0;
  if (st.firstDay() == 0 ||
      (st.firstDay() <= pMaster_->nDays() && st.lastDay() >= pMaster_->nDays()))
    return dualValues_[nurseNum];
  return 0;
}

// add a given constraint to the column
void RosterAssignmentConstraint::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Pattern &col) const {
  cons->push_back(assignmentCons_[col.nurseNum()]);
  coeffs->push_back(1);
}

std::string RosterAssignmentConstraint::toString() const {
  std::stringstream buff;
  for (const auto &pN : pMaster_->liveNurses())
    buff << "# Duals for nurse " << pN->num_ << ": "
         << -dualValues_[pN->num_] << std::endl;
  return buff.str();
}

std::string RosterAssignmentConstraint::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  buff << "#   | Assignment: "
       << -getDualCost(nurseNum, st, nullptr) << std::endl;
  return buff.str();
}

void RosterAssignmentConstraint::build() {
  char name[255];
  // build the roster assignment constraint for each nurse
  assignmentCons_.resize(pMaster_->nNurses());
  for (int i = 0; i < pMaster_->nNurses(); i++) {
    snprintf(name, sizeof(name), "feasibilityAssignmentVar_N%d", i);
    MyVar *feasibilityVar;
    pModel()->createPositiveFeasibilityVar(&feasibilityVar, name);
    feasibilityVars_.push_back(feasibilityVar);
    snprintf(name, sizeof(name), "assignmentCons_N%d", i);
    pModel()->createEQConsLinear(
        &assignmentCons_[i], name, 1, {feasibilityVar}, {1});
  }
}


RotationGraphConstraint::RotationGraphConstraint(
    MasterProblem *pMaster,
    vector2D<PBoundedResource> masterRotationGraphResources):
    ConstraintMP(pMaster, true),
    masterRotationGraphResources_(std::move(masterRotationGraphResources)),
    restsPerDay_(pScenario_->nNurses()),
    restVars_(pScenario_->nNurses()),
    initialStateRCSolutions_(pScenario_->nNurses()) {
  build();
}

// update the dual values of the constraints based on the current solution
void RotationGraphConstraint::updateDuals() {
  startWorkDualValues_.clear();
  endWorkDualValues_.clear();
  for (int n=0; n < pMaster_->nNurses(); n++) {
    const PRCGraph &pG = pRotationGraphs_[n];
    vector<double> startDualValues(pMaster_->nDays()),
        endDualValues(pMaster_->nDays());

    // START
    // get dual value associated to the origin node in the graph,
    // i.e. the node on the previous day : either the source or a maxRest node
    startDualValues[0] = pModel()->getDual(
        restCons_[n][pG->pSource(0)->id], true);
    // get dual values associated to the work flow constraints
    // don't take into account the last which is a sink
    for (int k = 1; k < pMaster_->nDays(); ++k)
      startDualValues[k] = pModel()->getDual(
          restCons_[n][maxRestNodes_[n][k - 1]->id], true);
    startWorkDualValues_.push_back(startDualValues);

    // END
    // get dual values associated to the work flow constraints,
    // i.e the firstRest node of the last day of the rotation
    // -1 corresponds to the coefficient
    for (int k = 0; k < pMaster_->nDays(); ++k)
      endDualValues[k] = -pModel()->getDual(
          restCons_[n][firstRestNodes_[n][k]->id], true);
    endWorkDualValues_.push_back(endDualValues);
  }
}

// update the dual values of the constraints randomly
void RotationGraphConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  startWorkDualValues_ = Tools::randomDoubleVector2D(
      pMaster_->nNurses(), pMaster_->nDays(),
      0, 7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
  endWorkDualValues_ = Tools::randomDoubleVector2D(
      pMaster_->nNurses(), pMaster_->nDays(),
      0, 7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
}

// return the dual cost of a stretch based on its consumption of
// the constraints
double RotationGraphConstraint::getDualCost(int nurseNum,
                   const Stretch &st,
                   const PAbstractShift &prevS) const {
  // if prevS is null, consider it's a rest
  bool work = prevS ? prevS->isWork() : false;
  double d = 0;
  int k = st.firstDay();
  if (k == 0) work = false;  // always count first worked day as a start
  for (const PShift &pS : st.pShifts()) {
    // start work
    if (!work && pS->isWork())
      d += startWorkDualValues_[nurseNum][k];
    // end work
    else if (work && pS->isRest())
      d += endWorkDualValues_[nurseNum][k-1];
    work = pS->isWork();
    k++;
  }
  // always count last work day as a end
  if (k >= pMaster_->nDays() && st.pShifts().back()->isWork())
    d += endWorkDualValues_[nurseNum][k-1];
  return d;
}

void RotationGraphConstraint::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Pattern &col) const {
  // add the arc corresponding to the rotation
  const PRCGraph &pG = pRotationGraphs_[col.nurseNum()];
  const PRCNode &pOrigin =
      (col.firstDay() == 0) ? pG->pSource(0) :
      maxRestNodes_[col.nurseNum()][col.firstDay()-1];
  const PRCNode &pTarget = firstRestNodes_[col.nurseNum()][col.lastDay()];
  pG->addSingleArc(pOrigin, pTarget, col, col.cost());
  // add the coefficient for both flow conservation constraints
  cons->push_back(restCons_[col.nurseNum()][pOrigin->id]);
  coeffs->push_back(1);
  cons->push_back(restCons_[col.nurseNum()][pTarget->id]);
  coeffs->push_back(-1);
}

std::string RotationGraphConstraint::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  buff << "#   | Rotation graph    : "
       << -getDualCost(nurseNum, st, nullptr) << std::endl;
  return buff.str();
}

double RotationGraphConstraint::getInitialStateCost(CostType costType) const {
  double cost = 0;
  for (int i = 0; i < pMaster_->nNurses(); i++)
    for (const auto &pArc : pRotationGraphs_[i]->pSource(0)->outArcs)
      if (pArc->stretch.pShifts().front()->isRest()) {
        double v = pModel()->getVarValue(restVars_[i][pArc->id]);
        if (v > pMaster_->epsilon())
          cost += v * initialStateRCSolutions_[i].costByType(costType);
      }
  return cost;
}

/*
 * Rotation constraints : build a rotation graph
 */
void RotationGraphConstraint::build() {
  pRotationGraphs_.clear();
  restVars_.clear();
  restsPerDay_.clear();
  restCons_.clear();
  // build the rotation network for each nurse
  for (const PLiveNurse &pN : pMaster_->liveNurses()) {
    // build a graph of 2 shifts: rest and work
    auto pG = std::make_shared<RCGraph>(pMaster_->nDays(), pMaster_->nShifts());
    createRotationNodes(pG, pN);
    createRotationArcs(pG, pN);
    createRotationArcsVars(pG, pN);
    createRotationNodesCons(pG, pN);
    pRotationGraphs_.push_back(pG);
  }
}

void RotationGraphConstraint::createRotationNodes(
    const PRCGraph &pG, const PLiveNurse &pN) {
  /* Create nodes: one source, nDays free rest and non free nodes */
  const PShift &pRest = pScenario_->pRestShift();
  PRCNode pSource =
      pG->addSingleNode(SOURCE_NODE, -1, pN->pStateIni_->pShift_);
  // the source is present
  vector<PRCNode> firstRestNodes, maxRestNodes;
  for (int k = 0; k < pMaster_->nDays(); k++) {
    NodeType nT = (k == pMaster_->nDays() - 1) ? SINK_NODE : PRINCIPAL_NETWORK;
    firstRestNodes.push_back(pG->addSingleNode(nT, k, pRest));
    maxRestNodes.push_back(pG->addSingleNode(nT, k, pRest));
  }
  firstRestNodes_.push_back(firstRestNodes);
  maxRestNodes_.push_back(maxRestNodes);
}

void RotationGraphConstraint::createRotationArcs(
    const PRCGraph &pG, const PLiveNurse &pN) {
  /* Create rest arcs */
  const PShift &pRest = pScenario_->pRestShift();
  vector<PRCNode> &firstRestNodes = firstRestNodes_[pN->num_],
      &maxRestNodes = maxRestNodes_[pN->num_];
  // if unlimited rest, take the min
  std::pair<int, double> minCons = minConsRest(pN), maxCons = maxConsRest(pN);
  for (int k = 0; k < pMaster_->nDays(); ++k) {
    // create all arcs from first rest node to max rest node
    PRCNode pOrigin = k == 0 ? pG->pSource(0) : firstRestNodes[k - 1];
    std::vector<PShift> pRestShifts;
    for (int l = minCons.first; l <= maxCons.first; l++) {
      if (k+l > pMaster_->nDays()) continue;
      pRestShifts.push_back(pRest);
      PRCNode pTarget = maxRestNodes[k - 1 + l];
      PRCArc pArc =
          pG->addSingleArc(pOrigin, pTarget, Stretch(k, pRestShifts), 0);
      addRotationRestCost(pArc, pN);
    }

    // create arc from max rest node to max rest node if
    // not first arc or not an hard UB
    if (k == 0 || maxCons.second == LARGE_SCORE) continue;
    pOrigin = maxRestNodes[k - 1];
    PRCNode pTarget = maxRestNodes[k];
    pG->addSingleArc(pOrigin, pTarget, Stretch(k, {pRest}), maxCons.second);
  }
}

std::pair<int, double> RotationGraphConstraint::minConsRest(
    const PLiveNurse &pN) {
  for (const auto &pR : masterRotationGraphResources_[pN->num_])
    if (pR->isHard()) {
      // Override LB if define a real hard bound
      if (pR->getLb() > 1) return {pR->getLb(), LARGE_SCORE};
    }
  // otherwise, just use 1 which is always the default LB
  // Do not use a soft LB, as the number of consecutive rest could be lower
  return {1, pScenario_->weights().WEIGHT_CONS_DAYS_OFF};
}

std::pair<int, double> RotationGraphConstraint::maxConsRest(
    const PLiveNurse &pN) {
  int maxR = 1;
  double cost = 0;
  for (const auto &pR : masterRotationGraphResources_[pN->num_]) {
    if (pR->isHard()) {
      // Override UB if define a real bound
      if (pR->getUb() < pMaster_->nDays()) return {pR->getUb(), LARGE_SCORE};
      if (pR->getLb() > maxR) {
        maxR = pR->getLb();
        cost = 0;
      }
    } else {
      // soft UB present
      if (pR->getUb() < pMaster_->nDays() && pR->getUb() > maxR) {
        maxR = pR->getUb();
        cost = pScenario_->weights().WEIGHT_CONS_DAYS_OFF;
      } else if (pR->getLb() > maxR) {
        maxR = pR->getLb();
        cost = 0;
      }
    }
  }
  return {maxR, cost};
}

void RotationGraphConstraint::addRotationRestCost(
    const PRCArc &pArc, const PLiveNurse &pN) {
#ifdef DBG
  // verify that it is indeed a rest arc
  for (const PShift &pS : pArc->stretch.pShifts())
    if (pS->isWork())
      Tools::throwError("addRotationRestCost works only for a rest arc.");
#endif
  // take a copy of the stretch and add a work shift at the end to ensure
  // to price correctly the rest costs LB if ends before the last day
  Stretch st = pArc->stretch;
  if (st.lastDay() < pMaster_->nDays() - 1)
    st.pushBack(Stretch(st.firstDay() + 1, {pScenario_->pAnyWorkShift()}));

  // build rotation
  RCSolution sol(st, 0, DBL_MAX);
  sol.resetCosts();

  // get the vector of resources of the subproblem
  auto pResources = pMaster_->getSPResources(pN);
  // create expander for stretch
  double c = pArc->cost;
  for (const auto &pR : pResources) {
    pR->initialize(*pArc->origin->pAShift, sol, pArc);
    sol.addCost(pArc->cost - c, pMaster_->resourceCostType(pR, pN));
    c = pArc->cost;
  }

  // create initial label
  PRCLabel pL;
  if (st.firstDay() == 0) {
    pL = std::make_shared<RCLabel>(pResources, *pN->pStateIni_);
  } else {
    const PShift &pS = pScenario_->pAnyWorkShift(pN->availableShifts_);
    State state(st.firstDay() - 1,  // day
                0,  // totalTimeWorked
                0,  // totalWeekendsWorked
                pN->minConsDaysWork(),  // consDaysWorked
                pScenario_->minConsShifts(pS->id),  // consShifts
                0,  // consDaysOff
                pS);  // pShift
    pL = std::make_shared<RCLabel>(pResources, state);
  }

  // expand
  c = pL->cost();
  for (const auto &pE : pArc->expanders) {
    ResourceValues &v = pL->getResourceValues(pE->resourceId);
    pE->expand(pL, &v);
    sol.addCost(pL->cost() - c,
        pMaster_->resourceCostType(pResources[pE->resourceId], pN));
    c = pL->cost();
  }

  // add the label cost to the arc
  pArc->addBaseCost(pL->cost());

  if (sol.firstDay() == 0 &&
      initialStateRCSolutions_[pN->num_].firstDay() == -1)
    initialStateRCSolutions_[pN->num_] = sol;
}

void RotationGraphConstraint::createRotationArcsVars(
    const PRCGraph &pG, const PLiveNurse &pN) {
  /* create variables for all arcs */
  char name[255];
  vector<MyVar*> vars;
  vector2D<MyVar*> varsByDay(pMaster_->nDays());
  MyVar *var;
  for (const PRCArc &pArc : pG->pArcs()) {
    // create a rotation from the stretch
    RotationPattern rot(
        RCSolution(pArc->stretch, pArc->cost, DBL_MAX), pN->num_);
    snprintf(name, sizeof(name), "RestingVars_N%d_%d_%d", pN->num_,
             pArc->stretch.firstDay(), pArc->stretch.lastDay());
    pModel()->createPositiveVar(
        &var, name, pArc->cost, rot.getCompactPattern());
    vars.push_back(var);
    for (int k = pArc->stretch.firstDay(); k <= pArc->stretch.lastDay(); k++)
      varsByDay[k].push_back(var);
  }
  restVars_.push_back(vars);
  restsPerDay_.push_back(varsByDay);
}

void RotationGraphConstraint::createRotationNodesCons(
    const PRCGraph &pG, const PLiveNurse &pN) {
  char name[255];
  const vector<MyVar*> restVars = restVars_[pN->num_];
  vector<MyCons*> cons;
  MyCons *con;
  // create a flow conservation constraint for each node
  for (const PRCNode &pNode : pG->pNodes()) {
    vector<MyVar*> vars;
    vector<double> coeffs;
    // flow in
    for (const PRCArc &pA : pNode->inArcs) {
      vars.push_back(restVars[pA->id]);
      coeffs.push_back(-1);
    }
    // flow out
    for (const PRCArc &pA : pNode->outArcs) {
      vars.push_back(restVars[pA->id]);
      coeffs.push_back(1);
    }
    // name
    snprintf(name, sizeof(name), "restNodes_N%d_%d", pN->num_, pNode->id);
    // Create flow constraints.
    // Out flow = 1 if source node, >=-1 if sink nodes, 0 otherwise
    if (pNode->type == SINK_NODE) {
      pModel()->createGEConsLinear(&con, name, -1, vars, coeffs);
    } else {
      int outFlow = 0;
      if (pNode->type == SOURCE_NODE) outFlow = 1;
      pModel()->createEQConsLinear(&con, name, outFlow, vars, coeffs);
    }
    cons.push_back(con);
  }
  restCons_.push_back(cons);
}