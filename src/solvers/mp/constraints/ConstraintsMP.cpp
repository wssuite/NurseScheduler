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

#include "ConstraintsMP.h"

#include <set>

#include "solvers/mp/MasterProblem.h"


ConstraintMP::ConstraintMP(MasterProblem *pMaster,
                           std::string _name,
                           bool impactColumns,
                           bool addConstraint):
    name(_name), pMaster_(pMaster), pScenario_(pMaster->pScenario()) {
  if (addConstraint) pMaster_->addConstraint(this, impactColumns);
}

Modeler * ConstraintMP::pModel() const { return pMaster_->pModel(); }

NursePositionCountConstraint::NursePositionCountConstraint(
    MasterProblem *pMaster) :
    ConstraintMP(pMaster, "Nurse position", true) {
  build();
}

void NursePositionCountConstraint::build() {
  // initialize vectors
  Tools::initVector3D<MyVar *>(&numberOfNursesByPositionVars_,
                               pMaster_->nDays(),
                               pScenario_->nShifts(),
                               pScenario_->nPositions(),
                               nullptr);
  Tools::initVector3D<MyCons *>(&numberOfNursesByPositionCons_,
                                pScenario_->nPositions(),
                                pMaster_->nDays(),
                                pScenario_->nShifts(),
                                nullptr);

  MyVar *var;
  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int p = 0; p < pScenario_->nPositions(); p++) {
        // creating variable
        char nameV[255];
        snprintf(nameV, sizeof(nameV), "nursesNumber_%d_%d_%d", k, s, p);
        pModel()->createPositiveVar(&var, nameV, 0);
        numberOfNursesByPositionVars_[k][s][p] = var;
        // creating constraint
        char nameC[255];
        snprintf(nameC, sizeof(nameC), "nursesNumberCons_%d_%d_%d", k, s, p);
        pModel()->createEQConsLinear(
            &numberOfNursesByPositionCons_[p][k][s], nameC, 0, {var}, {-1});
      }
    }
  }
}

// update the dual values of the constraints based on the current solution
void NursePositionCountConstraint::updateDuals() {
  dualValues_.clear();
  resetMaxDualValues();
  for (const PLiveNurse &pNurse : pMaster_->pLiveNurses()) {
    vector2D<double> duals(pMaster_->nDays());
    for (int k = 0; k < pMaster_->nDays(); ++k) {
      duals[k] = pModel()->getDuals(
          numberOfNursesByPositionCons_[pNurse->positionId()][k]);
      for (double d : duals[k])
        if (maxDualValues_[pNurse->num_] < d) maxDualValues_[pNurse->num_] = d;
    }
    dualValues_.push_back(duals);
  }
}


// update the dual values of the constraints randomly
void NursePositionCountConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  if (!pScenario_->isWeightsDefined())
    Tools::throwError("Cannot use random duals if weights is not defined.");
  // TODO(JO): below, every dual cost is initialized, so I fear that they
  //  will all be used to modify arcs costs, even in the roster-based
  //  decomposition
  dualValues_.clear();
  for (const PLiveNurse &pNurse : pMaster_->pLiveNurses()) {
    vector2D<double> duals;
    if (useInputData) {
      // This following 2D vector will contain all the dual costs corresponding
      // to each day and to each shift associated with this day.
      Tools::initVector2D(
          &duals, pMaster_->nDays(), pScenario_->nShifts(), 0.0);
      // All the values in the vector are initialized to 0 and only a few of
      // these will be randomly replaced by a value corresponding to the
      // optimal demand weight multiplied by a random coefficient
      // between 1 and 3.
      // The number of days-shifts whose dual cost will be changed is given by
      // the parameter 'nPerturbations'.
      std::set<int> daysToBeModified;
      if (nPerturbations > pMaster_->nDays())
        nPerturbations = pMaster_->nDays();
      while (daysToBeModified.size() < nPerturbations)
        daysToBeModified.insert(Tools::randomInt(0, pMaster_->nDays()-1));
      for (auto day : daysToBeModified) {
        int idShift = Tools::randomInt(1, pScenario_->nShifts()-1);
        int coeff = Tools::randomInt(1, 3);
        duals[day][idShift] = coeff*pScenario_->weights().underCoverage;
      }
    } else {
      duals = Tools::randomDoubleVector2D(
          pMaster_->nDays(), pScenario_->nShifts(),
          0, 3*pScenario_->weights().underCoverage);
      for (auto &dVect : duals) dVect[pScenario_->pRestShift()->id] = 0;
    }
    dualValues_.push_back(duals);
  }
}

// return the dual cost of a stretch based on its consumption of
// the constraints
double NursePositionCountConstraint::getDualCost(
    int nurseNum,
    const Stretch &st,
    const PAbstractShift &prevS) const {
  double d = 0;
  const auto &duals = dualValues_[nurseNum];
  for (int k = st.firstDayId(); k <= st.lastDayId(); k++)
    if (st.pShift(k)->isWork())
      d += duals[k % pMaster_->nDays()][st.pShift(k)->id];
  return d;
}

// add a given constraint to the column
void NursePositionCountConstraint::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Column &col) const {
  int p = pMaster_->pLiveNurses()[col.nurseNum()]->positionId();
  for (int k = col.firstDayId(); k <= col.lastDayId(); ++k) {
    int s = col.shift(k);
    if (s < 0) {
      for (const auto &c : numberOfNursesByPositionCons_[p][k])
        if (c) {
          cons->push_back(c);
          coeffs->push_back(1);
        }
    } else if (pScenario_->isWorkShift(s)) {  // if work
      cons->push_back(numberOfNursesByPositionCons_[p][k][s]);
      coeffs->push_back(1);
    }
  }
}

std::string NursePositionCountConstraint::toString() const {
  std::stringstream buff;
  for (const auto &pN : pMaster_->pLiveNurses()) {
    buff << "# Duals for nurse " << pN->num_ << ":" << std::endl;
    const auto &duals = dualValues_[pN->num_];
    for (int k=0; k < pMaster_->nDays(); ++k) {
      buff << "#   | Work day " << k << ":";
      for (const auto &pS : pScenario_->pShifts()) {
        if (pS->isRest()) continue;
        buff << " (" << pS->name << ", " << -duals[k][pS->id] << ")";
      }
      buff << std::endl;
    }
  }
  return buff.str();
}

std::string NursePositionCountConstraint::toString(
    int nurseNum, const Stretch &st) const {
  std::stringstream buff;
  const auto &duals = dualValues_[nurseNum];
  for (int k = st.firstDayId(); k <= st.lastDayId(); ++k) {
    const PShift &pS = st.pShift(k);
    if (pS->isWork())
      buff << "#   | Work day " << k << ": "
           << -duals[k][pS->id]
           << std::endl;
  }
  return buff.str();
}

AllocationConstraint::AllocationConstraint(MasterProblem *pMaster) :
    ConstraintMP(pMaster, "Skills allocation") {
  // build the constraints
  build();
}

void AllocationConstraint::build() {
  // initialize links between skills and positions
  vector2D<int> positionsPerSkill(pScenario_->nSkills()),
      positionsPerAltSkill(pScenario_->nSkills()),
      allSkillsPerPosition;
  int p = 0;
  for (const auto &position : pScenario_->pPositions()) {
    for (int sk : position->skills())
      positionsPerSkill[sk].push_back(p);
    for (int sk : position->altSkills())
      positionsPerAltSkill[sk].push_back(p);
    allSkillsPerPosition.push_back(position->allSkills());
    p++;
  }

  // initialize vectors
  Tools::initVector4D<MyVar *>(&skillsAllocVars_,
                               pMaster_->nDays(),
                               pScenario_->nShifts(),
                               pScenario_->nSkills(),
                               pScenario_->nPositions(),
                               nullptr);
  Tools::initVector3D<MyCons *>(&feasibleSkillsAllocCons_,
                                pMaster_->nDays(),
                                pScenario_->nShifts(),
                                pScenario_->nPositions(),
                                nullptr);

  const auto &numberOfNursesByPositionVars =
      pMaster_->getNursePositionCountVars();

  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
        if (!pS->hasSkill(sk)) continue;
        for (int p : positionsPerSkill[sk]) {
          char nameV[255];
          snprintf(nameV, sizeof(nameV),
                   "skillsAllocVar_%d_%d_%d_%d", k, s, sk, p);
          pModel()->createPositiveVar(
              &skillsAllocVars_[k][s][sk][p], nameV, 0);
        }
        for (int p : positionsPerAltSkill[sk]) {
          double c = pScenario_->pPosition(p)->skillCost(sk);
          if (isInfeasibleCost(c))
            continue;
          char nameV[255];
          snprintf(nameV, sizeof(nameV),
                   "altSkillsAllocVar_%d_%d_%d_%d", k, s, sk, p);
          pModel()->createPositiveVar(
              &skillsAllocVars_[k][s][sk][p], nameV, c);
        }
      }

      for (int p = 0; p < pScenario_->nPositions(); p++) {
        // adding variables and building skills allocation constraints
        vector<MyVar *> vars = {numberOfNursesByPositionVars[k][s][p]};
        vector<double> coeffs = {1};  // coeff of numberOfNursesByPositionVars
        for (int sk : allSkillsPerPosition[p])
          if (skillsAllocVars_[k][s][sk][p]) {
            vars.push_back(skillsAllocVars_[k][s][sk][p]);
            coeffs.push_back(-1);
          }
        char nameC[255];
        snprintf(nameC, sizeof(nameC),
                 "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
        pModel()->createEQConsLinear(
            &feasibleSkillsAllocCons_[k][s][p], nameC, 0, vars, coeffs);
      }
    }
  }
}

DemandConstraint::DemandConstraint(MasterProblem *pMaster, int demandIndex,
                                   const std::string &name, bool buildAll) :
    ConstraintMP(pMaster, "Demand "+name),
    pDemand_(pMaster->pDemands().at(demandIndex)),
    demandIndex_(demandIndex),
    buildAll_(buildAll),
    demandCons_(pMaster->nDays()) {
  build();
}

void DemandConstraint::update() {
  PDemand pOldDemand = pDemand_;
  pDemand_ = pMaster_->pDemands().at(demandIndex_);
  if (pOldDemand->isGE() && !pDemand_->isGE())
    Tools::throwError("New demand is GE while old one is not");
  if (pOldDemand->isLE() && !pDemand_->isLE())
    Tools::throwError("New demand is LE while old one is not");
  if (pOldDemand->isEQ() && !pDemand_->isEQ())
    Tools::throwError("New demand is EQ while old one is not");

  for (int k = 0; k < pMaster_->nDays(); k++)
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
        double c = pDemand_->cost(k, s, sk);
        bool newSoft = isSoftCost(c),
                oldSoft = isSoftCost(pOldDemand->cost(k, s, sk));
        if (newSoft != oldSoft)
          Tools::throwError("New demand is %s while old is {}.",
                            newSoft ? "soft" : "hard",
                            oldSoft ? "soft" : "hard");
        // update constraint
        demandCons_[k][s][sk]->setLhs(demand()[k][s][sk]);
        MyVar *v = underCovVars_[k][s][sk];
        if (v) v->setCost(c);
        v = overCovVars_[k][s][sk];
        if (v) v->setCost(c);
      }
    }
}


void DemandConstraint::build() {
  // initialize vectors
  Tools::initVector3D<MyVar*>(&underCovVars_, pMaster_->nDays(),
                              pScenario_->nShifts(),
                              pScenario_->nSkills(),
                              nullptr);
  Tools::initVector3D<MyVar*>(&overCovVars_, pMaster_->nDays(),
                              pScenario_->nShifts(),
                              pScenario_->nSkills(),
                              nullptr);
  Tools::initVector3D<MyCons *>(&demandCons_,
                                pMaster_->nDays(),
                                pScenario_->nShifts(),
                                pScenario_->nSkills(),
                                nullptr);
  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++)
        buildCons(k, pS->id, sk);
    }
  }
}

bool DemandConstraint::buildCons(int k, int s, int sk) {
  // check if need a constraint
  int d = demand()[k][s][sk];
  double c = pDemand_->cost(k, s, sk);
  if (!buildAll_ &&
      (c < pMaster_->epsilon() || (pDemand_->isGE() && d == 0)))
    return false;

  // create slack
  MyVar *underVar = nullptr, *overVar = nullptr;
  if (!pDemand_->isLE()) {
    char nameV[255];
    const char *baseV = isSoftCost(c) ? "underCov" : "underFeas";
    snprintf(nameV, sizeof(nameV),
             "%s%sVar_%d_%d_%d", baseV, prefix_.c_str(), k, s, sk);
    if (isSoftCost(c)) {
      pModel()->createPositiveVar(&underVar, nameV, c);
    } else {
      pModel()->createPositiveFeasibilityVar(&underVar, nameV);
    }
    underCovVars_[k][s][sk] = underVar;
  }

  if (!pDemand_->isGE()) {
    char nameV[255];
    const char *baseV = isSoftCost(c) ? "overCov" : "overFeas";
    snprintf(nameV, sizeof(nameV),
             "%s%sVar_%d_%d_%d", baseV, prefix_.c_str(), k, s, sk);
    if (isSoftCost(c)) {
      pModel()->createPositiveVar(&overVar, nameV, c);
    } else {
      pModel()->createPositiveFeasibilityVar(&overVar, nameV);
    }
    overCovVars_[k][s][sk] = overVar;
  }

  // adding variables and building demand constraints
  vector<MyVar *> vars;
  for (MyVar *v : pMaster_->getSkillsAllocVars()[k][s][sk])
    if (v) vars.push_back(v);
  vector<double> coeffs(vars.size(), 1);
  if (underVar) {
    vars.push_back(underVar);
    coeffs.push_back(1);
  }
  if (overVar) {
    vars.push_back(overVar);
    coeffs.push_back(-1);
  }

  char nameC[255];
  snprintf(nameC, sizeof(nameC),
           "%sCons_%d_%d_%d", prefix_.c_str(), k, s, sk);
  if (pDemand_->isLE())
    pModel()->createLEConsLinear(
            &demandCons_[k][s][sk], nameC, d, vars, coeffs);
  else if (pDemand_->isGE())
    pModel()->createGEConsLinear(
            &demandCons_[k][s][sk], nameC, d, vars, coeffs);
  else
    pModel()->createEQConsLinear(
            &demandCons_[k][s][sk], nameC, d, vars, coeffs);
  return true;
}

const vector3D<int>& DemandConstraint::demand() const {
  return pDemand_->demand();
}


