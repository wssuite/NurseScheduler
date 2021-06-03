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


ConstraintMP::ConstraintMP(MasterProblem *pMaster, bool impactColumns):
    pMaster_(pMaster),  pScenario_(pMaster->pScenario()) {
  if (impactColumns) pMaster_->addColumnConstraint(this);
}

Modeler * ConstraintMP::pModel() const { return pMaster_->pModel(); }

NursePositionCountConstraint::NursePositionCountConstraint(
    MasterProblem *pMaster) :
    ConstraintMP(pMaster) {
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

  char name[255];
  MyVar *var;
  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int p = 0; p < pScenario_->nPositions(); p++) {
        // creating variable
        snprintf(name, sizeof(name), "nursesNumber_%d_%d_%d", k, s, p);
        pModel()->createPositiveVar(&var, name, 0);
        numberOfNursesByPositionVars_[k][s][p] = var;
        // creating constraint
        snprintf(name, sizeof(name), "nursesNumberCons_%d_%d_%d", k, s, p);
        pModel()->createEQConsLinear(
            &numberOfNursesByPositionCons_[p][k][s], name, 0, {var}, {-1});
      }
    }
  }
}

// update the dual values of the constraints based on the current solution
void NursePositionCountConstraint::updateDuals() {
  dualValues_.clear();
  for (const PLiveNurse &pNurse : pMaster_->liveNurses()) {
    vector2D<double> duals(pMaster_->nDays());
    for (int k = 0; k < pMaster_->nDays(); ++k)
      duals[k] = pModel()->getDuals(
          numberOfNursesByPositionCons_[pNurse->pPosition_->id_][k]);
    dualValues_.push_back(duals);
  }
}


// update the dual values of the constraints randomly
void NursePositionCountConstraint::randomUpdateDuals(
    bool useInputData, int nPerturbations) {
  // TODO(JO): below, every dual cost is initialized, so I fear that they
  //  will all be used to modify arcs costs, even in the roster-based
  //  decomposition
  dualValues_.clear();
  for (const PLiveNurse &pNurse : pMaster_->liveNurses()) {
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
        duals[day][idShift] = coeff*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND;
      }
    } else {
      duals = Tools::randomDoubleVector2D(
          pMaster_->nDays(), pScenario_->nShifts(),
          0, 3*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND);
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
  for (int k = st.firstDay(); k <= st.lastDay(); k++)
    d += duals[k % pMaster_->nDays()][st.pShift(k)->id];
  return d;
}

// add a given constraint to the column
void NursePositionCountConstraint::addConsToCol(
    std::vector<MyCons *> *cons,
    std::vector<double> *coeffs,
    const Pattern &col) const {
  int p = pMaster_->liveNurses()[col.nurseNum()]->pPosition()->id_;
  for (int k = col.firstDay(); k <= col.lastDay(); ++k) {
    int s = col.shift(k);
    if (pScenario_->isAnyShift(s)) {
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
  for (const auto &pN : pMaster_->liveNurses()) {
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
  for (int k = st.firstDay(); k <= st.lastDay(); ++k) {
    const PShift &pS = st.pShift(k);
    if (pS->isWork())
      buff << "#   | Work day " << k << ": "
           << -duals[k][pS->id]
           << std::endl;
  }
  return buff.str();
}

AllocationConstraint::AllocationConstraint(
    MasterProblem *pMaster) :
    ConstraintMP(pMaster),
    positionsPerSkill_(pScenario_->nSkills()) {
  // initialize links between skills and positions
  int p = 0;
  for (const auto &position : pScenario_->pPositions()) {
    for (int sk : position->skills_)
      positionsPerSkill_[sk].push_back(p);
    skillsPerPosition_.push_back(position->skills_);
    p++;
  }
  // build the constraints
  build();
}

void AllocationConstraint::build() {
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

  char name[255];
  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++)
        for (int p : positionsPerSkill_[sk]) {
          snprintf(name, sizeof(name),
                   "skillsAllocVar_%d_%d_%d_%d", k, s, sk, p);
          pModel()->createPositiveVar(
              &skillsAllocVars_[k][s][sk][p], name, 0);
        }

      for (int p = 0; p < pScenario_->nPositions(); p++) {
        // adding variables and building skills allocation constraints
        vector<MyVar *> vars = {numberOfNursesByPositionVars[k][s][p]};
        for (int sk : skillsPerPosition_[p])
          vars.push_back(skillsAllocVars_[k][s][sk][p]);
        vector<double> coeffs(vars.size(), -1);
        coeffs[0] = 1;  // coeff of numberOfNursesByPositionVars
        snprintf(name, sizeof(name),
                 "feasibleSkillsAllocCons_%d_%d_%d", k, s, p);
        pModel()->createEQConsLinear(
            &feasibleSkillsAllocCons_[k][s][p], name, 0, vars, coeffs);
      }
    }
  }
}

DemandConstraint::DemandConstraint(
    MasterProblem *pMaster, bool minDemand, bool soft, double weight) :
    ConstraintMP(pMaster), minDemand_(minDemand),
    name_(minDemand ? "minDemand" : "optDemand"), soft_(soft), weight_(weight),
    demandCons_(pMaster->nDays()) {
  if (soft_) slackVars_.resize(pMaster->nDays());
  build();
}

void DemandConstraint::updateDemand() {
  for (int k = 0; k < pMaster_->nDays(); k++)
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++)
        demandCons_[k][s][sk]->setLhs(demand()[k][s][sk]);
    }
}


void DemandConstraint::build() {
  // initialize vectors
  Tools::initVector3D<MyVar*>(&slackVars_, pMaster_->nDays(),
                              pScenario_->nShifts(),
                              pScenario_->nSkills(),
                              nullptr);
  Tools::initVector3D<MyCons *>(&demandCons_,
                                pMaster_->nDays(),
                                pScenario_->nShifts(),
                                pScenario_->nSkills(),
                                nullptr);

  const auto &skillsAllocVars = pMaster_->getSkillsAllocVars();

  char name[255];
  for (int k = 0; k < pMaster_->nDays(); k++) {
    for (const PShift &pS : pScenario_->pShifts()) {
      // forget resting shift
      if (pS->isRest()) continue;
      int s = pS->id;
      for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
        // create slack
        MyVar *var;
        if (soft_) {
          snprintf(name, sizeof(name),
                   "slack%sVar_%d_%d_%d", name_.c_str(), k, s, sk);
          pModel()->createPositiveVar(&var, name, weight_);
        } else {
          snprintf(name, sizeof(name),
                   "feasibility%sVar_%d_%d_%d", name_.c_str(), k, s, sk);
          pModel()->createPositiveFeasibilityVar(&var, name);
        }
        slackVars_[k][s][sk] = var;

        // adding variables and building demand constraints
        vector<MyVar *> vars = {var};
        for (MyVar *v : skillsAllocVars[k][s][sk])
          if (v) vars.push_back(v);
        vector<double> coeffs(vars.size(), 1);

        snprintf(name, sizeof(name),
                 "%sCons_%d_%d_%d", name_.c_str(), k, s, sk);
        pModel()->createGEConsLinear(&demandCons_[k][s][sk],
                                     name,
                                     demand()[k][s][sk],
                                     vars,
                                     coeffs);
      }
    }
  }
}

const vector3D<int>& DemandConstraint::demand() const {
  if (minDemand_)
    return pMaster_->pDemand()->minDemand_;
  return pMaster_->pDemand()->optDemand_;
}


