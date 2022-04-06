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

#include "RotationMP.h"

#include <memory>
#include <algorithm>
#include <utility>

#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"

using std::vector;
using std::map;
using std::pair;
using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------

// when branching on this pattern, this method add the corresponding forbidden
// shifts to the set.
// It will forbid all the shifts that would be worked on a day that is already
// covered by this pattern.
// Moreover, there needs to be a resting day before and after each rotation,
// so the shifts can also be forbidden on these two days
// (if the rotation is not at an extremity of the horizon).
void RotationPattern::addForbiddenShifts(
    std::set<std::pair<int, int> > *forbiddenShifts,
    int nbShifts, PDemand pDemand) const {
  // from the previous day to the day after the end of the rotation,
  // forbid any work shifts
  for (int day = firstDay()-1; day <= lastDay()+1; day++) {
    if (day < pDemand->firstDay_) continue;
    if (day >= pDemand->firstDay_ + pDemand->nDays_) continue;
    for (int i = 1; i < nbShifts; ++i)
      forbiddenShifts->insert(pair<int, int>(day, i));
  }
}

void RotationPattern::checkReducedCost(const DualCosts &dualCosts,
                                       bool printBadPricing) const {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of a rotation:
   ************************************************/
  // add a rest shift at the end to mark the end of the rotation
  Stretch st(*this);
  PShift pRest = std::make_shared<Shift>();
  st.pushBack(pRest);
  double reducedCost = cost_ - dualCosts.getCost(nurseNum_, st, pRest);

  // Display: set to true if you want to display the details of the cost
  if (std::fabs(reducedCost_ - reducedCost) / (1 - reducedCost) > EPSILON) {
    // if do not print and not throwing an error
    if (!printBadPricing && reducedCost_ > reducedCost + EPSILON) return;

    cout << "# " << endl;
    cout << "# " << endl;
    cout << "# Bad dual cost: "
         << reducedCost_ << " != " << reducedCost << endl;
    cout << "# " << endl;
    cout << "#   | Base cost     : + " << cost_ << endl;
    cout << costsToString();
    cout << dualCosts.toString(nurseNum_, st);
    std::cout << toString() << "# " << endl;

    // throw an error only when a significant misprice
    // Indeed, if the real reduced cost (reducedCost) is greater than the one
    // found by the pricing, some columns with a positive reduced cost could be
    // generated.
    // The reason why the other situation can arise is that some path in the
    // subproblem could under estimate the real cost. These paths won't be
    // found when the subproblems are solved at optimality, but could  be
    // present when using heuristics.
    if (reducedCost_ < reducedCost + EPSILON) {
      dualCosts.getCost(nurseNum_, st, pRest);
      Tools::throwError("Invalid pricing of a rotation.");
    }
  }
}


//-----------------------------------------------------------------------------
//
//  C l a s s   R o t a t i o n M P
//
// Build and solve the master problem of the column generation scheme with
// rotations
//
//-----------------------------------------------------------------------------

RotationMP::RotationMP(const PScenario& pScenario,
                       SolverType solver) :
    MasterProblem(pScenario, solver) {}

RotationMP::~RotationMP() {
  delete rotationGraphConstraint_;
  delete dynamicConstraints_;
  delete totalShiftDurationConstraint_;
  delete totalWeekendConstraint_;
}

PPattern RotationMP::getPattern(MyVar *var) const {
  return std::make_shared<RotationPattern>(var, pScenario_);
}

// build the, possibly fractional, roster corresponding to the solution
// currently stored in the model
vector3D<double> RotationMP::fractionalRoster() const {
  vector3D<double> roster = MasterProblem::fractionalRoster();

  // add rest values
  for (auto &nurseRoster : roster)
    for (auto &dayRoster : nurseRoster) {
      // fill the rest shift
      double restValue = 1;
      for (int s = 1; s < nShifts(); ++s)
        restValue -= dayRoster[s];
      dayRoster[0] = restValue;
    }

  return roster;
}

// build the rostering problem
void RotationMP::build(const SolverParam &param) {
  /* build the base of the model */
  MasterProblem::build(param);

  /* Rotation constraints */
  rotationGraphConstraint_ =
      new RotationGraphConstraint(this, masterRotationGraphResources_);

  /* Min/Max constraints */
  buildResourceCons(param);

  /* Dynamic constraints */
  dynamicConstraints_ = new DynamicConstraints(
      this, totalShiftDurationConstraint_, totalWeekendConstraint_);
}

// Build the columns corresponding to the initial solution
void RotationMP::initializeSolution(const vector<Roster> &solution) {
  if (solution.empty()) return;

  // rotations are added for each nurse of the initial solution
  string baseName("initialRotation");
  // build the rotations of each nurse
  for (int i = 0; i < pScenario_->nNurses(); ++i) {
    // load the roster of nurse i
    Roster roster = solution[i];

    int firstDay;
    bool workedOnPreviousDay = false;
    std::vector<PShift> pShifts;
    // build all the successive rotation of this nurse
    for (int k = 0; k < pDemand_->nDays_; ++k) {
      // shift=0 => rest
      const PShift &pShift = roster.pShift(k);
      // if work, insert the shift in the map
      if (pShift->isWork()) {
        if (!workedOnPreviousDay) firstDay = k;
        pShifts.push_back(pShift);
        workedOnPreviousDay = true;
      } else if (workedOnPreviousDay) {
        // if stop to work, build the rotation
        RotationPattern rotation(firstDay, pShifts, i);
        computePatternCost(&rotation);
        pModel_->addInitialColumn(createColumn(rotation, baseName.c_str()));
        pShifts.clear();
        workedOnPreviousDay = false;
      }
    }
    // if work on the last day, build the rotation
    if (workedOnPreviousDay) {
      RotationPattern rotation(firstDay, pShifts, i);
      computePatternCost(&rotation);
      pModel_->addInitialColumn(createColumn(rotation, baseName.c_str()));
    }
  }
}

// separate the resources between the ones that will be managed by
// the master, the master rotation graph, and the sub problems
// must initialize spResources_
void RotationMP::splitPResources() {
  // put all the consecutive resources in the sub problem
  // all the consecutive rests are also in the rotation graph
  // the others are master constraints
  spResources_.clear();
  masterRotationGraphResources_.clear();
  masterConstraintResources_.clear();
  for (const auto &m : pResources_) {
    vector<PResource> pSPResources;
    vector<PBoundedResource> pRotPResources, pMPResources;
    for (const auto &p : m) {
      // cast to shared_ptr<BoundedResource>
      PBoundedResource pBR =
          std::dynamic_pointer_cast<BoundedResource>(p.first);
      if (pBR == nullptr)
        Tools::throwError("RotationMP cannot handle a resource "
                          "that does not derive from BoundedResource.");
      // cast hard constraint
      if (pBR->isHard()) {
        auto pHR = std::dynamic_pointer_cast<HardConsShiftResource>(pBR);
        // if hard cons -> add to sub problem and rotation if rest
        if (pHR) {
          pHR->setId(pSPResources.size());
          pSPResources.push_back(pHR);
          if (pHR->pShift()->isRest())
            pRotPResources.push_back(pHR);
        } else {
          // add to master
          pBR->setId(pMPResources.size());
          pMPResources.push_back(pBR);
        }
      } else {
        auto pSR = std::dynamic_pointer_cast<SoftConsShiftResource>(pBR);
        // if soft cons -> add to sub problem and rotation if rest
        if (pSR) {
          pSR->setId(pSPResources.size());
          pSPResources.push_back(pSR);
          if (pSR->pShift()->isRest())
            pRotPResources.push_back(pSR);
        } else {
          // add to master
          pBR->setId(pMPResources.size());
          pMPResources.push_back(pBR);
        }
      }
    }

    // masterRotationGraphResources_ can contains at max one soft and one hard
    // constraint on Rest
    if (pRotPResources.size() > 2)
      Tools::throwError("Cannot have more than two consecutive constraints "
                        "on rest shifts.");
    if (pRotPResources.size() == 2)
      if (pRotPResources.front()->isHard() ^ pRotPResources.back()->isHard())
        Tools::throwError("Cannot have two hard or two soft consecutive "
                          "constraints on rest shifts.");

    spResources_.push_back(pSPResources);
    masterRotationGraphResources_.push_back(pRotPResources);
    masterConstraintResources_.push_back(pMPResources);
  }
}

//------------------------------------------------------------------------------
// Build the variable of the rotation as well as all the affected constraints
// with their coefficients. if s=-1, the nurse i works on all shifts
//------------------------------------------------------------------------------
MyVar *RotationMP::addColumn(int nurseNum, const RCSolution &solution) {
  // Build rotation from RCSolution
  RotationPattern rotation(solution, nurseNum);
  rotation.treeLevel_ = pModel_->getCurrentTreeLevel();
#ifdef DBG
  computePatternCost(&rotation);
  DualCosts dualCosts(this);
  rotation.checkReducedCost(dualCosts, pPricer()->isLastRunOptimal());
  checkIfPatternAlreadyPresent(rotation);
#endif
  return createColumn(rotation, "rotation");
}

void RotationMP::buildResourceCons(const SolverParam &param) {
  totalShiftDurationConstraint_ = new TotalShiftDurationConstraint(this);
  totalWeekendConstraint_ = new TotalWeekendConstraint(this);

  for (const PLiveNurse &pN : theLiveNurses_) {
    for (const PBoundedResource &pR : masterConstraintResources_[pN->num_]) {
      if (pR->isHard()) {
//       auto pHR = std::dynamic_pointer_cast<HardConsWeekendShiftResource>(pR);
//        if (pHR) {
//          buildConsWeekendShiftResourceCons(pHR, nullptr, *pN);
//          continue;
//        }
        auto pHR2 =
            std::dynamic_pointer_cast<HardTotalShiftDurationResource>(pR);
        if (pHR2) {
          totalShiftDurationConstraint_->addConstraintFor(pHR2, nullptr, *pN);
          continue;
        }
        auto pHR3 = std::dynamic_pointer_cast<HardTotalWeekendsResource>(pR);
        if (pHR3) {
          totalWeekendConstraint_->addConstraintFor(pHR3, nullptr, *pN);
          continue;
        }
        Tools::throwError("Hard constraint not defined in "
                          "RotationMP::buildResourceCons");
      } else {
//       auto pSR = std::dynamic_pointer_cast<SoftConsWeekendShiftResource>(pR);
//        if (pSR) {
//          buildConsWeekendShiftResourceCons(nullptr, pSR.get(), *pN);
//          continue;
//        }
        auto pSR2 =
            std::dynamic_pointer_cast<SoftTotalShiftDurationResource>(pR);
        if (pSR2) {
          totalShiftDurationConstraint_->addConstraintFor(nullptr, pSR2, *pN);
          continue;
        }
        auto pSR3 = std::dynamic_pointer_cast<SoftTotalWeekendsResource>(pR);
        if (pSR3) {
          totalWeekendConstraint_->addConstraintFor(nullptr, pSR3, *pN);
          continue;
        }
        Tools::throwError("Soft constraint not defined in "
                          "RotationMP::buildResourceCons");
      }
    }
  }
}

double RotationMP::getColumnsCost(CostType costType) const {
  double cost = 0;
  if (costType == CONS_REST_COST)
    return pModel_->getTotalCost(rotationGraphConstraint_->getVariables())
        //        + pModel_->getTotalCost(longRestingVars_);
        // cost for empty rotation: rotation for initial state followed by
        // rest -> already included in longRestingVars_
        - rotationGraphConstraint_->getInitialStateCost(costType)
            // just initial rest costs;
        + MasterProblem::getColumnsCost(
            CONS_REST_COST, pModel_->getActiveColumns());

  cost = MasterProblem::getColumnsCost(costType, pModel_->getActiveColumns());
  if (costType == ROTATION_COST)  // add rest costs + historical costs
    cost += pModel_->getTotalCost(rotationGraphConstraint_->getVariables());
//        + pModel_->getTotalCost(longRestingVars_);
  else  // add historical non resting costs
    cost += rotationGraphConstraint_->getInitialStateCost(costType);
  return cost;
}

double RotationMP::getDaysCost() const {
  return pModel_->getTotalCost(totalShiftDurationConstraint_->getVariables());
}

double RotationMP::getWeekendCost() const {
  return pModel_->getTotalCost(totalWeekendConstraint_->getVariables());
}
