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
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/modeler/BcpModeler.h"
#include "solvers/mp/RCPricer.h"
#include "solvers/mp/TreeManager.h"
#include "solvers/mp/sp/rcspp/resources/AlternativeShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"

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

// when branching on this column, this method add the corresponding forbidden
// shifts to the set.
// It will forbid all the shifts that would be worked on a day that is already
// covered by this column.
// Moreover, there needs to be a resting day before and after each rotation,
// so the shifts can also be forbidden on these two days
// (if the rotation is not at an extremity of the horizon).
void RotationColumn::addForbiddenShifts(
    std::set<std::pair<int, int> > *forbiddenShifts,
    int nbShifts, PDemand pDemand) const {
  // from the previous day to the day after the end of the rotation,
  // forbid any work shifts
  for (int day = firstDayId()-1; day <= lastDayId()+1; day++) {
    if (day < pDemand->firstDayId_) continue;
    if (day >= pDemand->firstDayId_ + pDemand->nDays_) continue;
    for (int i = 1; i < nbShifts; ++i)
      forbiddenShifts->insert(pair<int, int>(day, i));
  }
}

void RotationColumn::checkReducedCost(const DualCosts &dualCosts,
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
    if (abs(reducedCost_ - reducedCost) > EPSILON) {
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
    MasterProblem(pScenario, solver),
    rotationGraphConstraint_(nullptr),
    dynamicConstraints_(nullptr),
    totalShiftDurationConstraint_(nullptr),
    totalWeekendConstraint_(nullptr) {}

RotationMP::RotationMP(const PScenario& pScenario,
                       SolverType solver,
                       const SolverParam &param):
    RotationMP(pScenario, solver) {
  RotationMP::build(param);
  setParameters(param);
}

RotationMP::~RotationMP() {
  delete rotationGraphConstraint_;
  delete dynamicConstraints_;
  delete totalShiftDurationConstraint_;
  delete totalWeekendConstraint_;
//  delete consWeekendConstraints_;
}

PColumn RotationMP::getPColumn(MyVar *var) const {
  return std::make_shared<RotationColumn>(var, pScenario_);
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

  /* Dynamic constraints (need min/max constraints) */
  if (dynamicWeights_.version() > 0)
    dynamicConstraints_ = new DynamicConstraints(
        this, totalShiftDurationConstraint_, totalWeekendConstraint_);
}

void RotationMP::update(const PDemand& pDemand) {
  if (dynamicConstraints_ == nullptr && dynamicWeights_.version() > 0)
    Tools::throwError("Cannot defined new dynamic weights if none were "
                      "defined at the initial build.");
  MasterProblem::update(pDemand);
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
        RotationColumn rotation(firstDay, pShifts, i);
        computeColumnCost(&rotation);
        pModel_->addInitialColumn(createColumn(rotation, baseName.c_str()));
        pShifts.clear();
        workedOnPreviousDay = false;
      }
    }
    // if work on the last day, build the rotation
    if (workedOnPreviousDay) {
      RotationColumn rotation(firstDay, pShifts, i);
      computeColumnCost(&rotation);
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
  for (const auto &vR : pResources_) {
    vector<PResource> pSPResources, pRotResources;
    vector<PBoundedResource> pRotBPResources, pMPResources;
    for (const auto &pR : vR) {
      // cast to shared_ptr<BoundedResource>
      PBoundedResource pBR =
          std::dynamic_pointer_cast<BoundedResource>(pR);
      // TODO(AL): this must be corrected, because we need to handle
      //  constraints that do not derive from BoundedResource: at this stage,
      //  if we deal with INRC2 instances, the corresponding constraints are
      //  the preferences and the complete weekends, which both belong to the
      //  subproblem, so the code below should be a simple fix
      if (pBR == nullptr) {
        if (pR->isInRotationMaster()) {
          std::cerr << "WARNING: Presently, the rotation decomposition cannot "
                       "handle non-bounded resources that belong to the "
                       "master problem. " << pR->name
                    << " is therefore ignored." << std::endl;
          continue;
        }
        // 1. Try and cast as alternative shift resource
        auto pRAlt = std::dynamic_pointer_cast<AlternativeShiftResource>(pR);
        if (pRAlt != nullptr) {
          pRAlt->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pRAlt);
          continue;
        }
        // 2. Try and cast as forbidden pattern resource
        auto pRForb = std::dynamic_pointer_cast<ForbiddenPatternResource>(pR);
        if (pRForb != nullptr) {
          pRForb->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pRForb);
          continue;
        }
        // 3. first try to cast as PreferenceResource
        auto pRPref = std::dynamic_pointer_cast<PreferenceResource>(pR);
        if (pRPref != nullptr) {
          pRPref->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pRPref);
          continue;
        }
        // 4. Try and cast as identical weekend resource
        auto pRIdent = std::dynamic_pointer_cast<IdentWeekendResource>(pR);
        if (pRIdent != nullptr) {
          pRIdent->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pRIdent);
          continue;
        }

        std::cerr << "WARNING: Presently, the rotation decomposition cannot "
                     "handle " << pR->name
                  << ". It is therefore ignored." << std::endl;
      } else if (pBR->isHard()) {
        // cast hard constraint
        auto pHR = std::dynamic_pointer_cast<HardConsShiftResource>(pBR);
        // if hard cons -> add to sub problem and rotation if rest
        if (pHR) {
          pHR->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pHR);
          if (pHR->pShift()->isRest())
            pRotBPResources.push_back(pHR);
        } else {
          // add to master
          pBR->setId(static_cast<int>(pMPResources.size()));
          pMPResources.push_back(pBR);
        }
      } else {
        auto pSR = std::dynamic_pointer_cast<SoftConsShiftResource>(pBR);
        // if soft cons -> add to sub problem and rotation if rest
        if (pSR) {
          pSR->setId(static_cast<int>(pSPResources.size()));
          pSPResources.push_back(pSR);
          if (pSR->pShift()->isRest())
            pRotBPResources.push_back(pSR);
        } else {
          // add to master
          pBR->setId(static_cast<int>(pMPResources.size()));
          pMPResources.push_back(pBR);
        }
      }
    }

    // masterRotationGraphResources_ can contains at max one soft and one hard
    // constraint on Rest
    if (pRotBPResources.size() > 2)
      Tools::throwError("Cannot have more than two consecutive constraints "
                        "on rest shifts.");
    if (pRotBPResources.size() == 2)
      if (pRotBPResources.front()->isHard() ^ pRotBPResources.back()->isHard())
        Tools::throwError("Cannot have two hard or two soft consecutive "
                          "constraints on rest shifts.");

    spResources_.push_back(pSPResources);
    masterRotationGraphResources_.push_back(pRotBPResources);
    masterConstraintResources_.push_back(pMPResources);
  }
}

//------------------------------------------------------------------------------
// Build the variable of the rotation as well as all the affected constraints
// with their coefficients. if s=-1, the nurse i works on all shifts
//------------------------------------------------------------------------------
MyVar *RotationMP::addColumn(int nurseNum, const RCSolution &solution) {
  // Build rotation from RCSolution
  RotationColumn rotation(solution, nurseNum);

#ifdef DBG
  computeColumnCost(&rotation);
  DualCosts dualCosts(this);
  rotation.checkReducedCost(dualCosts, pPricer()->isLastRunOptimal());
  checkIfColumnAlreadyPresent(rotation, true);
#endif

  return createColumn(rotation, "rotation");
}

void RotationMP::buildResourceCons(const SolverParam &param) {
  totalShiftDurationConstraint_ = new TotalShiftDurationConstraint(this);
  totalWeekendConstraint_ = new TotalWeekendConstraint(this);

  for (const PLiveNurse &pN : theLiveNurses_) {
    for (const PBoundedResource &pR : masterConstraintResources_[pN->num_]) {
      if (pR->isHard()) {
        auto pHR = std::dynamic_pointer_cast<HardConsWeekendShiftResource>(pR);
        if (pHR) {
//          if (consWeekendConstraints_ == nullptr)
//            consWeekendConstraints_ = new ConsWeekendConstraint(this);
//          consWeekendConstraints_->addConstraintFor(pHR, nullptr, *pN);
          continue;
        }
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
        std::cerr << "Hard constraint " << pR->name << " not defined in "
                  << "RotationMP::buildResourceCons()." << std::endl;
      } else {
        auto pSR = std::dynamic_pointer_cast<SoftConsWeekendShiftResource>(pR);
        if (pSR) {
//          if (consWeekendConstraints_ == nullptr)
//            consWeekendConstraints_ = new ConsWeekendConstraint(this);
//          consWeekendConstraints_->addConstraintFor(nullptr, pSR, *pN);
          continue;
        }
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

        std::cerr << "Soft constraint " << pR->name << " not defined in "
                  << "RotationMP::buildResourceCons()." << std::endl;
      }
    }
  }
}

// return the costs of all active columns associated to the type
// add the initial costs associated to the rest arcs
// add the cost of the rest arcs
std::map<CostType, double> RotationMP::getColumnsCosts() const {
  std::map<CostType, double> costs = MasterProblem::getColumnsCosts(),
      rotGraphCosts = rotationGraphConstraint_->getCosts();
  for (const auto &p : rotGraphCosts)
    costs[p.first] += p.second;
  return costs;
}
