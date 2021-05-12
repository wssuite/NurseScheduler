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

void RotationPattern::checkReducedCost(const PDualCosts &pCosts,
                                       bool printBadPricing) {
  // check if pNurse points to a nurse
  if (nurseNum_ == -1)
    Tools::throwError("LiveNurse = NULL");

  /************************************************
   * Compute all the dual costs of a rotation:
   ************************************************/

  double reducedCost(cost_);

  /* Working dual cost */
  for (int k = firstDay(); k <= lastDay(); ++k)
    reducedCost -= pCosts->workedDayShiftCost(k, shift(k));
  /* Start working dual cost */
  reducedCost -= pCosts->startWorkCost(firstDay());
  /* Stop working dual cost */
  reducedCost -= pCosts->endWorkCost(lastDay());
  /* Working on weekend */
  if (Tools::isSunday(firstDay()))
    reducedCost -= pCosts->workedWeekendCost();
  for (int k = firstDay(); k <= lastDay(); ++k)
    if (Tools::isSaturday(k))
      reducedCost -= pCosts->workedWeekendCost();


  // Display: set to true if you want to display the details of the cost
  if (std::fabs(reducedCost_ - reducedCost) / (1 - reducedCost) > 1e-3) {
    // if do not print and not throwing an error
    if (!printBadPricing && reducedCost_ > reducedCost + 1e-3) return;

    cout << "# " << endl;
    cout << "# " << endl;
    cout << "# Bad dual cost: "
         << reducedCost_ << " != " << reducedCost << endl;
    cout << "# " << endl;
    cout << "#   | Base cost     : + " << cost_ << endl;
    cout << costsToString();

    for (int k = firstDay(); k <= lastDay(); ++k)
      cout << "#   | Work day-shift " << k << ": - "
           << pCosts->workedDayShiftCost(k, shift(k)) << endl;
    cout << "#   | Start work    " << firstDay() << ": - "
         << pCosts->startWorkCost(firstDay()) << endl;
    cout << "#   | Finish Work   " << lastDay() << ": - "
         << pCosts->endWorkCost(lastDay()) << endl;
    if (Tools::isSunday(firstDay()))
      cout << "#   | Weekends      : - " << pCosts->workedWeekendCost() << endl;
    for (int k = firstDay(); k <= lastDay(); ++k)
      if (Tools::isSaturday(k))
        cout << "#   | Weekends      : - "
             << pCosts->workedWeekendCost() << endl;
    std::cout << toString() << "# " << endl;

    // throw an error only when a significant misprice
    // Indeed, if the real reduced cost (reducedCost) is greater than the one
    // found by the pricing, some columns with a positive reduced cost could be
    // generated.
    // The reason why the other situation can arise is that some path in the
    // subproblem could under estimate the real cost. These paths won't be
    // found when the subproblems are solved at optimality, but could  be
    // present when using heuristics.
    if (reducedCost_ < reducedCost + 1e-3)
      Tools::throwError("Invalid pricing of a rotation.");
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
    totalShiftDurationVars_(pScenario->nNurses()),
    totalWeekendVars_(pScenario->nNurses()),
    totalShiftDurationCons_(pScenario->nNurses()),
    totalWeekendCons_(pScenario->nNurses()),
    restsPerDay_(pScenario->nNurses()),
    restVars_(pScenario->nNurses()),
    initialStateRCSolutions_(pScenario->nNurses()),
    minWorkedDaysAvgVars_(pScenario->nNurses()),
    maxWorkedDaysAvgVars_(pScenario->nNurses()),
    maxWorkedWeekendAvgVars_(pScenario_->nNurses()),
    minWorkedDaysContractAvgVars_(pScenario->nContracts()),
    maxWorkedDaysContractAvgVars_(pScenario->nContracts()),
    maxWorkedWeekendContractAvgVars_(pScenario_->nContracts()),
    minWorkedDaysAvgCons_(pScenario->nNurses()),
    maxWorkedDaysAvgCons_(pScenario->nNurses()),
    maxWorkedWeekendAvgCons_(pScenario_->nNurses()),
    minWorkedDaysContractAvgCons_(pScenario->nContracts()),
    maxWorkedDaysContractAvgCons_(pScenario->nContracts()),
    maxWorkedWeekendContractAvgCons_(pScenario_->nContracts()) {
  // initialize the vectors indicating whether the min/max total constraints
  // with averaged bounds are considered
  Tools::initVector(&isMinWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedDaysAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMaxWorkedWeekendAvgCons_, pScenario_->nNurses(), false);
  Tools::initVector(&isMinWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedDaysContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
  Tools::initVector(&isMaxWorkedWeekendContractAvgCons_,
                    pScenario_->nContracts(),
                    false);
}

RotationMP::~RotationMP() = default;

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
  buildRotationGraph(param);

  /* Min/Max constraints */
  buildResourceCons(param);

  /* Dynamic constraints */
  buildDynamicCons(param);
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
    bool workedLastDay = false;
    std::vector<PShift> pShifts;
    // build all the successive rotation of this nurse
    for (int k = 0; k < pDemand_->nDays_; ++k) {
      // shift=0 => rest
      const PShift &pShift = roster.pShift(k);
      // if work, insert the shift in the map
      if (pShift->isWork()) {
        if (!workedLastDay) firstDay = k;
        pShifts[k] = pShift;
        workedLastDay = true;
      } else if (workedLastDay) {
        // if stop to work, build the rotation
        RotationPattern rotation(firstDay, pShifts, i);
        computePatternCost(&rotation);
        pModel_->addInitialColumn(addRotation(rotation, baseName.c_str()));
        pShifts.clear();
        workedLastDay = false;
      }
    }
    // if work on the last day, build the rotation
    if (workedLastDay) {
      RotationPattern rotation(firstDay, pShifts, i);
      computePatternCost(&rotation);
      pModel_->addInitialColumn(addRotation(rotation, baseName.c_str()));
    }
  }
}

PDualCosts RotationMP::buildDualCosts(PLiveNurse pNurse) const {
  return std::make_shared<RotationDualCosts>(
      getShiftsDualValues(pNurse),
      getStartWorkDualValues(pNurse),
      getEndWorkDualValues(pNurse),
      getWorkedWeekendDualValue(pNurse),
      getConstantDualvalue(pNurse));
}

vector2D<double> RotationMP::getShiftsDualValues(PLiveNurse pNurse) const {
  vector2D<double> dualValues = MasterProblem::getShiftsDualValues(pNurse);

  const int i = pNurse->num_;
  const int p = pNurse->pContract_->id_;

  /* Dynamic constraints */
  double minWorkedDaysAvg =
      isMinWorkedDaysAvgCons_[i] ?
      pModel_->getDual(minWorkedDaysAvgCons_[i], true) : 0.0;
  double maxWorkedDaysAvg =
      isMaxWorkedDaysAvgCons_[i] ?
      pModel_->getDual(maxWorkedDaysAvgCons_[i], true) : 0.0;

  double minWorkedDaysContractAvg =
      isMinWorkedDaysContractAvgCons_[p] ?
      pModel_->getDual(minWorkedDaysContractAvgCons_[p], true) : 0.0;
  double maxWorkedDaysContractAvg =
      isMaxWorkedDaysContractAvgCons_[p] ?
      pModel_->getDual(maxWorkedDaysContractAvgCons_[p], true) : 0.0;

  /* Min/Max constraints */
  double dynamicDual = minWorkedDaysAvg + minWorkedDaysContractAvg
      + maxWorkedDaysAvg + maxWorkedDaysContractAvg;


  for (int k = 0; k < nDays(); ++k) {
    vector<double> &dualValues2 = dualValues[k];
    for (const PShift &pS : pScenario_->pShifts()) {
      if (pS->isRest()) continue;
      // adjust the dual if the resource is active and
      // in function of the time duration of the shift
      double d = dynamicDual * pS->duration;
      for (const auto &p : totalShiftDurationCons_[pNurse->num_])
        if (p.first->pResource()->isActive(k, *pS))
          d += pModel_->getDual(p.second) * p.first->duration(pS);
      dualValues2[pS->id - 1] += d;
    }
  }

  return dualValues;
}

vector<double> RotationMP::getStartWorkDualValues(
    const PLiveNurse& pNurse) const {
  int i = pNurse->num_;
  const PRCGraph &pG = pRotationGraphs_[i];
  vector<double> dualValues(nDays());

  // get dual value associated to the origin node in the graph,
  // i.e. the node on the previous day : either the source or a maxRest node
  dualValues[0] = pModel_->getDual(restCons_[i][pG->pSource(0)->id], true);

  // get dual values associated to the work flow constraints
  // don't take into account the last which is a sink
  for (int k = 1; k < nDays(); ++k)
    dualValues[k] = pModel_->getDual(
        restCons_[i][maxRestNodes_[i][k - 1]->id], true);

  return dualValues;
}

vector<double> RotationMP::getEndWorkDualValues(const PLiveNurse& pNurse)
const {
  const int i = pNurse->num_;
  const PRCGraph &pG = pRotationGraphs_[i];
  vector<double> dualValues(nDays());

  // get dual values associated to the work flow constraints,
  // i.e the firstRest node of the last day of the rotation
  // -1 corresponds to the coefficient
  for (int k = 0; k < nDays(); ++k)
    dualValues[k] = -pModel_->getDual(
        restCons_[i][firstRestNodes_[i][k]->id], true);

  return dualValues;
}

double RotationMP::getWorkedWeekendDualValue(const PLiveNurse& pNurse) const {
  int id = pNurse->num_;
  double dualVal = 0;
  if (isMaxWorkedWeekendAvgCons_[id])
    dualVal += pModel_->getDual(maxWorkedWeekendAvgCons_[id], true);
  if (isMaxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_])
    dualVal += pModel_->getDual(
        maxWorkedWeekendContractAvgCons_[pNurse->pContract_->id_], true);
  for (const auto &p : totalWeekendCons_[pNurse->num_])
    dualVal += pModel_->getDual(p.second);
  return dualVal;
}

PDualCosts RotationMP::buildRandomDualCosts(bool optimalDemandConsidered,
                                            int NDaysShifts) const {
  return std::make_shared<RotationDualCosts>(
      getRandomWorkedDualCosts(optimalDemandConsidered, NDaysShifts),
      Tools::randomDoubleVector(
          pDemand_->nDays_,
          0, 7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND),
      Tools::randomDoubleVector(
          pDemand_->nDays_,
          0,
          7*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND),
      Tools::randomDouble(0, 2*pScenario_->weights().WEIGHT_TOTAL_WEEKENDS),
      Tools::randomDouble(-10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND,
                          10*pScenario_->weights().WEIGHT_OPTIMAL_DEMAND));
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
  PDualCosts costs = buildDualCosts(theLiveNurses_[nurseNum]);
  rotation.checkReducedCost(costs, pPricer_->isLastRunOptimal());
  checkIfPatternAlreadyPresent(rotation.getCompactPattern());
#endif
  return addRotation(rotation, "rotation", false);
}

MyVar *RotationMP::addRotation(const RotationPattern &rotation,
                               const char *baseName,
                               bool coreVar) {
  // nurse index
  const int nurseNum = rotation.nurseNum();

  // Column var, its name, and affected constraints with their coefficients
  MyVar *var;
  char name[255];
  vector<MyCons *> cons;
  vector<double> coeffs;

  /* resource constraints */
  addResourceConsToCol(rotation, &cons, &coeffs);

  /* Dynamic constraints */
  int nbWeekends = Tools::nWeekendsInInterval(rotation.firstDay(),
                                              rotation.lastDay());
  addDynamicConsToCol(
      &cons, &coeffs, nurseNum, rotation.duration(), nbWeekends);

  /* Skills coverage constraints */
  addSkillsCoverageConsToCol(&cons, &coeffs, rotation);

  snprintf(name, sizeof(name), "%s_N%d", baseName, nurseNum);
  if (coreVar) {
    pModel_->createPositiveVar(&var,
                               name,
                               rotation.cost(),
                               rotation.getCompactPattern());
    for (unsigned int i = 0; i < cons.size(); i++)
      pModel_->addCoefLinear(cons[i], var, coeffs[i]);
  } else {
    /* Rotation constraints
      They are added only for real rotations to be sure that
      the artificial variables can always be used to create a feasible solution
     */
    addRotationConsToCol(rotation, &cons, &coeffs);

    pModel_->createIntColumn(&var,
                             name,
                             rotation.cost(),
                             rotation.getCompactPattern(),
                             rotation.reducedCost(),
                             cons,
                             coeffs);
  }
  return var;
}

/*
 * Rotation constraints : build a rotation graph
 */
void RotationMP::buildRotationGraph(const SolverParam &param) {
  pRotationGraphs_.clear();
  restVars_.clear();
  restsPerDay_.clear();
  restCons_.clear();
  // build the rotation network for each nurse
  for (const PLiveNurse &pN : theLiveNurses_) {
    // build a graph of 2 shifts: rest and work
    auto pG = std::make_shared<RCGraph>(nDays(), nShifts());
    createRotationNodes(pG, pN);
    createRotationArcs(pG, pN);
    createRotationArcsVars(pG, pN);
    createRotationNodesCons(pG, pN);
    pRotationGraphs_.push_back(pG);
  }
}

void RotationMP::createRotationNodes(const PRCGraph &pG, const PLiveNurse &pN) {
  /* Create nodes: one source, nDays free rest and non free nodes */
  const PShift &pRest = pScenario_->pRestShift();
  PRCNode pSource =
      pG->addSingleNode(SOURCE_NODE, -1, pN->pStateIni_->pShift_);
  // the source is present
  vector<PRCNode> firstRestNodes, maxRestNodes;
  for (int k = 0; k < nDays(); k++) {
    NodeType nT = (k == nDays() - 1) ? SINK_NODE : PRINCIPAL_NETWORK;
    firstRestNodes.push_back(pG->addSingleNode(nT, k, pRest));
    maxRestNodes.push_back(pG->addSingleNode(nT, k, pRest));
  }
  firstRestNodes_.push_back(firstRestNodes);
  maxRestNodes_.push_back(maxRestNodes);
}

void RotationMP::createRotationArcs(const PRCGraph &pG, const PLiveNurse &pN) {
  /* Create rest arcs */
  const PShift &pRest = pScenario_->pRestShift();
  vector<PRCNode> &firstRestNodes = firstRestNodes_[pN->num_],
      &maxRestNodes = maxRestNodes_[pN->num_];
  // if unlimited rest, take the min
  std::pair<int, double> minCons = minConsRest(pN), maxCons = maxConsRest(pN);
  for (int k = 0; k < nDays(); ++k) {
    // create all arcs from first rest node to max rest node
    PRCNode pOrigin = k == 0 ? pG->pSource(0) : firstRestNodes[k - 1];
    std::vector<PShift> pRestShifts;
    for (int l = minCons.first; l <= maxCons.first; l++) {
      if (k+l > nDays()) continue;
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

std::pair<int, double> RotationMP::minConsRest(const PLiveNurse &pN) {
  for (const auto &pR : masterRotationGraphResources_[pN->num_])
    if (pR->isHard()) {
      // Override LB if define a real hard bound
      if (pR->getLb() > 1) return {pR->getLb(), LARGE_SCORE};
    }
  // otherwise, just use 1 which is always the default LB
  // Do not use a soft LB, as the number of consecutive rest could be lower
  return {1, pScenario_->weights().WEIGHT_CONS_DAYS_OFF};
}

std::pair<int, double> RotationMP::maxConsRest(const PLiveNurse &pN) {
  int maxR = 1;
  double cost = 0;
  for (const auto &pR : masterRotationGraphResources_[pN->num_]) {
    if (pR->isHard()) {
      // Override UB if define a real bound
      if (pR->getUb() < nDays()) return {pR->getUb(), LARGE_SCORE};
      if (pR->getLb() > maxR) {
        maxR = pR->getLb();
        cost = 0;
      }
    } else {
      // soft UB present
      if (pR->getUb() < nDays() && pR->getUb() > maxR) {
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

void RotationMP::addRotationRestCost(const PRCArc &pArc, const PLiveNurse &pN) {
#ifdef DBG
  // verify that it is indeed a rest arc
  for (const PShift &pS : pArc->stretch.pShifts())
    if (pS->isWork())
      Tools::throwError("addRotationRestCost works only for a rest arc.");
#endif
  // take a copy of the stretch and add a work shift at the end to ensure
  // to price correctly the rest costs LB if ends before the last day
  Stretch st = pArc->stretch;
  if (st.lastDay() < nDays() - 1)
    st.pushBack(Stretch(st.firstDay() + 1, {pScenario_->pAnyWorkShift()}));

  // build rotation
  RCSolution sol(st, 0, DBL_MAX);
  sol.resetCosts();

  // get the vector of resources of the subproblem
  auto pResources = getSPResources(pN);
  // create expander for stretch
  double c = pArc->cost;
  for (const auto &pR : pResources) {
    pR->initialize(*pArc->origin->pAShift, sol, pArc);
    sol.addCost(pArc->cost - c, resourceCostType(pR, pN));
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
    sol.addCost(
        pL->cost() - c, resourceCostType(pResources[pE->resourceId], pN));
    c = pL->cost();
  }

  // add the label cost to the arc
  pArc->addBaseCost(pL->cost());

  if (sol.firstDay() == 0 &&
      initialStateRCSolutions_[pN->num_].firstDay() == -1)
    initialStateRCSolutions_[pN->num_] = sol;
}

void RotationMP::createRotationArcsVars(
    const PRCGraph &pG, const PLiveNurse &pN) {
  /* create variables for all arcs */
  char name[255];
  vector<MyVar*> vars;
  vector2D<MyVar*> varsByDay(nDays());
  MyVar *var;
  for (const PRCArc &pArc : pG->pArcs()) {
    // create a rotation from the stretch
    RotationPattern rot(
        RCSolution(pArc->stretch, pArc->cost, DBL_MAX), pN->num_);
    snprintf(name, sizeof(name), "RestingVars_N%d_%d_%d", pN->num_,
        pArc->stretch.firstDay(), pArc->stretch.lastDay());
    pModel_->createPositiveVar(&var, name, pArc->cost, rot.getCompactPattern());
    vars.push_back(var);
    for (int k = pArc->stretch.firstDay(); k <= pArc->stretch.lastDay(); k++)
      varsByDay[k].push_back(var);
  }
  restVars_.push_back(vars);
  restsPerDay_.push_back(varsByDay);
}

void RotationMP::createRotationNodesCons(
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
      pModel_->createGEConsLinear(&con, name, -1, vars, coeffs);
    } else {
      int outFlow = 0;
      if (pNode->type == SOURCE_NODE) outFlow = 1;
      pModel_->createEQConsLinear(&con, name, outFlow, vars, coeffs);
    }
    cons.push_back(con);
  }
  restCons_.push_back(cons);
}

void RotationMP::addRotationConsToCol(const RotationPattern &pRot,
                                     std::vector<MyCons *> *cons,
                                     std::vector<double> *coeffs) {
  // add the arc corresponding to the rotation
  const PRCGraph &pG = pRotationGraphs_[pRot.nurseNum()];
  const PRCNode &pOrigin = (pRot.firstDay() == 0) ?
      pG->pSource(0) : maxRestNodes_[pRot.nurseNum()][pRot.firstDay()-1];
  const PRCNode &pTarget = firstRestNodes_[pRot.nurseNum()][pRot.lastDay()];
  pG->addSingleArc(pOrigin, pTarget, pRot, pRot.cost());
  // add the coefficient for both flow conservation constraints
  cons->push_back(restCons_[pRot.nurseNum()][pOrigin->id]);
  coeffs->push_back(1);
  cons->push_back(restCons_[pRot.nurseNum()][pTarget->id]);
  coeffs->push_back(-1);
}

void RotationMP::buildResourceCons(const SolverParam &param) {
  resourcesVars.resize(nNurses());
  for (const PLiveNurse &pN : theLiveNurses_) {
    for (const PBoundedResource &pR : masterConstraintResources_[pN->num_]) {
      if (pR->isHard()) {
        auto pHR = std::dynamic_pointer_cast<HardConsWeekendShiftResource>(pR);
        if (pHR) {
          buildConsWeekendShiftResourceCons(pHR.get(), nullptr, *pN);
          continue;
        }
        auto pHR2 =
            std::dynamic_pointer_cast<HardTotalShiftDurationResource>(pR);
        if (pHR2) {
          buildTotalShiftDurationResourceCons(pHR2.get(), nullptr, *pN);
          continue;
        }
        auto pHR3 = std::dynamic_pointer_cast<HardTotalWeekendsResource>(pR);
        if (pHR3) {
          buildTotalWeekendsResourceCons(pHR3.get(), nullptr, *pN);
          continue;
        }
        Tools::throwError("Hard constraint not defined in "
                          "RotationMP::buildResourceCons");
      } else {
        auto pSR = std::dynamic_pointer_cast<SoftConsWeekendShiftResource>(pR);
        if (pSR) {
          buildConsWeekendShiftResourceCons(nullptr, pSR.get(), *pN);
          continue;
        }
        auto pSR2 =
            std::dynamic_pointer_cast<SoftTotalShiftDurationResource>(pR);
        if (pSR2) {
          buildTotalShiftDurationResourceCons(nullptr, pSR2.get(), *pN);
          continue;
        }
        auto pSR3 = std::dynamic_pointer_cast<SoftTotalWeekendsResource>(pR);
        if (pSR3) {
          buildTotalWeekendsResourceCons(nullptr, pSR3.get(), *pN);
          continue;
        }
        Tools::throwError("Hard constraint not defined in "
                          "RotationMP::buildResourceCons");
      }
    }
  }
}

std::pair<vector<MyVar*>, vector<MyCons*>>
RotationMP::buildBoundedResourceCons(
    HardBoundedResource *pHR, SoftBoundedResource *pSR, const LiveNurse &pN) {
  char name[255];
  vector<double> coeffs;
  vector<MyVar*> vars;
  vector<MyCons*> cons;
  // compute bounds and slacks
  int lb = 0, ub = INT_MAX, lb_slack = 0, ub_slack = 0;
  bool lb_soft = false, ub_soft = false;
  // update bounds
  if (pHR) {
    lb = pHR->getLb();
    ub = pHR->getUb();
  }
  if (pSR) {
    if (pSR->getLb() > lb) {
      // update slack
      lb_soft = true;
      lb_slack = pSR->getLb() - lb;
      lb = pSR->getLb();
    }
    if (pSR->getUb() < ub) {
      // update slack
      ub_soft = true;
      // if no hard constraint on the UB, the slack is not bounded
      ub_slack = ub < INT_MAX ? ub - pSR->getUb() : INT_MAX;
      ub = pSR->getUb();
    }
  }

  /* create constraints and slack variables if soft constraint */
  const char* cName = pSR ? pSR->name.c_str() : pHR->name.c_str();
  MyCons *pC;
  // LB slack and constraint
  if (lb_soft) {
    vector<MyVar*> lbVars;
    vector<double>  lbCoeffs;
    if (pSR) {
      MyVar *vLB;
      snprintf(name, sizeof(name), "%s_lb_N%d_slack", cName, pN.num_);
      pModel_->createPositiveVar(&vLB, name, pSR->getLbCost(), {}, 0, lb_slack);
      lbVars = {vLB}; lbCoeffs = {1};
      vars.push_back(vLB);
    }
    snprintf(name, sizeof(name), "%s_lb_N%d", cName, pN.num_);
    pModel_->createGEConsLinear(&pC, name, lb, lbVars, lbCoeffs);
    cons.push_back(pC);
  }
  // UB slack and constraint
  if (ub_soft) {
    vector<MyVar*> ubVars;
    vector<double>  ubCoeffs;
    if (pSR) {
      MyVar *vUB;
      snprintf(name, sizeof(name), "%s_ub_N%d_slack", cName, pN.num_);
      pModel_->createPositiveVar(&vUB, name, pSR->getUbCost(), {}, 0, ub_slack);
      ubVars = {vUB}; ubCoeffs = {-1};
      vars.push_back(vUB);
    }
    snprintf(name, sizeof(name), "%s_ub_N%d", cName, pN.num_);
    pModel_->createLEConsLinear(&pC, name, ub, ubVars, ubCoeffs);
    cons.push_back(pC);
  }

  return {vars, cons};
}

void RotationMP::buildTotalShiftDurationResourceCons(
    HardTotalShiftDurationResource *pHR,
    SoftTotalShiftDurationResource *pSR,
    const LiveNurse &pN) {
  auto p = buildBoundedResourceCons(pHR, pSR, pN);
  TotalShiftDuration *pR = pHR;
  if (!pR) pR = pSR;
  totalShiftDurationVars_[pN.num_][pR] = p.first;
  totalShiftDurationCons_[pN.num_][pR] = p.second;
}

void RotationMP::buildTotalWeekendsResourceCons(
    HardTotalWeekendsResource *pHR,
    SoftTotalWeekendsResource *pSR,
    const LiveNurse &pN) {
  auto p = buildBoundedResourceCons(pHR, pSR, pN);
  TotalWeekend *pR = pHR;
  if (!pR) pR = pSR;
  totalWeekendVars_[pN.num_][pR] = p.first;
  totalWeekendCons_[pN.num_][pR] = p.second;
}

void RotationMP::addResourceConsToCol(const RotationPattern &rotation,
                                      std::vector<MyCons *> *cons,
                                      std::vector<double> *coeffs) {
  const PLiveNurse &pN = theLiveNurses_[rotation.nurseNum()];
  // if previous day is the initial state of the nurse
  PAbstractShift prevAShift;
  if (rotation.firstDay() == 0)
    prevAShift = pN->pStateIni_->pShift_;
  // otherwise, it's a rest shift
  else
    prevAShift =  pScenario_->pRestShift();

  for (const auto &p : totalShiftDurationCons_[rotation.nurseNum()]) {
    int c = p.first->computeConsumption(rotation);
    for (MyCons *pC : p.second) {
      cons->push_back(pC);
      coeffs->push_back(c);
    }
    assert(c == rotation.duration());
  }

  for (const auto &p : totalWeekendCons_[rotation.nurseNum()]) {
    PShift prevS = rotation.firstDay() == 0 ?
        pN->pStateIni_->pShift_ : pScenario_->pRestShift();
    bool ready = !p.first->pShift()->includes(*prevS);
    int c = p.first->computeConsumption(rotation, &ready);
    for (MyCons *pC : p.second) {
      cons->push_back(pC);
      coeffs->push_back(c);
    }
    assert(c == Tools::nWeekendsInInterval(rotation.firstDay(),
                                           rotation.lastDay()));
    assert(p.second.size() == 1);
  }
}

/*
 * Dynamic constraints
 */
void RotationMP::buildDynamicCons(const SolverParam &param) {
  char name[255];
  for (int i = 0; i < pScenario_->nNurses(); i++) {
    // add constraints on the total number of shifts to satisfy bounds that
    // correspond to the global bounds averaged over the weeks
    if (!minTotalShiftsAvg_.empty() && !maxTotalShiftsAvg_.empty()
        && !weightTotalShiftsAvg_.empty()) {
      // only add the constraint if is tighter than the already added constraint
      if (minTotalShiftsAvg_[i] > minTotalShifts_[i]) {
        snprintf(name, sizeof(name), "minWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&minWorkedDaysAvgVars_[i],
                                   name,
                                   weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "minWorkedDaysAvgCons_N%d", i);
        vector<MyVar*> vars = {minWorkedDaysAvgVars_[i]};
        vector<double> coeffs = {1};
        // if a LB constraint
        if (minTotalShifts_[i] > 0) {
          vars.push_back(totalShiftDurationVars_[i].begin()->second.at(0));
          coeffs.push_back(1);
        }
        pModel_->createGEConsLinear(
            &minWorkedDaysAvgCons_[i],
            name,
            minTotalShiftsAvg_[i],
            vars,
            coeffs);

        isMinWorkedDaysAvgCons_[i] = true;
      }

      if (maxTotalShiftsAvg_[i] < maxTotalShifts_[i]) {
        snprintf(name, sizeof(name), "maxWorkedDaysAvgVar_N%d", i);
        pModel_->createPositiveVar(&maxWorkedDaysAvgVars_[i],
                                   name,
                                   weightTotalShiftsAvg_[i]);

        snprintf(name, sizeof(name), "maxWorkedDaysAvgCons_N%d", i);
        MyVar *maxWorkedDaysVar =
            totalShiftDurationVars_[i].begin()->second.back();
        pModel_->createLEConsLinear(
            &maxWorkedDaysAvgCons_[i],
            name,
            maxTotalShiftsAvg_[i],
            {maxWorkedDaysVar, maxWorkedDaysAvgVars_[i]},
            {-1, -1});

        isMaxWorkedDaysAvgCons_[i] = true;
      }
    }

    // STAB: not implemented there yet
    if (!maxTotalWeekendsAvg_.empty() && !weightTotalWeekendsAvg_.empty()
        && maxTotalWeekendsAvg_[i] < theLiveNurses_[i]->maxTotalWeekends()
            - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_) {
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgVar_N%d", i);
      pModel_->createPositiveVar(&maxWorkedWeekendAvgVars_[i],
                                 name,
                                 weightTotalWeekendsAvg_[i]);
      snprintf(name, sizeof(name), "maxWorkedWeekendAvgCons_N%d", i);
      MyVar *maxWorkedWeekendVar = totalWeekendVars_[i].begin()->second.at(0);
      pModel_->createLEConsLinear(
          &maxWorkedWeekendAvgCons_[i],
          name,
          maxTotalWeekendsAvg_[i]
              - theLiveNurses_[i]->pStateIni_->totalWeekendsWorked_,
          {maxWorkedWeekendVar, maxWorkedWeekendAvgVars_[i]},
          {-1, -1});

      isMaxWorkedWeekendAvgCons_[i] = true;
    }
  }

  for (int p = 0; p < pScenario_->nContracts(); ++p) {
    if (!minTotalShiftsContractAvg_.empty()
        && !maxTotalShiftsContractAvg_.empty()
        && !weightTotalShiftsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "minWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&minWorkedDaysContractAvgVars_[p],
                                 name,
                                 weightTotalShiftsContractAvg_[p]);
      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedDaysContractAvgVars_[p],
                                 name,
                                 weightTotalShiftsContractAvg_[p]);

      snprintf(name, sizeof(name), "minWorkedDaysContractAvgCons_P%d", p);
      pModel_->createGEConsLinear(&minWorkedDaysContractAvgCons_[p],
                                  name,
                                  minTotalShiftsContractAvg_[p],
                                  {minWorkedDaysContractAvgVars_[p]},
                                  {1});

      snprintf(name, sizeof(name), "maxWorkedDaysContractAvgCons_P%d", p);
      pModel_->createLEConsLinear(&maxWorkedDaysContractAvgCons_[p],
                                  name,
                                  maxTotalShiftsContractAvg_[p],
                                  {maxWorkedDaysContractAvgVars_[p]},
                                  {-1});

      isMinWorkedDaysContractAvgCons_[p] = true;
      isMaxWorkedDaysContractAvgCons_[p] = true;
    }

    if (!maxTotalWeekendsContractAvg_.empty()
        && !weightTotalWeekendsContractAvg_.empty()) {
      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgVar_P%d", p);
      pModel_->createPositiveVar(&maxWorkedWeekendContractAvgVars_[p],
                                 name,
                                 weightTotalWeekendsContractAvg_[p]);

      snprintf(name, sizeof(name), "maxWorkedWeekendContractAvgCons_C%d", p);
      pModel_->createLEConsLinear(&maxWorkedWeekendContractAvgCons_[p],
                                  name,
                                  maxTotalWeekendsContractAvg_[p],
                                  {maxWorkedWeekendContractAvgVars_[p]},
                                  {-1});

      isMaxWorkedWeekendContractAvgCons_[p] = true;
    }
  }
}

int RotationMP::addDynamicConsToCol(vector<MyCons *> *cons,
                                    vector<double> *coeffs,
                                    int i,
                                    int nbDays,
                                    int nbWeekends) {
  int nbCons(0);
  int p = theLiveNurses_[i]->pContract_->id_;
  if (isMinWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons->push_back(minWorkedDaysAvgCons_[i]);
    coeffs->push_back(nbDays);
  }
  if (isMaxWorkedDaysAvgCons_[i]) {
    ++nbCons;
    cons->push_back(maxWorkedDaysAvgCons_[i]);
    coeffs->push_back(nbDays);
  }
  if (isMinWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons->push_back(minWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(nbDays);
  }
  if (isMaxWorkedDaysContractAvgCons_[p]) {
    ++nbCons;
    cons->push_back(maxWorkedDaysContractAvgCons_[p]);
    coeffs->push_back(nbDays);
  }

  if (nbWeekends) {
    if (isMaxWorkedWeekendAvgCons_[i]) {
      ++nbCons;
      cons->push_back(maxWorkedWeekendAvgCons_[i]);
      coeffs->push_back(nbWeekends);
    }

    if (isMaxWorkedWeekendContractAvgCons_[p]) {
      ++nbCons;
      cons->push_back(maxWorkedWeekendContractAvgCons_[p]);
      coeffs->push_back(nbWeekends);
    }
  }

  return nbCons;
}

double RotationMP::getColumnsCost(CostType costType) const {
  double cost = 0;
  if (costType == CONS_REST_COST)
    return pModel_->getTotalCost(restVars_)
        //        + pModel_->getTotalCost(longRestingVars_);
        // cost for empty rotation: rotation for initial state followed by
        // rest -> already included in longRestingVars_
        - getInitialStateCost(costType)
            // just initial rest costs;
        + MasterProblem::getColumnsCost(
            CONS_REST_COST, pModel_->getActiveColumns());

  cost = MasterProblem::getColumnsCost(costType, pModel_->getActiveColumns());
  if (costType == ROTATION_COST)  // add rest costs + historical costs
    cost += pModel_->getTotalCost(restVars_);
//        + pModel_->getTotalCost(longRestingVars_);
  else  // add historical non resting costs
    cost += getInitialStateCost(costType);
  return cost;
}

double RotationMP::getDaysCost() const {
  return pModel_->getTotalCost(totalShiftDurationVars_);
}

double RotationMP::getWeekendCost() const {
  return pModel_->getTotalCost(totalWeekendVars_);
}

double RotationMP::getInitialStateCost(CostType costType) const {
  double cost = 0;
  for (int i = 0; i < nNurses(); i++)
    for (const auto &pArc : pRotationGraphs_[i]->pSource(0)->outArcs)
      if (pArc->stretch.pShifts().front()->isRest()) {
        double value = pModel_->getVarValue(restVars_[i][pArc->id]);
        if (value > epsilon())
          cost += value * initialStateRCSolutions_[i].costByType(costType);
      }
  return cost;
}
