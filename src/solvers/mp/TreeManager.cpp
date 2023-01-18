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

#include "solvers/mp/TreeManager.h"

#include <list>

#include "solvers/mp/MasterProblem.h"

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::set;

//////////////////////////////////////////////////////////////
//
// B C P T R E E S T A T S   M E T H O D S
//
//////////////////////////////////////////////////////////////

RestTree::RestTree(PScenario pScenario,
                   PDemand pDemand,
                   double epsilon,
                   bool printCurrentNode) :
    MyTree(epsilon, printCurrentNode),
    pScenario_(pScenario),
    pDemand_(pDemand),
    statsRestByDay_(pDemand_->nDays_),
    statsWorkByDay_(pDemand_->nDays_),
    statsRestByNurse_(pScenario->nNurses()),
    statsWorkByNurse_(pScenario->nNurses()) {}

// Forbid the shifts (that are already fixed by the previous branching decision)
// in the generation  of new column.
// This method is useful for the RC pricer
void RestTree::addForbiddenShifts(PLiveNurse pNurse,
                                  set<pair<int, int> > *forbidenShifts) {
  MyNode *pNode = currentNode_.get();
  vector<MyVar *> arcs;
  while (pNode->pParent_) {
    MyNode *currentNode = pNode;
    pNode = pNode->pParent_;
    // If on a rest node, we forbid shifts only if rest is forced
    // In that case, no column can include a shift on the fixed resting day
    //
    auto *restNode = dynamic_cast<RestNode*>(currentNode);
    if (restNode) {
      if (restNode->pNurse_ == pNurse) {
        if (restNode->rest_) {
          for (int i = 1; i < pNurse->pScenario_->nShifts(); ++i)
            forbidenShifts->insert(pair<int, int>(restNode->day_, i));
        } else {
          forbidenShifts->insert(pair<int, int>(restNode->day_, 0));
        }
      }
      continue;
    }

    // If on a shift node, it is natural to forbid all the work shifts
    // in RC pricing.
    auto *shiftNode = dynamic_cast<ShiftNode*>(currentNode);
    if (shiftNode) {
      if (shiftNode->pNurse_ == pNurse)
        for (int s : shiftNode->forbiddenShifts_)
          forbidenShifts->insert(pair<int, int>(shiftNode->day_, s));
      continue;
    }

    // If in a column node , forbid the shifts depending on the pattern
    auto *columnsNode = dynamic_cast<ColumnsNode*>(currentNode);
    if (columnsNode) {
      for (PColumn pCol : columnsNode->columns_)
        if (pCol->nurseNum() == pNurse->num_)
          pCol->addForbiddenShifts(forbidenShifts,
                                  pScenario_->nShifts(),
                                  pDemand_);
      continue;
    }

    // If on a shift demand node or var node, nothing to forbid
    if (dynamic_cast<CoverageNode*>(currentNode) ||
        dynamic_cast<VarNode*>(currentNode))
      continue;

    Tools::throwError("Type of node not recognized.");
  }
}

void RestTree::updateStats(MyNode *node) {
  RestNode *restNode = dynamic_cast<RestNode *>(node);
  if (restNode != 0) {
    double increaseLB = node->getBestLB() - node->pParent_->getBestLB();
    if (restNode->rest_) {
      statsRestByDay_[restNode->day_].first++;
      statsRestByDay_[restNode->day_].second += increaseLB;
      statsRestByNurse_[restNode->pNurse_->num_].first++;
      statsRestByNurse_[restNode->pNurse_->num_].second += increaseLB;
    } else {
      statsWorkByDay_[restNode->day_].first++;
      statsWorkByDay_[restNode->day_].second += increaseLB;
      statsWorkByNurse_[restNode->pNurse_->num_].first++;
      statsWorkByNurse_[restNode->pNurse_->num_].second += increaseLB;
    }
  }

  ShiftNode *shiftNode = dynamic_cast<ShiftNode *>(node);
  if (shiftNode != 0) {
    double increaseLB = node->getBestLB() - node->pParent_->getBestLB();
    if (!shiftNode->work_) {
      statsRestByDay_[shiftNode->day_].first++;
      statsRestByDay_[shiftNode->day_].second += increaseLB;
      statsRestByNurse_[shiftNode->pNurse_->num_].first++;
      statsRestByNurse_[shiftNode->pNurse_->num_].second += increaseLB;
    } else {
      statsWorkByDay_[shiftNode->day_].first++;
      statsWorkByDay_[shiftNode->day_].second += increaseLB;
      statsWorkByNurse_[shiftNode->pNurse_->num_].first++;
      statsWorkByNurse_[shiftNode->pNurse_->num_].second += increaseLB;
    }
  }

  ColumnsNode *colsNode = dynamic_cast<ColumnsNode *>(node);
  if (colsNode != 0) {
    statsCols_.first++;
    statsCols_.second += node->getBestLB() - node->pParent_->getBestLB();
  }
}

string RestTree::writeBranchStats() const {
  std::stringstream rep;
  rep << "";

  int nbNurses = statsRestByNurse_.size();
  int firstDay = pDemand_->firstDayId_, nbDays = pDemand_->nDays_;

  rep << "Stats on Columns" << std::endl;
  char buffer0[100];
  snprintf(buffer0,
           sizeof(buffer0),
           "Has branched %3d times with an average increased "
           "of the lower bound of %.2f",
           statsCols_.first,
           statsCols_.second / statsCols_.first);
  rep << buffer0 << std::endl;
  rep << "-------------------------------------" << std::endl;

  rep << std::endl;
  rep << "Stats by Days" << std::endl;
  rep << "\t\t  ";
  for (int day = 0; day < nbDays; day++)
    rep << "|  " << Day::toDayOfWeekShortName(firstDay + day).at(0)
    << "  ";
  rep << "|" << std::endl;
  rep << writeOneStat("Rest", statsRestByDay_);
  rep << writeOneStat("Work", statsWorkByDay_);
  rep << "-------------------------------------" << std::endl;

  rep << std::endl;
  rep << "Stats by Nurse" << std::endl;
  rep << "\t\t  ";
  for (int n = 0; n < nbNurses; n++) {
    if (n < 10) rep << "|" << pScenario_->pNurse(n)->name_ << " ";
    else
      rep << "|" << pScenario_->pNurse(n)->name_;
  }
  rep << "|" << std::endl;
  rep << writeOneStat("Rest", statsRestByNurse_);
  rep << writeOneStat("Work", statsWorkByNurse_);
  rep << "-------------------------------------" << std::endl;

  return rep.str();
}

string RestTree::writeOneStat(string name,
                              const vector<pair<int, double>> &stats) const {
  std::stringstream rep;

  rep << name << "\t\t  ";
  for (unsigned int n = 0; n < stats.size(); n++) {
    const pair<int, double> &p = stats[n];
    if (p.first) {
      char buffer[100];
      snprintf(buffer, sizeof(buffer), "|%3.2f", p.second / p.first);
      rep << buffer;
    } else {
      rep << "| --- ";
    }
  }
  rep << "|" << std::endl;

  return rep.str();
}

void RestTree::reset() {
  MyTree::reset();
  statsRestByDay_.clear();
  statsRestByDay_.resize(pDemand_->nDays_);
  statsWorkByDay_.clear();
  statsWorkByDay_.resize(pDemand_->nDays_);
  statsRestByNurse_.clear();
  statsRestByNurse_.resize(pScenario_->nNurses());
  statsWorkByNurse_.clear();
  statsWorkByNurse_.resize(pScenario_->nNurses());
  statsCols_.first = 0;
  statsCols_.second = 0;
}

//////////////////////////////////////////////////////////////
//
// B R A N C H I N G  R U L E
//
//////////////////////////////////////////////////////////////


/*************************************************************
 * Diving branching rule: dive then close to .5
 *************************************************************/
Score* addScore(const Score &sc, Tools::FixedSizeList<Score> *bestScores) {
  for (auto it= bestScores->begin(); it != bestScores->end(); ++it)
    if (sc.score > it->score)
      return bestScores->insert(it, sc);
  // if there is still some place, insert score at the end
  if (bestScores->list().size() < bestScores->fixedSize_)
    return bestScores->push_back(sc);
  return nullptr;
}

double ScoreVarCloseHalf::score(PLiveNurse pNurse,
                                int day,
                                const std::vector<int> &shifts,
                                const std::vector<double> &values,
                                double baseScore) const {
  // choose the variable closest to .5
  // look for the value the closest to .5
  double value = 0;
  for (double v : values) value += v;
  baseScore += std::min(value - floor(value), ceil(value) - value);
  if (pNurse->pContract()->isWeekend(day)) baseScore += weekendAdvantage_;
  return baseScore;
}

double ScoreVarCloseHalfWeekendDecrement::score(
    PLiveNurse pNurse,
    int day,
    const std::vector<int> &shifts,
    const std::vector<double> &values,
    double baseScore) const {
  // choose the variable closest to .5
  // start first by weekend days starting from the end,
  // and then week days from the end
  double value = 0;
  for (double v : values) value += v;
  baseScore += std::min(value - floor(value), ceil(value) - value);
  baseScore -= day;
  // ensure that weekend are first
  if (!pNurse->pContract()->isWeekend(day)) baseScore -= pNurse->nbDays_;
  return baseScore;
}

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem *master,
                                     RestTree *tree,
                                     const char *name,
                                     bool randomSwapOfChilfren) :
    MyBranchingRule(name),
    pMaster_(master),
    pTree_(tree),
    pModel_(master->pModel()),
    randomSwapOfChilfren_(randomSwapOfChilfren) {
  // either ScoreVarCloseHalf or ScoreVarBestExpectedLBImprovement
  scoreFunc_ = std::unique_ptr<ScoreVar>(
      new ScoreVarCloseHalf(this, pMaster_->parameters().weekendAdvantage_));
//  scoreFunc_ = std::unique_ptr<ScoreVar>(
//      new ScoreVarCloseHalfWeekendDecrement(this));
}

// add all good candidates
bool DiveBranchingRule::column_node(
    std::vector<MyPBranchingCandidate> *candidates) {
  // retrieve maximum residual value available to branch on columns
  // The sum of the residual value (1 - column value) of the selected columns
  // must be lower than valueMax
  double valueMax = pModel_->getParameters().branchColumnUntilValue_;
  // if negative, use the default value
  if (valueMax < 0) valueMax = pMaster_->getBranchColumnValueMax();
  // if valueMax < epsilon, cannot branch
  if (valueMax < pModel_->epsilon())
    Tools::throwError("Cannot select any columns when branching on columns as "
                      "the maximum residual value (branchColumnUntilValue) is "
                      "too close to 0. You should either increase this value or"
                      "deactivate this feature by setting "
                      "nbDiveIfBranchOnColumns to 0.");

  // look for fractional columns
  // Fix all column above BRANCH_LB
  // search the good candidates
  vector<pair<MyVar *, double> > colCandidates;
  vector<MyVar *> integerFixingCandidates;
  std::vector<MyVar *> fixingCandidates;
  std::vector<MyVar *> otherFixingCandidates;
  vector<PColumn> columns;

  for (MyVar *var : pModel_->getActiveColumns()) {
//    // Several columns are kept -> so could remain fractional
//    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
//    if (pMaster_->isFixNurse(var->getNurseNum())) continue;

    // ROLLING/LNS: do not branch on relaxed columns
    // check if all days are relaxed
    bool relaxed = true;
    for (int k = Column::firstDay(var); k <= Column::lastDay(var); ++k)
      if (!pMaster_->isRelaxDay(k)) {
        relaxed = false;
        break;
      }
    if (relaxed)
      continue;

    // if we have already branched on this variable, continue
    if (var->getLB() > pModel_->epsilon()) continue;

    double value = pModel_->getVarValue(var);
    // if var is non null
    if (value > 1 - pModel_->epsilon()) {
      integerFixingCandidates.push_back(var);
      columns.push_back(pMaster_->getPColumn(var));
      continue;
    }
    if (value > pModel_->epsilon())
      colCandidates.emplace_back(var, 1 - value);
  }

  stable_sort(colCandidates.begin(),
              colCandidates.end(),
              Tools::compareObject<MyVar *>);

  if (valueMax < pModel_->epsilon())
    return false;
  fixingCandidates = chooseColumns(
      colCandidates,
      &columns,
      &valueMax,
      pModel_->getParameters().branchColumnDisjoint_ ?
      DayDisjointComparator() : ColumnsComparator());

  // if branching on columns, try disjoint shifts instead  of days
  if (pModel_->getParameters().branchColumnDisjoint_) {
    ShiftDisjointComparator comp2 = ShiftDisjointComparator();
    otherFixingCandidates = chooseColumns(colCandidates,
                                          &columns,
                                          &valueMax,
                                          comp2);
  }

  if (fixingCandidates.empty() && otherFixingCandidates.empty()
      && integerFixingCandidates.empty())
    return false;

  /* update candidate */
  // if no candidate, create one
  if (candidates->empty())
    candidates->push_back(
        std::make_shared<RCBranchingCandidate>(pTree_->getCurrentNode()));
  for (const MyPBranchingCandidate &candidate : *candidates) {
    int nodeIndex = 0;
    candidate->createNewChild(nodeIndex);
    const MyPBranchingNode &node = candidate->getChild(nodeIndex);
    for (MyVar *var : fixingCandidates) {
      int index = candidate->addBranchingVar(var);
      node->setLb(index, 1);
    }
    for (MyVar *var : otherFixingCandidates) {
      int index = candidate->addBranchingVar(var);
      node->setLb(index, 1);
    }
    for (MyVar *var : integerFixingCandidates) {
      int index = candidate->addBranchingVar(var);
      node->setLb(index, 1);
    }

    // Find the column to deactivate
    for (MyVar *var : pModel_->getActiveColumns()) {
      if (var->getUB() == 0)
        continue;
      for (PColumn nodeCol : columns) {
        // do not deactivate if activePat:
        // 1/ applied to a different nurse
        if (Column::nurseNum(var) != nodeCol->nurseNum())
          continue;
        PColumn pCol = pMaster_->getPColumn(var);
        // 2/ will be used (==nodeCol)
        // 3/ is disjoint with nodeCol
        if (pCol->equals(nodeCol) ||
            pCol->isDisjointWith(nodeCol, false))
          continue;

        // add the variable to the candidate
        int index = candidate->addBranchingVar(var);

        // check if the shift is present in shifts
        // set the UB to 0 for the non-possible rotations
        node->setUb(index, 0);

        break;
      }
    }

    /* update tree */
    // append fixing candidates
    auto *can = dynamic_cast<RCBranchingCandidate*>(candidate.get());
    can->addColumnsNode(nodeIndex, columns);
  }

  return true;
}

bool DiveBranchingRule::branching_candidates(
    int nCandidates,
    std::vector<MyPBranchingCandidate> *candidates) {
  // test each branching decision and stop if enough candidate has been found
  for (auto f : branchFunctions_) {
    // invoke the branch function
    ((*this).*f)(nCandidates - candidates->size(), candidates);
    if (candidates->size() >= nCandidates)
      return true;
  }
  return !candidates->empty();
}

//-----------------------------------------------------------------------------
// Branch on a number of nurses working on a shift with for given skill
//-----------------------------------------------------------------------------

void DiveBranchingRule::branchOnNumberNurses(
    int nCandidates, std::vector<MyPBranchingCandidate> *candidates) {
  // try to find a fractional coverage which is surrounded by
  // fractional coverage
  // count only to a certain ratio (discount) of the other variables
  // Need a certain score to branch:
  // if minScore > .5, the fractional coverage is surrounded by other fractional
  double discount = .7, minScore = pModel_->epsilon();
  const vector4D<MyVar*> &skillsAllocVars = pMaster_->getSkillsAllocVars();
  const vector3D<MyCons*> &optDemandCons = pMaster_->getOptDemandCons();

  // find best variable on which to branch: the closest to .5 and surrounded
  // by fractional on the same shift: the day before and after for any skill
  // (as most of the time a transfer can be operated)
  // 1 - compute daily demand
  int nbDays = pMaster_->nDays(), nbShifts = pMaster_->nShifts();
  vector2D<double> dailyShiftScores;
  Tools::initVector2D(&dailyShiftScores, nbShifts - 1, nbDays, .0);
  for (int k = 0; k < nbDays; ++k)
    for (int s = 1; s < nbShifts; ++s) {
      double &dailyScore = dailyShiftScores[s - 1][k];
      double v = pModel_->getVarValue(skillsAllocVars[k][s]);
      double frac = std::min(ceil(v) - v, v - floor(v));
      dailyScore += frac;
    }
  // 2 - find best one
  Tools::FixedSizeList<Score> bestScores(nCandidates);
  for (int s = 1; s < nbShifts; ++s) {
    const auto &dailyScores = dailyShiftScores[s - 1];
    double dailyScore = 0;  // score on the day
    for (int k = 0; k < nbDays; ++k) {
      // update prevScore with previous dailyScore
      double prevScore = dailyScore;
      // score on current day and next one
      dailyScore = dailyScores[k];
      // if 0 -> no fractional variable -> go to next day
      if (dailyScore < pModel_->epsilon()) continue;
      double nextScore = k < nbDays - 1 ? dailyScores[k + 1] : 0;
      Score sc0(dailyScore + discount * (prevScore + nextScore));
      Score *sc = addScore(sc0, &bestScores);
      if (sc) {
        sc->k = k;
        sc->s = s;
      }
    }
  }

#ifdef NS_DEBUG
  // if did not find enough fractional coverage
  if (bestScores.list().empty()) {
    std::cout << pMaster_->currentSolToString() << std::endl;
    std::cout << pMaster_->allocationToString() << std::endl;
    std::cout << pMaster_->coverageToString() << std::endl;
  } else {
    std::cout << pMaster_->coverageToString() << std::endl;
  }
#endif

  // build children
  for (const Score& sc : bestScores.list()) {
    auto candidate =
        std::make_shared<RCBranchingCandidate>(pTree_->getCurrentNode());
    int indexFloor = candidate->createNewChild(),
        indexCeil = candidate->createNewChild();
    const MyPBranchingNode &floorNode = candidate->getChild(indexFloor),
        &ceilNode = candidate->getChild(indexCeil);

    // add the cut
    vector<MyVar *> vars;
    const auto allocVars = skillsAllocVars[sc.k][sc.s];
    for (int sk = 0; sk < pMaster_->pScenario()->nSkills(); ++sk)
      for (MyVar *v : allocVars[sk])
        if (v) vars.push_back(v);
    vector<double> coeffs(vars.size(), 1);

    // create a new cons
    char name[50];
    snprintf(name, sizeof(name), "CoverageCons_%d_%d", sc.k, sc.s);
    MyCons *cut;
    pModel_->createCutLinear(
        &cut, name, 0, pMaster_->nNurses(), vars, coeffs);

    /* update candidate */
    int index = candidate->addNewBranchingCons(cut);
    double v = pModel_->getVarValue(skillsAllocVars[sc.k][sc.s]);
    floorNode->setRhs(index, floor(v));
    ceilNode->setLhs(index, ceil(v));

    // Push back new nodes in the tree
    candidate->addDemandNode(indexFloor, name, cut->getLhs(), floor(v));
    candidate->addDemandNode(indexCeil, name, ceil(v), cut->getRhs());

    // Here : random swap to decide the order of the siblings
    randomSwapLastChildrenIfEnable(candidate);

    candidates->push_back(std::move(candidate));
  }
}

//-----------------------------------------------------------------------------
// Branch on optimal demand variable
//-----------------------------------------------------------------------------
void DiveBranchingRule::branchOnOptDemand(
    int nCandidates, std::vector<MyPBranchingCandidate> *candidates) {
  // find the variable closest to .5
  Tools::FixedSizeList<Score> bestScores(nCandidates);
  const vector3D<MyVar *> optDemandVars = pMaster_->getOptDemandVars();

  int nbDays = pMaster_->nDays(), nbShifts = pMaster_->nShifts();
  for (int k = 0; k < nbDays; ++k)
    for (int s = 1; s < nbShifts; ++s) {
      for (int sk = 0; sk < pMaster_->pScenario()->nSkills(); ++sk) {
        double v = pModel_->getVarValue(optDemandVars[k][s][sk]);
        Score sc0(std::min(ceil(v) - v, v - floor(v)));
        Score *sc = addScore(sc0, &bestScores);
        if (sc)
          sc->var = optDemandVars[k][s][sk];
      }
    }

#ifdef NS_DEBUG
  if (bestScores.list().empty()) {
    std::cout << pMaster_->currentSolToString() << std::endl;
    std::cout << pMaster_->allocationToString() << std::endl;
    std::cout << pMaster_->coverageToString() << std::endl;
  }
#endif

  // build children
  for (const Score& sc : bestScores.list()) {
    auto candidate =
        std::make_shared<RCBranchingCandidate>(pTree_->getCurrentNode());
    int indexFloor = candidate->createNewChild(),
        indexCeil = candidate->createNewChild();
    const MyPBranchingNode &floorNode = candidate->getChild(indexFloor),
        &ceilNode = candidate->getChild(indexCeil);

    // add the variable to the candidate
    MyVar *bestVar = sc.var;
    int index = candidate->addBranchingVar(bestVar);
    double v = pModel_->getVarValue(bestVar);
    floorNode->setUb(index, floor(v));
    ceilNode->setLb(index, ceil(v));

    // Push back new nodes in the tree
    candidate->addVarNode(indexFloor, bestVar, bestVar->getLB(), floor(v));
    candidate->addVarNode(indexCeil, bestVar, ceil(v), bestVar->getUB());

    // Here : random swap to decide the order of the siblings
    randomSwapLastChildrenIfEnable(candidate);

    candidates->push_back(std::move(candidate));
  }
}

//-----------------------------------------------------------------------------
// Branch on a set of resting arcs
//-----------------------------------------------------------------------------
void DiveBranchingRule::branchOnRestDay(
    int nCandidates, std::vector<MyPBranchingCandidate> *candidates) {
  // find the variable closest to .5
  Tools::FixedSizeList<Score> bestScores(nCandidates);

  std::vector<double> baseScores(pMaster_->nNurses(), 0);
  if (pModel_->getParameters().branchBaseScore_ > 1e-9)
    baseScores =
        computeNurseBaseScore(pModel_->getParameters().branchBaseScore_);

  for (PLiveNurse pNurse : pMaster_->pLiveNurses()) {
//    // Several columns are kept -> so could remain fractional
//    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
//    if (pMaster_->isFixNurse(pNurse->num_)) continue;

    for (int k = 0; k < pMaster_->nDays(); ++k) {
      // ROLLING/LNS: do not branch on a relaxed day
      if (pMaster_->isRelaxDay(k))
        continue;

      double value =
          pModel_->getVarValue(pMaster_->getRestVarsPerDay(pNurse, k));

      // The value has to be not integer
      if (value < pModel_->epsilon() || value > 1 - pModel_->epsilon())
        continue;

      Score sc0(scoreFunc_->score(
          pNurse, k, {0}, {value}, baseScores[pNurse->num_]));
      Score *sc = addScore(sc0, &bestScores);
      if (sc) {
        sc->k = k;
        sc->pNurse = pNurse;
      }
    }
  }

  for (const Score& sc : bestScores.list()) {
    auto candidate =
        std::make_shared<RCBranchingCandidate>(pTree_->getCurrentNode());
    int index1 = candidate->createNewChild(),
        index2 = candidate->createNewChild();
    const MyPBranchingNode &restNode = candidate->getChild(index1),
        &workNode = candidate->getChild(index2);
    buildRestNodesCut(
        candidate, sc.pNurse, sc.k, true, restNode, workNode);
    deactivateColumns(
        candidate, sc.pNurse->num_, sc.k, {0}, workNode, restNode);

    // create potential nodes for the tree
    candidate->addRestNode(index1, sc.pNurse, sc.k, true);
    candidate->addRestNode(index2, sc.pNurse, sc.k, false);

    // Here : random swap to decide the order of the siblings
    randomSwapLastChildrenIfEnable(candidate);

    candidates->push_back(std::move(candidate));
  }
}


//-----------------------------------------------------------------------------
// Branch on a shift
//-----------------------------------------------------------------------------

void DiveBranchingRule::branchOnShifts(
    int nCandidates, std::vector<MyPBranchingCandidate> *candidates) {
  // find the variable closest to .5
  Tools::FixedSizeList<Score> bestScores(nCandidates);

  std::vector<double> baseScores(pMaster_->nNurses(), 0);
  if (pModel_->getParameters().branchBaseScore_ > 1e-9)
    baseScores =
        computeNurseBaseScore(pModel_->getParameters().branchBaseScore_);

  // compute the solution for each nurse, day, shift
  vector3D<double> fractionalRoster = pMaster_->fractionalRoster();

  // search for the best branching decision (set of shifts the closest to .5)
  for (PLiveNurse pNurse : pMaster_->pLiveNurses()) {
//    // Several columns are kept -> so could remain fractional
//    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
//    if (pMaster_->isFixNurse(pNurse->num_)) continue;

    for (int k = 0; k < pMaster_->nDays(); ++k) {
      // ROLLING/LNS: do not branch on a relaxed day
      if (pMaster_->isRelaxDay(k))
        continue;

      // Get the fraction that is spent on the resting shift and recover an
      // indexation of the working shifts from 1 to NbShifts_
      vector<std::pair<int, double>> fractionalNurseDay;
      bool isFractional = false;
      for (int s = 0; s < pMaster_->nShifts(); ++s) {
        double val = fractionalRoster[pNurse->num_][k][s];
        // check if solution is fractional
        if (pModel_->epsilon() < val && val < 1 - pModel_->epsilon())
          isFractional = true;
        fractionalNurseDay.emplace_back(s, -val);
      }
      // if not fractional, continue
      if (!isFractional)
        continue;

      // sort the scores
      sort(fractionalNurseDay.begin(),
           fractionalNurseDay.end(),
           Tools::compareObject<double>);

      // look for the balance between the shifts:
      // the sum of the weight in currentShifts should be as close as possible
      // to 0.5.
      // scoreNode1 is the sum of the weights of the shifts in currentShifts
      // nbShifts is the number of shifts in currentshifts
      // scoreNode2, nbSbhifts2 are the same for the shifts not in currentShifts
      double scoreNode1 = 0, scoreNode2 = 0;
      int nbshifts1 = 0, nbshifts2 = 0;
      vector<int> currentShifts;
      vector<double> values;
      for (const pair<int, double> &p : fractionalNurseDay) {
        if (p.second < pModel_->epsilon()) {
          if (nbshifts1 >= nbshifts2) {
            ++nbshifts2;
          } else {
            ++nbshifts1;
            values.push_back(-p.second);
            currentShifts.push_back(p.first);
          }
        } else if (scoreNode1 > scoreNode2) {
          ++nbshifts2;
          scoreNode2 += p.second;
        } else {
          ++nbshifts1;
          scoreNode1 += p.second;
          values.push_back(-p.second);
          currentShifts.push_back(p.first);
        }
      }

      // look for the value that is the most balanced scores
      Score sc0(scoreFunc_->score(
          pNurse, k, currentShifts, values, baseScores[pNurse->num_]));
      Score *sc = addScore(sc0, &bestScores);
      if (sc) {
        sc->k = k;
        sc->pNurse = pNurse;
        sc->forbiddenShifts = currentShifts;
      }
    }
  }

  for (const Score& sc : bestScores.list()) {
    // compute the forbidden shifts
    vector<int> complementaryShifts;
    bool canRest = false;
    for (int i = 0; i < pMaster_->nShifts(); ++i) {
      if (find(sc.forbiddenShifts.begin(), sc.forbiddenShifts.end(), i)
          == sc.forbiddenShifts.end()) {
        complementaryShifts.push_back(i);
        // if the forbidden shift 0 (rest) is in the complementary,
        // it is possible to rest for the first node
        if (i == 0) canRest = true;
      }
    }

    // create and build the nodes
    auto candidate =
        std::make_shared<RCBranchingCandidate>(pTree_->getCurrentNode());
    int index1 = candidate->createNewChild(),
        index2 = candidate->createNewChild();
    const MyPBranchingNode &restNode = candidate->getChild(index1),
        &workNode = candidate->getChild(index2);
    buildRestNodesCut(candidate,
                      sc.pNurse,
                      sc.k,
                      false,
                      restNode,
                      workNode);

    // Find the columns to deactivate
    deactivateColumns(candidate,
                      sc.pNurse->num_,
                      sc.k,
                      sc.forbiddenShifts,
                      // principal node:
                      // rest node if can rest, work node otherwise
                      canRest ? restNode : workNode,
                      canRest ? workNode : restNode);  // complementary node

    // add the node to the tree
    // rest node : some work shifts are forbidden
    // If canRest is true, just work shift are in the forbiddenShifts
    candidate->addShiftNode(index1, sc.pNurse, sc.k, false,
                            canRest ? sc.forbiddenShifts : complementaryShifts);
    // work node : rest shift is forbidden
    // If canRest is true, the rest shift is in the complementaryShifts
    candidate->addShiftNode(index2, sc.pNurse, sc.k, true,
                            canRest ? complementaryShifts : sc.forbiddenShifts);

    // Here : random swap to decide the order of the siblings
    randomSwapLastChildrenIfEnable(candidate);

    candidates->push_back(std::move(candidate));
  }
}

void DiveBranchingRule::buildRestNodesCut(
    const MyPBranchingCandidate &candidate,
    PLiveNurse pNurse,
    int day,
    bool forceRest,
    const MyPBranchingNode &restNode,
    const MyPBranchingNode &workNode) const {
  // creating the branching cut
  char name[50];
  snprintf(name, sizeof(name), "RestBranchingCons_N%d_%d", pNurse->num_, day);
  vector<MyVar *> restingArcs = pMaster_->getRestVarsPerDay(pNurse, day);
  vector<double> coeffs;
  Tools::initVector(&coeffs, restingArcs.size(), 1.0);
  // create a new cons
  MyCons *cons;
  pModel_->createCutLinear(&cons, name, 0, 1, restingArcs, coeffs);

  /* update candidate */
  int index = candidate->addNewBranchingCons(cons);
  if (forceRest) {
    restNode->setLhs(index, 1);
    restNode->setRhs(index, 1);
  }
  workNode->setLhs(index, 0);
  workNode->setRhs(index, 0);
}

void DiveBranchingRule::deactivateColumns(
    const MyPBranchingCandidate &candidate,
    int nurseNum,
    int day,
    std::vector<int> forbiddenShifts,
    const MyPBranchingNode &forbiddenNode,
    const MyPBranchingNode &complementaryNode) const {
  // Find the columns to deactivate
  for (MyVar *var : pModel_->getActiveColumns()) {
    if (var->getUB() < pModel_->epsilon())
      continue;

    if (Column::nurseNum(var) == nurseNum
        && Column::firstDay(var) <= day
        && day <= Column::lastDay(var)) {
      // add the variable to the candidate
      int index = candidate->addBranchingVar(var);
      PColumn pat = pMaster_->getPColumn(var);

      // check if the shift is present in shifts
      // remove the non-possible rotations
      if (find(forbiddenShifts.begin(),
               forbiddenShifts.end(),
               pat->shift(day)) != forbiddenShifts.end())
        forbiddenNode->setUb(index, 0);
      else
        complementaryNode->setUb(index, 0);
    }
  }
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

vector<MyVar *> DiveBranchingRule::chooseColumns(
    const vector<pair<MyVar *, double>> &candidates,
    vector<PColumn> *columns,
    double *maxValue,
    const ColumnsComparator &comparator) {
  vector<MyVar *> fixingCandidates;
  for (const pair<MyVar *, double> &p : candidates) {
    if (*maxValue < p.second) continue;

    PColumn pat1 = pMaster_->getPColumn(p.first);
    // check if this rotation is totally disjoint with all the others
    // if not should be disjoint for the shift and the nurse
    bool isDisjoint = true;
    for (PColumn pat2 : *columns)
      if (!comparator.is_disjoint(pat1, pat2)) {
        isDisjoint = false;
        break;
      }
    if (isDisjoint) {
      fixingCandidates.push_back(p.first);
      *maxValue -= p.second;
      columns->push_back(pat1);
    }
  }

  return fixingCandidates;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool DiveBranchingRule::compareColumnCloseToInt(
    const std::pair<MyVar *, double> &obj1,
    const std::pair<MyVar *, double> &obj2) {
  double frac1 = obj1.second - floor(obj1.second),
      frac2 = obj2.second - floor(obj2.second);
  double closeToInt1 = 1 - frac1, closeToInt2 = 1 - frac2;
  return (closeToInt1 < closeToInt2);
}


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool DiveBranchingRule::compareColumnCloseTo5(
    const std::pair<MyVar *, double> &obj1,
    const std::pair<MyVar *, double> &obj2) {
  double frac1 = obj1.second - floor(obj1.second),
      frac2 = obj2.second - floor(obj2.second);
  double closeTo5_1 = fabs(0.5 - frac1), closeTo5_2 = fabs(0.5 - frac2);
  return (closeTo5_1 < closeTo5_2);
}

MasterProblem *DiveBranchingRule::getMaster() const {
  return pMaster_;
}

Modeler *DiveBranchingRule::getModel() const {
  return pModel_;
}

// compute a base score for each nurse in order to advantage certain nurses
// in the branching selection process
std::vector<double> DiveBranchingRule::computeNurseBaseScore(double coeff) {
  // compute the cost per nurse as well as the number of active columns
  std::vector<double> costs(pMaster_->nNurses(), 0),
                      nCols(pMaster_->nNurses(), 0);
  for (MyVar *var : pModel_->getActiveColumns()) {
    double v = pModel_->getVarValue(var);
    if (v < pModel_->epsilon()) continue;

    int nurseId = Column::nurseNum(var);
    costs[nurseId] += var->getCost() * v;
    nCols[nurseId]++;
  }

  // compute the max of each vector
  int maxCost = 0;
  int maxCol = 0;
  for (int n=0; n < pMaster_->nNurses(); ++n) {
    if (costs[n] > maxCost) maxCost = costs[n];
    if (nCols[n] > maxCol) maxCol = nCols[n];
  }

  // compute a score based on the normalized value
  std::vector<double> scores;
  scores.reserve(pMaster_->nNurses());
  for (int n=0; n < pMaster_->nNurses(); ++n) {
    scores.push_back(
        coeff * ((1+ costs[n]) / (1 + maxCost)) * (nCols[n] / maxCol) / 2);
//    std::cout << "N" << n << ": " << scores.back() << "(" << costs[n]
//              << ", " << nCols[n] << ")" << std::endl;
  }
  return scores;
}
