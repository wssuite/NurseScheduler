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

RestTree::RestTree(PScenario pScenario, PDemand pDemand, double epsilon) :
    MyTree(epsilon),
    pScenario_(pScenario),
    pDemand_(pDemand),
    statsRestByDay_(pDemand_->nDays_),
    statsWorkByDay_(pDemand_->nDays_),
    statsRestByNurse_(pScenario->nbNurses()),
    statsWorkByNurse_(pScenario->nbNurses()) {}

// Forbid the shifts (that are already fixed by the previous branching decision)
// in the generation  of new column.
// This method is useful for the RC pricer
void RestTree::addForbiddenShifts(PLiveNurse pNurse,
                                  set<pair<int, int> > *forbidenShifts) {
  MyNode *node = currentNode_;
  vector<MyVar *> arcs;
  while (node->pParent_) {
    MyNode *currentNode = node;
    node = node->pParent_;
    // If on a rest node, we forbid shifts only if rest is forced
    // In that case, no column can include a shift on the fixed resting day
    //
    RestNode *restNode = dynamic_cast<RestNode *>(currentNode);
    if (restNode) {
      if (restNode->pNurse_ == pNurse) {
        if (restNode->rest_) {
          for (int i = 1; i < pNurse->pScenario_->nbShifts_; ++i)
            forbidenShifts->insert(pair<int, int>(restNode->day_, i));
        } else {
          forbidenShifts->insert(pair<int, int>(restNode->day_, 0));
        }
      }
      continue;
    }

    // If on a shift node, it is natural to forbid all the work shifts
    // in RC pricing.
    ShiftNode *shiftNode = dynamic_cast<ShiftNode *>(currentNode);
    if (shiftNode) {
      if (shiftNode->pNurse_ == pNurse)
        for (int s : shiftNode->forbiddenShifts_)
          forbidenShifts->insert(pair<int, int>(shiftNode->day_, s));
      continue;
    }

    // If in a column node , forbid the shifts depending on the pattern
    ColumnsNode *columnsNode = dynamic_cast<ColumnsNode *>(currentNode);
    if (columnsNode) {
      for (PPattern pat : columnsNode->patterns_)
        if (pat->nurseNum_ == pNurse->num_)
          pat->addForbiddenShifts(forbidenShifts,
                                  pScenario_->nbShifts(),
                                  pDemand_);
      continue;
    }

    // If on a shift demand node or var node, nothing to forbid
    if (dynamic_cast<CoverageNode *>(currentNode) ||
        dynamic_cast<VarNode *>(currentNode))
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
  int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nDays_;

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
    rep << "|  " << Tools::intToDay(
        firstDay + day).at(0) << "  ";
  rep << "|" << std::endl;
  rep << writeOneStat("Rest", statsRestByDay_);
  rep << writeOneStat("Work", statsWorkByDay_);
  rep << "-------------------------------------" << std::endl;

  rep << std::endl;
  rep << "Stats by Nurse" << std::endl;
  rep << "\t\t  ";
  for (int n = 0; n < nbNurses; n++) {
    if (n < 10) rep << "|" << pScenario_->theNurses_[n]->name_ << " ";
    else
      rep << "|" << pScenario_->theNurses_[n]->name_;
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
  statsWorkByDay_.clear();
  statsRestByNurse_.clear();
  statsWorkByNurse_.clear();
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

double ScoreVarCloseHalf::score(PLiveNurse pNurse,
                                int day,
                                const std::vector<int> &shifts,
                                const std::vector<double> &values) const {
  // choose the variable closest to .5
  // look for the value the closest to .5
  double value = 0;
  for (double v : values) value += v;
  double sc = std::min(value - floor(value),
                       ceil(value) - value);
  if (Tools::isWeekend(day)) sc += advantage;
  return sc;
}

double ScoreVarBestExpectedLBImprovement::score(
    PLiveNurse pNurse,
    int day,
    const std::vector<int> &shifts,
    const std::vector<double> &values) const {
  if (shifts.size() != values.size())
    Tools::throwError("Shifts and  Values vectors do not have the same size.");
  // choose the variable that yields the biggest expected improvement
  // compute reduced cost associated to the shift
  PDualCosts pDCosts = pRule->getPDualCosts()[pNurse->num_];
  // look at the worst dual costs for the shifts
  double ceilCost = DBL_MAX, floorCost = 0;
  for (int i = 0; i < shifts.size(); i++) {
    int s = shifts[i];
    double v = values[i];
    double dc = pDCosts->constant() + pDCosts->workedDayShiftCost(day, s);
    // add worked weekend dual cost if working, otherwise remove it
    if (Tools::isWeekend(day))
      dc += pDCosts->workedWeekendCost() * (s ? 1 : -1);
    // use absolute value
    dc = std::fabs(dc);
    // if use shift -> not using the other shifts -> worst case
    double ceilDC = dc * (ceil(v) - v);
    if (ceilDC < ceilCost) ceilCost = ceilDC;
    // if not using shift -> all the shifts won't be used -> just add them
    floorCost += dc * (v - floor(v));
  }

  // return the worst (closest to 0)
  return std::min(abs(floorCost), abs(ceilCost));
}

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem *master,
                                     RestTree *tree,
                                     const char *name) :
    MyBranchingRule(name),
    pMaster_(master),
    tree_(tree),
    pModel_(master->getModel()) {
  // either ScoreVarCloseHalf or ScoreVarBestExpectedLBImprovement
  scoreFunc_ = std::unique_ptr<ScoreVar>(
      new ScoreVarCloseHalf(this));
}

// add all good candidates
bool DiveBranchingRule::column_candidates(MyBranchingCandidate *candidate) {
  // look for fractional columns
  // Fix all column above BRANCH_LB
  // search the good candidates
  vector<pair<MyVar *, double> > candidates;
  vector<MyVar *> integerFixingCandidates;
  std::vector<MyVar *> fixingCandidates;
  std::vector<MyVar *> otherFixingCandidates;
  vector<PPattern> patterns;

  if ((pModel_->getParameters().branchColumnUntilValue_
      && pModel_->getParameters().branchColumnDisjoint_)
      || (!pModel_->getParameters().branchColumnUntilValue_
          && !pModel_->getParameters().branchColumnDisjoint_))
    Tools::throwError("DiveBranchingRule::branch_on_column: "
                      "the branching rule is not chosen properly!");

  for (MyVar *var : pModel_->getActiveColumns()) {
    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
    if (pMaster_->isFixNurse(var->getNurseNum())) continue;

    // ROLLING/LNS: do not branch on relaxed columns
    // check if all days are relaxed
    bool relaxed = true;
    for (int k = Pattern::firstDay(var); k <= Pattern::lastDay(var); ++k)
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
      patterns.emplace_back(pMaster_->getPattern(var));
      continue;
    }
    if (value > pModel_->epsilon())
      candidates.push_back(pair<MyVar *, double>(var, 1 - value));
  }

  stable_sort(
      candidates.begin(), candidates.end(), Tools::compareObject<MyVar *>);

  // Try the branching that mimics the heuristic
  if (pModel_->getParameters().branchColumnUntilValue_) {
    if (candidates.size() > 0) {
      double valueLeft = 1 - pModel_->epsilon();
      for (pair<MyVar *, double> &p : candidates) {
        if (p.second > valueLeft)
          break;
        else
          valueLeft -= p.second;
        fixingCandidates.push_back(p.first);
      }
    }
  }

  if (pModel_->getParameters().branchColumnDisjoint_) {
    double valueMax = pMaster_->getBranchColumnValueMax();
    DayDisjointComparator comp1 = DayDisjointComparator();
    fixingCandidates = chooseColumns(candidates,
                                     &patterns,
                                     &valueMax,
                                     comp1);

    ShiftDisjointComparator comp2 = ShiftDisjointComparator();
    otherFixingCandidates = chooseColumns(candidates,
                                          &patterns,
                                          &valueMax,
                                          comp2);
  }

  if (fixingCandidates.empty() && otherFixingCandidates.empty()
      && integerFixingCandidates.empty())
    return false;

  /* update candidate */
  MyBranchingNode &node = candidate->getChild(candidate->createNewChild());
  for (MyVar *var : fixingCandidates) {
    int index = candidate->addBranchingVar(var);
    node.setLb(index, 1);
  }
  for (MyVar *var : otherFixingCandidates) {
    int index = candidate->addBranchingVar(var);
    node.setLb(index, 1);
  }
  for (MyVar *var : integerFixingCandidates) {
    int index = candidate->addBranchingVar(var);
    node.setLb(index, 1);
  }

  // Find the column to deactivate
  for (MyVar *var : pModel_->getActiveColumns()) {
    if (var->getUB() == 0)
      continue;
    for (PPattern nodePat : patterns) {
      // do not deactivate if activePat:
      // 1/ applied to a different nurse
      if (Pattern::nurseNum(var) != nodePat->nurseNum_)
        continue;
      PPattern pat = pMaster_->getPattern(var);
      // 2/ will be used (==nodePat)
      // 3/ is disjoint with nodePat
      if (pat->equals(nodePat) ||
          pat->isDisjointWith(nodePat, false))
        continue;

      // add the variable to the candidate
      int index = candidate->addBranchingVar(var);

      // check if the shift is present in shifts
      // set the UB to 0 for the non-possible rotations
      node.setUb(index, 0);

      break;
    }
  }

  /* update tree */
  // append fixing candidates
  tree_->pushBackNewColumnsNode(patterns);

  return true;
}

bool DiveBranchingRule::branching_candidates(MyBranchingCandidate *candidate) {
  // store the dual costs, can be useful when deciding branching decisions
  pDualCosts_.clear();
  pDualCosts_.reserve(pMaster_->nNurses());
  for (PLiveNurse pNurse : pMaster_->liveNurses())
    pDualCosts_.emplace_back(pMaster_->buildDualCosts(pNurse));

  // test each branching decision and stop if a candidate has been found
  int nChildren = candidate->getChildren().size();
  for (auto f : branchFunctions_) {
    ((*this).*f)(candidate);  // invoke the branch function
    if (candidate->getChildren().size() > nChildren)
      return true;
  }
  return false;
}

//-----------------------------------------------------------------------------
// Branch on a number of nurses working on a shift with for given skill
//-----------------------------------------------------------------------------

void DiveBranchingRule::branchOnNumberNurses(MyBranchingCandidate *candidate) {
  // try to find a fractional coverage which is surrounded by
  // fractional coverage
  // count only to 90% the other variables
  // Need a certain score to branch:
  // if minScore > .5, the fractional coverage is surrounded by other fractional
  double discount = .9, minScore = .6;
  const vector4D<MyVar *> &skillsAllocVars = pMaster_->getSkillsAllocVars();

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
      for (int sk = 0; sk < pMaster_->pScenario()->nbSkills(); ++sk) {
        double v = pModel_->getVarValue(skillsAllocVars[k][s - 1][sk]);
        dailyScore += std::min(ceil(v) - v, v - floor(v));
      }
    }
  // 2 - find best one
  int bestDay = -1, bestShift = -1, bestSkill = -1;
  double bestScore = 0;
  for (int s = 1; s < nbShifts; ++s) {
    double dailyScore = 0;  // score on the day
    for (int k = 0; k < nbDays; ++k) {
      // update prevScore with previous dailyScore
      double prevScore = dailyScore;
      // score on current day and next one
      dailyScore = dailyShiftScores[s - 1][k];
      // if 0 -> no fractional variable -> go to next day
      if (dailyScore == 0) continue;
      double nextScore = k < nbDays - 1 ? dailyShiftScores[s - 1][k + 1] : 0;
      double score = discount * (prevScore + dailyScore + nextScore);
//      for (int sk = 0; sk < pMaster_->getScenario()->nbSkills(); ++sk) {
//        // count at 100% current day and skill
//        double v = pModel_->getVarValue(skillsAllocVars[k][s - 1][sk]);
//        double sc = (1 - discount) * std::min(ceil(v) - v, v - floor(v));
      // check if best
      if (score > bestScore) {
        bestScore = score;
        bestDay = k;
        bestShift = s;
//          bestSkill = sk;
//        }
      }
    }
  }

  // if did not find enough fractional coverage -> stop
  if (bestScore < minScore) return;

  std::cout << pMaster_->coverageToString() << std::endl;

  // otherwise, build children
  int indexFloor = candidate->createNewChild(),
      indexCeil = candidate->createNewChild();
  MyBranchingNode &floorNode = candidate->getChild(indexFloor),
      &ceilNode = candidate->getChild(indexCeil);

  // add the cut
  vector<MyVar *> vars;
  vector<double> coeffs;
  for (int sk = 0; sk < pMaster_->pScenario()->nbSkills(); ++sk)
    for (int p : pMaster_->getPositionsForSkill(sk)) {
      vars.push_back(skillsAllocVars[bestDay][bestShift - 1][sk][p]);
      coeffs.push_back(1);
    }

  // create a new cons
  char name[50];
  snprintf(name, sizeof(name), "CoverageCons_%d_%d_%d",
           bestDay, bestShift, bestSkill);
  MyCons *cut;
  pModel_->createCutLinear(
      &cut, name, 0, pMaster_->nNurses(), vars, coeffs);

  /* update candidate */
  int index = candidate->addNewBranchingCons(cut);
  double v = pModel_->getVarValue(
      skillsAllocVars[bestDay][bestShift - 1]);  // [bestSkill]
  floorNode.setRhs(index, floor(v));
  ceilNode.setLhs(index, ceil(v));

  // Push back new nodes in the tree
  tree_->pushBackNewDemandNode(name, cut->getLhs(), floor(v));
  tree_->pushBackNewDemandNode(name, ceil(v), cut->getRhs());

  // Here : random swap to decide the order of the siblings
  randomSwapLastChildren(candidate);
}

//-----------------------------------------------------------------------------
// Branch on optimal demand variable
//-----------------------------------------------------------------------------
void DiveBranchingRule::branchOnOptDemand(MyBranchingCandidate *candidate) {
  // find the variable closest to o.5
  MyVar *bestVar = nullptr;
  double bestScore = 0;
  const vector3D<MyVar *> optDemandVars = pMaster_->getOptDemandVars();

  int nbDays = pMaster_->nDays(), nbShifts = pMaster_->nShifts();
  for (int k = 0; k < nbDays; ++k)
    for (int s = 1; s < nbShifts; ++s) {
      for (int sk = 0; sk < pMaster_->pScenario()->nbSkills(); ++sk) {
        double v = pModel_->getVarValue(optDemandVars[k][s - 1][sk]);
        double score = std::min(ceil(v) - v, v - floor(v));
        if (score > bestScore + pModel_->epsilon()) {
          bestScore = score;
          bestVar = optDemandVars[k][s - 1][sk];
        }
      }
    }

  if (bestVar == nullptr) return;

  // otherwise, build children
  int indexFloor = candidate->createNewChild(),
      indexCeil = candidate->createNewChild();
  MyBranchingNode &floorNode = candidate->getChild(indexFloor),
      &ceilNode = candidate->getChild(indexCeil);

  // add the variable to the candidate
  int index = candidate->addBranchingVar(bestVar);
  double v = pModel_->getVarValue(bestVar);
  floorNode.setUb(index, floor(v));
  ceilNode.setLb(index, ceil(v));

  // Push back new nodes in the tree
  tree_->pushBackNewVarNode(bestVar, bestVar->getLB(), floor(v));
  tree_->pushBackNewVarNode(bestVar, ceil(v), bestVar->getUB());

  // Here : random swap to decide the order of the siblings
  randomSwapLastChildren(candidate);
}

//-----------------------------------------------------------------------------
// Branch on a set of resting arcs
//-----------------------------------------------------------------------------
void DiveBranchingRule::branchOnRestDay(MyBranchingCandidate *candidate) {
  // set of rest variable closest to .5
  int bestDay = -1;
  PLiveNurse pBestNurse(nullptr);
  double bestScore = -DBL_MAX;

  for (PLiveNurse pNurse : pMaster_->liveNurses()) {
    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
    if (pMaster_->isFixNurse(pNurse->num_)) continue;

    for (int k = 0; k < pMaster_->nDays(); ++k) {
      // ROLLING/LNS: do not branch on a relaxed day
      if (pMaster_->isRelaxDay(k))
        continue;

      double value =
          pModel_->getVarValue(pMaster_->getRestVarsPerDay(pNurse, k));

      // The value has to be not integer
      if (value < pModel_->epsilon() || value > 1 - pModel_->epsilon())
        continue;

      double sc = scoreFunc_->score(pNurse, k, {0}, {value});
      if (sc > bestScore) {
        bestDay = k;
        pBestNurse = pNurse;
        bestScore = sc;
      }
    }
  }

  if (pBestNurse == nullptr) return;

  int index1 = candidate->createNewChild(),
      index2 = candidate->createNewChild();
  MyBranchingNode &restNode = candidate->getChild(index1),
      &workNode = candidate->getChild(index2);
  buildRestNodesCut(
      candidate, pBestNurse, bestDay, true, &restNode, &workNode);
  deactivateColumns(
      candidate, pBestNurse->num_, bestDay, {0}, &workNode, &restNode);

  // Push back new nodes in the tree
  tree_->pushBackNewRestNode(pBestNurse, bestDay, true);
  tree_->pushBackNewRestNode(pBestNurse, bestDay, false);

  // Here : random swap to decide the order of the siblings
  randomSwapLastChildren(candidate);
}


//-----------------------------------------------------------------------------
// Branch on a set of resting arcs
//-----------------------------------------------------------------------------

void DiveBranchingRule::branchOnShifts(MyBranchingCandidate *candidate) {
  // set of rest variable closest to .5
  int bestDay = -1;
  PLiveNurse pBestNurse(nullptr);
  double bestScore = -DBL_MAX;
  vector<int> forbiddenShifts;

  // compute the solution for each nurse, day, shift
  vector3D<double> fractionalRoster = pMaster_->fractionalRoster();

  // search for the best branching decision (set of shifts the closest to .5)
  for (PLiveNurse pNurse : pMaster_->liveNurses()) {
    // ROLLING/LNS: do not branch on a nurse whose columns are fixed
    if (pMaster_->isFixNurse(pNurse->num_)) continue;

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
        fractionalNurseDay.push_back({s, -val});
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
            values.push_back(p.second);
            currentShifts.push_back(p.first);
          }
        } else if (scoreNode1 > scoreNode2) {
          ++nbshifts2;
          scoreNode2 += p.second;
        } else {
          ++nbshifts1;
          scoreNode1 += p.second;
          values.push_back(p.second);
          currentShifts.push_back(p.first);
        }
      }

      // look for the value that is the most balanced scores
      double sc = scoreFunc_->score(pNurse, k, currentShifts, values);
      if (sc > bestScore) {
        bestDay = k;
        pBestNurse = pNurse;
        bestScore = sc;
        forbiddenShifts = currentShifts;
      }
    }
  }

  if (pBestNurse != nullptr) {
    // compute the forbidden shifts
    vector<int> complementaryShifts;
    bool canRest = false;
    for (int i = 0; i < pMaster_->nShifts(); ++i) {
      if (find(forbiddenShifts.begin(), forbiddenShifts.end(), i)
          == forbiddenShifts.end()) {
        complementaryShifts.push_back(i);
        // if the forbidden shift 0 (rest) is in the complementary,
        // it is possible to rest for the first node
        if (i == 0) canRest = true;
      }
    }

    // create and build the nodes
    int index1 = candidate->createNewChild(),
        index2 = candidate->createNewChild();
    MyBranchingNode &node1 = candidate->getChild(index1),
        &node2 = candidate->getChild(index2);
    buildRestNodesCut(candidate,
                      pBestNurse,
                      bestDay,
                      false,
                      canRest ? &node1 : &node2,
                      canRest ? &node2 : &node1);

    // Find the columns to deactivate
    deactivateColumns(candidate,
                      pBestNurse->num_,
                      bestDay,
                      forbiddenShifts,
                      &node1,
                      &node2);

    // add the node to the tree
    tree_->pushBackNewShiftNode(pBestNurse, bestDay, !canRest, forbiddenShifts);
    tree_->pushBackNewShiftNode(pBestNurse,
                                bestDay,
                                canRest,
                                complementaryShifts);
  }
}

void DiveBranchingRule::buildRestNodesCut(MyBranchingCandidate *candidate,
                                          PLiveNurse pNurse,
                                          int day,
                                          bool forceRest,
                                          MyBranchingNode *restNode,
                                          MyBranchingNode *workNode) const {
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
    MyBranchingCandidate *candidate,
    int nurseNum,
    int day,
    std::vector<int> forbiddenShifts,
    MyBranchingNode *forbiddenNode,
    MyBranchingNode *complementaryNode) const {
  // Find the columns to deactivate
  for (MyVar *var : pModel_->getActiveColumns()) {
    if (var->getUB() == 0)
      continue;

    if (Pattern::nurseNum(var) == nurseNum
        && Pattern::firstDay(var) <= day
        && Pattern::lastDay(var) >= day) {
      // add the variable to the candidate
      int index = candidate->addBranchingVar(var);
      PPattern pat = pMaster_->getPattern(var);

      // check if the shift is present in shifts
      // set the UB to 0 for the non-possible rotations
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
    vector<PPattern> *patterns,
    double *maxValue,
    const ColumnsComparator &comparator) {
  vector<MyVar *> fixingCandidates;
  for (const pair<MyVar *, double> &p : candidates) {
    if (*maxValue < p.second) continue;

    PPattern pat1 = pMaster_->getPattern(p.first);
    // check if this rotation is totally disjoint with all the others
    // if not should be disjoint for the shift and the nurse
    bool isDisjoint = true;
    for (PPattern pat2 : *patterns)
      if (!comparator.is_disjoint(pat1, pat2)) {
        isDisjoint = false;
        break;
      }
    if (isDisjoint) {
      fixingCandidates.push_back(p.first);
      *maxValue -= p.second;
      patterns->push_back(pat1);
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

const vector<PDualCosts> &DiveBranchingRule::getPDualCosts() const {
  return pDualCosts_;
}
