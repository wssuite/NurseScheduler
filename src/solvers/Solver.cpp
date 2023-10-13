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

#include "Solver.h"

#include <cmath>
#include <algorithm>
#include <random>
#include <utility>

#include "solvers/mp/RotationMP.h"
#include "solvers/mp/RosterMP.h"
#include "solvers/mp/sp/rcspp/RCLabel.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

//-----------------------------------------------------------------------------
//  C l a s s   S t a t N u r s e C t
// The instances of this class gather the status of the constraints that relate
// to the nurses.
//-----------------------------------------------------------------------------

// Constructor and destructor
StatCtNurse::~StatCtNurse() {}

// initialize the statuses
void StatCtNurse::init(int nbDays) {
  nbDays_ = nbDays;
  costTotalDays_ = 0;
  costMissingDays_ = 0;
  costExceedingDays_ = 0;
  costTotalWeekEnds_ = 0;

  // initialize all the cost vectors
  Tools::initVector(&costConsShifts_, nbDays_, .0);
  Tools::initVector(&costConsDays_, nbDays_, .0);
  Tools::initVector(&costConsDaysOff_, nbDays_, .0);
  Tools::initVector(&costPref_, nbDays_, .0);
  Tools::initVector(&costWeekEnd_, nbDays_, .0);

  // initialize the violation vector
  Tools::initVector(&violSuccShifts_, nbDays_, false);
  Tools::initVector(&violSkill_, nbDays_, false);
}

//-----------------------------------------------------------------------------
//  C l a s s   L i v e N u r s e
// A live nurse is a nurse whose characteristics can evolve depending on
// the demand and on the planning that is being built
// They are needed in the solvers to duplicate the static nurses and define new
// attribute that can be modified.
//-----------------------------------------------------------------------------

// Constructor
LiveNurse::LiveNurse(const Nurse &nurse,
                     PScenario pScenario,
                     int nbDays,
                     int firstDay,
                     State *pStateIni,
                     PPreferences pPreferences,
                     int num) :
    Nurse(num, nurse),
    pScenario_(std::move(pScenario)),
    nbDays_(nbDays),
    firstDay_(firstDay),
    nurseNum_(nurse.num_),
    pStateIni_(pStateIni),
    pPreferences_(std::move(pPreferences)),
    pPosition_(nullptr),
    minWorkDaysNoPenaltyConsDays_(-1),
    maxWorkDaysNoPenaltyConsDays_(-1),
    minWorkDaysNoPenaltyTotalDays_(-1),
    maxWorkDaysNoPenaltyTotalDays_(-1),
    minAvgWorkDaysNoPenaltyTotalDays_(-1),
    maxAvgWorkDaysNoPenaltyTotalDays_(-1) {
  roster_.init(firstDay, nbDays, pScenario_->pRestShift());
  statCt_.init(nbDays);

  // check if initial state is valid
  if (pStateIni->pShift_->id >= 0 &&
      pStateIni->pShift_->id < pScenario_->nShifts()
      && isShiftNotAvailNorAlt(pStateIni->pShift_->id))
    Tools::throwError("Invalid initial shift %s, "
                      "as it is unavailable for nurse %s",
                      pStateIni->pShift_->name.c_str(),
                      name_.c_str());

  // initialize the states at each day
  states_.push_back(*pStateIni);
  for (int day = 0; day < nbDays_; day++) {
    State nextState;
    nextState.addDayToState(states_[day], pScenario_->pRestShift());
    states_.push_back(nextState);
  }

  // create new Resources from BaseResource
  for (const auto &pBR : pBaseResources_) {
    auto *pR = dynamic_cast<Resource *>(pBR->clone());
    pResources_.emplace_back(pR);
    auto pSDR = std::dynamic_pointer_cast<SoftTotalShiftDurationResource>(
        pResources_.back());
    if (pSDR) {
      totalShiftDurationResource_ = pSDR;
    } else {
      auto pWR = std::dynamic_pointer_cast<SoftTotalWeekendsResource>(
          pResources_.back());
      if (pWR)
        totalWeekendResource_ = pWR;
    }
  }

  // add preferences
  for (const auto &p : wishes())
    pResources_.push_back(std::make_shared<PreferenceResource>(
        std::make_shared<Day>(p.first), p.second,
        pScenario_->shiftsFactory().pEndShift()));
}

LiveNurse::~LiveNurse() = default;

const std::map<int, Wish> &LiveNurse::wishes() const {
  return pPreferences_->nurseWishes(nurseNum_);
}

// returns true if the nurse wishes the day-shift off
double LiveNurse::wishCostOfTheShift(int day, const PShift &pS) const {
  return pPreferences_->wishCostOfTheShift(nurseNum_, day, pS);
}

// returns true if the nurses reached the maximum number of consecutive worked
// days or is resting and did not reach the minimum number of resting days yet
// if consecutive number of shifts will only be reached by violating maximum
// number of worked days, go to rest only if consecutive working days penalty
// is the larger
bool LiveNurse::needRest(int day) {
  State state = states_[day];

  if (state.pShift_->isWork() && state.consDaysWorked_ >= maxConsDaysWork()) {
    return (state.consShifts_ >=
        pScenario_->minConsShiftsOfType(state.pShift_->type))
        || (pScenario_->weights().consShifts
            <= pScenario_->weights().consDaysWork);
  }
  return state.pShift_->isRest() && state.consDaysOff_ <= minConsDaysOff();
}

// returns true if the nurse needs to work one more day to reach the minimum
// number of consecutive working days or consecutive shifts
// if consecutive number of shifts will only be reached by violating maximum
// number of worked days, go to work only if consecutive shift penalty is
// the larger
bool LiveNurse::needWork(int day) {
  State state = states_[day];

  if (state.pShift_->isWork()) {
    if (state.consDaysWorked_ < minConsDaysWork()) {
      return true;
    } else if (state.consShifts_
        < pScenario_->minConsShiftsOfType(state.pShift_->type)) {
      return state.consDaysWorked_ < maxConsDaysWork() ||
          (pScenario_->weights().consShifts
              > pScenario_->weights().consDaysWork);
    }
  }

  return false;
}

// return true if the nurse is free to go to rest or work more without penalty
bool LiveNurse::isFreeToChoose(int day) {
  return !needWork(day) && !needRest(day);
}

// check the satisfaction of the hard constraints
// check the soft constraints and record the costs of the violations and the
// remaining margin for the satisfied ones.
void LiveNurse::checkConstraints(const Roster &roster,
                                 const vector<State> &states,
                                 int nbDays,
                                 bool payExcessImmediately,
                                 StatCtNurse *stat) {
  // check the satisfaction of the hard constraints and record the violations
  for (int day = 0; day < nbDays; day++) {
    // Check that the nurse has the assigned skill
    if (roster.pShift(day)->isWork()) {
      stat->violSkill_[day] = !hasSkill(roster.skill(day));
    } else {
      stat->violSkill_[day] = false;
    }

    // Check the forbidden successor constraint
    //
    stat->violSuccShifts_[day] =
        !roster.pShift(day)->canSucceed(*states[day].pShift_);
  }

  // check the soft constraints and record the costs of the violations and the
  // remaining margin for the satisfied ones.
  //
  for (int day = 1; day <= nbDays; day++) {
    // shift assigned on the previous day
    const PShift &prevPS = states[day - 1].pShift_, &pS = states[day].pShift_;

    // first look at consecutive working days or days off
    //
    int missingDays = 0, extraDays = 0;
    stat->costConsDays_[day - 1] = 0;
    stat->costConsDaysOff_[day - 1] = 0;

    // compute the violations of consecutive working days an
    if (pS->isWork()) {
      if (prevPS->isRest())
        missingDays = minConsDaysOff() - states[day - 1].consDaysOff_;

      stat->costConsDaysOff_[day - 1] +=
          (missingDays > 0) ? pScenario_->weights().consDaysOff
              * missingDays : 0;
      stat->costConsDays_[day - 1] +=
          (states[day].consDaysWorked_ > maxConsDaysWork())
          ? pScenario_->weights().consDaysWork : 0;
    } else {
      if (prevPS->isWork())
        missingDays = minConsDaysWork() - states[day - 1].consDaysWorked_;
      extraDays = states[day].consDaysOff_ - maxConsDaysOff();

      stat->costConsDays_[day - 1] +=
          (missingDays > 0) ? pScenario_->weights().consDaysWork
              * missingDays : 0;
      stat->costConsDaysOff_[day - 1] +=
          (extraDays > 0) ? pScenario_->weights().consDaysOff : 0;
    }

    // check the consecutive same shifts
    //
    stat->costConsShifts_[day - 1] = 0;
    int missingShifts = 0;

    // count the penalty for minimum consecutive shifts only for the previous
    // day when the new shift is different
    if (prevPS->isWork() && pS->type != prevPS->type) {
      missingShifts = pScenario_->minConsShiftsOfType(prevPS->type)
          - states[day - 1].consShifts_;
      stat->costConsShifts_[day - 1] +=
          (missingShifts > 0) ? pScenario_->weights().consShifts
              * missingShifts : 0;
    }

    // count the penalty for maximum consecutive shifts when the shift is worked
    // the last day will then be counted
    if (pS->isWork()) {
      stat->costConsShifts_[day - 1] +=
          (states[day].consShifts_ > pScenario_->maxConsShiftsOfType(pS->type))
          ? pScenario_->weights().consShifts : 0;
    }

    // check the preferences
    stat->costPref_[day - 1] = wishCostOfTheShift(day - 1, pS);

    // check the complete weekend (only if the nurse requires them)
    // this cost is only assigned to the sundays
    stat->costWeekEnd_[day - 1] = 0;
    if (pContract_->isLastWeekendDay(day - 1) && needCompleteWeekends()) {
      if (prevPS->isWork() ^ pS->isWork())
        stat->costWeekEnd_[day - 1] = pScenario_->weights().completeWeekend;
    }
  }  // end for day

  // get the costs due to total number of working days and weekends
  // pay for LB
  int n = std::max(0, (minTotalShifts() - states[nbDays].totalTimeWorked_));
  stat->costMissingDays_ = pScenario_->weights().totalShifts * n;

  // pay for UB
  n = std::max(0, states[nbDays].totalTimeWorked_ - maxTotalShifts());
  stat->costExceedingDays_ = pScenario_->weights().totalShifts * n;
  stat->costTotalDays_ = stat->costMissingDays_ + stat->costExceedingDays_;

  n = std::max(0, states[nbDays].totalWeekendsWorked_ - maxTotalWeekends());
  stat->costTotalWeekEnds_ = pScenario_->weights().totalWeekends * n;
}

// Build States from the roster
void LiveNurse::buildStates() {
  for (unsigned int k = 1; k < states_.size(); ++k)
    states_[k].addDayToState(states_[k - 1], roster_.pShift(k - 1));
}

//-----------------------------------------------------------------------------
// Compute the maximum and minimum number of working days from the input
// current state and in the next nbDays without getting any penalty for
// consecutive working days/days-off
//-----------------------------------------------------------------------------
std::pair<int, int> LiveNurse::computeMinMaxDaysNoPenaltyConsDay(
    State *pCurrentState, int nbDays) {
  int minWorkDaysNoPenaltyConsDays = 0, maxWorkDaysNoPenaltyConsDays = 0;

  // first treat the history
  //
  // remaining number of days when building a maximum and minimum working days
  // periods
  int nbDaysMin = nbDays, nbDaysMax = nbDays;
  if (pCurrentState->pShift_ && pCurrentState->pShift_->isWork()) {
    // if the last shift was not a rest
    // keep working until maximum consecutive number of shifts for maximum
    // working days or until minimum consecutive number of shifts for minimum
    // working days
    maxWorkDaysNoPenaltyConsDays =
        std::min(nbDays,
                 std::max(0,
                          maxConsDaysWork() - pCurrentState->consDaysWorked_));
    minWorkDaysNoPenaltyConsDays =
        std::min(nbDays,
                 std::max(0,
                          minConsDaysWork() - pCurrentState->consDaysWorked_));

    // perform minimum rest for maximum working days and maximum rest for
    // minimum working days
    nbDaysMax = nbDays
        - std::min(nbDays, maxWorkDaysNoPenaltyConsDays + minConsDaysOff());
    nbDaysMin = nbDays
        - std::min(nbDays, minWorkDaysNoPenaltyConsDays + maxConsDaysOff());
  } else if (pCurrentState->pShift_ && pCurrentState->pShift_->isRest()) {
    // if the last shift was a rest
    // keep resting until minimum consecutive number of rests for maximum
    // working days or until maximum consecutive number of rests for minimum
    // working days
    nbDaysMax = nbDays -
        std::min(nbDays,
                 std::max(0, minConsDaysOff() - pCurrentState->consDaysOff_));
    nbDaysMin = nbDays -
        std::min(nbDays,
                 std::max(0, maxConsDaysOff() - pCurrentState->consDaysOff_));
  }

  // lengths of the stints maximizing and minimizing the percentage of working
  // days respectively
  //
  int lengthStintMax = maxConsDaysWork() + minConsDaysOff();
  int lengthStintMin = minConsDaysWork() + maxConsDaysOff();

  // perform as many of these stints as possible
  //
  maxWorkDaysNoPenaltyConsDays +=
      (nbDaysMax / lengthStintMax) * maxConsDaysWork()
          + std::min(nbDaysMax % lengthStintMax, maxConsDaysWork());
  minWorkDaysNoPenaltyConsDays +=
      (nbDaysMin / lengthStintMin) * minConsDaysWork()
          + std::min(nbDaysMin % lengthStintMax, minConsDaysWork());

  return {minWorkDaysNoPenaltyConsDays, maxWorkDaysNoPenaltyConsDays};
}

// Print the contract type + preferences
void LiveNurse::printContractAndPreferences(const PScenario &pScenario) const {
  std::cout << "# Preferences:" << std::endl;
  Preferences::toString(wishes());
  std::cout << "# Contract :   ";
  std::cout << "Days [" << pContract_->minConsDaysOff_ << "<"
            << pContract_->maxConsDaysOff_ << "]   ";
  for (int s = 1; s < pScenario->nShiftTypes(); s++) {
    std::cout << pScenario->shiftType(s) << " ["
              << pScenario->minConsShiftsOfType(s) << "<"
              << pScenario->maxConsShiftsOfType(s) << "]   ";
  }
  std::cout << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# " << std::endl;
}

//-----------------------------------------------------------------------------
//  C l a s s   S o l v e r
//  Solves the offline problem
//  From a given problem (number of weeks, nurses, etc.),
//  can compute a solution.
//-----------------------------------------------------------------------------
vector<PLiveNurse> createLiveNurse(const PScenario &pScenario) {
  vector<PLiveNurse> liveNurses;

  for (int i = 0; i < pScenario->nNurses(); i++) {
    liveNurses.push_back(std::make_shared<LiveNurse>(
        *(pScenario->pNurse(i)),
        pScenario,
        pScenario->pDemand()->nDays_,
        pScenario->pDemand()->firstDayId_,
        &(*pScenario->pInitialState())[i],
        pScenario->pWeekPreferences(),
        i));
  }
  return liveNurses;
}
// Specific constructor
Solver::Solver(const PScenario &pScenario) :
    pScenario_(pScenario),
    pDemand_(pScenario->pDemand()),
    pPreferences_(pScenario->pWeekPreferences()),
    pInitState_(pScenario->pInitialState()),
    theLiveNurses_(createLiveNurse(pScenario)),
    dynamicWeights_(pScenario, createLiveNurse(pScenario)),
    // create the timer that records the lifetime of the solver and start it
    timerTotal_("Solver main", true),
    job_(),
    status_(UNSOLVED),
    totalCostUnderStaffing_(-1),
    maxTotalStaffNoPenalty_(-1),
    isPreprocessedSkills_(false),
    isPreprocessedNurses_(false) {
  if (pDemand_ == nullptr)
    Tools::throwError("Scenario cannot be initialized "
                      "with a nullptr for the demand.");
  if (pPreferences_ == nullptr)
    Tools::throwError("Scenario cannot be initialized "
                      "with a nullptr for the preferences.");
  if (pInitState_ == nullptr)
    Tools::throwError("Scenario cannot be initialized "
                      "with a nullptr for the initial states.");
  // initialize the preprocessed data of the skills
  for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
    maxStaffPerSkillNoPenalty_.push_back(-1.0);
    maxStaffPerSkillAvgWork_.push_back(-1.0);
    skillRarity_.push_back(1.0);
  }

  // modify some parameters depending on the benchmark
  if (pScenario_->isINRC_) {
    param_.optimalAbsoluteGap_ = 1.0;
    param_.absoluteGap_ = 1.0;
  }
}

Solver::~Solver() {}

Solver *Solver::newSolver(
    PScenario pScenario, Algorithm algo, SPType spType, SolverType sType) {
  switch (algo) {
    case GENCOL:
      switch (spType) {
        case LONG_ROTATION:
        case ALL_ROTATION:return new RotationMP(pScenario, sType);
          break;
        case ROSTER:return new RosterMP(pScenario, sType);
        default:Tools::throwError("The subproblem type is not handled yet");
      }
      break;
    default:Tools::throwError("The algorithm is not handled yet");
      break;
  }
  return nullptr;
}

// Load a solution in the solver
void Solver::loadSolution(const vector<Roster> &solution) {
  if (solution.size() != pScenario_->nNurses())
    Tools::throwError("Solver::loadSolution: there is not one roster "
                      "per nurse in the input solution!");

  solution_ = solution;
  for (int n = 0; n < pScenario_->nNurses(); n++) {
    PLiveNurse pNurse = theLiveNurses_[n];
    pNurse->roster(solution[n]);
    pNurse->buildStates();
  }
}

// copy the solver solution to this solver
void Solver::copySolution(Solver *pSolver) {
  loadSolution(pSolver->solution());
}

//------------------------------------------------------------------------
// Preprocess the nurses
// go through the nurses to collect data regarding the potential shift and
// skill coverage of the nurses
//-------------------------------------------------------------------------
void Solver::preprocessTheNurses() {
  // local variables for conciseness of the code
  //
  int nbDays = pDemand_->nDays_;

  this->specifyNursePositions();
  if (!isPreprocessedSkills_) this->preprocessTheSkills();
  this->computeMinMaxDaysNoPenaltyTotalDays();
  this->computeMinMaxDaysNoPenaltyConsDays();

  maxTotalStaffNoPenalty_ = 0;
  maxTotalStaffAvgWork_ = 0;
  for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
    maxStaffPerSkillNoPenalty_[sk] = 0;
    maxStaffPerSkillAvgWork_[sk] = 0;
  }
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    // add the working maximum number of working days to the maximum staffing
    //
    maxTotalStaffNoPenalty_ +=
        std::min(pNurse->maxWorkDaysNoPenaltyConsDays_,
                 pNurse->maxWorkDaysNoPenaltyTotalDays_);
    maxTotalStaffAvgWork_ +=
        static_cast<int>(pNurse->maxAvgWorkDaysNoPenaltyTotalDays_);

    // the staffing per skill is very rough here since Nurses can have
    // multiple skills
    // we thus weight the covering power of each nurse for each skill according
    // to the rarities of the skills
    double totalRarity = 0;
    for (int sk : pNurse->skills_)
      totalRarity += skillRarity_[sk];
    double maxWorkNoPenalty = std::min(pNurse->maxWorkDaysNoPenaltyConsDays_,
                                       pNurse->maxWorkDaysNoPenaltyTotalDays_);
    for (int sk : pNurse->skills_) {
      maxStaffPerSkillNoPenalty_[sk] +=
          (skillRarity_[sk] / totalRarity) * maxWorkNoPenalty;
      maxStaffPerSkillAvgWork_[sk] += (skillRarity_[sk] / totalRarity)
          * pNurse->maxAvgWorkDaysNoPenaltyTotalDays_;
    }
  }


  // initialize to zero the satisfied demand
  //
  Tools::initVector3D(&satisfiedDemand_,
                      nbDays,
                      pScenario_->nShifts(),
                      pScenario_->nSkills(),
                      0);

  isPreprocessedNurses_ = true;
}

// Find the position of each nurse
void Solver::specifyNursePositions() {
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    // the skills of the nurse need to be compared to the skills of each
    // existing position to determine the position of the nurse
    bool isPosition = true;
    for (const auto &pPosition : pScenario_->pPositions()) {
      isPosition = true;
      if (pPosition->nSkills() == pNurse->nSkills()) {
        for (int j = 0; j < pNurse->nSkills(); j++) {
          if (pNurse->skills_[j] != pPosition->skills_[j]) {
            isPosition = false;
            break;
          }
        }
      } else {
        isPosition = false;
      }

      if (isPosition) {
        pNurse->pPosition_ = pPosition;
        break;
      }
    }
    if (!isPosition) {
      Tools::throwError("The nurse has no position!");
    }
  }
}

// Compute the maximum and minimum number of working days in the period of
// the demand without getting any penalty for the total number of working days
void Solver::computeMinMaxDaysNoPenaltyTotalDays() {
  // number of days that will be covered after the current demand
  int nbDaysFuture =
      7 * (pScenario_->nWeeks() - pScenario_->thisWeek()) - pDemand_->nDays_;

  // For each contract, compute the maximum and minimum number of working days
  // that can be done after the current demand without ensuing penalties due to
  // the min/max numbers of working days and days off
  //
  vector<int> minWorkDaysFutureNoPenaltyConsDays(pScenario_->nContracts(), 0);
  vector<int> maxWorkDaysFutureNoPenaltyConsDays(pScenario_->nContracts(), 0);
  for (const auto &pContract : pScenario_->pContracts()) {
    // initialization
    minWorkDaysFutureNoPenaltyConsDays[pContract->id_] = 0;
    maxWorkDaysFutureNoPenaltyConsDays[pContract->id_] = 0;

    // we compute the minimum number of working days by assuming that the
    // first next week can start with a rest and we compute the maximum number
    // of working days by assuming that the first next week can start with work
    //

    // Compute the length of the stint that minimimizes/maximizes the number of
    // working days and  that respects the min/max number of consecutive working
    // days and days off
    // (a stint is a work rotation followed by a break)
    int lengthStintMin =
        pContract->minConsDaysWork_ + pContract->maxConsDaysOff_;
    int lengthStintMax =
        pContract->maxConsDaysWork_ + pContract->minConsDaysOff_;

    // In each case, the nurse would perform the min/max stint until there is
    // not enough remaining days to complete a stint. It will then start with
    // rest for the min stints or with work for the max stint
    minWorkDaysFutureNoPenaltyConsDays[pContract->id_] +=
        (nbDaysFuture / lengthStintMin) * pContract->minConsDaysWork_
            + std::max(0,
                       (nbDaysFuture % lengthStintMin)
                           - pContract->maxConsDaysOff_);
    maxWorkDaysFutureNoPenaltyConsDays[pContract->id_] +=
        (nbDaysFuture / lengthStintMax) * pContract->maxConsDaysWork_
            + std::min((nbDaysFuture % lengthStintMax),
                       pContract->maxConsDaysWork_);
  }

  // For each nurse, deduce the interval [dayMin,dayMax] defined as:
  // if the nurse works d days, with d in [dayMin,dayMax], then there exists a
  // sequence of future stints  such that the nurse respects the bounds on
  // the total number of working days without penalty
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    pNurse->minWorkDaysNoPenaltyTotalDays_ =
        std::max(0, pNurse->minTotalShifts()
            - pNurse->totalTimeWorked()
            - maxWorkDaysFutureNoPenaltyConsDays[pNurse->contractId()]);
    pNurse->maxWorkDaysNoPenaltyTotalDays_ =
        std::max(0, pNurse->maxTotalShifts()
            - pNurse->totalTimeWorked()
            - minWorkDaysFutureNoPenaltyConsDays[pNurse->contractId()]);
  }

  // Get the interval [avgMin,avgMax] defined as:
  // if the nurse works d days, with d in [avgMin,avgMax], then if it works
  // the same number of days every week, the nurse will respect the bounds on
  // the total number of working days without penalty
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    double demandMinAvgPerWeek =
        pNurse->minTotalShifts() * 1.0 / pScenario_->nWeeks();
    double demandMaxAvgPerWeek =
        pNurse->maxTotalShifts() * 1.0 / pScenario_->nWeeks();
    double demandMinAvgUntilThisDemand = demandMinAvgPerWeek
        * (pScenario_->thisWeek() + pDemand_->nDays_ / 7.0);
    double demandMaxAvgUntilThisDemand = demandMaxAvgPerWeek
        * (pScenario_->thisWeek() + pDemand_->nDays_ / 7.0);

    pNurse->minAvgWorkDaysNoPenaltyTotalDays_ =
        demandMinAvgUntilThisDemand - pNurse->totalTimeWorked();
    pNurse->maxAvgWorkDaysNoPenaltyTotalDays_ =
        std::max(0.0, demandMaxAvgUntilThisDemand - pNurse->totalTimeWorked());
  }
}

// compute the maximum and minimum number of working days in the period of
// the demand without getting any penalty for the number of consecutive
// shifts
// RqJO: this neglects the constraint of complete weekends and the
// preferences ; they should be added later
void Solver::computeMinMaxDaysNoPenaltyConsDays() {
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    auto p = pNurse->computeMinMaxDaysNoPenaltyConsDay(pNurse->pStateIni_,
                                                       pDemand_->nDays_);
    pNurse->minWorkDaysNoPenaltyConsDays_ = p.first;
    pNurse->maxWorkDaysNoPenaltyConsDays_ = p.second;
  }
}

//------------------------------------------------------------------------
// Preprocess the skills to get their rarity
// the value depends on the minimum demand for this skill, on the number
// of nurses that have the skill and on the number of skills per nurse
// that have the skill
//------------------------------------------------------------------------
void Solver::preprocessTheSkills() {
  // this vector will contain for each skill: the number of weighted nurses
  // that have the skill ; the weight of each nurse is its number of skills
  vector<double> nbNursesWeighted;

  for (int sk = 0; sk < pScenario_->nSkills(); sk++) {
    nbNursesWeighted.push_back(.0);
    for (const PLiveNurse &pNurse : theLiveNurses_)
      if (pNurse->hasSkill(sk))
        nbNursesWeighted[sk] +=
            pNurse->maxTotalShifts() * 1.0 / pow(pNurse->nSkills(), 2);

    // the skill rarity is the ratio of the the demand for the skill to the
    // weighted number of nurses that have the skill
    skillRarity_[sk] = pDemand_->minPerSkill_[sk] / nbNursesWeighted[sk];
  }

  // update the rarities of the skills in the scenario
  for (int p = 0; p < pScenario_->nPositions(); p++)
    pScenario_->pPosition(p)->updateRarities(skillRarity_);

  isPreprocessedSkills_ = true;
}

//------------------------------------------------------------------------------
// Create the vector of sorted nurses
// The nurses are ordered according to their position and the nurses that have
// the same position are shuffled
//------------------------------------------------------------------------------
void Solver::sortShuffleTheNurses() {
  vector<PPosition> positionsSorted;
  vector<vector<PLiveNurse>> nursePerPosition;

  // first, sort the position in the order in which the nurses should be treated
  for (int p = 0; p < pScenario_->nPositions(); p++) {
    positionsSorted.push_back(pScenario_->pPosition(p));
  }
  std::sort(positionsSorted.begin(), positionsSorted.end(), comparePositions);

  // then organize the nurses depending on their position
  for (int p = 0; p < pScenario_->nPositions(); p++) {
    vector<PLiveNurse> emptyVector;
    nursePerPosition.push_back(emptyVector);
  }
  for (int n = 0; n < pScenario_->nNurses(); n++) {
    PLiveNurse pNurse = theLiveNurses_[n];
    nursePerPosition[pNurse->pPosition()->id()].push_back(pNurse);
  }

  // shuffle the nurses that have the same position
  for (int p = 0; p < pScenario_->nPositions(); p++) {
    std::shuffle(nursePerPosition[p].begin(), nursePerPosition[p].end(),
                 Tools::getANewRandomGenerator());
  }

  // fill the sorted vector of live nurses
  theNursesSorted_.clear();
  for (int p = 0; p < pScenario_->nPositions(); p++) {
    int id = positionsSorted[p]->id();
    for (const auto &pN : nursePerPosition[id])
      theNursesSorted_.push_back(pN);
  }

  // todo: think about sorting the nurses according to their contract, we might
  // want to use the full-time nurses first
}

//------------------------------------------------------------------------------
// Initialize the greedy by preprocessing all the input attributes and sorting
// the shifts, skills, nurses
//------------------------------------------------------------------------------
void Solver::preprocessData() {
  // Preprocess the attributes of the greedy solver
  // the result of the preprocessing will be very useful to sort the attributes
  // before greedily covering the demand
  //
  if (!pDemand_->isPreprocessed_) pDemand_->preprocessMinDemand();
  if (!isPreprocessedSkills_) this->preprocessTheSkills();
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  // sort the nurses
  sortShuffleTheNurses();
}

//-----------------------------------------------------------------------------
// Compute the weights o the violation of the min/max number of working days
// For now, the update depends only on the initial states and on the contract
// of the nurses, on the number of days on the demand, on the number of weeks
// already treated and on the number of weeks left
// The required data on the nurses is mostly computed in preprocessTheNurses
//-----------------------------------------------------------------------------
void Solver::boundsAndWeights(WeightStrategy strategy) {
  // The nurses must be preprocessed to retrieve the information relative to the
  // past activity of the nurses and to their capacity to work more
  // in the future
  if (!isPreprocessedNurses_) this->preprocessTheNurses();

  dynamicWeights_.boundsAndWeights(strategy, pDemand_->nDays_);
}

// Compute min/max bounds as the ratio :
// number of days in demand / total number of remaining days
void Solver::computeBoundsAccordingToInitialState() {
  computeBoundsAccordingToRatio(1);
}

// Compute min/max bounds as the ratio :
// number of days in demand / total number of remaining days
void Solver::computeBoundsAccordingToDemandSize() {
  double ratio =
      nDays() / (7.0 * (pScenario_->nWeeks() - pScenario_->thisWeek()));
  computeBoundsAccordingToRatio(ratio);
}

void Solver::computeBoundsAccordingToRatio(double ratio) {
  // initialize the minimum and maximum number of total working days
  for (const auto &pN : theLiveNurses_) {
    const auto &pTotalDuration = pN->totalShiftDurationResource_;
    const auto &pTotalWeekend = pN->totalWeekendResource_;
    // default min and max
    int b = std::floor((pN->minTotalShifts() -
        pN->pStateIni_->totalTimeWorked_) * ratio + 1e-3);
    pTotalDuration->setLb(std::max(0, b));
    b = std::ceil((pN->maxTotalShifts() -
        pN->pStateIni_->totalTimeWorked_) * ratio - 1e-3);
    pTotalDuration->setUb(std::max(0, b));
    b = std::ceil((pN->maxTotalWeekends()
        - pN->pStateIni_->totalWeekendsWorked_) * ratio - 1e-3);
    pTotalWeekend->setUb(std::max(0, b));
  }
}

//------------------------------------------------------------------------
// Compare two nurses based on their position
// the function is used to sort the nurses in ascending rank of their
// position
// if their positions have the same rank, then the smaller nurse is found
// by a lexicographic comparison of the rarity of the skills of the nurses
//------------------------------------------------------------------------
bool compareNurses(const PLiveNurse &n1, const PLiveNurse &n2) {
  return comparePositions(n1->pPosition(), n2->pPosition());
}

//------------------------------------------------------------------------
// Compare two positions to sort them
// Three possible cases can happen
// 1) same positions
// 2) same rank: the first position to be treated is that with the rarest skill
// or the largest number of skills
// 3) the first position to be treated is that with the smaller rank
//------------------------------------------------------------------------
bool comparePositions(const PPosition &p1, const PPosition &p2) {
  if (p1->id() == p2->id()) {
    return false;
  } else if (p1->rank() == p2->rank()) {
    // the skillRarity vector is ALWAYS sorted in descending order, because the
    // updateRarities is the only setter for skillRarity and it sorts the vector
    for (int sk = 0; sk < std::min(p1->nSkills(), p2->nSkills()); sk++) {
      if (p1->skillRarity(sk) != p2->skillRarity(sk)) {
        return p1->skillRarity(sk) > p2->skillRarity(sk);
      }
      return p1->nSkills() > p2->nSkills();
    }
  } else {
    return p1->rank() < p2->rank();
  }
  return true;
}

//------------------------------------------------------------------------------
// Check the feasibility of the demand with these nurses
//-----------------------------------------------------------------------------

bool Solver::checkFeasibility() {
  Tools::throwError("Solver::checkFeasibility() is not implemented.");
  return true;
}

//------------------------------------------------------------------------------
// Count the fraction of current solution that is integer
//------------------------------------------------------------------------------
double Solver::computeFractionOfIntegerInCurrentSolution() const {
  vector3D<double> fracRoster = fractionalRoster();
  int nbFractional = 0;
  int nbDayNurse = 0;
  for (int nurse = 0; nurse < nNurses(); nurse++) {
    for (int day = 0; day < nDays(); day++) {
      nbDayNurse++;
      for (int shift = 1; shift < nShifts(); shift++) {
        double activity = fracRoster[nurse][day][shift];
        if ((activity < 1 - epsilon()) && (activity > epsilon())) {
          nbFractional++;
          break;
        }
      }
    }
  }

  return 1.0 - nbFractional * 1.0 / nbDayNurse;
}

////------------------------------------------------------------------------
//// Compute the fractional penalty due to weekend that should be paid if in a
//// model with complete plannings
////------------------------------------------------------------------------
// double Solver::computeFractionalWeekendPenalty() {
//  vector3D<double> fracRoster = fractionalRoster();
//  int nbWeeks = pScenario_->nWeeks() - pScenario_->thisWeek();
//  double fractionalWeekendPenalty = 0.0;
//  for (int nurse = 0; nurse < nNurses(); nurse++) {
//    int nbActiveWeekends = 0;
//    double penalizedFraction = 1.0;
//    for (int w = 0; w < nbWeeks; w++) {
//      double weekendActivity;
//      double activitySaturday = 0.0;
//      double activitySunday = 0.0;
//      for (int shift = 1; shift < nShifts(); shift++) {
//        activitySaturday += fracRoster[nurse][7 * w + 5][shift];
//        activitySunday += fracRoster[nurse][7 * w + 6][shift];
//      }
//      weekendActivity = std::max(activitySaturday, activitySunday);
//      if (weekendActivity > epsilon()) {
//        nbActiveWeekends++;
//        penalizedFraction = std::min(penalizedFraction, weekendActivity);
//      }
//    }
//    fractionalWeekendPenalty += penalizedFraction
//        * std::max(0.0, nbActiveWeekends - maxTotalWeekends_[nurse]);
//  }
//  return fractionalWeekendPenalty * pScenario_->weights().totalWeekends;
// }

// get the total cost of the current solution
// the solution is simply given by the roster of each nurse
double Solver::computeSolutionCost(int nbDays, bool payExcessImmediately) {
  double totalCost = 0;
  int nbNurses = pScenario_->nNurses();
  int nbShifts = pScenario_->nShifts(), nbSkills = pScenario_->nSkills();

  // reset the satisfied demand to compute it from scratch
  Tools::initVector3D(&satisfiedDemand_, nbDays, nbShifts, nbSkills, 0);

  // first add the individual cost of each nurse
  bool reachedHorizon =
      (pScenario_->thisWeek() + nbDays / 7 >= pScenario_->nWeeks());
  for (int n = 0; n < nbNurses; n++) {
    PLiveNurse pNurse = theLiveNurses_[n];
    StatCtNurse &stat = pNurse->statCt_;
    pNurse->checkConstraints(
        pNurse->roster_, pNurse->states_, nbDays, payExcessImmediately, &stat);

    for (int day = 0; day < nbDays; day++) {
      if (pNurse->roster_.pShift(day)->isWork()) {
        int s = pNurse->roster_.pShift(day)->id;
        satisfiedDemand_[day][s][pNurse->roster_.skill(day)]++;
      }
    }

    // add total costs -> for INRC2 use the original cost of the problem
    if (pScenario_->isINRC2_) {
      for (int day = 0; day < nbDays; day++) {
        totalCost += stat.costConsDays_[day] + stat.costConsDaysOff_[day] +
            stat.costConsShifts_[day] + stat.costPref_[day]
            + stat.costWeekEnd_[day];
      }
      if (reachedHorizon) {
        totalCost += stat.costTotalDays_ + stat.costTotalWeekEnds_;
      } else if (payExcessImmediately) {
        totalCost += stat.costExceedingDays_ + stat.costTotalWeekEnds_;
      }
    } else {
      totalCost += pNurse->roster_.cost();
    }
  }

  // add the cost of non-optimal demand
  for (int day = 0; day < nbDays; day++) {
    for (int sh = 1; sh < nbShifts; sh++) {
      for (int sk = 0; sk < nbSkills; sk++) {
        // check min demand
        int diff = satisfiedDemand_[day][sh][sk] -
            pDemand_->minDemand_[day][sh][sk];
        // hard if over-coverage forbidden
        if (pScenario_->weights().overCoverage < 0) {
          if (diff < 0) return XLARGE_SCORE;
        } else {
          if (diff < 0) {
            totalCost -= diff * pScenario_->weights().underCoverage;
          } else if (diff > 0) {
            totalCost += diff * pScenario_->weights().overCoverage;
          }
        }
        // check opt demand
        if (pDemand_->isOptDemand_) {
          int missingStaff;
          missingStaff = std::max(0, pDemand_->optDemand_[day][sh][sk]
              - satisfiedDemand_[day][sh][sk]);
          totalCost += pScenario_->weights().underCoverage * missingStaff;
        }
      }
    }
  }

  return totalCost;
}

// get aggregate information on the solution and write them in a string
string Solver::solutionStatisticsToString(int nbDays) {
  // start by computing every individual cost of the nurses
  computeSolutionCost(nbDays);

  // compute for each contract
  // 1) the number of penalized working days
  // 2) the number of extra days that could have been worked without penalty
  // 3) the total preferences cost
  // compute for each week and each contract
  // 1) the number of assignments
  // 2) the preferences cost
  // compute the same statistics for each skill and each week
  // compute the total number of assignments and the total minimum and
  // optimum demand
  vector<int> unusedDays(pScenario_->nContracts(), 0);
  vector<int> extraDays(pScenario_->nContracts(), 0);
  vector<double> costPrefPerWeek(pScenario_->nWeeks(), .0);
  vector<double> costPrefPerContract(pScenario_->nContracts(), .0);
  vector2D<int> assignmentsPerContractAndWeek;
  Tools::initVector2D(&assignmentsPerContractAndWeek,
                      pScenario_->nContracts(),
                      pScenario_->nWeeks(),
                      0);
  int totalAssignments = 0;
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    int c = pNurse->pContract_->id_;

    unusedDays[c] = std::max(0, pNurse->maxTotalShifts()
        - pNurse->state(nbDays).totalTimeWorked_);
    extraDays[c] = std::max(0, pNurse->state(nbDays).totalTimeWorked_
        - pNurse->maxTotalShifts());
    totalAssignments += pNurse->state(nbDays).totalTimeWorked_;

    for (int w = 0; w < pScenario_->nWeeks(); w++) {
      assignmentsPerContractAndWeek[c][w] +=
          pNurse->state(7 * (w + 1)).totalTimeWorked_
              - pNurse->state(7 * w).totalTimeWorked_;
    }

    for (int day = 0; day < nbDays; day++) {
      int week = (day / 7);
      costPrefPerWeek[week] += pNurse->statCt_.costPref_[day];
      costPrefPerContract[c] += pNurse->statCt_.costPref_[day];
    }
  }

  string strStats;

  return strStats;
}

//------------------------------------------------
// Display functions
//------------------------------------------------

// Return the solution at day k
vector<Roster> Solver::solutionAtDay(int k) {
  vector<Roster> ans;
  ans.reserve(solution_.size());
  int nbDays = k + 1;
  for (const Roster &r : solution_) {
    // if we do not need to cut the current solution, just add it directly
    if (r.nDays() == nbDays) {
      ans.push_back(r);
      continue;
    }
    // otherwise, build a new roster from the current solution
    int firstDay = r.firstDayId();
    vector<PShift> pShifts = r.pShifts();
    pShifts.resize(nbDays);
    vector<int> skills = r.skills();
    skills.resize(nbDays);
    Roster s(firstDay, pShifts, skills);
    ans.push_back(s);
  }
  return ans;
}

// Extend the rosters in the solution with the days covered by the input
// solution
void Solver::extendSolution(vector<Roster> solutionExtension) {
  // detect an error if the sizes of the extension is not the same as that
  // of the solution
  if (solutionExtension.size() != solution_.size()) {
    Tools::throwError("Solver::extendSolution: "
                      "the extension has not the same size as the solution");
  }

  // Extend each roster in the solution
  for (unsigned int i = 0; i < solution_.size(); i++) {
    solution_[i].pushBack(solutionExtension[i]);
  }
}

// return the final states of the nurses
vector<State> Solver::finalStates() {
  vector<State> pFinalStates;
  int nbDays = pDemand_->nDays_;
  for (const PLiveNurse &pNurse : theLiveNurses_) {
    pFinalStates.push_back(pNurse->state(nbDays));
  }

  return pFinalStates;
}

// return the states of the nurses at day k
vector<State> Solver::statesOfDay(int k) {
  vector<State> pStatesOfDayK;
  for (const PLiveNurse &pNurse : theLiveNurses_)
    pStatesOfDayK.push_back(pNurse->state(k));
  return pStatesOfDayK;
}

// display the whole solution
string Solver::solutionToString() {
  return solutionToString(pDemand_->firstDayId_,
                          pDemand_->nDays_,
                          pScenario_->thisWeek());
}

// display the whole solution week by week for nbWeeks weeks.
vector<string> Solver::solutionToString(int nbWeeks) {
  vector<string> solutions;

  // build the solution for each week
  int firstDay = pDemand_->firstDayId_;
  for (int w = 0; w < nbWeeks; ++w) {
    solutions.push_back(solutionToString(
        firstDay, 7, pScenario_->thisWeek() + w));
    firstDay += 7;
  }

  return solutions;
}

// display the solution between firstDay and firstDay+nbDays in the required
// format
string Solver::solutionToString(int firstDay, int nbDays, int firstWeek) {
  std::stringstream rep;
  int nbNurses = pScenario_->nNurses();

  // write to stringstream that can then be printed in any output file
  // follow the template described by the competition
  rep << "SOLUTION" << std::endl;
  rep << firstWeek << " " << pScenario_->name() << std::endl;
  rep << std::endl;

  // compute the total number of assignments that are not rests
  // if no shift is assigned to a nurse on given, it still counts
  int nbAssignments = 0;
  std::vector<int> endDays;
  for (const auto &pN : theLiveNurses_) {
    int endDay = firstDay + nbDays;
    if (pN->roster_.nDays() < endDay) {
      std::cerr << "Roster for nurse " << pN->name_ << " has only "
                << pN->roster_.nDays() << " days (less than the requested "
                << endDay << ")" << std::endl;
      endDay = pN->roster_.nDays();
    }
    endDays.push_back(endDay);
    for (int day = firstDay; day < endDay; day++) {
      if (pN->roster_.pShift(day)->isWork()) nbAssignments++;
    }
  }

  rep << "ASSIGNMENTS = " << nbAssignments << std::endl;
  for (int n = 0; n < nbNurses; n++) {
    for (int day = firstDay; day < endDays[n]; day++) {
      const PShift &pS = theLiveNurses_[n]->roster_.pShift(day);
      if (pS->isWork()) {
        rep << theLiveNurses_[n]->name_ << " "
            << Day::toDayOfWeekShortName(day) << " ";
        int skill = theLiveNurses_[n]->roster_.skill(day);
        rep << pS->name << " " << pScenario_->skillName(skill) << std::endl;
      }
    }
  }
  return rep.str();
}

void Solver::solutionToTxt(string outdir) {
  Tools::LogOutput solS(outdir, false);
  solS << solutionToString();
}

void Solver::solutionToUI(string outdir) {
  string fileName = outdir + pScenario_->name() + ".txt";
  Tools::LogOutput uiStream(fileName, false);
  uiStream << pScenario_->getHeader();
  uiStream << "ASSIGNMENTS" << std::endl;
  for (const auto &pN : theLiveNurses_) {
    for (int day = 0; day < pScenario_->nDays(); day++) {
      const PShift &pS = pN->roster_.pShift(day);
      if (pS->isWork()) {
        uiStream << pN->nurseNum_ << ",";
        tm *date = Tools::dateForDay(&pScenario_->startDate(), day);
        uiStream << date->tm_year + 1900 <<
                 "-" << date->tm_mon + 1 << "-" << date->tm_mday << ",";
        uiStream << pS->name << ",";
        int skillId = pN->roster_.skill(day);
        uiStream << pN->pScenario_->skillName(skillId) << std::endl;
      }
    }
  }

  // Add stats to Json
  uiStream << "REPORT" << std::endl << "{";
  map<string, double> costs = costsConstraintsByName();
  for (auto const &p : costs)
    uiStream << "\"" << p.first << "\":" << p.second << "," << std::endl;
  uiStream << R"("Total cost":)" << objValue_ << "," << std::endl;
  uiStream << R"("Optimization status":")" << statusToString.at(status_)
           << "\"," << std::endl;
  uiStream << "}" << std::endl;
}

void Solver::solutionToXmlINRC(string outdir) {
  if (outdir.empty()) outdir = param_.outdir_;
  string outXml = outdir + pScenario_->name() + ".xml";
  Tools::LogOutput xmlStream(outXml, false);

  xmlStream << "<Solution>\n";
  xmlStream << "    <SchedulingPeriodID>" << pScenario_->name()
            << "</SchedulingPeriodID>\n";
  xmlStream << "    <Competitor>Legrain, Omer</Competitor>\n";
  for (const auto &pN : theLiveNurses_) {
    for (int day = 0; day < pScenario_->nDays(); day++) {
      const PShift &pS = pN->roster_.pShift(day);
      if (pS->isWork()) {
        xmlStream << "    <Assignment>\n";
        tm *date = Tools::dateForDay(&pScenario_->startDate(), day);
        xmlStream << "        <Date>" << date->tm_year + 1900 <<
                  "-" << date->tm_mon + 1 << "-" << date->tm_mday
                  << "</Date>\n";
        xmlStream << "        <Employee>" << pN->num_ << "</Employee>\n";
        xmlStream << "        <ShiftType>" << pS->name << "</ShiftType>\n";
        xmlStream << "    </Assignment>\n";
      }
    }
  }
  xmlStream << "</Solution>\n";
}

std::string Solver::solutionToSolINRC() {
  std::stringstream rep;

  rep << "Obj. = 100";

  for (const auto &pN : theLiveNurses_) {
    rep << "\n";
    for (int day = 0; day < pScenario_->nDays(); day++) {
      const PShift &pS = pN->roster_.pShift(day);
      rep << pS->id - 1 << " ";
    }
  }
  rep << "\nTime = 100.0s";
  return rep.str();
}

// display the whole solution in a more readable format and append advanced
// information on the solution quality
string Solver::writeResourceCosts() {
  if (pScenario_->isINRC2_)
    return writeResourceCostsINRC2();
  return "";
}

string Solver::writeResourceCostsPerNurse(bool resetInitialState) {
  std::stringstream rep;
  int name_w = 10, shift_w = 2;

  double totalCost = 0, totalS = 0, totalWE = 0;
  for (const auto &pN : theLiveNurses_) {
    totalCost += pN->roster_.cost();
    rep << "================================" << std::endl;
    rep << "Employee: " << pN->name_
        << ", Contract: " << pN->pContract_->name_ << std::endl;
    rep << "================================" << std::endl;
    if (pN->totalShiftDurationResource_) {
      int nS = pN->states_[nDays()].totalTimeWorked_;
      if (resetInitialState) nS -= pN->states_[0].totalTimeWorked_;
      rep << "\tTotal Shifts: " << nS << " for ["
          << pN->totalShiftDurationResource_->getLb() << "("
          << pN->totalShiftDurationResource_->getLbCost() << ") , "
          << pN->totalShiftDurationResource_->getUb() << "("
          << pN->totalShiftDurationResource_->getUbCost() << ")]";
      double c = pN->totalShiftDurationResource_->getCost(nS);
      totalS += c;
      if (abs(c) > 1e-3) rep << " -> " << c;
      rep << std::endl;
    }
    if (pN->totalWeekendResource_) {
      int nWE = pN->states_[nDays()].totalWeekendsWorked_;
      if (resetInitialState) nWE -= pN->states_[0].totalWeekendsWorked_;
      rep << "\tTotal Weekends: " << nWE << " for ["
          << pN->totalWeekendResource_->getLb() << "("
          << pN->totalWeekendResource_->getLbCost() << ") , "
          << pN->totalWeekendResource_->getUb() << "("
          << pN->totalWeekendResource_->getUbCost() << ")]";
      double c = pN->totalWeekendResource_->getCost(nWE);
      totalWE += c;
      if (abs(c) > 1e-3) rep << " -> " << c;
      rep << std::endl;
    }
    for (int t = CONS_SHIFTS_COST; t <= TOTAL_WEEKEND_COST; t++) {
      auto ct = (CostType) t;
      if (pN->roster_.costs().at(ct) > 0) {
        rep << "\t" << namesByCostType.at(ct)
            << ": " << pN->roster_.costs().at(ct) << std::endl;
      }
    }
  }
  rep << std::endl;
  rep << "================================" << std::endl;
  rep << "Total shifts costs = " << totalS << std::endl;
  rep << "Total weekend costs = " << totalWE << std::endl;
  rep << "Total individual costs = " << totalCost << std::endl;
  rep << "================================" << std::endl;

  return rep.str();
}

string Solver::writeResourceCostsINRC2() {
  std::stringstream rep;
  int nbNurses = pScenario_->nNurses(), nbShifts = pScenario_->nShifts();
  int nbSkills = pScenario_->nSkills();
  int firstDay = pDemand_->firstDayId_, nbDays = pDemand_->nDays_;
  int name_w = 10, shift_w = 2;

  // if not using, default resources, the following costs doesn't mean anything
  if (!useDefaultResources())
    return rep.str();

  // compute the total cost and in the mean time update the structures of each
  // live nurse that contains all the required information on soft and hard
  // constraints satisfaction
  //
  double totalCost = computeSolutionCost(nbDays);

  // store temporarily the data that is about to be written
  //
  int violMinCover = 0, violReqSkill = 0, violForbiddenSucc = 0;
  double costOptCover = 0, costTotalDays = 0, costTotalWeekEnds = 0;
  double costConsDays = 0, costConsDaysOff = 0, costConsShifts = 0,
      costPref = 0, costWeekEnds = 0;

  // constraints related to the demand
  for (int day = firstDay; day < firstDay + nbDays; day++)
    for (int sh = 1; sh < nbShifts; sh++)
      for (int sk = 0; sk < nbSkills; sk++) {
        violMinCover += std::max(0, pDemand_->minDemand_[day][sh][sk]
            - satisfiedDemand_[day][sh][sk]);
        costOptCover += pScenario_->weights().underCoverage
            * std::max(0, pDemand_->optDemand_[day][sh][sk]
                - satisfiedDemand_[day][sh][sk]);
      }

  bool reachedHorizon =
      (pScenario_->thisWeek() + nbDays / 7 >= pScenario_->nWeeks());
  for (int n = 0; n < nbNurses; n++) {
    PLiveNurse pNurse = theLiveNurses_[n];

    costTotalDays += reachedHorizon ?
                     pNurse->statCt_.costTotalDays_
                                    : pNurse->statCt_.costExceedingDays_;
    costTotalWeekEnds += pNurse->statCt_.costTotalWeekEnds_;

    for (int day = firstDay; day < firstDay + nbDays; day++) {
      // record the violations
      int skill = pNurse->roster_.skill(day);
      const PShift &pS = pNurse->roster_.pShift(day);
      if (pS->isWork()) violReqSkill += pNurse->hasSkill(skill) ? 0 : 1;
      if (!pS->canSucceed(*pNurse->states_[day].pShift_))
        violForbiddenSucc += 1;

      // the other costs per soft constraint can be read from the stat structure
      costConsDays += pNurse->statCt_.costConsDays_[day];
      costConsDaysOff += pNurse->statCt_.costConsDaysOff_[day];
      costConsShifts += pNurse->statCt_.costConsShifts_[day];
      costPref += pNurse->statCt_.costPref_[day];
      costWeekEnds += pNurse->statCt_.costWeekEnd_[day];
    }
  }

  // write the status of hard and soft constraints
  //
  rep << "Hard constraints violations\n";
  rep << "---------------------------\n";
  rep << "Minimal coverage constraints: " << violMinCover << std::endl;
  rep << "Required skill constraints: " << violReqSkill << std::endl;
  rep << "Illegal shift type succession constraints: " << violForbiddenSucc
      << std::endl;
  rep << "Single assignment per day: 0" << std::endl;

  rep << "\nCost per constraint type\n";
  rep << "------------------------\n";
  rep << "Total assignment constraints: " << costTotalDays << std::endl;
  rep << "Consecutive working days constraints: " << costConsDays << std::endl;
  rep << "Consecutive days off constraints: " << costConsDaysOff << std::endl;
  rep << "Consecutive shifts constraints: " << costConsShifts << std::endl;
  rep << "Preferences: " << costPref << std::endl;
  rep << "Max working weekend: " << costTotalWeekEnds << std::endl;
  rep << "Complete weekends: " << costWeekEnds << std::endl;
  rep << "Optimal coverage constraints: " << costOptCover << std::endl;

  rep << "\n---------------------------\n";
  rep << "\nTotal cost: " << totalCost << std::endl;

  return rep.str();
}

string Solver::solutionToLogString() {
  std::stringstream rep;
  int nbNurses = pScenario_->nNurses(), nbShifts = pScenario_->nShifts();
  int nbSkills = pScenario_->nSkills();
  int firstDay = pDemand_->firstDayId_, nbDays = pDemand_->nDays_;
  int name_w = 10, shift_w = 2;
  rep << "Complete shift schedule" << std::endl << std::endl;
  rep << std::setw(name_w) << "";
  for (int day = firstDay; day < firstDay + nbDays; day++) {
    if (day != firstDay && Day::isFirstDayOfWeek(day))
      rep << " |";
    rep << " | " << std::setw(shift_w)
        << Day::toDayOfWeekShortName(day).at(0);
  }
  rep << " |" << std::endl;

  int sep_w = name_w + 1 + (shift_w + 3) * nbDays + 2 * (nbDays / 7) + 1;
  string sep(sep_w, '=');
  rep << sep << std::endl;

  for (int n = 0; n < nbNurses; n++) {
    PLiveNurse pNurse = theLiveNurses_[n];
    rep << std::setw(name_w) << pNurse->name_;
    for (int day = firstDay; day < firstDay + nbDays; day++) {
      if (day != firstDay && Day::isFirstDayOfWeek(day))
        rep << " |";
      const PShift &pS = pNurse->roster_.pShift(day);
      if (pS->isWork())
        rep << " | " << std::setw(shift_w) << pS->name.substr(0, shift_w);
      else
        rep << " | " << string(shift_w, '-');
    }
    rep << " |" << std::endl;
  }
  rep << std::endl;

  return rep.str();
}

// When a solution of multiple consecutive weeks is available,
// display the complete solution in the log and write the solution of
// the weeks separately
bool Solver::displaySolutionMultipleWeeks(bool addNumNursesInFileName,
                                          std::string outDir) {
  if (outDir.empty()) outDir = param_.outdir_;
  // treat the case where the solver was unable to find a feasible solution
  if (status_ == INFEASIBLE) {
    std::cerr
        << "Solver::displaySolutionMultipleWeeks: "
           "The solver was not able to find a solution"
        << std::endl;
    return false;
  }

  // write separately the solutions of each week in the required output format
  int nWeeks = pScenario_->nDays() / 7, firstWeek = pScenario_->thisWeek(),
      lastWeek = firstWeek + nWeeks;
  // just write the week for which the indices were provided.
  // if empty, write all of them
  vector<int> weekIndices = param_.weekIndices_;
  if (weekIndices.empty()) {
    for (int w = pScenario_->thisWeek(); w < lastWeek; ++w)
      weekIndices.push_back(w);
  } else if (weekIndices.size() <= nWeeks) {
    nWeeks = weekIndices.size();
  } else {
    std::cerr << "Too many weeks (" << weekIndices.size()
              << ") have been requested through weekIndices_." << std::endl;
  }
  // write the weeks within weekIndices
  vector<string> sols = solutionToString(nWeeks);
  for (int w = 0; w < nWeeks; ++w) {
    string solutionFile = outDir + "sol";
    if (addNumNursesInFileName)
      solutionFile += "-n" + std::to_string(nNurses());
    solutionFile += "-week" + std::to_string(weekIndices[w]) + ".txt";
    Tools::LogOutput solutionStream(solutionFile, false);
    solutionStream << sols[w];
  }

  return true;
}
