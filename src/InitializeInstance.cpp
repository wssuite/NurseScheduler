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

#include "InitializeInstance.h"

#include <dirent.h>
#include <algorithm>
#include <memory>
#include <set>
#include <utility>

#include "ReadWrite.h"
#include "solvers/mp/RotationMP.h"
#include "tools/Tools.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

#include "ReadWrite.h"

PScenario buildInstance(const InputPaths &inputPaths) {
  std::cout << "Origin guessed: " << inputPaths.origin() << std::endl;

//  if (REST_SHIFT_ID != 0)
//    Tools::throwError("Scenario reader works only with REST_SHIFT_ID = 0.");
  // UI file
  if (inputPaths.ui()) {
    // Read the scenario
    return ReadWrite::readScenarioUI(inputPaths.scenario(),
                                     inputPaths.instance());
  } else if (inputPaths.inrc2()) {
    // Read the scenario
    PScenario pScenario = ReadWrite::readScenarioINRC2(inputPaths.scenario());
    // Read the demand and preferences and link them with the scenario
    ReadWrite::readINRC2Weeks(inputPaths.weeks(), pScenario);
    // Read the history
    ReadWrite::readHistoryINRC2(inputPaths.history(), pScenario);

    // Initialize the resources
    initializeResourcesINRC2(pScenario);
    return pScenario;
  } else if (inputPaths.inrc()) {
    return ReadWrite::readINRCInstance(inputPaths.scenario());
  } else if (inputPaths.nrp()) {
    return ReadWrite::readNRPInstance(inputPaths.scenario());
  } else {
    Tools::throwError("The instance file is not recognized.");
    return nullptr;
  }
}

void initializeResourcesINRC2(const PScenario &pScenario) {
// initialize all the resources of the nurses
  int nbDays = pScenario->nDays();
  vector<State> *pInitialState = pScenario->pInitialState();
  const Weights &weights = pScenario->weights();

  for (const auto &pN : pScenario->pNurses()) {
    const auto &stateIni = pInitialState->at(pN->num_);
    // a. total assignments
    pN->addBaseResource(
        std::make_shared<SoftTotalShiftDurationResource>(
            std::max(0, pN->minTotalShifts() - stateIni.totalTimeWorked_),
            std::max(0, pN->maxTotalShifts() - stateIni.totalTimeWorked_),
            weights.totalShifts,
            weights.totalShifts,
            pScenario->shiftsFactory().pAnyWorkShift(),
            nbDays,
            pScenario->maxDuration()));

    // b. consecutive assignments
    pN->addBaseResource(std::make_shared<SoftConsShiftResource>(
        pN->minConsDaysWork(),
        pN->maxConsDaysWork(),
        weights.consDaysWork,
        weights.consDaysWork,
        pScenario->shiftsFactory().pAnyWorkShift(),
        CONS_WORK_COST,
        nbDays,
        stateIni.consDaysWorked_,
        false));

    // c. total weekend resource; in INRC2 weekends are always on saturday and
    // sunday
    pN->addBaseResource(std::make_shared<SoftTotalWeekendsResource>(
        std::max(0, pN->maxTotalWeekends() - stateIni.totalWeekendsWorked_),
        weights.totalWeekends,
        pScenario->shiftsFactory().pAnyWorkShift(),
        nbDays));

    // d. complete weekend resource
    if (pN->needCompleteWeekends())
      pN->addBaseResource(std::make_shared<SoftIdentWeekendResource>(
          std::make_shared<ShiftWorkComparator>(),
          weights.completeWeekend));

    // e. consecutive shift resources (rest included)
    // initialize resources on the number of consecutive shifts of each type
    for (int st = 0; st < pScenario->nShiftTypes(); st++) {
      shared_ptr<AbstractShift> pAShift =
          pScenario->shiftsFactory().pAnyTypeShift(st);
      if (pAShift->isWork()) {
        int consShiftsInitial =
            pAShift->includes(*stateIni.pShift_) ?
            stateIni.consShifts_ : 0;
        pN->addBaseResource(std::make_shared<SoftConsShiftResource>(
            pScenario->minConsShiftsOfType(st),
            pScenario->maxConsShiftsOfType(st),
            weights.consShifts,
            weights.consShifts,
            pAShift,
            CONS_SHIFTS_COST,
            nbDays,
            consShiftsInitial));
      } else if (pAShift->isRest()) {
        // needed to evaluate the historical cost
        pN->addBaseResource(std::make_shared<SoftConsShiftResource>(
            pN->minConsDaysOff(),
            pN->maxConsDaysOff(),
            weights.consDaysOff,
            weights.consDaysOff,
            pAShift,
            CONS_REST_COST,
            nbDays,
            stateIni.consDaysOff_));
      }
    }
  }
}

/*****************************************************************************
* Separate the scenario into multiple scenarios that only focus on the nurses
* whose positions are in the same connected component of positions
******************************************************************************/

vector<PScenario> divideScenarioIntoConnectedPositions(
    const PScenario &pScenario) {
  vector<PScenario> scenariosPerComponent;

  // First, identify the connected components of the rcspp of positions
  pScenario->computeConnectedPositions();

  for (int c = 0; c < pScenario->nbOfConnectedComponentsOfPositions(); c++) {
    const vector<PPosition>
        &positionsInTheComponent = pScenario->componentOfConnectedPositions(c);
    const vector<PNurse> &nursesInTheComponent =
        pScenario->nursesInConnectedComponentOfPositions(c);

    // retrieve a vector containing the indices of the skills contained in the
    // component (decreasing order).
    // use a set first, because it manages duplicate skills automatically.
    std::set<int> skillsInTheComponent;
    for (const PPosition &pPosition : positionsInTheComponent)
      for (int skill : pPosition->skills())
        skillsInTheComponent.insert(skill);

    vector<int>
        skillsVector(skillsInTheComponent.begin(), skillsInTheComponent.end());
    std::stable_sort(skillsVector.begin(),
                     skillsVector.end(),
                     Tools::compareDecreasing);

    // build a vector containing the skills that need to be removed from the
    // scenario.
    vector<int> skillsToRemove;
    skillsToRemove.reserve(pScenario->nSkills());
    for (int skill = 0; skill < pScenario->nSkills(); skill++)
      skillsToRemove.push_back(skill);
    for (int skill : skillsVector)
      skillsToRemove.erase(skillsToRemove.begin() + skill);
    std::stable_sort(skillsToRemove.begin(),
                     skillsToRemove.end(),
                     Tools::compareDecreasing);

    // shorten the vectors intToSkill and skillToInt to match the list of skills
    vector<string> intToSkill(pScenario->nSkills());
    for (int sk : skillsToRemove) {
      intToSkill.erase(intToSkill.begin() + sk);
    }

    // create the demand that relates only to input skills
    PDemand pDemand = pScenario->pDemand();

    // erase the skills to remove from the minimum and optimal demands
    vector3D<int> minDemand = pDemand->minDemand_;
    for (int day = 0; day < pDemand->nDays_; day++)
      for (int shift = 0; shift < pDemand->nShifts_; shift++)
        for (int skill : skillsToRemove) {
          minDemand[day][shift][skill] = 0;
        }
    PDemand pDemandInTheComponent;
    if (pDemand->isOptDemand_) {
      vector3D<int> optDemand = pDemand->optDemand_;
      for (int day = 0; day < pDemand->nDays_; day++)
        for (int shift = 0; shift < pDemand->nShifts_; shift++)
          for (int skill : skillsToRemove) {
            optDemand[day][shift][skill] = 0;
          }
      pDemandInTheComponent =
          std::make_shared<Demand>(pDemand->nDays_,
                                   pDemand->firstDayId_,
                                   pDemand->nShifts_,
                                   pDemand->nSkills_,
                                   pDemand->name_,
                                   minDemand,
                                   optDemand);
    } else {
      pDemandInTheComponent =
          std::make_shared<Demand>(pDemand->nDays_,
                                   pDemand->firstDayId_,
                                   pDemand->nShifts_,
                                   pDemand->nSkills_,
                                   pDemand->name_,
                                   minDemand);
    }

    // create the preferences that relate only to the nurses of the connected
    // component.
    PPreferences pPreferencesInTheComponent =
        pScenario->pWeekPreferences()->keep(nursesInTheComponent);

    // Create the new scenario
    PScenario pScenarioInTheConnectedComponent =
        std::make_shared<Scenario>(pScenario,
                                   nursesInTheComponent,
                                   pDemandInTheComponent,
                                   pPreferencesInTheComponent);

    // create the initial states that relate only to the nurses of the connected
    // component.
    vector<State> intialStatesInTheComponent;
    vector<State> *pInitialState = pScenario->pInitialState();
    intialStatesInTheComponent.reserve(nursesInTheComponent.size());
    for (const PNurse &pNurse : nursesInTheComponent) {
      intialStatesInTheComponent.push_back(pInitialState->at(pNurse->num_));
    }
    pScenarioInTheConnectedComponent->setInitialState(
        intialStatesInTheComponent);

    // Push back in the vector of scenarios
    //
    scenariosPerComponent.push_back(pScenarioInTheConnectedComponent);
  }

  return scenariosPerComponent;
}

/******************************************************************************
* When a solution of multiple consecutive weeks is available, load it in a
* solver for all the weeks and  display the results
******************************************************************************/
void displaySolutionMultipleWeeks(const string &dataDir,
                                  const string &instanceName,
                                  int historyIndex,
                                  const vector<int> &weekIndices,
                                  const vector<Roster> &solution,
                                  Status status,
                                  const string &outDir) {
  if (outDir.empty()) return;

  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = static_cast<int>(weekIndices.size());
  string catWeeks;
  for (int w = 0; w < nbWeeks; w++) catWeeks += std::to_string(weekIndices[w]);
  string logPath = outDir + "Log-" + catWeeks + ".txt";
  Tools::LogOutput outStream(logPath, false);

  // treat the case where the solver was unable to find a feasible solution
  if (status == INFEASIBLE) {
    outStream << "The solver was not able to find a solution\n";
    return;
  }

  // load the solution in a new solver
  InputPaths inputPaths(dataDir, instanceName, historyIndex, weekIndices);
  PScenario pScen = buildInstance(inputPaths);
  auto pSolver = std::make_unique<Solver>(pScen);
  pSolver->loadSolution(solution);

  // write the log file for all the weeks
  outStream << pSolver->writeResourceCosts();

  // write separately the solutions of each week in the required output format
  vector<string> solutions = pSolver->solutionToString(nbWeeks);
  for (int w = 0; w < nbWeeks; ++w) {
    string solutionFile = outDir;
    solutionFile += "Sol-";
    solutionFile += instanceName;
    solutionFile += "-";
    solutionFile += catWeeks;
    solutionFile += "-";
    solutionFile += std::to_string(weekIndices[w]);
    solutionFile += "-";
    solutionFile += std::to_string(w);
    solutionFile += ".txt";
    Tools::LogOutput solutionStream(solutionFile, false);
    solutionStream << solutions[w];
  }
}

void displaySolutionMultipleWeeks(const InputPaths &inputPaths,
                                  const vector<Roster> &solution,
                                  Status status,
                                  const string &outDir) {
  if (outDir.empty()) return;

  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = inputPaths.nbWeeks();
  string catWeeks;
  for (int w = 0; w < nbWeeks; w++) catWeeks += inputPaths.week(w);
  string logPath = outDir + "Log-" + catWeeks + ".txt";
  Tools::LogOutput outStream(logPath, false);

  // treat the case where the solver was unable to find a feasible solution
  if (status == INFEASIBLE) {
    outStream << "The solver was not able to find a solution\n";
    return;
  }

  // load the solution in a new solver
  PScenario pScen = buildInstance(inputPaths);
  auto pSolver = std::make_unique<Solver>(pScen);
  pSolver->loadSolution(solution);

  // write the log file for all the weeks
  outStream << pSolver->writeResourceCosts();

  // write separately the solutions of each week in the required output format
  vector<string> solutions = pSolver->solutionToString(nbWeeks);
  for (int w = 0; w < nbWeeks; ++w) {
    string solutionFile =
        outDir + "Sol-" + inputPaths.instance() + "-" + catWeeks + "-"
            + inputPaths.week(w) + "-" + std::to_string(w) + ".txt";
    Tools::LogOutput solutionStream(solutionFile, false);
    solutionStream << solutions[w];
  }
}

/******************************************************************************
* Compute and record stats on all the demand files of all the instances in the
* input directory
******************************************************************************/

void computeStatsOnTheDemandsOfAllInstances(const string &inputDir) {
  struct dirent *dirp;

  // Open the input directory
  DIR *dp = opendir(inputDir.c_str());
  if (dp == nullptr) {
    Tools::throwError("Error while opening ");
  } else {
    std::cout << "Reading from directory " << inputDir << std::endl;
  }
  while ((dirp = readdir(dp))) {
    std::string filename(dirp->d_name);

    // The instance names start with "WD"
    if (filename[0] != 'n') continue;
    ReadWrite::compareDemands((string) (inputDir + filename),
                              (string) ("outfiles/statDemands/" + filename
                                  + ".txt"));
  }
}