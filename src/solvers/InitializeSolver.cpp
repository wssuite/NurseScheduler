/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 *  full license detail.
 */

#include "InitializeSolver.h"

#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <memory>
#include <set>

#include "tools/ReadWrite.h"
#include "solvers/mp/RotationMP.h"
#include "tools/MyTools.h"

// initialize the counter of objects
unsigned int MyObject::s_count = 0;
unsigned int Pattern::s_count = 0;

using std::string;
using std::vector;
using std::map;
using std::pair;

/******************************************************************************
* Read the arguments in non compact format
******************************************************************************/
InputPaths *readNonCompactArguments(int argc, char **argv) {
  InputPaths *pInputPaths = nullptr;

  // Read the arguments and store them in inputPaths
  //
  int narg = 1;
  while (narg < argc) {
    std::cout << "arg = " << argv[narg] << " " << argv[narg + 1] << std::endl;
    // Attention usine a gaz: problem with the java simulator that add quote
    // marks in the arguments, which messes with the open file methods
    // the code below is here to remove these quote marks
    //
    string str(argv[narg + 1]);
    std::size_t found = str.find("\"");
    while (found != std::string::npos) {
      str.erase(found, 1);
      found = str.find("\"");
    }

    if (!strcmp(argv[narg], "--sce")) {
      pInputPaths->scenario(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--his")) {
      pInputPaths->history(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--week")) {
      pInputPaths->addWeek(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--sol")) {
      pInputPaths->solutionPath(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--param")) {
      pInputPaths->paramFile(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--timeout")) {
      pInputPaths->timeOut(std::stod(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--rand")) {
      pInputPaths->randSeed(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-type")) {
      pInputPaths->SPType(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--n-threads")) {
      pInputPaths->nThreads(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-strategy")) {
      pInputPaths->SPStrategy(std::stoi(str));
      narg += 2;
    } else {
      Tools::throwError("main: the argument does not match the expected list!");
    }
  }
  // Throw an error if a necessary input file is missing
  if (pInputPaths->scenario().empty() || pInputPaths->history().empty()
      || pInputPaths->weeks().empty()) {
    throw Tools::myException(
        "readNonCompactArguments: A necessary file name is missing!",
        __LINE__);
  }

  // Non compact format is only for release versions, so no log file is required
  pInputPaths->logPath("");

  return pInputPaths;
}

/******************************************************************************
* Read the arguments in compact format
******************************************************************************/

InputPaths *readCompactArguments(int argc, char **argv) {
  // Default arguments are set to enable simple call to the function without
  // argument.
  std::string dataDir = "", instanceName = "", solutionPath = "", logPath = "",
      paramFile = "";
  std::string SPType = "LONG";
  int historyIndex = 0, randSeed = 0, nTreads = 1, SPStrategy = 0;
  std::vector<int> weekIndices;
  double timeOut = LARGE_TIME;

  // Read the arguments and store them in inputPaths
  //
  int narg = 1;
  while (narg < argc) {
    const char *arg = argv[narg];
    std::cout << "arg = " << arg << " " << argv[narg + 1] << std::endl;
    // WARNING: problem with the java simulator that add quote
    // marks in the arguments, which messes with the open file methods
    // the code below is here to remove these quote marks
    //
    string str(argv[narg + 1]);
    std::size_t found = str.find("\"");
    while (found != std::string::npos) {
      str.erase(found, 1);
      found = str.find("\"");
    }

    if (!strcmp(arg, "--dir")) {
      dataDir = str;
      narg += 2;
    } else if (!strcmp(arg, "--instance")) {
      instanceName = str;
      narg += 2;
    } else if (!strcmp(arg, "--his")) {
      historyIndex = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--weeks")) {
      weekIndices = Tools::parseList(str, '-');
      narg += 2;
    } else if (!strcmp(arg, "--sol")) {
      solutionPath = str;
      narg += 2;
    } else if (!strcmp(arg, "--log")) {
      logPath = str;
      narg += 2;
    } else if (!strcmp(arg, "--param")) {
      paramFile = str;
      narg += 2;
    } else if (!strcmp(arg, "--timeout")) {
      timeOut = std::stod(str);
      narg += 2;
    } else if (!strcmp(arg, "--rand")) {
      randSeed = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--sp-type")) {
      SPType = str;
      narg += 2;
    } else if (!strcmp(argv[narg], "--n-threads")) {
      nTreads = std::stoi(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-strategy")) {
      SPStrategy = std::stoi(str);
      narg += 2;
    } else {
      std::stringstream err_buff;
      err_buff << "main: the argument (" << arg
               << ") does not match the expected list!";
      Tools::throwError(err_buff.str());
    }
  }

  // Initialize the input paths
  //
  InputPaths *pInputPaths =
      new InputPaths(dataDir,
                     instanceName,
                     historyIndex,
                     weekIndices,
                     solutionPath,
                     logPath,
                     paramFile,
                     timeOut,
                     randSeed,
                     SPType,
                     SPStrategy,
                     nTreads);

  return pInputPaths;
}

/************************************************************************
* Initialize the week scenario by reading the input files
*************************************************************************/

PScenario initializeScenario(string scenPath,
                             string demandPath,
                             string historyPath,
                             string logPath) {
  // Initialize demand and preferences
  PDemand pDemand(nullptr);
  PPreferences pPref(nullptr);

  // Read the scenario
  PScenario pScen = ReadWrite::readScenario(scenPath);

  // Read the demand and preferences and link them with the scenario
  ReadWrite::readWeek(demandPath, pScen, &pDemand, &pPref);
  pScen->linkWithDemand(pDemand);
  pScen->linkWithPreferences(pPref);

  // Read the history
  ReadWrite::readHistory(historyPath, pScen);

  // Check that the scenario was read properly if logfile specified in input
  if (!logPath.empty()) {
    Tools::LogOutput logStream(logPath);
    logStream << pScen->toString() << std::endl;
    logStream << pScen->pWeekDemand()->toString(true) << std::endl;
  }

  return pScen;
}

PScenario initializeScenario(const InputPaths &inputPaths, string logPath) {
  // Initialize demand and preferences
  PDemand pDemand(nullptr);
  PPreferences pPref(nullptr);

  // Read the scenario
  PScenario pScen = ReadWrite::readScenario(inputPaths.scenario());

  // Read the demand and preferences and link them with the scenario
  ReadWrite::readWeek(inputPaths.week(0), pScen, &pDemand, &pPref);
  pScen->linkWithDemand(pDemand);
  pScen->linkWithPreferences(pPref);

  // Read the history
  ReadWrite::readHistory(inputPaths.history(), pScen);

  // Check that the scenario was read properly if logfile specified in input
  if (!logPath.empty()) {
    Tools::LogOutput logStream(logPath);
    logStream << pScen->toString() << std::endl;
    logStream << pScen->pWeekDemand()->toString(true) << std::endl;
  }

  return pScen;
}

/*****************************************************************************
* Initialize the scenario for multiple weeks
* When calling this function, the intent is to solve all the weeks at once
******************************************************************************/

PScenario initializeMultipleWeeks(string dataDir,
                                  string instanceName,
                                  int historyIndex,
                                  vector<int> weekIndices,
                                  string logPath) {
  // build the paths of the input files
  InputPaths inputPaths(dataDir, instanceName, historyIndex, weekIndices);

  // Read the scenario
  PScenario pScen = ReadWrite::readScenario(inputPaths.scenario());

  // Read the demand and preferences and link them with the scenario
  ReadWrite::readWeeks(inputPaths.weeks(), pScen);

  // Read the history
  ReadWrite::readHistory(inputPaths.history(), pScen);

  // Check that the scenario was read properly if logfile specified in input
  if (!logPath.empty()) {
    Tools::LogOutput logStream(logPath);
    logStream << pScen->toString() << std::endl;
    logStream << pScen->pWeekDemand()->toString(true) << std::endl;
  }

  return pScen;
}

PScenario initializeMultipleWeeks(const InputPaths &inputPaths,
                                  string logPath) {
  // Read the scenario
  PScenario pScen = ReadWrite::readScenario(inputPaths.scenario());

  // Read the demand and preferences and link them with the scenario
  ReadWrite::readWeeks(inputPaths.weeks(), pScen);

  // Read the history
  ReadWrite::readHistory(inputPaths.history(), pScen);

  // Check that the scenario was read properly if logfile specified in input
  if (!logPath.empty()) {
    Tools::LogOutput logStream(logPath);
    logStream << pScen->toString() << std::endl;
    logStream << pScen->pWeekDemand()->toString(true) << std::endl;
  }

  return pScen;
}

/*****************************************************************************
* Separate the scenario into multiple scenarios that only focus on the nurses
* whose positions are in the same connected component of positions
******************************************************************************/

vector<PScenario> divideScenarioIntoConnectedPositions(PScenario pScenario) {
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
    for (PPosition pPosition : positionsInTheComponent)
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
    for (int skill = 0; skill < pScenario->nbSkills_; skill++)
      skillsToRemove.push_back(skill);
    for (int skill : skillsVector)
      skillsToRemove.erase(skillsToRemove.begin() + skill);
    std::stable_sort(skillsToRemove.begin(),
                     skillsToRemove.end(),
                     Tools::compareDecreasing);

    // shorten the vectors intToSkill and skillToInt to match the list of skills
    vector<string> intToSkill(pScenario->intToSkill_);
    map<string, int> skillToInt(pScenario->skillToInt_);

    for (int skill : skillsToRemove) {
      skillToInt.erase(pScenario->intToSkill_[skill]);
      intToSkill.erase(intToSkill.begin() + skill);
    }

    // create the demand that relates only to input skills
    //
    PDemand pDemand = pScenario->pWeekDemand();

    // erase the skills to remove from the minimum and optimal demands
    vector3D<int> minDemand = pDemand->minDemand_;
    vector3D<int> optDemand = pDemand->optDemand_;
    for (int day = 0; day < pDemand->nbDays_; day++)
      for (int shift = 0; shift < pDemand->nbShifts_; shift++)
        for (int skill : skillsToRemove) {
          minDemand[day][shift][skill] = 0;
          optDemand[day][shift][skill] = 0;
        }
    PDemand pDemandInTheComponent =
        std::make_shared<Demand>(pDemand->nbDays_,
                                 pDemand->firstDay_,
                                 pDemand->nbShifts_,
                                 pDemand->nbSkills_,
                                 pDemand->name_,
                                 minDemand,
                                 optDemand);

    // create the preferences that relate only to the nurses of the connected
    // component.
    PPreferences pPreferencesInTheComponent =
        std::make_shared<Preferences>(nursesInTheComponent,
                                      pDemand->nbDays_,
                                      pScenario->nbShifts_);
    PPreferences pPreferences = pScenario->pWeekPreferences();

    // only keep the demand of the nurses in the component
    for (PNurse pNurse : nursesInTheComponent)
      for (const auto &itDay : pPreferences->nurseWishesOff(pNurse->id_))
        for (const auto &itShift : itDay.second)
          pPreferencesInTheComponent->addShiftOff(pNurse->id_,
                                                  itDay.first,
                                                  itShift.shift,
                                                  itShift.level);

    // only keep the demand of the nurses in the component
    for (PNurse pNurse : nursesInTheComponent)
      for (const auto &itDay : pPreferences->nurseWishesOn(pNurse->id_))
        for (const auto &itShift : itDay.second)
          pPreferencesInTheComponent->addShiftOn(pNurse->id_,
                                                 itDay.first,
                                                 itShift.shift,
                                                 itShift.level);

    // Create the new scenario
    //
    PScenario pScenarioInTheConnectedComponent =
        std::make_shared<Scenario>(pScenario,
                                   nursesInTheComponent,
                                   pDemandInTheComponent,
                                   pPreferencesInTheComponent);

    // create the initial states that relate only to the nurses of the connected
    // component.
    vector<State> intialStatesInTheComponent;
    vector<State> *pInitialState = pScenario->pInitialState();

    for (PNurse pNurse : nursesInTheComponent) {
      intialStatesInTheComponent.push_back(pInitialState->at(pNurse->id_));
    }
    pScenarioInTheConnectedComponent->setInitialState(
        intialStatesInTheComponent);

    // Push back in the vector of scenarios
    //
    scenariosPerComponent.push_back(pScenarioInTheConnectedComponent);
  }

  return scenariosPerComponent;
}

/*****************************************************************************
* Create a solver of the class specified by the input algorithm type
******************************************************************************/
Solver *setSolverWithInputAlgorithm(PScenario pScen, Algorithm algorithm) {
  Solver *pSolver = nullptr;
  switch (algorithm) {
    case GENCOL:
      // DBG: add solver type as option: S_CLP, S_GUROBI ...
      pSolver = new RotationMP(pScen,
                               pScen->pWeekDemand(),
                               pScen->pWeekPreferences(),
                               pScen->pInitialState(),
                               S_CLP);
      break;
    default: Tools::throwError("The algorithm is not handled yet");
      break;
  }
  return pSolver;
}

/******************************************************************************
* When a solution of multiple consecutive weeks is available, load it in a
* solver for all the weeks and  display the results
******************************************************************************/
void displaySolutionMultipleWeeks(string dataDir,
                                  string instanceName,
                                  int historyIndex,
                                  const vector<int> &weekIndices,
                                  const vector<Roster> &solution,
                                  Status status,
                                  string outDir) {
  if (outDir.empty()) return;

  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = weekIndices.size();
  string catWeeks;
  for (int w = 0; w < nbWeeks; w++) catWeeks += std::to_string(weekIndices[w]);
  string logPath = outDir + "Log-" + catWeeks + ".txt";
  Tools::LogOutput outStream(logPath);

  // treat the case where the solver was unable to find a feasible solution
  if (status == INFEASIBLE) {
    outStream << "The solver was not able to find a solution\n";
    return;
  }

  // load the solution in a new solver
  PScenario pScen =
      initializeMultipleWeeks(dataDir, instanceName, historyIndex, weekIndices);
  Solver *pSolver = new Solver(pScen,
                               pScen->pWeekDemand(),
                               pScen->pWeekPreferences(),
                               pScen->pInitialState());
  pSolver->loadSolution(solution);

  // write the log file for all the weeks
  outStream << pSolver->solutionToLogString();

  // write separately the solutions of each week in the required output format
  vector<string> solutions = pSolver->solutionToString(nbWeeks);
  for (int w = 0; w < nbWeeks; ++w) {
    string solutionFile = outDir + "Sol-" + instanceName + "-" + catWeeks + "-"
        + std::to_string(weekIndices[w]) + "-" + std::to_string(w) + ".txt";
    Tools::LogOutput solutionStream(solutionFile);
    solutionStream << solutions[w];
  }
  delete pSolver;
}

void displaySolutionMultipleWeeks(InputPaths inputPaths,
                                  const vector<Roster> &solution,
                                  Status status,
                                  string outDir) {
  if (outDir.empty()) return;

  // initialize the log stream
  // first, concatenate the week numbers
  int nbWeeks = inputPaths.nbWeeks();
  string catWeeks;
  for (int w = 0; w < nbWeeks; w++) catWeeks += inputPaths.week(w);
  string logPath = outDir + "Log-" + catWeeks + ".txt";
  Tools::LogOutput outStream(logPath);

  // treat the case where the solver was unable to find a feasible solution
  if (status == INFEASIBLE) {
    outStream << "The solver was not able to find a solution\n";
    return;
  }

  // load the solution in a new solver
  PScenario pScen = initializeMultipleWeeks(inputPaths);
  Solver *pSolver = new Solver(pScen,
                               pScen->pWeekDemand(),
                               pScen->pWeekPreferences(),
                               pScen->pInitialState());
  pSolver->loadSolution(solution);

  // write the log file for all the weeks
  outStream << pSolver->solutionToLogString();

  // write separately the solutions of each week in the required output format
  vector<string> solutions = pSolver->solutionToString(nbWeeks);
  for (int w = 0; w < nbWeeks; ++w) {
    string solutionFile =
        outDir + "Sol-" + inputPaths.instance() + "-" + catWeeks + "-"
            + inputPaths.week(w) + "-" + std::to_string(w) + ".txt";
    Tools::LogOutput solutionStream(solutionFile);
    solutionStream << solutions[w];
  }

  delete pSolver;
}

/******************************************************************************
* Compute and record stats on all the demand files of all the instances in the
* input directory
******************************************************************************/

void computeStatsOnTheDemandsOfAllInstances(string inputDir) {
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
