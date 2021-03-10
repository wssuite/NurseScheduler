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

#include <exception>

#include "ParseArguments.h"
#include "tools/ReadWrite.h"
#include "tools/DemandGenerator.h"
#include "solvers/StochasticSolver.h"
// #include "CbcModeler.h"
#include "tools/Tools.h"
#include "solvers/InitializeSolver.h"

using std::string;
using std::vector;
using std::pair;

/******************************************************************************
* Solve one week inside the stochastic process
******************************************************************************/
void solveOneWeek(InputPaths *pInputPaths) {
  const std::string &solPath = pInputPaths->solutionPath();
  string::size_type found = solPath.find_last_of(".");
  string logPathIni = solPath.substr(0, found),
      logPath = logPathIni + "Log.txt";
  Tools::LogOutput logStream(logPath);


  // set the scenario
  logStream << "# Initialize the scenario" << std::endl;
  PScenario pScen = initializeScenario(*pInputPaths);

  // set the options of the stochastic solver
  // (the corresponding method needs to be change manually for tests)
  logStream << "# Set the options" << std::endl;
  StochasticSolverOptions options;
  setStochasticSolverOptions(
      &options, pScen, solPath, logPathIni, pInputPaths->timeOut());

  // check if a parameters file in the the in the log path ini
  // This files are hard coded as there are no options to give them
  // to the simulator of the INRCII
  string pathIni = "";
  found = solPath.find_last_of("/");
  if (found != string::npos)
    pathIni = solPath.substr(0, found + 1);
  string stochasticOptions = pathIni + "solverOptions.txt",
      generationOptions = pathIni + "generationOptions.txt",
      evaluationOptions = pathIni + "evaluationOptions.txt";
  try {
    logStream << "Stochastic options:" << std::endl <<
              ReadWrite::readStochasticSolverOptions(stochasticOptions,
                                                     &options)
              << std::endl;
  } catch (const char *ex) {}
  try {
    logStream << "Generation options:" << std::endl <<
              ReadWrite::readSolverOptions(generationOptions,
                                           &options.generationParameters_)
              << std::endl;
  } catch (const char *ex) {}
  try {
    logStream << "Evaluation options:" << std::endl <<
              ReadWrite::readSolverOptions(evaluationOptions,
                                           &options.evaluationParameters_)
              << std::endl;
  } catch (const char *ex) {}

  // get history demands by reading the custom file
  //
  vector<PDemand> demandHistory;
  demandHistory.push_back(pScen->pWeekDemand());
  if (!pInputPaths->customInputFile().empty())
    ReadWrite::readCustom(
        pInputPaths->customInputFile(), pScen, &demandHistory);

  Solver *pSolver = new StochasticSolver(pScen, options, demandHistory);

  logStream << "# Solve the week" << std::endl;
  pSolver->solve();
  int solutionStatus = pSolver->status();
  logStream << "# Solution status = " << solutionStatus << std::endl;

  Tools::LogOutput solStream(solPath);
  solStream << pSolver->solutionToString() << std::endl;

  //  release memory
  delete pSolver;
  logStream.close();

  // Write the solution in the required output format
  //
  if (!pInputPaths->customOutputFile().empty()) {
    ReadWrite::writeCustom(pInputPaths->customOutputFile(),
                           pInputPaths->week(0),
                           pInputPaths->customInputFile());
  }
  std::cout << "Custom output file : "
            << pInputPaths->customOutputFile() << std::endl;
  // Todo: the method that writes the history file corresponding to the
  // solution
  // string outputHistoryFile("history-week");
  // outputHistoryFile += std::to_string(pScen->thisWeek()) + ".txt";
  // std::cout << "Output history file: " << outputHistoryFile << std::endl;
}

/******************************************************************************
* Test a solution on multiple weeks
* In this method, the weeks are solved sequentially without knowledge of future
* demand
******************************************************************************/

pair<double, int> testMultipleWeeksStochastic(string dataDir,
                                              string instanceName,
                                              int historyIndex,
                                              vector<int> weekIndices,
                                              StochasticSolverOptions
                                              stochasticSolverOptions,
                                              string outdir,
                                              std::vector<int> seeds) {
  // build the paths of the input files
  InputPaths inputPaths(dataDir, instanceName, historyIndex, weekIndices);

  // initialize the scenario object of the first week
  PScenario pScen = initializeScenario(inputPaths);

  // solve the problem for each week and store the solution in the vector below
  vector<Roster> solution;
  int nbWeeks = weekIndices.size();
  Status solutionStatus;

  vector<PDemand> demandHistory;
  double currentCost = 0;
  int nbSched = 0;

  for (int week = 0; week < nbWeeks; week++) {
    if (week >= seeds.size())
      seeds.push_back(std::rand());
    std::srand(seeds[week]);

    demandHistory.push_back(pScen->pWeekDemand());

    StochasticSolver *pSolver = new StochasticSolver(pScen,
                                                     stochasticSolverOptions,
                                                     demandHistory,
                                                     currentCost);

    currentCost += pSolver->solve();
    nbSched += pSolver->getNbSchedules();
    printf("Current cost = %.2f \n", currentCost);
    solutionStatus = pSolver->status();
    if (solutionStatus == INFEASIBLE) {
      delete pSolver;
      break;
    }

    // update the overall solution with the solution of the week that was just
    // treated
    // warning: we must make sure that the main solution in pSolver applies only
    // to the demand of one week
    vector<Roster> weekSolution = pSolver->solutionAtDay(6);
    if (solution.empty()) {
      solution = weekSolution;
    } else {
      for (int n = 0; n < pScen->nNurses(); n++)
        solution[n].push_back(weekSolution[n]);
    }

    std::cout << pSolver->solutionToLogString() << std::endl;

    // prepare the scenario for next week if we did not reach the last week yet
    if (week < nbWeeks - 1) {
      // Initialize demand and preferences
      PDemand pDemand(nullptr);
      PPreferences pPref(nullptr);

      // Read the demand and preferences and link them with the scenario
      ReadWrite::readWeek(inputPaths.week(week + 1), pScen, &pDemand, &pPref);

      // read the initial state of the new week from the last state of the
      // last week
      // modify the dayId_ to show that this is the first day of the new week
      vector<State> initialStates = pSolver->statesOfDay(6);
      for (int i = 0; i < pScen->nNurses(); i++) {
        initialStates[i].dayId_ = 0;
      }

      // update the scenario to treat next week
      pScen->updateNewWeek(pDemand, pPref, initialStates);
    }

    delete pSolver;
  }

  printf("Total cost = %.2f \n", currentCost);

  string seedOutfile = outdir + "seeds.txt";
  Tools::LogOutput seedStream(seedOutfile, true);
  char str[50];
  snprintf(str, sizeof(str),
           "Cost %.2f; NbGene %d; Seeds", currentCost, nbSched);
  seedStream << str;
  for (int s : seeds)
    seedStream << " " << s;
  seedStream << std::endl;

  // Display the solution
  // displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex,
  // weekIndices, solution, solutionStatus, outDir);

  return {currentCost, nbSched};
}

/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char **argv) {
  std::cout << "Number of arguments= " << argc << std::endl;

  // Detect errors in the number of arguments
  //
  if (argc % 2 != 1) {
    Tools::throwError("main: There should be an even number of arguments!");
  } else if (argc > 1 && (argc < 9 || argc > 17)) {
    Tools::throwError(
        "main: There is either too many or not enough arguments!");
  }

  // Simulate default behavior for a test instance
  //
  if (argc == 1) {
    std::cout << "Running the default method..." << std::endl;
    string dataDir = "datasets/";
    string instanceName = "n030w4";
    int historyIndex = 1;
    vector<int> weekIndices = {6, 7, 5, 3};
    StochasticSolverOptions stochasticSolverOptions;
    stochasticSolverOptions.totalTimeLimitSeconds_ = 40;
    string outdir = "outfiles/" + instanceName + "/";
    std::vector<int> seeds = {50, 35, 70, 80};
    testMultipleWeeksStochastic(dataDir,
                                instanceName,
                                historyIndex,
                                weekIndices,
                                stochasticSolverOptions,
                                outdir,
                                seeds);
  } else {
    // Nominal behavior of the executable, as required by INRCII
    // Retrieve the file names in arguments
    //
    InputPaths *pInputPaths = readArguments(argc, argv);

    // Initialize the random seed
    //
    srand(pInputPaths->randSeed());

    // Solve the week
    solveOneWeek(pInputPaths);

    delete pInputPaths;
  }
}
