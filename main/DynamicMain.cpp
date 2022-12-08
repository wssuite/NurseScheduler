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

#include "ParseArguments.h"
#include "ReadWrite.h"
#include "tools/DemandGenerator.h"
#include "solvers/StochasticSolver.h"
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
  string::size_type found = solPath.find_last_of('.');
  string logPathIni = solPath.substr(0, found),
      logPath = logPathIni + "Log.txt";
  Tools::LogOutput logStream(logPath, false, true);


  // set the scenario
  logStream << "# Initialize the scenario" << std::endl;
  PScenario pScen = initializeScenarioINRC2(*pInputPaths);

  // set the options of the stochastic solver
  // (the corresponding method needs to be change manually for tests)
  logStream << "# Set the options" << std::endl;
  StochasticSolverOptions options;
  options.setStochasticSolverOptions(
      pScen, solPath, logPathIni, pInputPaths->timeOut());

  // check if a parameters file in the the in the log path ini
  // This files are hard coded as there are no options to give them
  // to the simulator of the INRCII
  string pathIni;
  found = solPath.find_last_of('/');
  if (found != string::npos)
    pathIni = solPath.substr(0, found + 1);
  string stochasticOptions = pathIni + "solverOptions.txt",
      generationOptions = pathIni + "generationOptions.txt",
      evaluationOptions = pathIni + "evaluationOptions.txt";
  try {
    logStream << "Stochastic options:" << std::endl;
    options.read(stochasticOptions);
  } catch (const char *ex) {}
  try {
    logStream << "Generation options:" << std::endl;
    options.generationParameters_.read(generationOptions);
  } catch (const char *ex) {}
  try {
    logStream << "Evaluation options:" << std::endl;
    options.evaluationParameters_.read(evaluationOptions);
  } catch (const char *ex) {}

  // get history demands by reading the custom file
  //
  vector<PDemand> demandHistory;
  demandHistory.push_back(pScen->pDemand());
  if (!pInputPaths->customInputFile().empty())
    ReadWrite::readCustom(
        pInputPaths->customInputFile(), pScen, &demandHistory);

  Solver *pSolver = new StochasticSolver(pScen, options, demandHistory);

  logStream << "# Solve the week" << std::endl;
  pSolver->solve();
  int solutionStatus = pSolver->status();
  logStream << "# Solution status = " << solutionStatus << std::endl;
//  pSolver->solutionToTxt(solPath);

  logStream << pSolver->writeResourceCostsPerNurse() << std::endl;

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
}

/******************************************************************************
* Test a solution on multiple weeks
* In this method, the weeks are solved sequentially without knowledge of future
* demand
******************************************************************************/

pair<double, int> testMultipleWeeksStochastic(
    const string& dataDir,
    const string& instanceName,
    int historyIndex,
    const vector<int>& weekIndices,
    StochasticSolverOptions stochasticSolverOptions,
    const string& outdir,
    std::vector<int> seeds) {
  // build the paths of the input files
  InputPaths inputPaths(dataDir, instanceName, historyIndex, weekIndices);

  // initialize the scenario object of the first week
  PScenario pScen = initializeScenarioINRC2(inputPaths);
  stochasticSolverOptions.setStochasticSolverOptions(
      pScen, outdir, "", stochasticSolverOptions.totalTimeLimitSeconds_, true);

  // solve the problem for each week and store the solution in the vector below
  vector<Roster> solution;
  int nbWeeks = weekIndices.size();
  Status solutionStatus;

  // whole scenario for the whole horizon
  PScenario pWholeScen = initializeScenarioINRC2(inputPaths);

  vector<PDemand> demandHistory;
  double partialCost = 0, totalCost = 0;
  int nbSched = 0;

  for (int week = 0; week < nbWeeks; week++) {
    auto rdm = Tools::getANewRandomGenerator();
    if (week >= seeds.size())
      seeds.emplace_back(rdm());
    Tools::initializeRandomGenerator(seeds[week]);

    demandHistory.push_back(pScen->pDemand());

    std::cout << pScen->toStringINRC2() << std::endl;

    auto *pSolver = new StochasticSolver(pScen,
                                         stochasticSolverOptions,
                                         demandHistory,
                                         partialCost);

    totalCost = pSolver->solve();
    nbSched += pSolver->getNGeneratedSolutions();
    partialCost += pSolver->computeSolutionCost(false);

    printf("Current cost = %.2f (partial cost = %.2f) \n\n",
           totalCost, partialCost);

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
        solution[n].pushBack(weekSolution[n]);
    }

    std::cout << "============== Current solver ===============" << std::endl;
    std::cout << pSolver->writeResourceCosts() << std::endl;

    // compute whole cost and compare to current cost
    std::cout << "============== Whole solver ===============" << std::endl;
    Solver wholeSolver(pWholeScen);
    wholeSolver.loadSolution(solution);
    double wholeCost = wholeSolver.computeSolutionCost();
    std::cout << wholeSolver.writeResourceCosts() << std::endl;

    if (abs(wholeCost - totalCost) >= 1) {
      std::cout << "Current solver:" << std::endl;
      std::cout << pSolver->writeResourceCostsINRC2() << std::endl;
      std::cout << "============================================" << std::endl;
      std::cout << "Whole solver:" << std::endl;
      std::cout << wholeSolver.writeResourceCostsINRC2() << std::endl;
      std::cout << wholeSolver.solutionToLogString() << std::endl;
      Tools::throwError("Issue with current cost (%.1f) that is different of "
                        "the real cost (%.1f).", totalCost, wholeCost);
    }

    // prepare the scenario for next week if we did not reach the last week yet
    if (week < nbWeeks - 1) {
      // Initialize demand and preferences
      PDemand pDemand(nullptr);
      PPreferences pPref(nullptr);

      // Read the demand and preferences and link them with the scenario
      ReadWrite::readWeekINRC2(inputPaths.week(week + 1),
                               pScen,
                               &pDemand,
                               &pPref);

      // read the initial state of the new week from the last state of the
      // last week
      // modify the dayId_ to show that this is the first day of the new week
      vector<State> initialStates = pSolver->statesOfDay(7);
      for (int i = 0; i < pScen->nNurses(); i++) {
        initialStates[i].dayId_ = 0;
      }

      // update the scenario to treat next week
      pScen->updateNewWeek(pDemand, pPref, initialStates);

      // append demand and preferences
      pWholeScen->pushBack(pDemand, pPref);
    } else {
      std::cout << wholeSolver.writeResourceCostsPerNurse() << std::endl;
    }

    delete pSolver;
  }

  printf("Total cost = %.2f \n", totalCost);

  string seedOutfile = outdir + "seeds.txt";
  Tools::LogOutput seedStream(seedOutfile, true);
  char str[50];
  snprintf(str, sizeof(str),
           "Cost %.2f; NbGene %d; Seeds", totalCost, nbSched);
  seedStream << str;
  for (int s : seeds)
    seedStream << " " << s;
  seedStream << std::endl;

  // Display the solution
  // displaySolutionMultipleWeeks(dataDir, instanceName, historyIndex,
  // weekIndices, solution, solutionStatus, outDir);

  return {totalCost, nbSched};
}

/******************************************************************************
* Main method
******************************************************************************/

int main(int argc, char **argv) {
  std::cout << "Number of arguments= " << argc << std::endl;

  // Detect errors in the number of arguments
  if (argc % 2 != 1) {
    Tools::throwError("main: There should be an even number of arguments!");
  } else if (argc > 1 && (argc < 9 || argc > 17)) {
    Tools::throwError(
        "main: There is either too many or not enough arguments!");
  }

  // Simulate default behavior for a test instance
  if (argc == 1) {
    std::cout << "Running the default method..." << std::endl;
    string dataDir = "datasets/INRC2/";
    string instanceName = "n030w4_1_2-7-0-9";
//    instanceName = "n060w4_1_6-1-1-5";
//    instanceName = "n030w8_1_2-7-0-9-3-6-0-6";
//    instanceName = "n110w8_0_2-1-1-7-2-6-4-7";
    std::vector<string> inst = Tools::tokenize<string>(instanceName, '_');
    int historyIndex = std::stoi(inst[1]);
    vector<int> weekIndices = Tools::tokenize<int>(inst[2], '-');
    instanceName = inst[0];
    std::vector<int> seeds = {68, 54, 78, 98, 68, 54, 78, 98};
    Tools::ThreadsPool::setMaxGlobalThreadsToMax();
    StochasticSolverOptions stochasticSolverOptions;
    stochasticSolverOptions.totalTimeLimitSeconds_ = 20;
    string outdir = "outfiles/" + instanceName + "/";
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
    InputPaths *pInputPaths = readArguments(argc, argv);

    // Initialize the random seed
    Tools::initializeRandomGenerator(pInputPaths->randSeed());

    // Solve the week
    solveOneWeek(pInputPaths);

    delete pInputPaths;
  }
}
