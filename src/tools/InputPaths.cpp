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

#include "tools/InputPaths.h"

/******************************************************************************
* The instances of InputPaths contain the paths of the input files of the
* problem
*******************************************************************************/
InputPaths::InputPaths() {
  instance_ = "";
  scenario_ = "";
  history_ = "";
}

InputPaths::InputPaths(std::string dataDir,
                       std::string instanceName,
                       int historyIndex,
                       std::vector<int> weekIndices,
                       std::string solutionPath,
                       std::string logPath,
                       std::string paramFile,
                       double timeOut,
                       int verbose,
                       int randSeed,
                       std::string SPType,
                       int SPStrategy,
                       std::string RCSPPType,
                       int nThreads) :
    instance_(instanceName),
    historyIndex_(historyIndex),
    weekIndices_(weekIndices),
    solutionPath_(solutionPath),
    logPath_(logPath),
    paramFile_(paramFile),
    verbose_(verbose),
    randSeed_(randSeed),
    timeOut_(timeOut),
    SPType_(SPType),
    SPStrategy_(SPStrategy),
    RCSPPType_(RCSPPType),
    nThreads_(nThreads) {
  std::string instanceDir = dataDir + instanceName + "/";
  // initialize the scenario and history file names
  scenario_ = instanceDir + "Sc-" + instanceName + ".txt";
  history_ = instanceDir + "H0" + "-" + instanceName + "-"
      + std::to_string(historyIndex) + ".txt";

  // initialize the file names for each week demand
  for (int week : weekIndices)
    weeks_.emplace_back(instanceDir + "WD-" + instanceName + "-" +
        std::to_string(week) + ".txt");
}
