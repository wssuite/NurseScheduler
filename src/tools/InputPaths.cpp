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
#include <iostream>
#include <string>
#include "tools/Tools.h"

/******************************************************************************
* The instances of InputPaths contain the paths of the input files of the
* problem
*******************************************************************************/
InputPaths::InputPaths() {
  instance_ = "";
  scenario_ = "";
  history_ = "";
}

InputPaths::InputPaths(const std::string &dataDir,
                       const std::string &instanceName,
                       int historyIndex,
                       std::vector<int> weekIndices,
                       const std::string &solutionPath,
                       const std::string &logPath,
                       const std::string &paramFile,
                       int timeOut,
                       int verbose,
                       int randSeed,
                       const std::string &SPType,
                       int SPStrategy,
                       const std::string &RCSPPType,
                       int nThreads,
                       int nCandidates,
                       const std::string &origin) :
    instance_(instanceName),
    scenario_(dataDir + instanceName),
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
    nThreads_(nThreads),
    nCandidates_(nCandidates) {
  if (origin.empty()) {
    guessOrigin();
  } else {
    this->origin(origin);
  }
  // if inrc2 format (should have a vector of weeks)
  if (inrc2()) {
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
}

InputPaths::InputPaths(const std::string &dataDir,
                       const std::string &instanceName,
                       std::vector<int> weekIndices,
                       const std::string &solutionPath,
                       const std::string &logPath,
                       const std::string &paramFile,
                       int timeOut,
                       int verbose,
                       int randSeed,
                       const std::string &SPType,
                       int SPStrategy,
                       const std::string &RCSPPType,
                       int nThreads,
                       int nCandidates,
                       const std::string &origin) :
    instance_(instanceName),
    scenario_(dataDir + instanceName),
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
    nThreads_(nThreads),
    nCandidates_(nCandidates) {
  if (origin.empty()) {
    guessOrigin();
  } else {
    this->origin(origin);
  }
  // if inrc2 format (should have a vector of weeks)
  if (inrc2()) {
    std::string instanceDir = dataDir + instanceName + "/";
    // initialize the scenario and history file names
    scenario_ = instanceDir + "Sc-" + instanceName + ".txt";

    // initialize the file names for each week demand
    for (int week : weekIndices)
      weeks_.emplace_back(instanceDir + "WD-" + instanceName + "-" +
          std::to_string(week) + ".txt");
  }
}

void InputPaths::guessOrigin() {
  origin_ = guessOrigin(scenario_);
}

InstanceOrigin InputPaths::guessOrigin(const std::string &path) const {
  std::size_t found = path.find_last_of('.');
  if (found == std::string::npos) {
    return INRCII;
  }

  std::fstream file;
  file.open(path.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read the file " << path << std::endl;
    std::cout << "The input file was not opened properly!" << std::endl;
    throw Tools::myException("The input file was not opened properly!",
                             __LINE__);
  }
  std::string strTmp;
  std::string l;
  // remove comments and empty lines
  // analyse the first meaningful line
  while (std::getline(file, l) && file.good()) {
    if (l.empty() || l[0] == '/' || l[0] == '#' || l[0] == '\n') {
      continue;
    } else {
      if (l.find("SCENARIO") != std::string::npos)
        return INRCII;
      if (l.find("SCHEDULING_PERIOD;") != std::string::npos)
        return INRC;
      if (l.find("SECTION_HORIZON") != std::string::npos)
        return NRP;
      if (l.find("HEADERS") != std::string::npos)
        return UI;
    }
  }

  Tools::throwException("The origin of the instance has not been recognized.");
  return UI;
}
