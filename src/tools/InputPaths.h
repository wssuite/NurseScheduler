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

#ifndef SRC_TOOLS_INPUTPATHS_H_
#define SRC_TOOLS_INPUTPATHS_H_

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <string>

#include "tools/Tools.h"

// Origin of the instance: it will fix the format of the input file
// NA stands for Not Available. It will be guessed based on the inputs
enum InstanceOrigin { NA, INRC, INRCII, NRP, UI };

static std::map<std::string, InstanceOrigin> instanceOriginsByName = {
        {"INRC", INRC},
        {"INRCII", INRCII},
        {"NRP", NRP},
        {"UI", UI},
        {"", NA}
};
static std::map<InstanceOrigin, std::string> namesByInstanceOrigin =
        Tools::buildNamesByType(instanceOriginsByName);

// The instances of this class contain the paths of the input paths of the
// problem
class InputPaths {
 public:
  InputPaths() {}

  InputPaths(const std::string &dataDir,
             const std::string &instanceName,
             int historyIndex = -1,
             std::vector<int> weekIndices = {},
             const std::string &solutionPath = "",
             const std::string &logPath = "",
             const std::string &paramFile = "",
             int timeOut = -1,
             int verbose = -1,
             int randSeed = 0,
             const std::string &SPType = "",
             int SPStrategy = -1,
             const std::string &RCSPPType = "",
             int nThreads = -1,
             int nCandidates = -1,
             const std::string &origin = "");

  InputPaths(const std::string &dataDir,
             const std::string &instanceName,
             std::vector<int> weekIndices,
             const std::string &solutionPath = "",
             const std::string &logPath = "",
             const std::string &paramFile = "",
             int timeOut = -1,
             int verbose = -1,
             int randSeed = 0,
             const std::string &SPType = "",
             int SPStrategy = -1,
             const std::string &RCSPPType = "",
             int nThreads = -1,
             int nCandidates = -1,
             const std::string &origin = "");

 protected:
  std::string instance_;
  InstanceOrigin origin_ = NA;
  std::string scenario_;
  int historyIndex_;
  std::string history_;
  std::vector<std::string> weeks_;
  std::vector<int> weekIndices_;
  std::string solutionPath_;
  std::string logPath_;
  std::string customInputFile_;
  std::string customOutputFile_;
  std::string paramFile_;
  int verbose_ = -1;
  int randSeed_ = 0;
  double timeOut_ = 3600;
  std::string SPType_ = "ROSTER";
  int SPStrategy_ = 0;
  std::string RCSPPType_ = "DEFAULT";
  int nThreads_ = 1;
  int nCandidates_ = -1;

 public:
  // get/set attributes
  const std::string &instance() const { return instance_; }
  const std::string &scenario() const { return scenario_; }
  void scenario(const std::string &scenario) { scenario_ = scenario; }

  void origin(const std::string &origin) {
    origin_ = instanceOriginsByName.at(Tools::toUpperCase(origin));
  }
  void origin(InstanceOrigin origin) { origin_ = origin; }
  bool inrc() const { return origin_ == INRC; }
  bool inrc2() const { return origin_ == INRCII; }
  bool nrp() const { return origin_ == NRP; }
  bool ui() const { return origin_ == UI; }
  const std::string &origin() const {
    return namesByInstanceOrigin.at(origin_);
  }

  void guessOrigin();
  InstanceOrigin guessOrigin(const std::string &path) const;

  int historyIndex() const { return historyIndex_; }
  const std::string &history() const { return history_; }
  void history(const std::string &history) { history_ = history; }

  const std::vector<std::string> &weeks() const { return weeks_; }
  const std::string &week(int w) const { return weeks_[w]; }
  int weekIndex(int w) const { return weekIndices_[w]; }
  int nbWeeks() const { return weeks_.size(); }
  void addWeek(const std::string &week) { weeks_.push_back(week); }

  const std::string &paramFile() const { return paramFile_; }
  void paramFile(const std::string &file) { paramFile_ = file; }
  const std::string &solutionPath() const { return solutionPath_; }
  void solutionPath(const std::string &path) { solutionPath_ = path; }
  const std::string &logPath() const { return logPath_; }
  void logPath(const std::string &path) { logPath_ = path; }
  const std::string &customInputFile() const { return customInputFile_; }
  void customInputFile(const std::string &file) { customInputFile_ = file; }
  const std::string &customOutputFile() const { return customOutputFile_; }
  void customOutputFile(const std::string &file) { customOutputFile_ = file; }

  int randSeed() const { return randSeed_; }
  void randSeed(int seed) { randSeed_ = seed; }
  int timeOut() const { return timeOut_; }
  void timeOut(int t) { timeOut_ = t; }

  int verbose() const { return verbose_; }
  void verbose(int verbose) { verbose_ = verbose; }

  const std::string &SPType() const { return SPType_; }
  void SPType(const std::string &t) { SPType_ = t; }

  int SPStrategy() const { return SPStrategy_; }
  void SPStrategy(int t) { SPStrategy_ = t; }

  const std::string &RCSPPType() const { return RCSPPType_; }
  void RCSPPType(const std::string &t) { RCSPPType_ = t; }

  int nThreads() const { return nThreads_; }
  void nThreads(int n) { nThreads_ = n; }

  int nCandidates() const { return nCandidates_; }
  void nCandidates(int n) { nCandidates_ = n; }
};

#endif  // SRC_TOOLS_INPUTPATHS_H_
