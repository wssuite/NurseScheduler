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

#ifndef SRC_TOOLS_INPUTPATHS_H_
#define SRC_TOOLS_INPUTPATHS_H_

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// The instances of this class contain the paths of the input paths of the
// problem
class InputPaths {
 public:
  InputPaths();
  InputPaths(std::string dataDir,
             std::string instanceName,
             int historyIndex,
             std::vector<int> weekIndices,
             std::string solutionPath = "",
             std::string logPath = "",
             std::string paramFile = "",
             double timeOut = 3600.0,
             int randSeed = 0,
             std::string SPType = "LONG",
             int SPStrategy = 1,
             int nThreads = 1);

 protected:
  std::string instance_;
  std::string scenario_;
  int historyIndex_;
  std::string history_;
  std::vector<std::string> weeks_;
  std::vector<int> weekIndices_;
  std::string solutionPath_ = "";
  std::string logPath_ = "";
  std::string paramFile_ = "";
  int randSeed_ = 0;
  double timeOut_ = 3600;
  std::string SPType_ = "LONG";
  int SPStrategy_ = 0;
  int nThreads_ = 1;

 public:
  // get/set attributes
  const std::string& instance() const { return instance_; }
  const std::string& scenario() const { return scenario_; }
  void scenario(std::string scenario) { scenario_ = scenario; }

  int historyIndex() const { return historyIndex_; }
  const std::string& history() const { return history_; }
  void history(std::string history) { history_ = history; }

  const std::vector<std::string>& weeks() const { return weeks_; }
  const std::string& week(int w) const { return weeks_[w]; }
  int weekIndex(int w) const { return weekIndices_[w]; }
  int nbWeeks() const { return weeks_.size(); }
  void addWeek(std::string week) { weeks_.push_back(week); }

  const std::string& paramFile() const { return paramFile_; }
  void paramFile(std::string file) { paramFile_ = file; }
  const std::string& solutionPath() const { return solutionPath_; }
  void solutionPath(std::string path) { solutionPath_ = path; }
  const std::string& logPath() const { return logPath_; }
  void logPath(std::string path) { logPath_ = path; }

  int randSeed() const { return randSeed_; }
  void randSeed(int seed) { randSeed_ = seed; }
  double timeOut() const { return timeOut_; }
  void timeOut(double t) { timeOut_ = t; }

  const std::string& SPType() const { return SPType_; }
  void SPType(std::string t) { SPType_ = t; }

  int SPStrategy() const { return SPStrategy_; }
  void SPStrategy(int t) { SPStrategy_ = t; }

  int nThreads() const { return nThreads_; }
  void nThreads(int n) { nThreads_ = n; }
};

#endif  // SRC_TOOLS_INPUTPATHS_H_
