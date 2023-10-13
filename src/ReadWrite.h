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

#ifndef SRC_READWRITE_H_
#define SRC_READWRITE_H_

#include <ctime>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>
#include <map>
#include <vector>

#include "data/Scenario.h"

#ifdef WIN32
#ifndef NAN
const unsigned int64 nan[2]={0xffffffff, 0x7fffffff};
#define NAN (*(const double *) nan)
#endif
#endif


//--------------------------------------------------------------------------
//
//  C l a s s   R e a d W r i t e
//
//  Contains the (static) functions to read the input and write the output
//
//--------------------------------------------------------------------------

class ReadWrite {
// All functions in this class shall be public
 public:
  //--------------------------------------------------------------------------
  // Methods that read all the input files and store the content in the
  // input scenario instance
  //
  // Read the scenario file and store the content in a Scenario instance
  //
  static PScenario readScenarioINRC2(const std::string &fileName);
  static PScenario readScenarioUI(const std::string &fileName,
                                  const string &name);
  // Read several week files and
  // store the content in one demand and one preference
  //
  static PDemand readINRC2Weeks(const std::vector<std::string> &strWeekFiles,
                                const PScenario &pScenario);
  // Read the Week file and store the content in a Scenario instance
  //
  static void readWeekINRC2(const std::string &strWeekFile,
                            const PScenario &pScenario,
                            PDemand *pDemand,
                            PPreferences *pPref);

  // Read the history file
  //
  static void readHistoryINRC2(const std::string &strHistoryFile,
                               const PScenario &pScenario);

  // Read the input custom file
  // Store the result in a vector of historical demands and
  // return the number of treated weeks
  //
  static int readCustom(const std::string &strCustomInputFile,
                        const PScenario &pScenario,
                        std::vector<PDemand> *demandHistory);

  // Print the main characteristics of all the demands of an input directory
  // This is done to find some invariant properties among demands
  static void compareDemands(const std::string &inputDir, std::string logFile);

  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  // Methods that write the ouputs of the solver

  // Write the solution file for the current week
  //
  // void writeSolution(std::string strCustomOutputFile, Solution* pSolution);
  //

  // Write the output custom file from values in the scenario and the solution
  // instances
  //
  static void writeUI(const string &path, const PScenario &pScenario);
  static void writeCustom(std::string stdCustomOutputFile,
                          const std::string &strWeekFile,
                          const std::string &strCustomInputFile = "");

  //--------------------------------------------------------------------------
  static PScenario readINRCInstance(const string &fileName);

  // Read the preferences from INRC file
  static void readINRCPreferences(std::fstream *pFile,
                                  const PScenario &pScenario,
                                  const vector<PNurse> &theNurses,
                                  const tm &tmStart,
                                  const PPreferences &pPref);

  // Read the forbidden patterns and store the corresponding resources in the
  // contracts
  static void readINRCPatterns(const vector<PContract> &pContracts,
                               std::fstream *pFile,
                               const ShiftsFactory &shiftsFactory);

  static void readINRCDemand(int nbDays,
                             int nbWeeks,
                             int nbSkills,
                             int nbShifts,
                             const tm &tmStart,
                             const vector<PShift> &pShifts,
                             std::fstream *pFile,
                             PDemand *pDemand);

  static PScenario readNRPInstance(const string &fileName);
};
#endif  // SRC_READWRITE_H_
