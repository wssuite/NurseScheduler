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

#ifndef SRC_PARSING_PARSEINRC2_H_
#define SRC_PARSING_PARSEINRC2_H_

#include <vector>
#include <string>

#include "tools/Tools.h"
#include "data/Scenario.h"

//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//
// Read the scenario file and store the content in a Scenario instance
//
PScenario readScenarioINRC2(const std::string &fileName);

// Read several week files and
// store the content one preference and return demands
vector<PDemand> readINRC2Weeks(
        const std::vector<std::string> &strWeekFiles,
        const PScenario &pScenario);

// Read the Week file and store the content in a Scenario instance
vector<PDemand> readWeekINRC2(
        const std::string &strWeekFile,
        const PScenario &pScenario,
        PPreferences *pPref);

// Read the history file
//
void readHistoryINRC2(const std::string &strHistoryFile,
                             const PScenario &pScenario);

// Read the input custom file
// Store the result in a vector of historical demands and
// return the number of treated weeks
//
int readCustom(const std::string &strCustomInputFile,
                      const PScenario &pScenario,
                      vector2D<PDemand> *pDemandsHistory);

// Print the main characteristics of all the demands of an input directory
// This is done to find some invariant properties among demands
void compareDemands(const std::string &inputDir, std::string logFile);

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
void writeCustom(std::string stdCustomOutputFile,
                        const std::string &strWeekFile,
                        const std::string &strCustomInputFile = "");


#endif  // SRC_PARSING_PARSEINRC2_H_
