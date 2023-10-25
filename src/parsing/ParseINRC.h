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

#ifndef SRC_PARSING_PARSEINRC_H_
#define SRC_PARSING_PARSEINRC_H_

#include <string>
#include <vector>

#include "tools/Tools.h"
#include "data/Scenario.h"


PScenario readINRCInstance(const string &fileName);

// Read the preferences from INRC file
void readINRCPreferences(std::fstream *pFile,
                         const PScenario &pScenario,
                         const vector<PNurse> &theNurses,
                         const tm &tmStart,
                         const PPreferences &pPref);

// Read the forbidden patterns and store the corresponding resources in the
// contracts
void readINRCPatterns(const vector<PContract> &pContracts,
                      std::fstream *pFile,
                      const ShiftsFactory &shiftsFactory);

PDemand readINRCDemand(int nbDays,
                       int nbWeeks,
                       int nbSkills,
                       int nbShifts,
                       const tm &tmStart,
                       const vector<PShift> &pShifts,
                       std::fstream *pFile);

#endif  // SRC_PARSING_PARSEINRC_H_
