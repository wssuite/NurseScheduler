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

#include "ParseNRP.h"

#include <algorithm>
#include <utility>
#include <memory>
#include <string>
#include <vector>

#include "solvers/mp/sp/rcspp/resources/UnwantedShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

//--------------------------------------------------------------------------
// Method that read an NRP input files and stores the content in the
// output scenario instance
//
PScenario readNRPInstance(const string &fileName) {
  std::fstream file;
  Tools::openFile(fileName, &file);
  string strTmp;
  int intTmp;

  // declare the attributes that will initialize the Scenario instance
  string instName, name;
  int nDays = 0, nWeeks = -1, nSkills = 1, nShifts = -1, nNurses = 0;
  vector<string> intToSkill = {"None"}, intToShift;
  map<string, int> skillToInt = {{"None", 0}}, shiftToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType, shiftIDToShiftTypeID;
  vector<int> shiftDurations;
  vector2D<int> forbiddenShiftSuccessors;
  vector2D<string> forbiddenShiftSuccessorsID;
  vector<int> firstShiftOfType;
  vector<PShift> pShifts;
  vector<PContract> pContracts;
  vector<PNurse> pNurses;

  // fetch instance name
  instName = Tools::tokenize<string>(fileName, '/').back();
  instName = Tools::tokenize<string>(instName, '.').front();

  // Read the horizon
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_HORIZON")) {
    Tools::throwError("The beginning of the NRP file is not as expected");
  }
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  nDays = std::stoi(strTmp);
  nWeeks = nDays / 7;

  // Read shifts
  // There is a total correspondence between shift and shift type
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_SHIFTS")) {
    Tools::throwError("The NRP file is not as expected");
  }
  // for the rest shift
  nShifts = 1;
  shiftToInt[REST_SHIFT] = REST_SHIFT_ID;
  intToShift.emplace_back(REST_SHIFT);
  shiftDurations.push_back(0);
  shiftIDToShiftTypeID.push_back(0);
  firstShiftOfType.push_back(0);
  forbiddenShiftSuccessors.emplace_back();

  // the other shifts
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // ShiftID, Length in min, Shifts which cannot follow this shift | separated
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    // ShiftID
    shiftToInt[tokens[0]] = nShifts;
    intToShift.push_back(tokens[0]);
    // Length in min
    shiftDurations.push_back(std::stoi(tokens[1]));
    // Forbidden successors
    forbiddenShiftSuccessorsID.push_back(
            Tools::tokenize<string>(tokens[2], '|'));

    // read next line
    nShifts++;
    Tools::readLine(&file, &strTmp);
  }

  // set forbidden successors and shift types
  for (int i = 1; i < nShifts; i++) {
    vector<int> forbiddenSuccessors;
    for (const string &s : forbiddenShiftSuccessorsID[i - 1])
      forbiddenSuccessors.push_back(shiftToInt[s]);
    std::sort(forbiddenSuccessors.begin(), forbiddenSuccessors.end());
    forbiddenShiftSuccessors.push_back(forbiddenSuccessors);

    // set shift types
    int type = 1;
    for (; type < firstShiftOfType.size(); type++) {
      int s = firstShiftOfType[type];
      // check if same forbidden successors
      const auto &fSucc = forbiddenShiftSuccessors.at(s);
      if (fSucc.size() != forbiddenSuccessors.size()) continue;
      // find each forbidden shift
      bool found = true;
      for (int s1 : forbiddenSuccessors) {
        // if s1 not found, break
        if (std::find(fSucc.begin(), fSucc.end(), s1) == fSucc.end()) {
          found = false;
          break;
        }
      }
      // check if all found
      if (found)
        break;
    }
    shiftIDToShiftTypeID.push_back(type);
    if (type == firstShiftOfType.size())
      firstShiftOfType.push_back(i);
  }

  // create PShifts
  for (int i = 0; i < nShifts; i++) {
    vector<int> successorList;
    auto it = forbiddenShiftSuccessors[i].begin();
    for (int j = 0; j < nShifts; j++) {
      if (it != forbiddenShiftSuccessors[i].end() && j == *it) {
        it++;
      } else {
        successorList.push_back(j);
      }
    }
    pShifts.push_back(std::make_shared<Shift>(
            intToShift[i], i, shiftIDToShiftTypeID[i],
            shiftDurations[i], successorList));
  }
  ShiftsFactory shiftFactory(pShifts);
  int maxDuration =
          *std::max_element(shiftDurations.begin(), shiftDurations.end());

  // Read nurses
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_STAFF")) {
    Tools::throwError("The NRP file is not as expected");
  }
  vector<State> initialState;
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // ID, MaxShifts, MaxTotalMinutes, MinTotalMinutes, MaxConsecutiveShifts,
    // MinConsecutiveShifts, MinConsecutiveDaysOff, MaxWeekends
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    // ID
    string nurseName = tokens[0];
    nurseNameToInt[nurseName] = nNurses;
    vector<PBaseResource> pResources;
    vector<int> availShifts;
    for (int i = 0; i < nShifts; i++) availShifts.push_back(i);
    // MaxConsecutiveShifts, MinConsecutiveShifts
    int maxCons = std::stoi(tokens[4]), minCons = std::stoi(tokens[5]);
    // DO NOT ENFORCE LB AT THE BEGINNING AND END
    pResources.push_back(std::make_shared<HardConsShiftResource>(
            minCons, maxCons, shiftFactory.pAnyWorkShift(),
            nDays, 0, false, false, false));
    // MinConsecutiveDaysOff
    int minConsDaysOff = std::stoi(tokens[6]);
    pResources.push_back(std::make_shared<HardConsShiftResource>(
            minConsDaysOff, nDays, shiftFactory.pAnyRestShift(),
            nDays, minConsDaysOff));
    // MaxShifts
    std::vector<std::pair<PShift, int>> maxShifts;
    for (const string &s : Tools::tokenize<string>(tokens[1], '|')) {
      vector<string> p = Tools::tokenize<string>(s, '=');
      PShift pS = pShifts[shiftToInt[p[0]]];
      int ub = std::stoi(p[1]);
      // just forbid the shift if ub = 0
      if (ub == 0) {
        auto it = std::find(availShifts.begin(), availShifts.end(), pS->id);
        availShifts.erase(it);
      } else if (ub < nDays) {
        auto pR = std::make_shared<HardTotalShiftDurationResource>(
                0, ub, pS, nDays, true);
        pR->maxDurationForRemainingDays(maxCons, minConsDaysOff);
        pResources.push_back(pR);
      }
    }
    // MaxTotalMinutes, MinTotalMinutes
    int maxTot = std::stoi(tokens[2]), minTot = std::stoi(tokens[3]);
    auto pR = std::make_shared<HardTotalShiftDurationResource>(
            minTot, maxTot, shiftFactory.pAnyWorkShift(),
            nDays, false, maxDuration);
    pR->maxDurationForRemainingDays(maxCons, minConsDaysOff);
    pResources.push_back(pR);
    // set min/max cons shifts
    minConsShiftType.push_back(minConsDaysOff);
    maxConsShiftType.push_back(nDays);
    minConsShiftType.resize(nShifts, minCons);
    maxConsShiftType.resize(nShifts, maxCons);
    // MaxWeekends
    int maxWE = std::stoi(tokens[7]);
    pResources.push_back(std::make_shared<HardTotalWeekendsResource>(
            0, maxWE, shiftFactory.pAnyWorkShift(), nDays));

    auto pC = std::make_shared<Contract>(
            nNurses, nurseName, minTot, maxTot, minCons, maxCons,
            minConsDaysOff, nDays, maxWE, false, pResources);
    pContracts.push_back(pC);

    vector<int> skills = {0};
    pNurses.push_back(std::make_shared<Nurse>(
            nNurses, nurseName, nSkills, skills, availShifts, pC));
    nNurses++;

    // Initialize the history of every nurse to empty state
    State nurseState(pShifts[0], 0, 0, 0, 0, 0, minConsDaysOff, minConsDaysOff);
    initialState.push_back(nurseState);

    // read next line
    Tools::readLine(&file, &strTmp);
  }

  // Read the references
  // Read hard days off
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_DAYS_OFF")) {
    Tools::throwError("The NRP file is not as expected");
  }
  PPreferences pPref = std::make_shared<Preferences>(nNurses, nDays, nShifts);
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // ID, Day
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    int nurseId = nurseNameToInt[tokens[0]];
    // add each day off
    for (int i = 1; i < tokens.size(); i++) {
      int day = std::stoi(tokens[i]);
      pPref->addShiftOff(nurseId, day, shiftFactory.pAnyWorkShift());
    }
    // read next line
    Tools::readLine(&file, &strTmp);
  }

  // Read soft shift on
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_SHIFT_ON_REQUESTS")) {
    Tools::throwError("The NRP file is not as expected");
  }
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // EmployeeID, Day, ShiftID, Weight
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    int nurseId = nurseNameToInt[tokens[0]], day = std::stoi(tokens[1]),
            s = shiftToInt[tokens[2]], w = std::stoi(tokens[3]);
    pPref->addShiftOn(nurseId, day, pShifts[s], w);
    // read next line
    Tools::readLine(&file, &strTmp);
  }

  // Read soft shift off
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_SHIFT_OFF_REQUESTS")) {
    Tools::throwError("The NRP file is not as expected");
  }
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // EmployeeID, Day, ShiftID, Weight
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    int nurseId = nurseNameToInt[tokens[0]], day = std::stoi(tokens[1]),
            s = shiftToInt[tokens[2]], w = std::stoi(tokens[3]);
    pPref->addShiftOff(nurseId, day, pShifts[s], w);
    // read next line
    Tools::readLine(&file, &strTmp);
  }

  // Read the demand
  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  if (!Tools::strStartsWith(strTmp, "SECTION_COVER")) {
    Tools::throwError("The NRP file is not as expected");
  }
  vector3D<int> demand;
  Tools::initVector3D(&demand, nDays, nShifts, nSkills);
  int underCov = -1, overCov = -1;

  Tools::readLinesUntilNotAComment(&file, &strTmp, false);
  while (!strTmp.empty()) {
    // Day, ShiftID, Requirement, Weight for under, Weight for over
    vector<string> tokens = Tools::tokenize<string>(strTmp, ',');
    int day = std::stoi(tokens[0]), s = shiftToInt[tokens[1]],
            req = std::stoi(tokens[2]),
            under = std::stoi(tokens[3]), over = std::stoi(tokens[4]);
    demand[day][s][0] = req;
    if (underCov < 0) {
      underCov = under;
      overCov = over;
    } else if (underCov != under) {
      Tools::throwError(
              "Cannot handle different under coverage weights: %d != %d.",
              underCov, under);
    } else if (overCov != over) {
      Tools::throwError(
              "Cannot handle different over coverage weights: %d != %d.",
              overCov, over);
    }
    // read next line
    Tools::readLine(&file, &strTmp);
  }
  vector<PDemand> pDemands = {
          std::make_shared<Demand>(nDays, 0, nShifts, 1, "under_requirements",
                                   demand, D_GE, underCov),
          std::make_shared<Demand>(nDays, 0, nShifts, 1, "over_requirement",
                                   demand, D_LE, overCov)
  };

  PScenario pScenario = std::make_shared<Scenario>(instName,
                                                   nWeeks,
                                                   intToSkill,
                                                   skillToInt,
                                                   pShifts,
                                                   vector<std::string>(),
                                                   pContracts,
                                                   pNurses);

  // link the scenario to the demand, preferences and history
  pScenario->linkWithDemand(pDemands);
  pScenario->linkWithPreferences(pPref);
  pScenario->setInitialState(initialState);

  return pScenario;
}
