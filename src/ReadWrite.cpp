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

#include "ReadWrite.h"

#include <dirent.h>
#include <cmath>
#include <ctime>
#include <cstring>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <regex>// NOLINT [build/c++11]

#include "tools/Tools.h"
#include "tools/parseToolsUI.h"
#include <boost/algorithm/string.hpp>
#include "data/Scenario.h"
#include "solvers/Solver.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "solvers/mp/sp/rcspp/resources/AlternativeShiftResource.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

//--------------------------------------------------------------------------
// Method that read an NRP input files and stores the content in the
// output scenario instance
//
PScenario ReadWrite::readNRPInstance(const string &fileName) {
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
  intToShift.push_back(REST_SHIFT);
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
            0, ub, pS, nDays, 1, 1);
        pR->maxDurationForRemainingDays(maxCons, minConsDaysOff);
        pResources.push_back(pR);
      }
    }
    // MaxTotalMinutes, MinTotalMinutes
    int maxTot = std::stoi(tokens[2]), minTot = std::stoi(tokens[3]);
    auto pR = std::make_shared<HardTotalShiftDurationResource>(
        minTot, maxTot, shiftFactory.pAnyWorkShift(), nDays, maxDuration);
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
        maxWE, shiftFactory.pAnyWorkShift(), nDays));

    auto pC = std::make_shared<Contract>(
        nNurses, nurseName, minTot, maxTot, minCons, maxCons,
        minConsDaysOff, nDays, maxWE, false, pResources);
    pContracts.push_back(pC);

    vector<int> skills = {0};
    pNurses.push_back(std::make_shared<Nurse>(
        nNurses, nurseName, nShifts, nSkills, skills, availShifts, pC));
    nNurses++;

    // Initialize the history of every nurse to empty state
    State nurseState(0, 0, 0, 0, minConsDaysOff, minConsDaysOff, pShifts[0]);
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
  PDemand pDemand =
      std::make_shared<Demand>(nDays, 0, nShifts, 1, "requirement", demand);

  auto pWeights = std::shared_ptr<Weights>(new Weights(
      underCov, 0, 0, 0, 0, {}, 0, 0, 0, overCov));

  PScenario pScenario = std::make_shared<Scenario>(instName,
                                                   nWeeks,
                                                   nSkills,
                                                   intToSkill,
                                                   skillToInt,
                                                   pShifts,
                                                   vector<std::string>(),
                                                   minConsShiftType,
                                                   maxConsShiftType,
                                                   nNurses,
                                                   pContracts,
                                                   nNurses,
                                                   pNurses,
                                                   nurseNameToInt,
                                                   pWeights,
                                                   "",
                                                   true,
                                                   false);

  // link the scenario to the demand, preferences and history
  pScenario->linkWithDemand(pDemand);
  pScenario->linkWithPreferences(pPref);
  pScenario->setInitialState(initialState);

  return pScenario;
}

//--------------------------------------------------------------------------
// Method that read an INRC input files and stores the content in the
// output scenario instance
//
PScenario ReadWrite::readINRCInstance(const string &fileName) {
  std::fstream file;
  Tools::openFile(fileName, &file);
  string strTmp;
  int intTmp;

  // declare the attributes that will initialize the Scenario instance
  string instName, name;
  int nWeeks = -1, nSkills = -1, nShifts = -1, nShiftTypes = -1,
      nContracts = -1, nNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType;
  map<string, int> skillToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType, shiftIDToShiftTypeID;
  vector2D<int> forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  vector<PContract> pContracts;
  vector<PNurse> pNurses;

  // Read the instName of the instance, and the scheduling horizon dates
  Tools::readUntilChar(&file, ';', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SCHEDULING_PERIOD")) {
    Tools::throwError("The beginning of the INRC file is not as expected");
  }

  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  file >> instName;
  instName.pop_back();
  Tools::readUntilChar(&file, ',', &strTmp);
  const std::tm tmStart(*Tools::readDateFromStr(strTmp));
  Tools::readUntilChar(&file, ';', &strTmp);
  const std::tm tmEnd(*Tools::readDateFromStr(strTmp));

  // set the number of days and weeks and the first week day in the horizon
  int nDays(0);
  nDays = tmEnd.tm_yday - tmStart.tm_yday + 1;
  nWeeks = std::ceil(nDays / 7);

  Day::setFirstDayOfWeek(
      (DayOfWeek) (tmStart.tm_wday >= 0 ? tmStart.tm_wday - 1 : 6));
  std::cout << "Starting week day = " << dayOfWeekToName(Day::firstDayOfWeek())
            << " ; number of days = " << nDays
            << " ; number of weeks = " << nWeeks;

  if (nDays < 0)
    Tools::throwError("The start of the horizon should be before its end");
  if (tmEnd.tm_year != tmStart.tm_year)
    Tools::throwError("The beginning and end of the scheduling horizon should"
                      " be on the same year");

  // set the list of skills
  Tools::readUntilChar(&file, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SKILLS "))
    Tools::throwError("The INRC file is not as expected");

  file >> nSkills;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  for (int i = 0; i < nSkills; i++) {
    file >> strTmp;
    strTmp.pop_back();
    intToSkill.push_back(strTmp);
    skillToInt[strTmp] = i;
  }

  // set the list of shift types: here there is major conversion to make,
  // since in INRC, a shift type also specifies a list of skills that can be
  // assigned to it
  //
  Tools::readUntilChar(&file, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SHIFT_TYPES "))
    Tools::throwError("The INRC file is not as expected");

  file >> nShifts;
  nShifts += 1;  // add the rest shift

  // in INRC, any shift can come after any other shift (possibly at the
  // cost of a penalty)
  vector<int> successorList;
  successorList.reserve(nShifts);
  for (int i = 0; i < nShifts; i++)
    successorList.push_back(i);

  // IMPORTANT : INSERT REST SHIFT !!!!!!
  intToShiftType.push_back(REST_SHIFT);
  // initialize the shift
  pShifts.push_back(std::make_shared<Shift>(
      REST_SHIFT, 0, 0, 0, successorList));
  nShiftTypes = 1;
  shiftIDToShiftTypeID.push_back(0);

  // read the name, type, time interval and skills of the shifts
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  string shName, typeName, skName;
  vector<int> skillIds;
  for (int i = 1; i < nShifts; i++) {
    // get the name and type name of the shift: the type name is cut to
    // keep only the first day
    file >> shName;
    shName.pop_back();
    Tools::readUntilChar(&file, ',', &strTmp);
    std::stringstream sstream(strTmp);
    sstream >> typeName;

    // get the time interval and duration of the shift
    Tools::readUntilChar(&file, ',', &strTmp);
    Tools::Time timeStart(Tools::readHourFromStr(strTmp));
    Tools::readUntilChar(&file, ',', &strTmp);
    Tools::Time timeEnd(Tools::readHourFromStr(strTmp));
    // look at previous shift types to see if it has already been created
    int type = -1;
//    for (int j = 1; j < i; j++) {
//      if (timeStart.equals(pShifts[j]->timeStart) &&
//          timeEnd.equals(pShifts[j]->timeEnd)) {
//        type = pShifts[j]->type;
//        break;
//      }
//    }
    if (type == -1) {
      type = nShiftTypes++;
      intToShiftType.push_back(typeName);
      shiftTypeToInt[typeName] = type;
    }
    shiftIDToShiftTypeID.push_back(type);


    // get the skills of the shift
    int nSk;
    skillIds.clear();
    file >> nSk;
    Tools::readUntilChar(&file, ',', &strTmp);
    for (int j = 0; j < nSk; j++) {
      file >> skName;
      skName.pop_back();
      skillIds.push_back(skillToInt.at(skName));
    }

    // initialize the shift
    pShifts.push_back(std::make_shared<Shift>(
        shName, i, type, 1, successorList, skillIds, timeStart, timeEnd));
  }
  forbiddenShiftSuccessors.resize(nShifts);

  ShiftsFactory shiftsFactory(pShifts);

  // Read the constraints of each contract
  Tools::readUntilChar(&file, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "CONTRACTS "))
    Tools::throwError("The INRC file is not as expected");

  file >> nContracts;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  // Read each contract type
  for (int i = 0; i < nContracts; i++) {
    string contractName;
    bool isOn;
    int lbOn, lbCost, lb, ubOn, ubCost, ub;
    file >> intTmp;
    if (intTmp != i)
      Tools::throwError("Contract ids are not ordered");
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> contractName;
    contractName.pop_back();

    // first field is for single assignment which is a hard constraint
    // according to documentation, we skip it
    Tools::readUntilChar(&file, ')', &strTmp);

    // Initialize all the resources of the contract
    // some characteristics of each contract will need to be specified at the
    // nurse level
    //
    vector<PBaseResource> pResources;

    // following pair of fields is for total assignments
    Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                               &ubOn, &ub, &ubCost);
    if (lbOn || ubOn) {
      pResources.push_back(std::make_shared<SoftTotalShiftDurationResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nDays,
          lbCost,
          ubCost,
          shiftsFactory.pAnyWorkShift(),
          nDays,
          1));
    }

    // consecutive working days
    Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                               &ubOn, &ub, &ubCost);
    if (lbOn || ubOn) {
      pResources.push_back(std::make_shared<SoftConsShiftResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nDays,
          lbCost,
          ubCost,
          shiftsFactory.pAnyWorkShift(),
          CONS_WORK_COST,
          nDays,
          0,
          true));
    }


    // consecutive free days
    if (lbOn || ubOn) {
      Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                                 &ubOn, &ub, &ubCost);
      pResources.push_back(std::make_shared<SoftConsShiftResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nDays,
          lbCost,
          ubCost,
          shiftsFactory.pAnyRestShift(),
          CONS_REST_COST,
          nDays,
          0,
          true));
    }

    // Below, we must record the properties of two weekend resources before
    // getting the corresponding definition of a weekend
    //
    // consecutive working weekends
    Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                               &ubOn, &ub, &ubCost);

    // maximum total number of weekends
    int totalWkndOn, totalWkndUb, totalWkndCost;
    Tools::readUbResource(&file, &totalWkndOn, &totalWkndUb, &totalWkndCost);

    // description of a weekend
    Tools::readUntilChar(&file, ' ', &strTmp);
    Tools::readUntilChar(&file, ',', &strTmp);
    vector<string> weekendDays;
    string weekDay;
    weekDay.push_back(strTmp.front());
    for (int j = 1; j < strTmp.size(); j++) {
      if (isupper(strTmp[j])) {
        weekendDays.push_back(weekDay);
        weekDay.clear();
      }
      weekDay.push_back(strTmp[j]);
    }
    weekendDays.push_back(weekDay);
    DayOfWeek firstWeekendDay = nameToDayOfWeek(weekendDays.front());
    DayOfWeek lastWeekendDay = nameToDayOfWeek(weekendDays.back());

    // initialize the consecutive and total weekend resources
    if (lbOn || ubOn) {
      pResources.push_back(std::make_shared<SoftConsWeekendShiftResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nWeeks,
          lbCost,
          ubCost,
          shiftsFactory.pAnyWorkShift(),
          nDays,
          firstWeekendDay,
          lastWeekendDay,
          true));
    }

    if (totalWkndOn) {
      pResources.push_back(std::make_shared<SoftTotalWeekendsResource>(
          totalWkndUb,
          totalWkndCost,
          shiftsFactory.pAnyWorkShift(),
          nDays,
          firstWeekendDay,
          lastWeekendDay));
    }

    // complete weekend constraint
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isOn;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> ubCost;
    if (isOn) {
      pResources.push_back(std::make_shared<SoftIdentWeekendResource>(
          std::make_shared<ShiftWorkComparator>(), ubCost,
          firstWeekendDay, lastWeekendDay));
    }

    // identical shift types during weekend
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isOn;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> ubCost;
    if (isOn) {
      pResources.push_back(std::make_shared<SoftIdentWeekendResource>(
          std::make_shared<ShiftTypeComparator>(), ubCost,
          firstWeekendDay, lastWeekendDay));
    }

    // no night shift before free weekends : this will be treated as a pattern
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isOn;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> ubCost;
    bool nightFound = false;
    for (const auto &p : shiftTypeToInt)
      if (p.first == "Night") {
        nightFound = true;
        break;
      }
    if (isOn && nightFound) {
      // create the forbidden pattern: one night on the day before weekend and
      vector<PAbstractShift> pAShifts;
      vector<PAbstractDay> pADays;
      if (firstWeekendDay == MONDAY)
        pADays.push_back(std::make_shared<WeekDay>(SUNDAY));
      else
        pADays.push_back(
            std::make_shared<WeekDay>((DayOfWeek) (firstWeekendDay - 1)));
      if (lastWeekendDay >= firstWeekendDay) {
        for (int d = firstWeekendDay; d <= lastWeekendDay; d++)
          pADays.push_back(std::make_shared<WeekDay>((DayOfWeek) (d)));
      } else {
        for (int d = firstWeekendDay; d <= SUNDAY; d++)
          pADays.push_back(std::make_shared<WeekDay>((DayOfWeek) (d)));
        for (int d = MONDAY; d <= lastWeekendDay; d++)
          pADays.push_back(std::make_shared<WeekDay>((DayOfWeek) (d)));
      }
      pAShifts.push_back(
          shiftsFactory.pAnyTypeShift(shiftTypeToInt.at("Night")));
      for (int d = 1; d < pADays.size(); d++)
        pAShifts.push_back(shiftsFactory.pAnyRestShift());
      pResources.push_back(std::make_shared<SoftForbiddenPatternResource>(
          Pattern(pAShifts, pADays),
          ubCost));
    }

    // two free days after night shift
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isOn;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> ubCost;
    if (isOn && nightFound) {
      PAbstractShift pAShift =
          shiftsFactory.pAnyTypeShift(shiftTypeToInt.at("Night"));
      pResources.push_back(std::make_shared<SoftFreeDaysAfterShiftResource>(
          pAShift, 2, ubCost));
    }

    // Alternative skill category
    // Use alternative shift and not alternative skills as shifts
    // include skills in INRC
    bool isAlternativeShift;
    double alternativeShiftCost;
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isAlternativeShift;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> alternativeShiftCost;
    // when alternative is false, any skill can be used by any nurse
    if (!isAlternativeShift) {
      isAlternativeShift = true;
      alternativeShiftCost = 0;
    }

    // Forbidden pattern ids
    int nPatterns;
    vector<int> patternIds;
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> nPatterns;
    Tools::readUntilChar(&file, ',', &strTmp);
    for (int j = 0; j < nPatterns; j++) {
      file >> intTmp;
      patternIds.push_back(intTmp);
    }

    pContracts.push_back(std::make_shared<Contract>(
        i,
        contractName,
        pResources,
        firstWeekendDay,
        lastWeekendDay,
        patternIds,
        false,
        0,
        isAlternativeShift,
        alternativeShiftCost));
    Tools::readUntilChar(&file, ';', &strTmp);
  }

  // Read the forbidden patterns and store the corresponding resources in the
  // contracts
  readINRCPatterns(pContracts, &file, shiftsFactory);

  // Read the list of nurses
  // BEWARE: there is no equal sign before the number of employees!
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  Tools::readUntilChar(&file, ' ', &strTmp);
  if (!Tools::strEndsWith(strTmp, "EMPLOYEES"))
    Tools::throwError("The INRC file is not as expected");

  file >> nNurses;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  for (int i = 0; i < nNurses; i++) {
    // BEWARE: nurse ids are in general not ordered in the INRC files!
    int nurseId;
    file >> nurseId;
    Tools::readUntilChar(&file, ',', &strTmp);
    string nurseName;
    Tools::readUntilChar(&file, ',', &nurseName);
    file >> intTmp;
    PContract pContract = pContracts[intTmp];
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> intTmp;
    Tools::readUntilChar(&file, ' ', &strTmp);
    vector<int> skills, availableShifts;
    // read the skills
    for (int j = 0; j < intTmp; j++) {
      Tools::readUntilOneOfTwoChar(&file, ' ', ';', &strTmp);
      skills.push_back(skillToInt.at(strTmp));
    }
    // sort the skill indices before initializing the nurse
    std::sort(skills.begin(), skills.end());

    // add the rest shift and sort the available shifts
    for (const auto &pS : pShifts)
      if (std::any_of(skills.begin(), skills.end(),
                      [pS](int sk) { return pS->hasSkill(sk); }))
        availableShifts.push_back(pS->id);
    std::sort(availableShifts.begin(), availableShifts.end());
    // check if rest is available
    if (availableShifts.front() != 0)
      Tools::throwError("Rest shift is not available for nurse %d", nurseId);

    // as the skills are included in the shift,
    // give all the skills to each nurses
    skills.clear();
    for (int i = 0; i < nSkills; ++i) skills.push_back(i);
    pNurses.push_back(std::make_shared<Nurse>(
        nurseId, nurseName, nShifts, nSkills, skills, availableShifts,
        pContract));
    nurseNameToInt[nurseName] = i;

    if (pContract->alternativeShift_)
      pNurses.back()->addBaseResource(std::make_shared<
          AlternativeShiftResource>(
          pNurses.back(), nShifts));
  }
  // sort the nurses per id to avoid errors
  std::stable_sort(pNurses.begin(),
                   pNurses.end(),
                   compareNursesById);

  // Check that all fields were initialized before initializing the scenario
  //
  if (nWeeks == -1 || nSkills == -1 || nShifts == -1 || nContracts == -1
      || nNurses == -1) {
    Tools::throwError("In readINRCInstances: missing fields in the "
                      "initialization");
  }

  PScenario pScenario = std::make_shared<Scenario>(instName,
                                                   nWeeks,
                                                   nSkills,
                                                   intToSkill,
                                                   skillToInt,
                                                   pShifts,
                                                   intToShiftType,
                                                   minConsShiftType,
                                                   maxConsShiftType,
                                                   nContracts,
                                                   pContracts,
                                                   nNurses,
                                                   pNurses,
                                                   nurseNameToInt,
                                                   std::make_shared<Weights>(),
                                                   "",
                                                   true,
                                                   false);
  pScenario->setStartDate(tmStart);


  // Read the demand
  PDemand pDemand;
  readINRCDemand(nDays,
                 nWeeks,
                 nSkills,
                 nShifts,
                 tmStart,
                 pShifts,
                 &file,
                 &pDemand);

  // Read preferences
  PPreferences pPref = std::make_shared<Preferences>(nNurses,
                                                     nDays,
                                                     nShifts);
  readINRCPreferences(&file,
                      pScenario,
                      pNurses,
                      tmStart,
                      pPref);

  // Initialize the history of every nurse to empty state
  vector<State> initialState;
  // Add a fictitious shift just for the initial state
  const PShift &pNoneShift =
      shiftsFactory.pNoneShift()->pIncludedShifts().front();
  for (int n = 0; n < nNurses; n++) {
    State nurseState(0, 0, 0, 0, 0, 0, pNoneShift);
    initialState.push_back(nurseState);
  }

  // link the scenario to the demand, preferences and history
  pScenario->linkWithDemand(pDemand);
  pScenario->linkWithPreferences(pPref);
  pScenario->setInitialState(initialState);

  return pScenario;
}

void ReadWrite::readINRCDemand(int nDays,
                               int nWeeks,
                               int nSkills,
                               int nShifts,
                               const tm &tmStart,
                               const vector<PShift> &pShifts,
                               std::fstream *pFile,
                               PDemand *pDemand) {
  // 1. read weekly demands
  string strTmp;
  int intTmp;
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DAY_OF_WEEK_COVER "))
    Tools::throwError("The INRC file is not as expected: "
                      "DAY_OF_WEEK_COVER not found");

  int nCovers;
  *pFile >> nCovers;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);

  // init the demand vectors
  map<string, int> shiftToInt;
  for (const auto &pS : pShifts) shiftToInt[pS->name] = pS->id;
  vector3D<int> minDemands;
  Tools::initVector3D(&minDemands,
                      nDays,
                      nShifts,
                      nSkills,
                      0);
  for (int i = 0; i < nCovers; i++) {
    *pFile >> strTmp;
    strTmp.pop_back();
    DayOfWeek dayOfWeek = nameToDayOfWeek(strTmp);
    *pFile >> strTmp;
    strTmp.pop_back();
    int shiftId = shiftToInt.at(strTmp);
    *pFile >> intTmp;
    if (pShifts[shiftId]->skills.size() != 1)
      Tools::throwError("A shift can include only one skill.");
    int skillId = pShifts[shiftId]->skills.front();
    for (int w = 0; w < nWeeks; w++) {
      int dayId = Day::getDayId(dayOfWeek, w);
      minDemands[dayId][shiftId][skillId] = intTmp;
    }
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // 2. read specific demands
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DATE_SPECIFIC_COVER "))
    Tools::throwError("The INRC file is not as expected: "
                      "DATE_SPECIFIC_COVER not found.");

  *pFile >> nCovers;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);

  for (int i = 0; i < nCovers; i++) {
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmDemand(*Tools::readDateFromStr(strTmp));
    int dayId(tmDemand.tm_yday - tmStart.tm_yday);
    *pFile >> strTmp;
    strTmp.pop_back();
    int shiftId = shiftToInt.at(strTmp);
    *pFile >> intTmp;
    if (pShifts[shiftId]->skills.size() != 1)
      Tools::throwError("A shift can include only one skill.");
    int skillId = pShifts[shiftId]->skills.front();
    minDemands[dayId][shiftId][skillId] = intTmp;
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // initialize pDemand
  *pDemand = std::make_shared<Demand>(nDays,
                                      0,
                                      nShifts,
                                      nSkills,
                                      "allWeeks",
                                      minDemands);
}

void ReadWrite::readINRCPatterns(const vector<PContract> &pContracts,
                                 std::fstream *pFile,
                                 const ShiftsFactory &shiftsFactory) {
  vector<SoftForbiddenPatternResource> patternResources;

  string strTmp;
  int intTmp;
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "PATTERNS "))
    Tools::throwError("The INRC file is not as expected: "
                      "PATTERNS not found.");

  *pFile >> intTmp;
  int nPatterns = intTmp;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  // Read the patterns and create one resources for each pattern
  map<string, int> shiftToInt;
  for (const auto &pS : shiftsFactory.pShifts())
    shiftToInt[pS->name] = pS->id;
  for (int i = 0; i < nPatterns; i++) {
    *pFile >> intTmp;
    if (intTmp != i)
      Tools::throwError("Pattern ids are not ordered");
    Tools::readUntilChar(pFile, ',', &strTmp);
    double cost;
    *pFile >> cost;
    Tools::readUntilChar(pFile, ',', &strTmp);
    *pFile >> intTmp;
    std::vector<PAbstractShift> pAShifts;
    std::vector<PAbstractDay> pADays;
    for (int j = 0; j < intTmp; j++) {
      Tools::readUntilChar(pFile, '(', &strTmp);
      Tools::readUntilChar(pFile, '|', &strTmp);
      if (strTmp == "Any") {
        pAShifts.push_back(shiftsFactory.pAnyWorkShift());
      } else if (strTmp == "None") {
        pAShifts.push_back(shiftsFactory.pAnyRestShift());
      } else {
        int shiftType = shiftToInt.at(strTmp);
        pAShifts.push_back(shiftsFactory.pAnyTypeShift(shiftType));
      }

      Tools::readUntilChar(pFile, ')', &strTmp);
      if (strTmp == "Any") {
        pADays.push_back(std::make_shared<AnyDay>());
      } else {
        DayOfWeek dayOfWeek = nameToDayOfWeek(strTmp);
        pADays.push_back(std::make_shared<WeekDay>(dayOfWeek));
      }
    }

    patternResources.emplace_back(
        SoftForbiddenPatternResource(
            Pattern(pAShifts, pADays), cost));
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // below, we create different resource pointers for each contract to make
  // sure that resources are not shared among nurses with different contracts
  for (const auto &pC : pContracts) {
    for (int patternId : pC->forbiddenPatternIds_) {
      SoftForbiddenPatternResource r = patternResources[patternId];
      pC->addResource(std::make_shared<SoftForbiddenPatternResource>(
          r.getPattern(), r.getCost()));
    }
  }
}

// Read the preferences of all the nurses in an INRC input file
void ReadWrite::readINRCPreferences(std::fstream *pFile,
                                    const PScenario &pScenario,
                                    const vector<PNurse> &pNurses,
                                    const tm &tmStart,
                                    const PPreferences &pPref) {
  // 1. ready day on preferences
  string strTmp;
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DAY_OFF_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "DAY_OFF_REQUESTS not found.");

  int nRequests;
  *pFile >> nRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nRequests; i++) {
    int nurseId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> cost;

    const Wish &wish = pPref->addShiftOff(
        nurseId, dayId, pScenario->shiftsFactory().pAnyWorkShift(), cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }


  // 2. read day on preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DAY_ON_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "DAY_ON_REQUESTS not found.");

  *pFile >> nRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nRequests; i++) {
    int nurseId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> cost;

    pPref->addShiftOn(nurseId, dayId,
                      pScenario->shiftsFactory().pAnyWorkShift(), cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // 3. read shift off preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SHIFT_OFF_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "SHIFT_OFF_REQUESTS not found.");

  *pFile >> nRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nRequests; i++) {
    int nurseId, shiftId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> strTmp;
    strTmp.pop_back();
    const PAbstractShift &pAS = pScenario->pShift(strTmp);
    *pFile >> cost;

    pPref->addShiftOff(nurseId, dayId, pAS, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // 4. read shift on preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SHIFT_ON_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "SHIFT_ON_REQUESTS not found.");

  *pFile >> nRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nRequests; i++) {
    int nurseId, shiftId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    Tools::readUntilChar(pFile, ',', &strTmp);
    *pFile >> strTmp;
    strTmp.pop_back();
    const PAbstractShift &pAS = pScenario->pShift(strTmp);
    *pFile >> cost;

    pPref->addShiftOn(nurseId, dayId, pAS, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }
}

//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
PScenario ReadWrite::readScenarioINRC2(const string &fileName) {
  std::fstream file;
  Tools::openFile(fileName, &file);
  string title;
  string strTmp;
  int intTmp;
  // declare the attributes that will initialize the Scenario instance
  //
  string name;
  int nDays = -1, nWeeks = -1, nSkills = -1, nShifts = -1,
      nShiftsType = -1, nContracts = -1, nNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType;
  map<string, int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType, shiftIDToShiftTypeID;
  vector<double> shiftDurations;
  vector2D<int> shiftTypeIDToShiftID,
      forbiddenShiftTypeSuccessors,
      forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  vector<PContract> pContracts;
  std::map<string, PContract> pContractsByName;
  vector<PNurse> pNurses;

  // default weights
  PWeights pWeights = std::make_shared<Weights>();

  // Fill the attributes of the scenario structure
  //

  // 1. read the name of the scenario
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "SCENARIO "))
    Tools::throwError("The INRC2 file is not as expected: "
                      "%s found instead of SCENARIO.", title.c_str());

  file >> name;

  // 2. read the number of weeks in scenario
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "WEEKS "))
    Tools::throwError("The INRC2 file is not as expected: "
                      "%s found instead of WEEKS.", title.c_str());

  file >> nWeeks;
  nDays = 7 * nWeeks;

  // 3. read the skills of the scenario
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "SKILLS "))
    Tools::throwError("The INRC2 file is not as expected: "
                      "%s found instead of SKILLS.", title.c_str());

  file >> nSkills;
  for (int i = 0; i < nSkills; i++) {
    file >> strTmp;
    intToSkill.push_back(strTmp);
    skillToInt.insert(pair<string, int>(strTmp, i));
  }

  // 4. read the different shift types and forbidden successions. In INRC2,
  // there is a bijection between shift types and shifts
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "SHIFT_TYPES "))
    Tools::throwError("The INRC2 file is not as expected:"
                      "%s found instead of SHIFT_TYPES.", title.c_str());

  //
  // Number of shifts : Given number + REST_SHIFT
  file >> intTmp;
  nShiftsType = intTmp + 1;

  // IMPORTANT : INSERT REST SHIFT !!!!!!
  // It is given 0 and 99 as bounds so that they never perturbate the cost
  intToShiftType.push_back(REST_SHIFT);
  intToShift.push_back(REST_SHIFT);
  shiftTypeToInt.insert(pair<string, int>(REST_SHIFT, 0));
  shiftToInt.insert(pair<string, int>(REST_SHIFT, 0));
  shiftTypeIDToShiftID.resize(nShiftsType);
  shiftIDToShiftTypeID.push_back(0);
  shiftTypeIDToShiftID[0].push_back(0);
  shiftDurations.push_back(0.0);
  minConsShiftType.push_back(0);
  maxConsShiftType.push_back(99);

  // Other shift types
  //
  for (int i = 1; i < nShiftsType; i++) {
    // Name
    file >> strTmp;
    intToShiftType.push_back(strTmp);
    shiftTypeToInt.insert(pair<string, int>(strTmp, i));
    Tools::readUntilChar(&file, '(', &strTmp);
    // Min consecutive
    file >> intTmp;
    minConsShiftType.push_back(intTmp);
    Tools::readUntilChar(&file, ',', &strTmp);
    // Max consecutive
    file >> intTmp;
    maxConsShiftType.push_back(intTmp);
    Tools::readLine(&file, &strTmp);
  }

  // 4.b. Forbidden successions for the shift type
  while (!file.eof() &&
      !Tools::strEndsWith(title, "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
    file >> title;

  if (file.eof())
    Tools::throwError("The INRC2 file is not as expected: "
                      "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS not found");

  forbiddenShiftTypeSuccessors.resize(nShiftsType);
  // Reading all lines
  for (int i = 1; i < nShiftsType; i++) {
    // Which current shift type ?
    string currentShiftType;
    file >> currentShiftType;
    int currentShiftTypeId = shiftTypeToInt.at(currentShiftType);
    // How many forbidden after it ?
    file >> intTmp;
    // Which ones are forbidden ?
    for (int j = 0; j < intTmp; j++) {
      file >> strTmp;
      forbiddenShiftTypeSuccessors[currentShiftTypeId].push_back(
          shiftTypeToInt.at(strTmp));
    }
    // make sure the forbidden successors are sorted in increasing order
    std::sort(forbiddenShiftTypeSuccessors[currentShiftTypeId].begin(),
              forbiddenShiftTypeSuccessors[currentShiftTypeId].end());
    // read end of line
    Tools::readLine(&file, &strTmp);
  }

  // 4.c. Check if SHIFTS is defined.
  // It is not a field present for the INRC2 format,
  // but it can be handled by our solver.
  // In case if it's not present (i.e., the default behavior for the INRC2),
  // we create one unique SHIFT for each SHITF_TYPE,
  // otherwise the shifts are read from the SHIFTS block in the file
  // Look at new1 for an example
  Tools::readUntilChar(&file, '=', &title);
  if (Tools::strEndsWith(title, "SHIFTS ")) {
    // Number of shifts
    file >> intTmp;
    nShifts = intTmp + 1;  // +1 for the REST_SHIFT

    // Shifts
    for (int i = 1; i < nShifts; i++) {
      // Name
      file >> strTmp;
      intToShift.push_back(strTmp);
      shiftToInt.insert(pair<string, int>(strTmp, i));
      // Duration (same unit than the min and max duration/days for the horizon)
      double duration;
      file >> duration;
      shiftDurations.push_back(duration);
      // shift type of the shift
      string currentShiftType;
      file >> currentShiftType;

      int currentShiftTypeId = shiftTypeToInt.at(currentShiftType);
      shiftIDToShiftTypeID.push_back(currentShiftTypeId);
      shiftTypeIDToShiftID[currentShiftTypeId].push_back(i);

      // read next line
      Tools::readLine(&file, &strTmp);
    }
    // read for the next field
    Tools::readUntilChar(&file, '=', &title);
  } else {
    // Create one shift for each shift type
    nShifts = nShiftsType;
    for (int i = 1; i < nShifts; i++) {
      std::string shiftType = intToShiftType.at(i);
      intToShift.push_back(shiftType);
      shiftToInt.insert(pair<string, int>(shiftType, i));
      shiftIDToShiftTypeID.push_back(i);
      shiftTypeIDToShiftID[i].push_back(i);
      shiftDurations.push_back(1.0);
    }
  }

  // 4.d. Build forbiddenShiftSuccessors by deducing it
  // from FORBIDDEN_SHIFT_TYPES_SUCCESSIONS
  forbiddenShiftSuccessors.resize(nShifts);
  for (int i = 0; i < nShiftsType; i++)
    for (int j : forbiddenShiftTypeSuccessors.at(i))
      for (int s1 : shiftTypeIDToShiftID[i])
        for (int s2 : shiftTypeIDToShiftID[j])
          forbiddenShiftSuccessors[s1].push_back(s2);

  // Then, Check if SHIFTS is FORBIDDEN_SHIFT_SUCCESSIONS.
  // It is not a field present for the INRC2 format,
  // but it can be handled by our solver.
  // It is organized as FORBIDDEN_SHIFT_TYPES_SUCCESSIONS except there is the
  // number of rows that are defined: FORBIDDEN_SHIFT_SUCCESSIONS = X
  // Look at new1 for an example
  if (Tools::strEndsWith(title, "FORBIDDEN_SHIFT_SUCCESSIONS ")) {
    // Reading all lines
    file >> intTmp;
    int nRows = intTmp;
    for (int i = 0; i < nRows; i++) {
      // Which current shift type ?
      string currentShift;
      file >> currentShift;
      int currentShiftId = shiftToInt.at(currentShift);
      // How many forbidden after it ?
      file >> intTmp;
      // Which ones are forbidden ?
      auto &f = forbiddenShiftSuccessors[currentShiftId];
      for (int j = 0; j < intTmp; j++) {
        file >> strTmp;
        int k = shiftToInt.at(strTmp);
        if (find(f.begin(), f.end(), k) == f.end())
          f.push_back(k);
      }
      // read end of line
      Tools::readLine(&file, &strTmp);
    }
    Tools::readUntilChar(&file, '=', &title);
  }

  // 4.e. Create shift structures
  for (int i = 0; i < nShifts; i++) {
    // make sure the forbidden successors are sorted in increasing order
    vector<int> &f = forbiddenShiftSuccessors[i];
    std::sort(f.begin(), f.end());
    // number of forbidden successors
    std::vector<int> successorList;
    for (int s = 0; s < nShifts; ++s) {
      if (find(f.begin(), f.end(), s) == f.end())
        successorList.push_back(s);
    }
    int type = shiftIDToShiftTypeID[i];
    pShifts.push_back(std::make_shared<Shift>(intToShift[i],
                                              i,
                                              type,
                                              shiftDurations[i],
                                              successorList));
  }

  // 5. read the contract types
  if (!Tools::strEndsWith(title, "CONTRACTS "))
    Tools::throwError("The INRC2 file is not as expected %s found "
                      "instead of CONTRACTS.", title.c_str());
  file >> intTmp;
  nContracts = intTmp;
  // read each contract type
  for (int i = 0; i < nContracts; i++) {
    string contractName;
    int minDays, maxDays, minConsWork, maxConsWork, minConsRest,
        maxConsRest, maxWeekends, isTotalWeekend;
    file >> contractName;
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> minDays;
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> maxDays;
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> minConsWork;
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> maxConsWork;
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> minConsRest;
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> maxConsRest;
    Tools::readUntilChar(&file, ' ', &strTmp);
    file >> maxWeekends;
    Tools::readUntilChar(&file, ' ', &strTmp);
    file >> isTotalWeekend;
    Tools::readLine(&file, &strTmp);

    pContracts.push_back(
        std::make_shared<Contract>(i,
                                   contractName,
                                   minDays,
                                   maxDays,
                                   minConsWork,
                                   maxConsWork,
                                   minConsRest,
                                   maxConsRest,
                                   maxWeekends,
                                   isTotalWeekend));
    pContractsByName[contractName] = pContracts.back();
  }

  // 6. read the list of nurses
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "NURSES "))
    Tools::throwError("The INRC2 file is not as expected: "
                      "%s found instead of NURSES.", title.c_str());

  file >> nNurses;
  for (int i = 0; i < nNurses; i++) {
    string nurseName, contractName;
    int n;
    vector<int> skills, availableShifts;
    // read everything on the line
    file >> nurseName;
    file >> contractName;
    file >> n;
    // read the skills of the nurse
    for (int j = 0; j < n; j++) {
      file >> strTmp;
      // if a skill
      auto it = skillToInt.find(strTmp);
      if (it != skillToInt.end()) {
        skills.push_back(it->second);
        continue;
      }
      // TODO(AL): the code below is not for INRC2 files either
      /*// if a shift
      it = shiftToInt.find(strTmp);
      if (it != shiftToInt.end()) {
        availableShifts.push_back(it->second);
        continue;
      }
      // if a shift type
      it = shiftTypeToInt.find(strTmp);
      if (it != shiftTypeToInt.end()) {
        for (int s : shiftTypeIDToShiftID[it->second])
          availableShifts.push_back(s);
        continue;
      }*/
    }
    // check if skills empty, otherwise sort the skill indices
    // before initializing the nurse
    if (skills.empty()) {
      if (nSkills > 1)
        Tools::throwError("Several skills have been defined, "
                          "but none of them were attributed to nurses.");
      skills = {0};  // add the only skill
    } else {
      std::sort(skills.begin(), skills.end());
    }

    // TODO(AL): should revise below to be specific to INRC2
    // if no shifts in the list of available shifts, put all of them, otherwise
    // sort them
    if (availableShifts.empty()) {
      for (int s = 0; s < nShifts; s++)
        availableShifts.push_back(s);
    } else {
      availableShifts.push_back(0);  // add the rest shift
      std::sort(availableShifts.begin(), availableShifts.end());
    }

    pNurses.push_back(std::make_shared<Nurse>(
        i, nurseName, nShifts, nSkills, skills, availableShifts,
        pContractsByName.at(contractName)));
    nurseNameToInt.insert(pair<string, int>(nurseName, i));
  }

  // Check that all fields were initialized before initializing the scenario
  //
  if (nWeeks == -1 || nSkills == -1 || nShifts == -1 || nContracts == -1
      || nNurses == -1) {
    Tools::throwError("In readScenarioINRC2: missing fields in the "
                      "initialization");
  }

  return std::make_shared<Scenario>(name,
                                    nWeeks,
                                    nSkills,
                                    intToSkill,
                                    skillToInt,
                                    pShifts,
                                    intToShiftType,
                                    minConsShiftType,
                                    maxConsShiftType,
                                    nContracts,
                                    pContracts,
                                    nNurses,
                                    pNurses,
                                    nurseNameToInt,
                                    pWeights,
                                    "",
                                    false,
                                    true);
}

PDemand ReadWrite::readINRC2Weeks(const std::vector<std::string> &strWeekFiles,
                                  const PScenario &pScenario) {
  // initialize pDemand
  PDemand pDemand;
  PPreferences pPref;

  for (const string &strWeekFile : strWeekFiles)
    if (!pDemand) {
      ReadWrite::readWeekINRC2(strWeekFile, pScenario, &pDemand, &pPref);
    } else {
      // load the next week
      PDemand nextDemand;
      PPreferences nextPref;
      ReadWrite::readWeekINRC2(strWeekFile, pScenario, &nextDemand, &nextPref);
      // update the current weeks
      pDemand->pushBack(nextDemand);
      pPref->pushBack(nextPref);
    }

  // link the scenario to the current demand and preferences
  pScenario->linkWithDemand(pDemand);
  pScenario->linkWithPreferences(pPref);

  return pDemand;
}

// Read the Week file and store the content in a Scenario instance
//
void ReadWrite::readWeekINRC2(const std::string &strWeekFile,
                              const PScenario &pScenario,
                              PDemand *pDemand, PPreferences *pPref) {
  // open the file
  std::fstream file;
  Tools::openFile(strWeekFile, &file);

  string title;
  string strTmp;
  int intTmp;

  // declare the attributes to be updated in the PScenario
  //
  string weekName;
  vector3D<int> minWeekDemand;
  vector3D<int> optWeekDemand;

  // fill the attributes when reading the week file
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    // Read the name of the week
    //
    if (Tools::strEndsWith(title, "WEEK_DATA")) {
      file >> weekName;
    } else if (Tools::strEndsWith(title, "REQUIREMENTS")) {
      // Read the requirements
      //
      string shiftName, skillName;
      int shiftId, skillId;
      // init the vectors
      Tools::initVector3D(&minWeekDemand,
                          7,
                          pScenario->nShifts(),
                          pScenario->nSkills(),
                          0);
      Tools::initVector3D(&optWeekDemand,
                          7,
                          pScenario->nShifts(),
                          pScenario->nSkills(),
                          0);

      // Do not take the rest shift into account here
      // (by initialization, requirements already at 0
      char c = '0';
      while (std::isalnum(c)) {  // check if character is alphanumeric
        // Read shift and skill
        file >> shiftName;
        file >> skillName;
        shiftId = pScenario->shift(shiftName);
        skillId = pScenario->skillId(skillName);
        // For every day in the week, read min and opt values
        for (int day = 0; day < 7; day++) {
          Tools::readUntilChar(&file, '(', &strTmp);
          file >> intTmp;
          minWeekDemand[day][shiftId][skillId] = intTmp;
          Tools::readUntilChar(&file, ',', &strTmp);
          file >> intTmp;
          optWeekDemand[day][shiftId][skillId] = intTmp;
        }
        Tools::readLine(&file, &strTmp);
        c = file.peek();
      }
    } else if (Tools::strEndsWith(title, "SHIFT_OFF_REQUESTS ")) {
      // Read the shift off requests
      if (!*pPref)
        *pPref = std::make_shared<Preferences>(pScenario->nNurses(),
                                               7,
                                               pScenario->nShifts());

      // Temporary vars
      const std::vector<double> &prefCosts = pScenario->weights().preferences;
      string nurseName, shift, day, strLevel;
      int nShifts, nurseNum;
      DayOfWeek dayOfWeek;
      PREF_LEVEL level = WEAK;
      file >> nShifts;
      for (int i = 0; i < nShifts; i++) {
        if (nurseName.empty())
          file >> nurseName;
        nurseNum = pScenario->nurse(nurseName);
        nurseName.clear();
        file >> shift;
        file >> day;
        dayOfWeek = shortNameToDayOfWeek(day);
        // in case there is no level defined for the preferences
        file >> strLevel;
        try {
          level = (PREF_LEVEL) std::stoi(strLevel);
        } catch (const std::invalid_argument &ia) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addShiftOff(nurseNum, dayOfWeek,
                                pScenario->shiftsFactory().pAnyWorkShift(),
                                prefCosts.at(level));
        } else {
          (*pPref)->addShiftOff(nurseNum, dayOfWeek, pScenario->pShift(shift),
                                prefCosts.at(level));
        }
      }
    } else if (Tools::strEndsWith(title, "SHIFT_ON_REQUESTS ")) {
      // Read the shift on requests
      if (!*pPref)
        *pPref = std::make_shared<Preferences>(pScenario->pNurses(),
                                               7,
                                               pScenario->nShifts());

      // Temporary vars
      const std::vector<double> &prefCosts = pScenario->weights().preferences;
      string nurseName, shift, day, strLevel;
      int nShifts, nurseNum;
      DayOfWeek dayOfWeek;
      PREF_LEVEL level = WEAK;
      file >> nShifts;
      for (int i = 0; i < nShifts; i++) {
        if (nurseName.empty())
          file >> nurseName;
        nurseNum = pScenario->nurse(nurseName);
        nurseName.clear();
        file >> shift;
        file >> day;
        dayOfWeek = shortNameToDayOfWeek(day);
        // in case there is no level defined for the preferences
        file >> strLevel;
        try {
          level = (PREF_LEVEL) std::stoi(strLevel);
        } catch (const std::invalid_argument &ia) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addShiftOn(nurseNum, dayOfWeek,
                               pScenario->shiftsFactory().pAnyWorkShift(),
                               prefCosts.at(level));
        } else {
          (*pPref)->addShiftOn(nurseNum, dayOfWeek, pScenario->pShift(shift),
                               prefCosts.at(level));
        }
      }
    }
  }

  // Define a new instance of demand
  *pDemand = std::make_shared<Demand>(7,
                                      0,
                                      pScenario->nShifts(),
                                      pScenario->nSkills(),
                                      weekName,
                                      minWeekDemand,
                                      optWeekDemand);
#ifdef NS_DEBUG
  std::cout << "Demand created" << std::endl;
#endif
}

// Read the history file
//
void ReadWrite::readHistoryINRC2(const std::string &strHistoryFile,
                                 const PScenario &pScenario) {
  if (strHistoryFile.empty()) {
    std::cout << "Cyclic option is enable, so no history file is loaded."
              << std::endl;
    pScenario->enableCyclic();
    return;
  }
  // open the file
  std::fstream file;
  Tools::openFile(strHistoryFile, &file);

  string title;
  string strTmp;

  // declare the attributes to be updated in the PScenario
  //
  int thisWeek;
  string weekName;
  vector<State> initialState;


  // fill the attributes of the week structure
  //
  while (file.good()) {
    Tools::readLine(&file, &title);

    // Read the index and name of the week
    //
    if (!strcmp(title.c_str(), "HISTORY")) {
      file >> thisWeek;
      file >> weekName;
      // Raise exception if it does not match the week previously read !
      if (strcmp(weekName.c_str(), (pScenario->weekName()).c_str()) != 0) {
        std::cout << "The given history file requires week " << weekName
                  << std::endl;
        std::cout << " but a different one (" << pScenario->weekName()
                  << ") has been given!" << std::endl;
        Tools::throwError("History file and week data file do not match!");
      }
    } else if (Tools::strEndsWith(title, "NURSE_HISTORY")) {
      // Read each nurse's initial state
      //
      for (int n = 0; n < pScenario->nNurses(); n++) {
        string nurseName, shiftName;
        // int nurseNum;
        int shiftId, totalTimeWorked, totalWeekendsWorked, consDaysWorked,
            consShiftWorked, consRest, consShifts;
        file >> nurseName;
        // nurseNum = pScenario->nurseNameToInt_.at(nurseName);
        file >> totalTimeWorked;
        file >> totalWeekendsWorked;
        file >> shiftName;
        if (shiftName == "None") shiftName = REST_SHIFT;
        shiftId = pScenario->shift(shiftName);
        file >> consShiftWorked;
        file >> consDaysWorked;
        file >> consRest;

        if (consRest == 0 && consDaysWorked == 0)
          Tools::throwError("History of nurse %s is invalid as one must "
                            "either work or rest.", nurseName.c_str());

        consShifts = (shiftId == 0) ? consRest : consShiftWorked;
        State nurseState(0,
                         totalTimeWorked,
                         totalWeekendsWorked,
                         consDaysWorked,
                         consShifts,
                         consRest,
                         pScenario->pShift(shiftId));
        initialState.push_back(nurseState);
      }
    }
  }
  pScenario->setThisWeek(thisWeek);
  pScenario->setInitialState(initialState);
}

// Read the input custom file
// Store the result in a vector of historical demands and
// return the number of treated weeks
//
int ReadWrite::readCustom(const string &strCustomInputFile,
                          const PScenario &pScenario,
                          vector<PDemand> *demandHistory) {
  // open the file
  std::fstream file;
  Tools::openFile(strCustomInputFile, &file);

  string title;
  int nWeeks;

  // get the custom information
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    // Read the file names of the past demand
    //
    if (!strcmp(title.c_str(), "PAST_DEMAND_FILES")) {
      file >> nWeeks;
      if (!nWeeks) continue;

      string strDemandFile;
      for (int i = 0; i < nWeeks; i++) {
        file >> strDemandFile;
        PDemand pDemand;
        PPreferences pPref;
        readWeekINRC2(strDemandFile, pScenario, &pDemand, &pPref);
        demandHistory->push_back(pDemand);
      }
    }
  }
  return nWeeks;
}

void ReadWrite::writeCustom(string strCustomOutputFile,
                            const string &strWeekFile,
                            const string &strCustomInputFile) {
  Tools::LogOutput outStream(std::move(strCustomOutputFile));

  // if there is no custom input file, this is the first week
  if (strCustomInputFile.empty()) {
    outStream << "PAST_DEMAND_FILES= " << 1 << std::endl;
    outStream << strWeekFile << std::endl;
    return;
  }

  // open the custom input file
  // we want the content of the input custom file in the custom output file
  std::fstream file;
  Tools::openFile(strCustomInputFile, &file);

  string title;
  int nWeeks;

  // fill the attributes of the week structure
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    // Read the file names of the past demand
    //
    if (!strcmp(title.c_str(), "PAST_DEMAND_FILES")) {
      file >> nWeeks;
      outStream << "PAST_DEMAND_FILES= " << nWeeks + 1 << std::endl;
      if (!nWeeks) continue;

      string strDemandFile;
      for (int i = 0; i < nWeeks; i++) {
        file >> strDemandFile;
        outStream << strDemandFile << std::endl;
      }
      outStream << strWeekFile << std::endl;
    }
  }
}

void writeUI(const string &path, const PScenario &pScenario) {
  string filename = path + "_ui.txt";
  Tools::LogOutput outStream(filename);

  // open the custom input file
  // we want the content of the input custom file in the custom output file
  std::fstream file;
  Tools::openFile(filename, &file);
  file << "SCHEDULING_PERIOD\n";
  file << pScenario->name() << ",1,";
  file << std::put_time(&pScenario->startDate(), "%Y-%m-%d") << ",";
  auto end_date = pScenario->startDate();
  end_date.tm_yday += pScenario->nDays();
  mktime(&end_date);
  file << std::put_time(&end_date, "%Y-%m-%d") << "\nEND\nSKILLS\n";
  for (int i; i < pScenario->nSkills(); i++)
    file << pScenario->skillName(i) << "\n";
  file << "END\nSHIFTS\n";
  for (int i; i < pScenario->nShifts(); i++)
    file << pScenario->shiftName(i) << ",00:00,00:00" << std::endl;
  file << "END\nSHIFTS_TYPES\n";
  for (int i = 0; i < pScenario->nShiftTypes(); i++) {
    file << pScenario->shiftType(i) << ", "
         << (pScenario->pShiftsOfType(i)).size();
    for (const auto &s : pScenario->pShiftsOfType(i))
      file << "," << s->name;
    file << std::endl;
  }
  file << "\nEND\nCONTRACTS\n{";
}

PScenario ReadWrite::readScenarioUI(
    const std::string &fileName, const string &name) {
  // declare the attributes that will initialize the Scenario instance
  //
  string header;
  bool contracts_processed = false;
  int nDays = -1, nWeeks = -1, nSkills = -1, nShifts = -1,
      nShiftsType = -1, nContracts = -1, nNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType;
  map<string, int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType, shiftIDToShiftTypeID;
  vector<double> shiftDurations;
  vector2D<int> shiftTypeIDToShiftID,
      forbiddenShiftTypeSuccessors,
      forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  // "input" -> 0 if input is a shift, 1 if input is a shift type,
  // 2 if input if a shift group
  map<string, int> inputToShiftGenre;
  vector<PContract> pContracts;
  std::map<string, PContract> pContractsByName;
  vector<PNurse> pNurses;
  map<std::string, vector<PBaseResource>> ruleSets;
  // default weights
  vector<double> pref_scale{10, 50, 200, 1000};
  PWeights pWeights = std::make_shared<Weights>(
      10, 10, 30, 30, 10, pref_scale, 1000, 500, 500, -1);
  PDemand pDemand;
  PPreferences pPref;
  vector<State> initialState;
  PScenario pScenario;
  std::tm startDay;
  struct HistoryPeriod hist;
  // open the file
  std::fstream file;
  std::cout << "Reading " << fileName << std::endl;
  file.open(fileName.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read the file " << fileName << std::endl;
    std::cout << "The input file was not opened properly!" << std::endl;

    throw Tools::myException("The input file was not opened properly!",
                             __LINE__);
  }
  string strTmp;
  std::string l;
  std::string buffer;
  // remove comments and empty lines
  while (std::getline(file, l) && file.good()) {
    if (l[0] == '#' || l[0] == '\n') {
      continue;
    } else {
      buffer += l + "\n";
    }
  }
  // Parse the buffer in meaningful blocks
  std::regex reg("END\n");
  std::sregex_token_iterator iter(buffer.begin(), buffer.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);

  for (const string &s : tokens) {
    std::stringstream X(s);
    std::string line;
    std::getline(X, line);
    boost::trim(line);
    if (line == "HEADERS") {
      header = s;
    } else if (line == "SCHEDULING_PERIOD") {
      const struct Sched_Period t = parse_scheduling_period(s);
      nWeeks = t.nWeeks;
      nDays = t.nDays;
      startDay = t.startDay;
      header.append(s);
      Day::setFirstDayOfWeek(
              (DayOfWeek) (startDay.tm_wday >= 0 ? startDay.tm_wday - 1 : 6));
    } else if (line == "SKILLS") {
      const struct Skills_Parsed t = parse_skills(s);
      nSkills = t.nbSkills;
      skillToInt = t.skillToInt;
      intToSkill = t.intToSkill;
    } else if (line == "SHIFTS") {
      const struct Shifts_Parsed t = parse_shifts(s);
      nShifts = t.nbShifts;
      shiftDurations = t.hoursInShift;
      intToShift = t.intToShift;
      shiftToInt = t.shiftToInt;
    } else if (line == "SHIFT_TYPES") {
      const struct ShiftTypes_Parsed
          t = parse_shiftsType(s, shiftToInt, intToShift, shiftDurations);
      inputToShiftGenre = t.inputToShiftGenre;
      pShifts = t.pShifts;
      nShiftsType = t.nbShiftsType;
      intToShiftType = t.intToShiftType;
      shiftTypeToInt = t.shiftTypeToInt;
      shiftTypeIDToShiftID = t.shiftTypeIDToShiftID;
      shiftIDToShiftTypeID = t.shiftIDToShiftTypeID;
    } else if (line == "CONTRACTS") {
      ruleSets = parse_contracts(s,
                                 ShiftsFactory(pShifts),
                                 inputToShiftGenre,
                                 shiftToInt,
                                 skillToInt,
                                 shiftTypeToInt,
                                 pWeights,
                                 nDays);
      contracts_processed = ruleSets.size();
    } else if (line == "CONTRACT_GROUPS" && contracts_processed) {
      ruleSets = parse_group_contracts(s, ruleSets);
    } else if (line == "EMPLOYEES") {
      // As I go thought the employee,
      // I craft them contracts based on the rule Sets as needed
      vector<int> allSkills, allShifts;
      for (int i = 0; i < nSkills; i++) {
        allSkills.push_back(i);
      }
      for (int i = 0; i < nShifts; i++) {
        allShifts.push_back(i);
      }
      struct Nurses_Parsed n = parse_nurses(s,
                                            pContractsByName,
                                            pContracts,
                                            ruleSets,
                                            pNurses,
                                            allSkills,
                                            allShifts);
      pContractsByName = n.pContractsByName;
      pContracts = n.pContracts;
      pNurses = n.pNurses;
      nNurses = pNurses.size();
      for (int n = 0; n < nNurses; n++) {
        nurseNameToInt[pNurses[n]->name_] = n;
      }
      nContracts = pContracts.size();

      pScenario = std::make_shared<Scenario>(name,
                                             nWeeks,
                                             nSkills,
                                             intToSkill,
                                             skillToInt,
                                             pShifts,
                                             intToShift,
                                             vector<int>(2, nDays),
                                             vector<int>(2, nDays),
                                             nContracts,
                                             pContracts,
                                             nNurses,
                                             pNurses,
                                             nurseNameToInt,
                                             pWeights,
                                             header,
                                             false,
                                             false);

    } else if (line == "HOSPITAL_DEMAND") {
      pScenario->setStartDate(startDay);
      pDemand = parse_demand(s, shiftToInt, skillToInt, nDays, startDay);
      std::cout << pDemand->toString(true) << std::endl;
      pScenario->linkWithDemand(pDemand);
    } else if (line == "PREFERENCES") {
      // parse s - line
      parse_preferences(s, pScenario);
    } else if (line == "HISTORY") {
      initialState = parse_history(s, pScenario, startDay);
    }
  }
  // Initialize the history of every nurse to empty state
  // Add a fictitious shift just for the initial state if no history given
  if (initialState.empty()) {
    const PShift &pNoneShift =
        pScenario->shiftsFactory().pNoneShift()->pIncludedShifts().front();
    for (int n = 0; n < nNurses; n++) {
      State nurseState(-1, 0, 0, 0, 0, 0, pNoneShift);
      initialState.push_back(nurseState);
    }
  }
  pScenario->setInitialState(initialState);
  return pScenario;
}

/************************************************************************
* Print the main characteristics of all the demands of an input directory
* This is done to find some invariant properties among demands
*************************************************************************/
void ReadWrite::compareDemands(const string &inputDir, string logFile) {
  struct dirent *dirp;
  Tools::LogOutput logStream(std::move(logFile), 8);

  vector2D<int> minPerShift, optPerShift, minPerSkill, optPerSkill;
  vector2D<int> minHighestPerSkill, optHighestPerSkill;
  vector<int> minTotal, optTotal;

  // Open the input directory
  DIR *dp = opendir(inputDir.c_str());
  if (dp == nullptr) {
    Tools::throwError("Error while opening ");
  } else {
    std::cout << "Reading from directory " << inputDir << std::endl;
  }

  // Read the scenario that appears in the directory
  unsigned found = inputDir.find_last_of('/');
  string instanceName = inputDir.substr(found + 1);
  string scenFile = inputDir + "/Sc-" + instanceName + ".txt";
  string historyFile = inputDir + "/H0-" + instanceName + "-0.txt";

  PScenario pScen = ReadWrite::readScenarioINRC2(scenFile);

  // Go through all the demand files of the directory
  int coDemand = 0;
  while ((dirp = readdir(dp))) {
    std::string filename(dirp->d_name);

    // The file names of week demands start with "WD"
    found = filename.find("WD");
    if (found > 0) continue;

    string filepath = inputDir;
    filepath += "/";
    filepath += filename;

    PDemand pDemand;
    PPreferences pPref;
    ReadWrite::readWeekINRC2(filepath, pScen, &pDemand, &pPref);

    logStream << "#####################################\n";
    logStream << "# DEMAND FILE: " << filepath << std::endl;
    logStream << "#####################################\n\n";
    logStream << pDemand->toString(true) << std::endl << std::endl;

    // record the advanced data on the demand
    minTotal.push_back(pDemand->minTotal_);
    optTotal.push_back(pDemand->optTotal_);
    minPerShift.push_back(pDemand->minPerShift_);
    optPerShift.push_back(pDemand->optPerShift_);
    minPerSkill.push_back(pDemand->minPerSkill_);
    optPerSkill.push_back(pDemand->optPerSkill_);
    minHighestPerSkill.push_back(pDemand->minHighestPerSkill_);
    optHighestPerSkill.push_back(pDemand->optHighestPerSkill_);

    // link the scenario with the first demand and preferences to be able to
    // retrieve information about the nurses
    if (!coDemand) {
      pScen->linkWithDemand(pDemand);
      pScen->linkWithPreferences(pPref);
    } else {
      coDemand++;
    }
  }

  // Also presolve the nurses to get statistics on the capacity of the nurses
  // to cover the demand
  //
  ReadWrite::readHistoryINRC2(historyFile, pScen);
  auto *pSolver = new Solver(pScen);
  pSolver->preprocessTheNurses();

  // Write a summary  of the advanced data computed for the demands
  logStream << "#####################################" << std::endl;
  logStream << "# SUMMARY OF THE STATISTICS" << std::endl;
  logStream << "#####################################";
  logStream << std::endl << std::endl;

  logStream
      << "# Total capacity of the nurses "
         "(without unavoidable penalty/with average work): ";
  logStream << Tools::itoa(pSolver->maxTotalStaffNoPenalty_) + "/"
      + Tools::itoa(pSolver->maxTotalStaffAvgWork_);
  logStream.endl();
  logStream.endl();

  logStream
      << "# Total capacity of the nurses per skill "
         "(without unavoidable penalty/with average work):" << std::endl;
  for (int i = 0; i < pScen->nSkills(); i++) {
    logStream << "# SK" + Tools::itoa(i) + ":";
    logStream << Tools::itoa(pSolver->maxStaffPerSkillNoPenalty_[i]) + "/" +
        Tools::itoa(pSolver->maxStaffPerSkillAvgWork_[i]);
    logStream.endl();
  }
  logStream.endl();

  logStream << "# Total demand (min/opt): \n";
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream.endl();
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream << Tools::itoa(minTotal[d]) + "/" + Tools::itoa(optTotal[d]);
  }
  logStream.endl();
  logStream.endl();

  logStream << "# Demand per shift\n";
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream.endl();
  for (int i = 1; i < pScen->nShifts(); i++) {
    logStream << "# SH" + Tools::itoa(i) + ":";
    for (unsigned int d = 0; d < minPerShift.size(); d++) {
      logStream << Tools::itoa(minPerShift[d][i]) + "/"
          + Tools::itoa(optPerShift[d][i]);
    }
    logStream.endl();
  }
  logStream.endl();

  logStream << "# Demand per skill\n";
  logStream << "#";
  for (unsigned int d = 0; d < minPerSkill.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream.endl();
  for (int i = 0; i < pScen->nSkills(); i++) {
    logStream.setWidth(12);
    logStream << "# " + pScen->skillName(i) + ":";
    logStream.setWidth(10);
    for (unsigned int d = 0; d < minPerShift.size(); d++) {
      logStream << Tools::itoa(minPerSkill[d][i]) + "/"
          + Tools::itoa(optPerSkill[d][i]);
    }
    logStream.endl();
  }
  logStream.endl();

  logStream << "# Highest demand per skill for one shift" << std::endl;
  logStream << "#";
  for (unsigned int d = 0; d < minPerSkill.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream.endl();
  for (int i = 0; i < pScen->nSkills(); i++) {
    logStream << "# SK" + Tools::itoa(i) + ":";
    for (unsigned int d = 0; d < minPerShift.size(); d++) {
      logStream << Tools::itoa(minHighestPerSkill[d][i]) + "/"
          + Tools::itoa(optHighestPerSkill[d][i]);
    }
    logStream.endl();
  }
  logStream.endl();

  logStream.setWidth(12);
  logStream << "# Agregate indicators" << std::endl;
  logStream << "# Total understaffing without unavoidable penalty" << std::endl;
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream << "Average" << "Std dev" << std::endl;
  logStream << "#";
  double averageMin = 0, averageOpt = 0, stdDevMin = 0, stdDevOpt = 0;
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream
        << Tools::itoa(minTotal[d] - pSolver->maxTotalStaffNoPenalty_) + "/" +
            Tools::itoa(optTotal[d] - pSolver->maxTotalStaffNoPenalty_);
    averageMin += static_cast<double>(minTotal[d] -
        pSolver->maxTotalStaffNoPenalty_)
        / static_cast<double>(minPerShift.size());
    averageOpt += static_cast<double>(optTotal[d] -
        pSolver->maxTotalStaffNoPenalty_)
        / static_cast<double>(minPerShift.size());
    stdDevMin += pow(minTotal[d] - pSolver->maxTotalStaffNoPenalty_, 2);
    stdDevOpt += pow(optTotal[d] - pSolver->maxTotalStaffNoPenalty_, 2);
  }
  stdDevMin =
      sqrt(stdDevMin / minPerShift.size() - pow(averageMin, 2));
  stdDevOpt =
      sqrt(stdDevOpt / minPerShift.size() - pow(averageOpt, 2));
  logStream << Tools::itoa(averageMin) + "/" + Tools::itoa(averageOpt);
  logStream << Tools::itoa(stdDevMin) + "/" + Tools::itoa(stdDevOpt);
  logStream.endl();
  logStream.endl();

  logStream << "# Total understaffing with average work" << std::endl;
  logStream << "#";
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream << "Average" << "Std dev" << std::endl;
  logStream << "#";
  averageMin = 0, averageOpt = 0, stdDevMin = 0, stdDevOpt = 0;
  for (unsigned int d = 0; d < minPerShift.size(); d++) {
    logStream
        << Tools::itoa(minTotal[d] - pSolver->maxTotalStaffAvgWork_) + "/" +
            Tools::itoa(optTotal[d] - pSolver->maxTotalStaffAvgWork_);
    averageMin += static_cast<double>(minTotal[d] -
        pSolver->maxTotalStaffAvgWork_)
        / static_cast<double>(minPerShift.size());
    averageOpt += static_cast<double>(optTotal[d] -
        pSolver->maxTotalStaffAvgWork_)
        / static_cast<double>(minPerShift.size());
    stdDevMin += pow(minTotal[d] - pSolver->maxTotalStaffAvgWork_, 2);
    stdDevOpt += pow(optTotal[d] - pSolver->maxTotalStaffAvgWork_, 2);
  }
  stdDevMin =
      sqrt(stdDevMin / minPerShift.size() - pow(averageMin, 2));
  stdDevOpt =
      sqrt(stdDevOpt / minPerShift.size() - pow(averageOpt, 2));
  logStream << Tools::itoa(averageMin) + "/" + Tools::itoa(averageOpt);
  logStream << Tools::itoa(stdDevMin) + "/" + Tools::itoa(stdDevOpt);
  logStream.endl();
  logStream.endl();

  logStream.setPrecision(1);
  logStream << "# Understaffing per skill without unavoidable penalty"
            << std::endl;
  logStream << "#";
  for (unsigned int d = 0; d < minPerSkill.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream << "Average" << "Std dev" << std::endl;
  for (int i = 0; i < pScen->nSkills(); i++) {
    logStream << "# " + pScen->skillName(i) + ":";
    averageMin = 0, averageOpt = 0, stdDevMin = 0, stdDevOpt = 0;
    for (unsigned int d = 0; d < minPerShift.size(); d++) {
      logStream << Tools::itoa(
          minPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i]) + "/"
          + Tools::itoa(
              optPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i]);
      averageMin += (minPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i])
          / minPerShift.size();
      averageOpt += (optPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i])
          / minPerShift.size();
      stdDevMin +=
          pow(minPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i], 2);
      stdDevOpt +=
          pow(optPerSkill[d][i] - pSolver->maxStaffPerSkillNoPenalty_[i], 2);
    }
    stdDevMin =
        sqrt(stdDevMin / minPerShift.size() - pow(averageMin, 2));
    stdDevOpt =
        sqrt(stdDevOpt / minPerShift.size() - pow(averageOpt, 2));
    logStream << Tools::itoa(averageMin) + "/" + Tools::itoa(averageOpt);
    logStream << Tools::itoa(stdDevMin) + "/" + Tools::itoa(stdDevOpt);
    logStream.endl();
  }
  logStream.endl();

  logStream << "# Understaffing per skill with average work" << std::endl;
  logStream << "#";
  for (unsigned int d = 0; d < minPerSkill.size(); d++) {
    logStream << "WD" + Tools::itoa(d);
  }
  logStream << "Average" << "Std dev" << std::endl;
  for (int i = 0; i < pScen->nSkills(); i++) {
    logStream << "# " + pScen->skillName(i) + ":";
    averageMin = 0, averageOpt = 0, stdDevMin = 0, stdDevOpt = 0;
    for (unsigned int d = 0; d < minPerShift.size(); d++) {
      logStream << Tools::itoa(
          minPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i]) + "/"
          + Tools::itoa(
              optPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i]);
      averageMin += (minPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i])
          / minPerShift.size();
      averageOpt += (optPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i])
          / minPerShift.size();
      stdDevMin +=
          pow(minPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i], 2);
      stdDevOpt +=
          pow(optPerSkill[d][i] - pSolver->maxStaffPerSkillAvgWork_[i], 2);
    }
    stdDevMin =
        sqrt(stdDevMin / minPerShift.size() - pow(averageMin, 2));
    stdDevOpt =
        sqrt(stdDevOpt / minPerShift.size() - pow(averageOpt, 2));
    logStream << Tools::itoa(averageMin) + "/" + Tools::itoa(averageOpt);
    logStream << Tools::itoa(stdDevMin) + "/" + Tools::itoa(stdDevOpt);
    logStream.endl();
  }
  logStream.endl();
}
