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

#include "tools/Tools.h"
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
// Method that read an INRC input files and stores the content in the
// output scenario instance
//
PScenario ReadWrite::readINRCInstance(const string& fileName) {
  std::fstream file;
  Tools::openFile(fileName, &file);
  string strTmp;
  int intTmp;

  // declare the attributes that will initialize the Scenario instance
  string instName, name;
  int nbWeeks = -1, nbSkills = -1, nbShifts = -1, nbShiftTypes = -1,
      nbContracts = -1, nbNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType;
  map<string, int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType,
      shiftIDToShiftTypeID, nbForbidShiftSucc;
  vector<double> shiftDurations;
  vector2D<int> shiftTypeIDToShiftID,
      forbiddenShiftTypeSuccessors,
      forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  vector<PContract> pContracts;
  vector<PNurse> theNurses;

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
  int nbDays(0);
  nbDays = tmEnd.tm_yday - tmStart.tm_yday + 1;
  nbWeeks = std::ceil(nbDays / 7);

  Day::setFirstDayOfWeek(
      (DayOfWeek) (tmStart.tm_wday >= 0 ? tmStart.tm_wday - 1 : 6));
  std::cout << "Starting week day = " << dayOfWeekToName(Day::firstDayOfWeek())
            << " ; number of days = " << nbDays
            << " ; number of weeks = " << nbWeeks;

  if (nbDays < 0)
    Tools::throwError("The start of the horizon should be before its end");
  if (tmEnd.tm_year != tmStart.tm_year)
    Tools::throwError("The beginning and end of the scheduling horizon should"
                      " be on the same year");

  // set the list of skills
  Tools::readUntilChar(&file, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SKILLS "))
    Tools::throwError("The INRC file is not as expected");

  file >> nbSkills;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  for (int i = 0; i < nbSkills; i++) {
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

  file >> nbShifts;
  nbShifts += 1;  // add the rest shift

  // in INRC, any shift can come after any other shift (possibly at the
  // cost of a penalty)
  vector<int> successorList;
  successorList.reserve(nbShifts);
  for (int i = 0; i < nbShifts; i++)
    successorList.push_back(i);

  // IMPORTANT : INSERT REST SHIFT !!!!!!
  intToShiftType.push_back(REST_SHIFT);
  shiftToInt[REST_SHIFT] = 0;
  shiftDurations.push_back(0.0);
  // initialize the shift
  pShifts.push_back(std::make_shared<Shift>(REST_SHIFT,
                                            0,
                                            0,
                                            shiftDurations[0],
                                            successorList));
  nbShiftTypes = 1;
  shiftTypeIDToShiftID.push_back(vector<int>());
  shiftIDToShiftTypeID.push_back(0);
  shiftTypeIDToShiftID[0].push_back(0);

  // read the name, type, time interval and skills of the shifts
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  string shName, typeName, skName;
  vector<int> skillIds;
  for (int i = 1; i < nbShifts; i++) {
    // get the name and type name of the shift: the type name is cut to
    // keep only the first day
    file >> shName;
    shName.pop_back();
    shiftToInt[shName] = i;
    Tools::readUntilChar(&file, ',', &strTmp);
    std::stringstream sstream(strTmp);
    sstream >> typeName;

    // get the time interval and duration of the shift
    Tools::readUntilChar(&file, ',', &strTmp);
    Tools::Time timeStart(Tools::readHourFromStr(strTmp));
    Tools::readUntilChar(&file, ',', &strTmp);
    Tools::Time timeEnd(Tools::readHourFromStr(strTmp));
    shiftDurations.push_back(1.0);
    // look at previous shift types to see if it has already been created
    int type = -1;
    for (int j = 1; j < i; j++) {
      if (timeStart.equals(pShifts[j]->timeStart) &&
          timeEnd.equals(pShifts[j]->timeEnd)) {
        type = pShifts[j]->type;
        break;
      }
    }
    if (type == -1) {
      type = nbShiftTypes++;
      intToShiftType.push_back(typeName);
      shiftTypeToInt[typeName] = type;
      shiftTypeIDToShiftID.push_back(vector<int>());
    }
    shiftIDToShiftTypeID.push_back(type);
    shiftTypeIDToShiftID[type].push_back(i);


    // get the skills of the shift
    int nbSk;
    skillIds.clear();
    file >> nbSk;
    Tools::readUntilChar(&file, ',', &strTmp);
    for (int j = 0; j < nbSk; j++) {
      file >> skName;
      skName.pop_back();
      skillIds.push_back(skillToInt.at(skName));
    }

    // initialize the shift
    pShifts.push_back(std::make_shared<Shift>(
        shName, i, type, shiftDurations[i],
        successorList, skillIds,
        timeStart, timeEnd));
  }
  forbiddenShiftTypeSuccessors.resize(nbShiftTypes);
  forbiddenShiftSuccessors.resize(nbShifts);
  nbForbidShiftSucc = vector<int>(nbShifts, 0);

  // Read the constraints of each contract
  Tools::readUntilChar(&file, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "CONTRACTS "))
    Tools::throwError("The INRC file is not as expected");

  file >> nbContracts;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  // Read each contract type
  for (int i = 0; i < nbContracts; i++) {
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
    PAbstractShift pWork = std::make_shared<AnyWorkShift>();
    if (lbOn || ubOn) {
      pResources.push_back(std::make_shared<SoftTotalShiftDurationResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nbDays,
          lbCost,
          ubCost,
          pWork,
          nbDays,
          *max_element(shiftDurations.begin(), shiftDurations.end())));
    }

    // consecutive working days
    Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                               &ubOn, &ub, &ubCost);
    if (lbOn || ubOn) {
      pResources.push_back(std::make_shared<SoftConsShiftResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nbDays,
          lbCost,
          ubCost,
          pWork,
          CONS_WORK_COST,
          nbDays,
          0,
          true));
    }


    // consecutive free days
    if (lbOn || ubOn) {
      Tools::readBoundedResource(&file, &lbOn, &lb, &lbCost,
                                 &ubOn, &ub, &ubCost);
      pResources.push_back(std::make_shared<SoftConsShiftResource>(
          lbOn ? lb : 0,
          ubOn ? ub : nbDays,
          lbCost,
          ubCost,
          pShifts[0],
          CONS_REST_COST,
          nbDays,
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
          ubOn ? ub : nbWeeks,
          lbCost,
          ubCost,
          pWork,
          nbDays,
          firstWeekendDay,
          lastWeekendDay,
          true));
    }

    if (totalWkndOn) {
      pResources.push_back(std::make_shared<SoftTotalWeekendsResource>(
          totalWkndUb,
          totalWkndCost,
          pWork,
          nbDays,
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
          std::make_shared<ShiftComparator>(), ubCost,
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
          std::make_shared<AnyOfTypeShift>(shiftTypeToInt.at("Night")));
      for (int d = 1; d < pADays.size(); d++)
        pAShifts.push_back(std::make_shared<AnyRestShift>());
      pResources.push_back(std::make_shared<SoftForbiddenPatternResource>(
          Pattern(pAShifts, pADays),
          ubCost));
    }

    // two free days after night shift
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isOn;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> ubCost;
    if (isOn) {
      PAbstractShift pAShift;
      bool isFound = false;
      for (const auto &pS : pShifts) {
        if (intToShiftType[pS->type] == "Night") {
          pAShift = pS;
          isFound = true;
          break;
        }
      }
      if (!isFound) {
        Tools::throwError("No shift with name Night was found");
      }
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
    int nbPatterns;
    vector<int> patternIds;
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> nbPatterns;
    Tools::readUntilChar(&file, ',', &strTmp);
    for (int j = 0; j < nbPatterns; j++) {
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
  readINRCPatterns(shiftToInt, pContracts, &file, pShifts);

  // Read the list of nurses
  // BEWARE: there is no equal sign before the number of employees!
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  Tools::readUntilChar(&file, ' ', &strTmp);
  if (!Tools::strEndsWith(strTmp, "EMPLOYEES"))
    Tools::throwError("The INRC file is not as expected");

  file >> nbNurses;
  Tools::readUntilAndWhileChar(&file, '/', &strTmp);
  for (int i = 0; i < nbNurses; i++) {
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
    for (int i = 0; i < nbSkills; ++i) skills.push_back(i);
    theNurses.push_back(std::make_shared<Nurse>(
        nurseId, nurseName, nbShifts, nbSkills, skills, availableShifts,
        pContract));
    nurseNameToInt[nurseName] = i;

    if (pContract->alternativeShift_)
      theNurses.back()->addBaseResource(std::make_shared<
          AlternativeShiftResource>(
          theNurses.back(), nbShifts));
  }
  // sort the nurses per id to avoid errors
  std::stable_sort(theNurses.begin(),
                   theNurses.end(),
                   compareNursesById);

  // Check that all fields were initialized before initializing the scenario
  //
  if (nbWeeks == -1 || nbSkills == -1 || nbShifts == -1 || nbContracts == -1
      || nbNurses == -1) {
    Tools::throwError("In readINRCInstances: missing fields in the "
                      "initialization");
  }

  PScenario pScenario = std::make_shared<Scenario>(instName,
                                                   nbWeeks,
                                                   nbSkills,
                                                   intToSkill,
                                                   skillToInt,
                                                   nbShifts,
                                                   shiftToInt,
                                                   shiftDurations,
                                                   shiftIDToShiftTypeID,
                                                   nbShiftTypes,
                                                   intToShiftType,
                                                   shiftTypeToInt,
                                                   shiftTypeIDToShiftID,
                                                   minConsShiftType,
                                                   maxConsShiftType,
                                                   nbForbidShiftSucc,
                                                   forbiddenShiftSuccessors,
                                                   pShifts,
                                                   nbContracts,
                                                   pContracts,
                                                   nbNurses,
                                                   theNurses,
                                                   nurseNameToInt,
                                                   std::make_shared<Weights>(),
                                                   true,
                                                   false);
  pScenario->setStartDate(tmStart);


  // Read the demand
  PDemand pDemand;
  readINRCDemand(nbDays,
                 nbWeeks,
                 nbSkills,
                 nbShifts,
                 tmStart,
                 shiftToInt,
                 pShifts,
                 &file,
                 &pDemand);

  // Read preferences
  PPreferences pPref = std::make_shared<Preferences>(nbNurses,
                                                     nbDays,
                                                     nbShifts);
  readINRCPreferences(&file,
                      shiftToInt,
                      pShifts,
                      theNurses,
                      tmStart,
                      pPref);

  // Initialize the history of every nurse to empty state
  vector<State> initialState;
  // Add a fictitious shift just for the initial state
  PShift noneShift = std::make_shared<Shift>(
      "NoneShift", nbShifts, nbShiftTypes, 0, successorList, false, false);
  for (int n = 0; n < nbNurses; n++) {
    // TODO(AL): we possibly have an issue here, because the initial history
    //  is empty in INRC, so I added a shift that is not rest nor work. Not
    //  sure it is 100% sufficient
    /* if (consRest == 0 && consDaysWorked == 0)
       Tools::throwError("History of nurse %s is invalid as one must "
                        "either work or rest.", nurseName.c_str()); */
    State nurseState(0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     noneShift);
    initialState.push_back(nurseState);
  }

  // link the scenario to the demand, preferences and history
  pScenario->linkWithDemand(pDemand);
  pScenario->linkWithPreferences(pPref);
  pScenario->setInitialState(initialState);

  return pScenario;
}

void ReadWrite::readINRCDemand(int nbDays,
                               int nbWeeks,
                               int nbSkills,
                               int nbShifts,
                               const tm &tmStart,
                               const map<string, int> &shiftToInt,
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

  int nbCovers;
  *pFile >> nbCovers;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);

  // init the demand vectors
  vector3D<int> minDemands;
  Tools::initVector3D(&minDemands,
                      nbDays,
                      nbShifts,
                      nbSkills,
                      0);
  for (int i = 0; i < nbCovers; i++) {
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
    for (int w = 0; w < nbWeeks; w++) {
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

  *pFile >> nbCovers;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);

  for (int i = 0; i < nbCovers; i++) {
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
  *pDemand = std::make_shared<Demand>(nbDays,
                                      0,
                                      nbShifts,
                                      nbSkills,
                                      "allWeeks",
                                      minDemands);
}

void ReadWrite::readINRCPatterns(const map<string, int> &shiftToInt,
                                 const vector<PContract> &pContracts,
                                 std::fstream *pFile,
                                 const vector<PShift> &pShifts) {
  vector<SoftForbiddenPatternResource> patternResources;

  string strTmp;
  int intTmp;
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "PATTERNS "))
    Tools::throwError("The INRC file is not as expected: "
                      "PATTERNS not found.");

  *pFile >> intTmp;
  int nbPatterns = intTmp;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  // Read the patterns and create one resources for each pattern
  for (int i = 0; i < nbPatterns; i++) {
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
        pAShifts.push_back(std::make_shared<AnyWorkShift>());
      } else {
        int shiftId = shiftToInt.at(strTmp);
        pAShifts.push_back(pShifts[shiftId]);
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
  for (const auto& pC : pContracts) {
    for (int patternId : pC->forbiddenPatternIds_) {
      SoftForbiddenPatternResource r = patternResources[patternId];
      pC->addResource(std::make_shared<SoftForbiddenPatternResource>(
          r.getPattern(), r.getCost()));
    }
  }
}

// Read the preferences of all the nurses in an INRC input file
void ReadWrite::readINRCPreferences(std::fstream *pFile,
                                    const map<string, int> &shiftToInt,
                                    const vector<PShift> &pShifts,
                                    const vector<PNurse> &theNurses,
                                    const tm &tmStart,
                                    const PPreferences &pPref) {
  // 1. ready day on preferences
  string strTmp;
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DAY_OFF_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "DAY_OFF_REQUESTS not found.");

  int nbRequests;
  *pFile >> nbRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nbRequests; i++) {
    int nurseId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> cost;

    const Wish &wish = pPref->addDayOff(nurseId, dayId, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }


  // 2. read day on preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "DAY_ON_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "DAY_ON_REQUESTS not found.");

  *pFile >> nbRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nbRequests; i++) {
    int nurseId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> cost;

    pPref->addDayOn(nurseId, dayId, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // 3. read shift off preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SHIFT_OFF_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "SHIFT_OFF_REQUESTS not found.");

  *pFile >> nbRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nbRequests; i++) {
    int nurseId, shiftId;
    double cost;
    *pFile >> nurseId;
    Tools::readUntilChar(pFile, ',', &strTmp);
    Tools::readUntilChar(pFile, ',', &strTmp);
    const tm tmRequest(*Tools::readDateFromStr(strTmp));
    int dayId(tmRequest.tm_yday - tmStart.tm_yday);
    *pFile >> strTmp;
    strTmp.pop_back();
    shiftId = shiftToInt.at(strTmp);
    *pFile >> cost;

    const PAbstractShift &pAS = pShifts[shiftId];
    pPref->addShiftOff(nurseId, dayId, pAS, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // 4. read shift on preferences
  Tools::readUntilChar(pFile, '=', &strTmp);
  if (!Tools::strEndsWith(strTmp, "SHIFT_ON_REQUESTS "))
    Tools::throwError("The INRC file is not as expected: "
                      "SHIFT_ON_REQUESTS not found.");

  *pFile >> nbRequests;
  Tools::readUntilAndWhileChar(pFile, '/', &strTmp);
  for (int i = 0; i < nbRequests; i++) {
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
    shiftId = shiftToInt.at(strTmp);
    *pFile >> cost;

    const PAbstractShift &pAS = pShifts[shiftId];
    pPref->addShiftOn(nurseId, dayId, pAS, cost);
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // create resources
  for (const PNurse &pN : theNurses)
    for (const auto &p : pPref->nurseWishes(pN->num_)) {
      pN->addBaseResource(
          std::make_shared<SoftPreferenceResource>(
              std::make_shared<Day>(p.first),
              p.second));
  }
}

//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
PScenario ReadWrite::readScenarioINRC2(const string& fileName) {
  std::fstream file;
  Tools::openFile(fileName, &file);
  string title;
  string strTmp;
  int intTmp;
  // declare the attributes that will initialize the Scenario instance
  //
  string name;
  int nbDays = -1, nbWeeks = -1, nbSkills = -1, nbShifts = -1,
      nbShiftsType = -1, nbContracts = -1, nbNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType;
  map<string, int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType,
      shiftIDToShiftTypeID, nbForbidShiftSucc;
  vector<double> shiftDurations;
  vector2D<int> shiftTypeIDToShiftID,
      forbiddenShiftTypeSuccessors,
      forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  vector<PContract> pContracts;
  std::map<string, PContract> pContractsByName;
  vector<PNurse> theNurses;

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

  file >> nbWeeks;
  nbDays = 7 * nbWeeks;

  // 3. read the skills of the scenario
  Tools::readUntilChar(&file, '=', &title);
  if (!Tools::strEndsWith(title, "SKILLS "))
    Tools::throwError("The INRC2 file is not as expected: "
                      "%s found instead of SKILLS.", title.c_str());

  file >> nbSkills;
  for (int i = 0; i < nbSkills; i++) {
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
  nbShiftsType = intTmp + 1;

  // IMPORTANT : INSERT REST SHIFT !!!!!!
  // It is given 0 and 99 as bounds so that they never perturbate the cost
  //
  if (REST_SHIFT_ID != 0)
    Tools::throwError("Scenario reader works only with REST_SHIFT_ID = 0.");
  intToShiftType.push_back(REST_SHIFT);
  intToShift.push_back(REST_SHIFT);
  shiftTypeToInt.insert(pair<string, int>(REST_SHIFT, 0));
  shiftToInt.insert(pair<string, int>(REST_SHIFT, 0));
  shiftTypeIDToShiftID.resize(nbShiftsType);
  shiftIDToShiftTypeID.push_back(0);
  shiftTypeIDToShiftID[0].push_back(0);
  shiftDurations.push_back(0.0);
  minConsShiftType.push_back(0);
  maxConsShiftType.push_back(99);

  // Other shift types
  //
  for (int i = 1; i < nbShiftsType; i++) {
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
    Tools::readUntilChar(&file, '\n', &strTmp);
  }

  // 4.b. Forbidden successions for the shift type
  while (!file.eof() &&
         !Tools::strEndsWith(title, "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
    file >> title;

  if (file.eof())
    Tools::throwError("The INRC2 file is not as expected: "
                      "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS not found");

  forbiddenShiftTypeSuccessors.resize(nbShiftsType);
  // Reading all lines
  for (int i = 1; i < nbShiftsType; i++) {
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
    Tools::readUntilChar(&file, '\n', &strTmp);
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
    nbShifts = intTmp + 1;  // +1 for the REST_SHIFT

    // Shifts
    for (int i = 1; i < nbShifts; i++) {
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
      Tools::readUntilChar(&file, '\n', &strTmp);
    }
    // read for the next field
    Tools::readUntilChar(&file, '=', &title);
  } else {
    // Create one shift for each shift type
    nbShifts = nbShiftsType;
    for (int i = 1; i < nbShifts; i++) {
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
  forbiddenShiftSuccessors.resize(nbShifts);
  for (int i = 0; i < nbShiftsType; i++)
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
      Tools::readUntilChar(&file, '\n', &strTmp);
    }
    Tools::readUntilChar(&file, '=', &title);
  }

  // 4.e. Create shift structures
  nbForbidShiftSucc.resize(nbShifts);
  for (int i = 0; i < nbShifts; i++) {
    // make sure the forbidden successors are sorted in increasing order
    vector<int> &f = forbiddenShiftSuccessors[i];
    std::sort(f.begin(), f.end());
    // number of forbidden successors
    nbForbidShiftSucc[i] = f.size();
    std::vector<int> successorList;
    for (int s = 0; s < nbShifts; ++s) {
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
  nbContracts = intTmp;
  // read each contract type
  for (int i = 0; i < nbContracts; i++) {
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
    Tools::readUntilChar(&file, '\n', &strTmp);

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

  file >> nbNurses;
  for (int i = 0; i < nbNurses; i++) {
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
      if (nbSkills > 1)
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
      for (int s = 0; s < nbShifts; s++)
        availableShifts.push_back(s);
    } else {
      availableShifts.push_back(0);  // add the rest shift
      std::sort(availableShifts.begin(), availableShifts.end());
    }

    theNurses.push_back(std::make_shared<Nurse>(
        i, nurseName, nbShifts, nbSkills, skills, availableShifts,
        pContractsByName.at(contractName)));
    nurseNameToInt.insert(pair<string, int>(nurseName, i));
  }

  // Check that all fields were initialized before initializing the scenario
  //
  if (nbWeeks == -1 || nbSkills == -1 || nbShifts == -1 || nbContracts == -1
      || nbNurses == -1) {
    Tools::throwError("In readScenarioINRC2: missing fields in the "
                      "initialization");
  }

  return std::make_shared<Scenario>(name,
                                    nbWeeks,
                                    nbSkills,
                                    intToSkill,
                                    skillToInt,
                                    nbShifts,
                                    shiftToInt,
                                    shiftDurations,
                                    shiftIDToShiftTypeID,
                                    nbShiftsType,
                                    intToShiftType,
                                    shiftTypeToInt,
                                    shiftTypeIDToShiftID,
                                    minConsShiftType,
                                    maxConsShiftType,
                                    nbForbidShiftSucc,
                                    forbiddenShiftSuccessors,
                                    pShifts,
                                    nbContracts,
                                    pContracts,
                                    nbNurses,
                                    theNurses,
                                    nurseNameToInt,
                                    pWeights,
                                    false,
                                    true);
}

PDemand ReadWrite::readINRC2Weeks(const std::vector<std::string>& strWeekFiles,
                                  const PScenario& pScenario) {
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
      pScenario->addAWeek();
    }

  // link the scenario to the current demand and preferences
  pScenario->linkWithDemand(pDemand);
  pScenario->linkWithPreferences(pPref);

  return pDemand;
}

// Read the Week file and store the content in a Scenario instance
//
void ReadWrite::readWeekINRC2(const std::string& strWeekFile,
                              const PScenario& pScenario,
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
        Tools::readUntilChar(&file, '\n', &strTmp);
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
      int nbShifts, nurseNum;
      DayOfWeek dayOfWeek;
      PREF_LEVEL level = WEAK;
      file >> nbShifts;
      for (int i = 0; i < nbShifts; i++) {
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
        } catch (const std::invalid_argument& ia) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addDayOff(nurseNum, dayOfWeek, prefCosts.at(level));
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
      int nbShifts, nurseNum;
      DayOfWeek dayOfWeek;
      PREF_LEVEL level = WEAK;
      file >> nbShifts;
      for (int i = 0; i < nbShifts; i++) {
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
        } catch (const std::invalid_argument& ia) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addDayOn(nurseNum, dayOfWeek, prefCosts.at(level));
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
#ifdef DBG
  std::cout << "Demand created" << std::endl;
#endif
}

// Read the history file
//
void ReadWrite::readHistoryINRC2(const std::string& strHistoryFile,
                                 const PScenario& pScenario) {
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
    Tools::readUntilChar(&file, '\n', &title);

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
int ReadWrite::readCustom(const string& strCustomInputFile,
                          const PScenario& pScenario,
                          vector<PDemand> *demandHistory) {
  // open the file
  std::fstream file;
  Tools::openFile(strCustomInputFile, &file);


  string title;
  int nbWeeks;

  // get the custom information
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    // Read the file names of the past demand
    //
    if (!strcmp(title.c_str(), "PAST_DEMAND_FILES")) {
      file >> nbWeeks;
      if (!nbWeeks) continue;

      string strDemandFile;
      for (int i = 0; i < nbWeeks; i++) {
        file >> strDemandFile;
        PDemand pDemand;
        PPreferences pPref;
        readWeekINRC2(strDemandFile, pScenario, &pDemand, &pPref);
        demandHistory->push_back(pDemand);
      }
    }
  }
  return nbWeeks;
}

void ReadWrite::writeCustom(string strCustomOutputFile,
                            const string& strWeekFile,
                            const string& strCustomInputFile) {
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
  int nbWeeks;

  // fill the attributes of the week structure
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    // Read the file names of the past demand
    //
    if (!strcmp(title.c_str(), "PAST_DEMAND_FILES")) {
      file >> nbWeeks;
      outStream << "PAST_DEMAND_FILES= " << nbWeeks + 1 << std::endl;
      if (!nbWeeks) continue;

      string strDemandFile;
      for (int i = 0; i < nbWeeks; i++) {
        file >> strDemandFile;
        outStream << strDemandFile << std::endl;
      }
      outStream << strWeekFile << std::endl;
    }
  }
}

/************************************************************************
* Print the main characteristics of all the demands of an input directory
* This is done to find some invariant properties among demands
*************************************************************************/
void ReadWrite::compareDemands(const string& inputDir, string logFile) {
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
