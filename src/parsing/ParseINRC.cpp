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

#include "ParseINRC.h"

#include <algorithm>
#include <memory>
#include <vector>
#include <string>

#include "solvers/mp/sp/rcspp/resources/UnwantedShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

//--------------------------------------------------------------------------
// Method that read an INRC input files and stores the content in the
// output scenario instance
//
PScenario readINRCInstance(const string &fileName) {
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
  map<string, double> altShiftsPerContract;
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

    // create type
    int type = type = nShiftTypes++;
    intToShiftType.push_back(typeName);
    shiftTypeToInt[typeName] = type;
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
              true));
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
              0, totalWkndUb,
              0, totalWkndCost,
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
    bool isUnwantedShift;
    double alternativeShiftCost;
    Tools::readUntilChar(&file, '(', &strTmp);
    file >> isUnwantedShift;
    Tools::readUntilChar(&file, '|', &strTmp);
    file >> alternativeShiftCost;
    // when alternative is false, any skill can be used by any nurse
    if (!isUnwantedShift)
      alternativeShiftCost = 0;
    altShiftsPerContract[contractName] = alternativeShiftCost;


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
            vector<int>(),
            vector<double>(),
            firstWeekendDay,
            lastWeekendDay,
            patternIds));
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
    try {
      // if name is an integer, add nurse in front
      std::stoi(nurseName);
      nurseName = "Nurse " + nurseName;
    } catch(...) {}
    file >> intTmp;
    PContract pContract = pContracts[intTmp];
    Tools::readUntilChar(&file, ',', &strTmp);
    file >> intTmp;
    Tools::readUntilChar(&file, ' ', &strTmp);
    vector<int> skills;
    // read the skills
    for (int j = 0; j < intTmp; j++) {
      Tools::readUntilOneOfTwoChar(&file, ' ', ';', &strTmp);
      skills.push_back(skillToInt.at(strTmp));
    }
    // sort the skill indices before initializing the nurse
    std::sort(skills.begin(), skills.end());

    // add the rest shift and sort the available shifts
    vector<int> availableShifts;
    // add all shifts (default for INRCI),
    // and make the alternative shift to be paid with UnwantedShiftResource
    vector<PShift> unwantedShifts;
    for (const auto &pS : pShifts) {
      availableShifts.push_back(pS->id);
      if (!std::any_of(skills.begin(), skills.end(),
                       [pS](int sk) { return pS->hasSkill(sk); }))
        unwantedShifts.push_back(pS);
    }
    std::sort(availableShifts.begin(), availableShifts.end());
    // check if rest is available
    if (availableShifts.front() != 0)
      Tools::throwError("Rest shift is not available for nurse %d", nurseId);

    // as the skills are included in the shift,
    // give all the skills to each nurses
    skills.clear();
    for (int sk = 0; sk < nSkills; ++sk) skills.push_back(sk);
    pNurses.push_back(std::make_shared<Nurse>(
            nurseId, nurseName, nSkills, skills, availableShifts,
            pContract));
    nurseNameToInt[nurseName] = i;

    if (altShiftsPerContract[pContract->name_] > 1e-3) {
      pNurses.back()->addBaseResource(std::make_shared<UnwantedShiftResource>(
              nShifts, unwantedShifts, altShiftsPerContract[pContract->name_]));
    }
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
                                                   intToSkill,
                                                   skillToInt,
                                                   pShifts,
                                                   intToShiftType,
                                                   pContracts,
                                                   pNurses);
  pScenario->setStartDate(tmStart);


  // Read the demand
  PDemand pDemand = readINRCDemand(
          nDays, nWeeks, nSkills, nShifts, tmStart, pShifts, &file);

  // Read preferences
  PPreferences pPref = std::make_shared<Preferences>(nNurses, nDays, nShifts);
  readINRCPreferences(&file, pScenario, pNurses, tmStart, pPref);

  // Initialize the history of every nurse to empty state
  vector<State> initialState;
  // Add a fictitious shift just for the initial state
  for (int n = 0; n < nNurses; n++)
    initialState.push_back(State::noneState(shiftsFactory));

  // link the scenario to the demand, preferences and history
  pScenario->linkWithDemand({pDemand});
  pScenario->linkWithPreferences(pPref);
  pScenario->setInitialState(initialState);

  return pScenario;
}

PDemand readINRCDemand(
        int nDays, int nWeeks, int nSkills, int nShifts, const tm &tmStart,
        const vector<PShift> &pShifts, std::fstream *pFile) {
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
  vector3D<int> demands;
  Tools::initVector3D(&demands, nDays, nShifts, nSkills, 0);
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
      demands[dayId][shiftId][skillId] = intTmp;
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
    demands[dayId][shiftId][skillId] = intTmp;
    Tools::readUntilChar(pFile, ';', &strTmp);
  }

  // return hard pDemand (hard by default)
  return std::make_shared<Demand>(
          nDays, 0, nShifts, nSkills, "allWeeks", demands, D_EQ);
}

void readINRCPatterns(const vector<PContract> &pContracts,
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
void readINRCPreferences(std::fstream *pFile,
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
