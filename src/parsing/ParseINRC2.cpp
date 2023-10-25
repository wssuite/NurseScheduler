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

#include "ParseINRC2.h"

#include <dirent.h>
#include <algorithm>
#include <map>
#include <memory>
#include <utility>
#include <string>
#include <vector>

#include "solvers/Solver.h"
#include "solvers/mp/sp/rcspp/resources/UnwantedShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"


using std::string;
using std::vector;
using std::map;
using std::pair;

//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
PScenario readScenarioINRC2(const string &fileName) {
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
            i, nurseName, nSkills, skills, availableShifts,
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
                                    intToSkill,
                                    skillToInt,
                                    pShifts,
                                    intToShiftType,
                                    minConsShiftType,
                                    maxConsShiftType,
                                    pContracts,
                                    pNurses,
                                    pWeights);
}

vector<PDemand> readINRC2Weeks(
        const std::vector<std::string> &strWeekFiles,
        const PScenario &pScenario) {
  // initialize pDemands
  vector<PDemand> pDemands;
  PPreferences pPref;
  for (const string &strWeekFile : strWeekFiles)
    if (pDemands.empty()) {
      pDemands = readWeekINRC2(strWeekFile, pScenario, &pPref);
    } else {
      // load the next week
      PPreferences nextPref;
      vector<PDemand> pNextD =
              readWeekINRC2(strWeekFile, pScenario, &nextPref);
      // update the current weeks
      for (int i = 0; i < pDemands.size(); i++)
        pDemands[i]->pushBack(pNextD[i]);
      pPref->pushBack(nextPref);
    }

  // link the scenario to the current demand and preferences
  pScenario->linkWithDemand(pDemands);
  pScenario->linkWithPreferences(pPref);

  return pDemands;
}

// Read the Week file and store the content in a Scenario instance
//
vector<PDemand> readWeekINRC2(
        const std::string &strWeekFile,
        const PScenario &pScenario,
        PPreferences *pPref) {
  // open the file
  std::fstream file;
  Tools::openFile(strWeekFile, &file);

  string title;
  string strTmp;
  int intTmp;

  // declare the attributes to be updated in the PScenario
  int nDays = 7;  // read always week by week
  string weekName;
  vector3D<int> minWeekDemand;
  vector3D<int> optWeekDemand;
  const Weights &weights = pScenario->weights();

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
                          nDays,
                          pScenario->nShifts(),
                          pScenario->nSkills(),
                          0);
      Tools::initVector3D(&optWeekDemand,
                          nDays,
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
        for (int day = 0; day < nDays; day++) {
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
        *pPref = std::make_shared<Preferences>(
                pScenario->nNurses(), nDays, pScenario->nShifts());

      // Temporary vars
      const std::vector<double> &prefCosts = weights.preferences;
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
        // check if a number
        Tools::trim(&strLevel);
        if (!strLevel.empty() && isdigit(strLevel[0]))
          level = (PREF_LEVEL) std::stoi(strLevel);
        // else, it's the next line: strLevel contains the next nurse name
        else
          nurseName = strLevel;

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
        *pPref = std::make_shared<Preferences>(
                pScenario->pNurses(), nDays, pScenario->nShifts());

      // Temporary vars
      const std::vector<double> &prefCosts = weights.preferences;
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

#ifdef NS_DEBUG
  std::cout << "Demand created" << std::endl;
#endif

  return {
          std::make_shared<Demand>(
                  nDays, 0, pScenario->nShifts(), pScenario->nSkills(),
                  weekName, minWeekDemand, D_GE),
          std::make_shared<Demand>(
                  nDays, 0, pScenario->nShifts(), pScenario->nSkills(),
                  weekName, optWeekDemand, D_GE, weights.underCoverage)
  };
}

// Read the history file
//
void readHistoryINRC2(const std::string &strHistoryFile,
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
  int thisWeek;
  string weekName;
  vector<State> initialState;


  // fill with the attributes of the week structure
  //
  while (file.good()) {
    Tools::readLine(&file, &title);

    // Read the index and name of the week
    //
    if (!strcmp(title.c_str(), "HISTORY")) {
      file >> thisWeek;
      file >> weekName;
      // Raise exception if it does not match the week previously read !
      if (strcmp(weekName.c_str(), (pScenario->demandName()).c_str()) != 0) {
        std::cout << "The given history file requires week " << weekName
                  << std::endl;
        std::cout << " but a different one (" << pScenario->demandName()
                  << ") has been given!" << std::endl;
        Tools::throwError("History file and week data file do not match!");
      }
    } else if (Tools::strEndsWith(title, "NURSE_HISTORY")) {
      // Read each nurse's initial state
      //
      for (int n = 0; n < pScenario->nNurses(); n++) {
        string nurseName, shiftName;
        // int nurseNum;
        int shiftId, totalTimeWorked, totalDaysWorked, totalWeekendsWorked,
                consDaysWorked, consShiftWorked, consRest, consShifts;
        file >> nurseName;
        // nurseNum = pScenario->nurseNameToInt_.at(nurseName);
        file >> totalDaysWorked;
        totalTimeWorked = totalDaysWorked;
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
        State nurseState(pScenario->pShift(shiftId),
                         0,
                         totalTimeWorked,
                         totalDaysWorked,
                         totalWeekendsWorked,
                         consDaysWorked,
                         consShifts,
                         consRest);
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
int readCustom(const string &strCustomInputFile,
               const PScenario &pScenario,
               vector2D<PDemand> *pDemandsHistory) {
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
        PPreferences pPref;
        vector<PDemand> pDemands =
                readWeekINRC2(strDemandFile, pScenario, &pPref);
        pDemandsHistory->push_back(pDemands);
      }
    }
  }
  return nWeeks;
}

void writeCustom(string strCustomOutputFile,
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

/************************************************************************
* Print the main characteristics of all the demands of an input directory
* This is done to find some invariant properties among demands
*************************************************************************/
void compareDemands(const string &inputDir, string logFile) {
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

  PScenario pScen = readScenarioINRC2(scenFile);

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

    PPreferences pPref;
    vector<PDemand> pDemands = readWeekINRC2(filepath, pScen, &pPref);

    logStream << "#####################################\n";
    logStream << "# DEMAND FILE: " << filepath << std::endl;
    logStream << "#####################################\n\n";
    for (const auto &pD : pDemands)
      logStream << pD->toString(true) << std::endl << std::endl;

    // record the advanced data on the demand
    minTotal.push_back(pDemands[0]->demandTotal_);
    optTotal.push_back(pDemands[1]->demandTotal_);
    minPerShift.push_back(pDemands[0]->demandPerShift_);
    optPerShift.push_back(pDemands[1]->demandPerShift_);
    minPerSkill.push_back(pDemands[0]->demandPerSkill_);
    optPerSkill.push_back(pDemands[1]->demandPerSkill_);
    minHighestPerSkill.push_back(pDemands[0]->demandHighestPerSkill_);
    optHighestPerSkill.push_back(pDemands[1]->demandHighestPerSkill_);

    // link the scenario with the first demand and preferences to be able to
    // retrieve information about the nurses
    if (!coDemand) {
      pScen->linkWithDemand(pDemands);
      pScen->linkWithPreferences(pPref);
    } else {
      coDemand++;
    }
  }

  // Also presolve the nurses to get statistics on the capacity of the nurses
  // to cover the demand
  readHistoryINRC2(historyFile, pScen);
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
    logStream << Tools::itoa(minTotal[d] - pSolver->maxTotalStaffNoPenalty_)
          + "/" + Tools::itoa(optTotal[d] - pSolver->maxTotalStaffNoPenalty_);
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
      stdDevMin += pow(minPerSkill[d][i] -
              pSolver->maxStaffPerSkillNoPenalty_[i], 2);
      stdDevOpt += pow(optPerSkill[d][i] -
              pSolver->maxStaffPerSkillNoPenalty_[i], 2);
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
