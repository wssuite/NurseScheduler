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
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <streambuf>
#include <string>
#include <utility>

#include "tools/Tools.h"
#include "data/Scenario.h"
#include "solvers/Solver.h"
#include "solvers/StochasticSolver.h"

#include <boost/assign/list_of.hpp>

using std::string;
using std::vector;
using std::map;
using std::pair;

std::map<std::string, WeightStrategy> stringToWeightStrategy =
    boost::assign::map_list_of("MAX", MAX)("MEAN", MEAN)("RANDOMMEANMAX",
                                                         RANDOMMEANMAX)(
        "BOUNDRATIO",
        BOUNDRATIO)("NO_STRAT", NO_STRAT);
std::map<std::string, RankingStrategy> stringToRankingStrategy =
    boost::assign::map_list_of("SCORE", RK_SCORE)("MEAN", RK_MEAN);


//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scenario file and store the content in a Scenario instance
//
PScenario ReadWrite::readScenario(string fileName) {
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

  string title;
  string strTmp;
  int intTmp;
  // declare the attributes that will initialize the Scenario instance
  //
  string name;
  int nbWeeks = -1, nSkills = -1, nShifts = -1, nbShiftsType = -1,
      nbContracts = -1, nbNurses = -1;
  vector<string> intToSkill, intToShift, intToShiftType, intToContract;
  map<string, int> skillToInt, shiftToInt, shiftTypeToInt, nurseNameToInt;
  vector<int> minConsShiftType, maxConsShiftType, hoursInShift,
      shiftIDToShiftTypeID, nbForbidShiftSucc;
  vector2D<int> shiftTypeIDToShiftID,
      forbiddenShiftTypeSuccessors,
      forbiddenShiftSuccessors;
  vector<PShift> pShifts;
  map<string, PConstContract> contracts;
  vector<PNurse> theNurses;

  bool foundShift = false;

  // default weights
  PWeights pWeights = std::make_shared<Weights>();

  // fill the attributes of the scenario structure
  //
  while (file.good()) {
    Tools::readUntilChar(&file, '=', &title);

    // Read the name of the scenario
    //
    if (Tools::strEndsWith(title, "SCENARIO ")) {
      file >> name;
    } else if (Tools::strEndsWith(title, "WEEKS ")) {
      // Read the number of weeks in scenario
      //
      file >> nbWeeks;
    } else if (Tools::strEndsWith(title, "SKILLS ")) {
      // Read the number of weeks in scenario
      //
      file >> nSkills;
      for (int i = 0; i < nSkills; i++) {
        file >> strTmp;
        intToSkill.push_back(strTmp);
        skillToInt.insert(pair<string, int>(strTmp, i));
      }
    } else if (Tools::strEndsWith(title, "SHIFT_TYPES ")) {
      // Read the different shift types and forbidden successions
      //
      // Number of shifts : Given number + REST_SHIFT
      file >> intTmp;
      nbShiftsType = intTmp + 1;

      // IMPORTANT : INSERT REST SHIFT !!!!!!
      // It is given 0 and 99 as bounds so that they never perturbate the cost
      //
      intToShiftType.push_back(REST_SHIFT);
      shiftTypeToInt.insert(pair<string, int>(REST_SHIFT, 0));
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

      while (!Tools::strEndsWith(title, "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"))
        file >> title;

      // Forbidden successions
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
    } else if (Tools::strEndsWith(title, "SHIFTS ")) {
      // Read the different shifts
      //
      foundShift = true;

      // IMPORTANT : INSERT REST SHIFT !!!!!!
      //
      intToShift.push_back(REST_SHIFT);
      shiftToInt.insert(pair<string, int>(REST_SHIFT, 0));
      hoursInShift.push_back(0);
      int currentShiftTypeId = shiftTypeToInt.at(REST_SHIFT);
      shiftIDToShiftTypeID.push_back(currentShiftTypeId);

      shiftTypeIDToShiftID.resize(nbShiftsType);
      shiftTypeIDToShiftID[0].push_back(0);

      // Number of shifts
      file >> intTmp;
      nShifts = intTmp + 1;

      // Shifts
      //
      for (int i = 1; i < nShifts; i++) {
        // Name
        file >> strTmp;
        intToShift.push_back(strTmp);
        shiftToInt.insert(pair<string, int>(strTmp, i));

        int hours;
        file >> hours;
        hoursInShift.push_back(hours);

        string currentShiftType;
        file >> currentShiftType;

        int currentShiftTypeId = shiftTypeToInt.at(currentShiftType);
        shiftIDToShiftTypeID.push_back(currentShiftTypeId);
        shiftTypeIDToShiftID[currentShiftTypeId].push_back(i);
      }
    }  else if (Tools::strEndsWith(title, "FORBIDDEN_SHIFT_SUCCESSIONS ")) {
      // Forbidden successions
      //
      forbiddenShiftSuccessors.resize(nShifts);
      // Reading all lines
      int n;
      file >> n;
      for (int i = 0; i < n; i++) {
        // Which current shift type ?
        string currentShift;
        file >> currentShift;
        int currentShiftId = shiftToInt.at(currentShift);
        // How many forbidden after it ?
        file >> intTmp;
        // Which ones are forbidden ?
        for (int j = 0; j < intTmp; j++) {
          file >> strTmp;
          forbiddenShiftSuccessors[currentShiftId].push_back(
              shiftToInt.at(strTmp));
        }
        // read end of line
        Tools::readUntilChar(&file, '\n', &strTmp);
      }
    } else if (Tools::strEndsWith(title, "CONTRACTS ")) {
      // Read the different contracts type
      //
      file >> intTmp;
      nbContracts = intTmp;
      // Read each contract type
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

        PConstContract pContract =
            std::make_shared<Contract>(i,
                                       contractName,
                                       minDays,
                                       maxDays,
                                       minConsWork,
                                       maxConsWork,
                                       minConsRest,
                                       maxConsRest,
                                       maxWeekends,
                                       isTotalWeekend,
                                       pWeights);
        contracts[contractName] = pContract;
        intToContract.push_back(contractName);
      }
    } else if (Tools::strEndsWith(title, "NURSES ")) {
      // Define shifts by default
      //  to be backward compatible with old style of input file
      //  (without the SHIFTS section)
      //  set default shifts
      if (!foundShift) {
        nShifts = nbShiftsType;
        intToShift = intToShiftType;
        shiftToInt = shiftTypeToInt;
        shiftTypeIDToShiftID.resize(nbShiftsType);
        for (int i = 0; i < nbShiftsType; i++) {
          // 1 as default (could be days, hours, ...) for non rest shift (>0)
          hoursInShift.push_back(i > 0);
          shiftIDToShiftTypeID.push_back(i);
          shiftTypeIDToShiftID[i].push_back(i);
        }
      }
      // Read all nurses
      //
      file >> nbNurses;
      for (int i = 0; i < nbNurses; i++) {
        string nurseName, contractName;
        int n;
        vector<int> skills, availableShifts;
        // Read everything on the line
        file >> nurseName;
        file >> contractName;
        file >> n;
        // either a skill, a shift or a shiftType
        for (int j = 0; j < n; j++) {
          file >> strTmp;
          // if a skill
          auto it = skillToInt.find(strTmp);
          if (it != skillToInt.end()) {
            skills.push_back(it->second);
            continue;
          }
          // if a shift
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

        // if no shifts, put all of them, otherwise sort them
        if (availableShifts.empty()) {
          for (int s = 0; s < nShifts; s++)
            availableShifts.push_back(s);
        } else {
          availableShifts.push_back(0);  // add the rest shift
          std::sort(availableShifts.begin(), availableShifts.end());
        }

        theNurses.emplace_back(std::make_shared<Nurse>(
            i, nurseName, nShifts, nSkills, skills, availableShifts,
            contracts.at(contractName)));
        nurseNameToInt.insert(pair<string, int>(nurseName, i));
      }
    } else if (title[title.size()-1] != '\n') {
      // if not the end of the file
      Tools::throwError("Field %s not recognized.", title.c_str());
    }
  }

  // Merge forbidden successors from shift type and shift
  nbForbidShiftSucc.resize(nShifts);
  forbiddenShiftSuccessors.resize(nShifts);
  for (int i = 1; i < nShifts; i++) {
    // add all the shifts corresponding to a forbidden type successor
    vector<int> &succs = forbiddenShiftSuccessors[i];
    for (int st : forbiddenShiftTypeSuccessors[shiftIDToShiftTypeID[i]])
      succs.insert(succs.end(),
          shiftTypeIDToShiftID[st].begin(), shiftTypeIDToShiftID[st].end());
    // make sure the forbidden successors are sorted in increasing order
    std::sort(succs.begin(), succs.end());
    // remove duplicate if any
    succs.erase(unique(succs.begin(), succs.end()), succs.end());
    // set nbForbidShiftSucc
    nbForbidShiftSucc[i] = succs.size();
  }

  // Create shift structures
  //
  for (int i = 0; i < nShifts; i++) {
    std::vector<int> successorList;
    vector<int> &f = forbiddenShiftSuccessors[i];
    for (int s = 0; s < nShifts; ++s) {
      if (find(f.begin(), f.end(), s) == f.end())
        successorList.push_back(s);
    }
    int type = shiftIDToShiftTypeID[i];
    pShifts.emplace_back(std::make_shared<Shift>(intToShift[i],
                                                 i,
                                                 type,
                                                 hoursInShift[i],
                                                 successorList,
                                                 minConsShiftType[type],
                                                 maxConsShiftType[type]));
  }

  // Check that all fields were initialized before initializing the scenario
  //
  if (nbWeeks == -1 || nSkills == -1 || nShifts == -1 || nbContracts == -1
      || nbNurses == -1) {
    Tools::throwError("In readScenario: missing fields in the initialization");
  }

  return std::make_shared<Scenario>(name,
                                    nbWeeks,
                                    nSkills,
                                    intToSkill,
                                    skillToInt,
                                    nShifts,
                                    intToShift,
                                    shiftToInt,
                                    hoursInShift,
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
                                    intToContract,
                                    contracts,
                                    nbNurses,
                                    theNurses,
                                    nurseNameToInt,
                                    pWeights);
}

PDemand ReadWrite::readWeeks(std::vector<std::string> strWeekFiles,
                             PScenario pScenario) {
  // initialize pDemand
  PDemand pDemand;
  PPreferences pPref;

  for (const string &strWeekFile : strWeekFiles)
    if (!pDemand) {
      ReadWrite::readWeek(strWeekFile, pScenario, &pDemand, &pPref);
    } else {
      // load the next week
      PDemand nextDemand;
      PPreferences nextPref;
      ReadWrite::readWeek(strWeekFile, pScenario, &nextDemand, &nextPref);
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
void ReadWrite::readWeek(std::string strWeekFile, PScenario pScenario,
                         PDemand *pDemand, PPreferences *pPref) {
  // open the file
  std::fstream file;
  std::cout << "Reading " << strWeekFile << std::endl;
  file.open(strWeekFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read the file " << strWeekFile << std::endl;
    std::cout << "The input file was not opened properly!" << std::endl;
    throw Tools::myException("The input file was not opened properly!",
                             __LINE__);
  }

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
        skillId = pScenario->skill(skillName);
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
      //
      if (!*pPref)
        *pPref = std::make_shared<Preferences>(pScenario->nNurses(),
                                               7,
                                               pScenario->nShifts());
      // Temporary vars
      string nurseName, shift, day, strLevel;
      int nbShifts, nurseNum, dayId;
      PREF_LEVEL level = WEAK;
      file >> nbShifts;
      for (int i = 0; i < nbShifts; i++) {
        if (nurseName.empty())
          file >> nurseName;
        nurseNum = pScenario->nurse(nurseName);
        nurseName.clear();
        file >> shift;
        file >> day;
        dayId = Tools::dayToInt(day);
        // in case there is no level defined for the preferences
        file >> strLevel;
        try {
          level = (PREF_LEVEL) std::stoi(strLevel);
        } catch (std::invalid_argument) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addDayOff(nurseNum, dayId, level);
        } else {
          // shiftId = pScenario->shiftTypeToInt_.at(shift);
          int shiftId = pScenario->shift(shift);
          (*pPref)->addShiftOff(nurseNum, dayId, shiftId, level);
        }
      }
    } else if (Tools::strEndsWith(title, "SHIFT_ON_REQUESTS ")) {
      // Read the shift on requests
      //
      if (!*pPref)
        *pPref = std::make_shared<Preferences>(pScenario->pNurses(),
                                               7,
                                               pScenario->nShifts());
      // Temporary vars
      string nurseName, shift, day, strLevel;
      int nbShifts, nurseNum, dayId;
      PREF_LEVEL level = WEAK;
      file >> nbShifts;
      for (int i = 0; i < nbShifts; i++) {
        if (nurseName.empty())
          file >> nurseName;
        nurseNum = pScenario->nurse(nurseName);
        nurseName.clear();
        file >> shift;
        file >> day;
        dayId = Tools::dayToInt(day);
        // in case there is no level defined for the preferences
        file >> strLevel;
        try {
          level = (PREF_LEVEL) std::stoi(strLevel);
        } catch (std::invalid_argument) {
          // has read the next line: strLevel contains the next nurse name
          nurseName = strLevel;
        }

        if (shift == "Any") {
          (*pPref)->addDayOn(nurseNum, dayId, level);
        } else {
          // shiftId = pScenario->shiftTypeToInt_.at(shift);
          int shiftId = pScenario->shift(shift);
          (*pPref)->addShiftOn(nurseNum, dayId, shiftId, level);
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
  std::cout << "Demand created" << std::endl;
}

// Read the history file
//
void ReadWrite::readHistory(std::string strHistoryFile, PScenario pScenario) {
  if (strHistoryFile.empty()) {
    std::cout << "Cyclic option is enable, so no history file is loaded."
              << std::endl;
    pScenario->enableCyclic();
    return;
  }
  // open the file
  std::fstream file;
  std::cout << "Reading " << strHistoryFile << std::endl;
  file.open(strHistoryFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read " << strHistoryFile << std::endl;
    Tools::throwError("The input file was not opened properly!");
  }

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
      if (strcmp(weekName.c_str(), (pScenario->weekName()).c_str())) {
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
                         pScenario->shiftIDToShiftTypeID(shiftId),
                         shiftId);
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
int ReadWrite::readCustom(string strCustomInputFile,
                          PScenario pScenario,
                          vector<PDemand> *demandHistory) {
  // open the file
  std::fstream file;
  std::cout << "Reading " << strCustomInputFile << std::endl;
  file.open(strCustomInputFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read " << strCustomInputFile << std::endl;
    Tools::throwError("The input file was not opened properly!");
  }

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
        readWeek(strDemandFile, pScenario, &pDemand, &pPref);
        demandHistory->push_back(pDemand);
      }
    }
  }
  return nbWeeks;
}

void ReadWrite::writeCustom(string strCustomOutputFile,
                            string strWeekFile,
                            string strCustomInputFile) {
  Tools::LogOutput outStream(strCustomOutputFile);

  // if there is no custom input file, this is the first week
  if (strCustomInputFile.empty()) {
    outStream << "PAST_DEMAND_FILES= " << 1 << std::endl;
    outStream << strWeekFile << std::endl;
    return;
  }

  // open the custom input file
  // we want the content of the input custom file in the custom output file
  std::fstream file;
  std::cout << "Reading " << strCustomInputFile << std::endl;
  file.open(strCustomInputFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read " << strCustomInputFile << std::endl;
    Tools::throwError("The input file was not opened properly!");
  }

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
* Read the options of the stochastic and ot the other solvers
*************************************************************************/
std::string ReadWrite::readStochasticSolverOptions(
    string strOptionFile, StochasticSolverOptions *options) {

  // open the file
  std::fstream file;
  std::cout << "Reading " << strOptionFile << std::endl;
  file.open(strOptionFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read " << strOptionFile << std::endl;
    Tools::throwError("The input file was not opened properly!");
  }

  string title;

  // fill the attributes of the options structure
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    if (!strcmp(title.c_str(), "withEvaluation")) {
      file >> options->withEvaluation_;
    }
    if (!strcmp(title.c_str(), "withIterativeDemandIncrease")) {
      file >> options->withIterativeDemandIncrease_;
    }
    if (!strcmp(title.c_str(), "generationCostPerturbation")) {
      file >> options->generationCostPerturbation_;
    }
    if (!strcmp(title.c_str(), "evaluationCostPerturbation")) {
      file >> options->evaluationCostPerturbation_;
    }
    if (!strcmp(title.c_str(), "generationAlgorithm")) {
      string strtmp;
      file >> strtmp;
      options->generationAlgorithm_ = AlgorithmsByName.at(strtmp);
    }
    if (!strcmp(title.c_str(), "evaluationAlgorithm")) {
      string strtmp;
      file >> strtmp;
      options->evaluationAlgorithm_ = AlgorithmsByName.at(strtmp);
    }
    if (!strcmp(title.c_str(), "rankingStrategy")) {
      string strtmp;
      file >> strtmp;
      options->rankingStrategy_ = stringToRankingStrategy[strtmp];
    }
    if (!strcmp(title.c_str(), "nExtraDaysGenerationDemands")) {
      file >> options->nExtraDaysGenerationDemands_;
    }
    if (!strcmp(title.c_str(), "nEvaluationDemands")) {
      file >> options->nEvaluationDemands_;
    }
    if (!strcmp(title.c_str(), "nDaysEvaluation")) {
      file >> options->nDaysEvaluation_;
    }
    if (!strcmp(title.c_str(), "nGenerationDemandsMax")) {
      file >> options->nGenerationDemandsMax_;
    }
  }

  std::ifstream fin(strOptionFile.c_str());
  std::ostringstream sout;
  while (fin.good())
    copy(std::istreambuf_iterator<char>(fin),
         std::istreambuf_iterator<char>(),
         std::ostreambuf_iterator<char>(sout));
  return sout.str();
}

std::string ReadWrite::readSolverOptions(string strOptionFile,
                                         SolverParam *options) {
  // open the file
  std::fstream file;
  std::cout << "Reading " << strOptionFile << std::endl;
  file.open(strOptionFile.c_str(), std::fstream::in);
  if (!file.is_open()) {
    std::cout << "While trying to read " << strOptionFile << std::endl;
    Tools::throwError("The input file was not opened properly!");
  }

  string title;

  // fill the attributes of the options structure
  //
  while (file.good()) {
    Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);

    if (!strcmp(title.c_str(), "maxSolvingTimeSeconds")) {
      file >> options->maxSolvingTimeSeconds_;
    }
    if (!strcmp(title.c_str(), "printEverySolution")) {
      file >> options->printEverySolution_;
    }
    if (!strcmp(title.c_str(), "absoluteGap")) {
      file >> options->maxAbsoluteGap_;
    }
    if (!strcmp(title.c_str(), "minRelativeGap")) {
      file >> options->minRelativeGap_;
    }
    if (!strcmp(title.c_str(), "relativeGap")) {
      file >> options->relativeGap_;
    }
    if (!strcmp(title.c_str(), "nbDiveIfMinGap")) {
      file >> options->nbDiveIfMinGap_;
    }
    if (!strcmp(title.c_str(), "nbDiveIfRelGap")) {
      file >> options->nbDiveIfRelGap_;
    }
    if (!strcmp(title.c_str(), "solveToOptimality")) {
      file >> options->solveToOptimality_;
    }
    if (!strcmp(title.c_str(), "weightStrategy")) {
      string strtmp;
      file >> strtmp;
      options->weightStrategy_ = stringToWeightStrategy[strtmp];
    }
    if (!strcmp(title.c_str(), "stopAfterXSolution")) {
      file >> options->stopAfterXSolution_;
    }
  }

  std::ifstream fin(strOptionFile.c_str());
  std::ostringstream sout;
  while (fin.good())
    copy(std::istreambuf_iterator<char>(fin),
         std::istreambuf_iterator<char>(),
         std::ostreambuf_iterator<char>(sout));
  return sout.str();
}

// Read the solution from multiple week solution files
//
vector<Roster> ReadWrite::readSolutionMultipleWeeks(
    vector<std::string> strWeekSolFiles, PScenario pScenario) {
  int nbWeeks = strWeekSolFiles.size();
  string strNurse, strDay, strShift, strSkill;
  int nurse, day, shift, skill;


  // initialize the solution
  vector<Roster> solution;
  for (int n = 0; n < pScenario->nNurses(); n++) {
    vector<int> shifts(7 * nbWeeks, 0);
    vector<int> skills(7 * nbWeeks, -1);

    Roster roster(7 * nbWeeks, 0, shifts, skills);
    solution.push_back(roster);
  }

  // read the solution of each week
  for (int w = 0; w < nbWeeks; w++) {
    std::string strSolFile = strWeekSolFiles[w];
    int firstDay = 7 * w;

    // open the week solution file
    std::fstream file;
    std::cout << "Reading " << strSolFile << std::endl;
    file.open(strSolFile.c_str(), std::fstream::in);
    if (!file.is_open()) {
      std::cout << "While trying to read " << strSolFile << std::endl;
      Tools::throwError("The input file was not opened properly!");
    }
    string title;

    // parse the file until reaching the number of assignments
    std::size_t found = title.find("ASSIGNMENTS");
    while (found != std::string::npos) {
      Tools::readUntilOneOfTwoChar(&file, '\n', '=', &title);
      found = title.find("ASSIGNMENTS");
      if (!file.good()) {
        Tools::throwError("The solution file is not formatted as it should!");
      }
    }
    Tools::readUntilOneOfTwoChar(&file, '\n', '\n', &title);

    // parse the assignments
    while (file.good()) {
      file >> strNurse >> strDay >> strShift >> strSkill;
      nurse = pScenario->nurse(strNurse);
      day = firstDay + Tools::dayToInt(strDay);
      shift = pScenario->shift(strShift);
      skill = pScenario->skill(strSkill);

      solution[nurse].assignTask(day, shift, skill);
    }
  }

  return solution;
}

/************************************************************************
* Print the main characteristics of all the demands of an input directory
* This is done to find some invariant properties among demands
*************************************************************************/
void ReadWrite::compareDemands(string inputDir, string logFile) {
  struct dirent *dirp;
  Tools::LogOutput logStream(logFile, 8);

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
  unsigned found = inputDir.find_last_of("/");
  string instanceName = inputDir.substr(found + 1);
  string scenFile = inputDir + "/Sc-" + instanceName + ".txt";
  string historyFile = inputDir + "/H0-" + instanceName + "-0.txt";

  PScenario pScen = ReadWrite::readScenario(scenFile);

  // Go through all the demand files of the directory
  int coDemand = 0;
  while ((dirp = readdir(dp))) {
    std::string filename(dirp->d_name);

    // The file names of week demands start with "WD"
    std::size_t found = filename.find("WD");
    if (found > 0) continue;

    string filepath = inputDir + "/" + filename;

    PDemand pDemand;
    PPreferences pPref;
    ReadWrite::readWeek(filepath, pScen, &pDemand, &pPref);

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

  // Also preprocess the nurses to get statistics on the capacity of the nurses
  // to cover the demand
  //
  ReadWrite::readHistory(historyFile, pScen);
  Solver *pSolver = new Solver(pScen);
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
    logStream << "# " + pScen->skill(i) + ":";
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
    averageMin += (minTotal[d] - pSolver->maxTotalStaffNoPenalty_)
        / minPerShift.size();
    averageOpt += (optTotal[d] - pSolver->maxTotalStaffNoPenalty_)
        / minPerShift.size();
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
    averageMin += (minTotal[d] - pSolver->maxTotalStaffAvgWork_)
        / minPerShift.size();
    averageOpt += (optTotal[d] - pSolver->maxTotalStaffAvgWork_)
        / minPerShift.size();
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
    logStream << "# " + pScen->skill(i) + ":";
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
    logStream << "# " + pScen->skill(i) + ":";
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
