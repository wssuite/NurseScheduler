//
// Copyright 2023 Flore Caye
//


#include <chrono>// NOLINT [build/c++11]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <regex>// NOLINT [build/c++11]
#include <string>
#include <memory>

#include <boost/algorithm/string.hpp>
#include "solvers/mp/sp/rcspp/resources/AlternativeShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "parseToolsUI.h"
#include "Tools.h"
#include "data/Nurse.h"
#include "data/Shift.h"

using std::string;
using std::vector;

struct Sched_Period parse_scheduling_period(const std::string &sch_per) {
  int nbWeeks, nbDays;
  std::tm firstDay;
  std::stringstream X(sch_per);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  std::regex reg(",");
  std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  string start_time;

  start_time = tokens[2];
  firstDay = *Tools::readDateFromStr(tokens[2]);
  nbDays = Tools::readDateFromStr(tokens[3])->tm_yday - firstDay.tm_yday + 1;
  nbWeeks = ceil(nbDays / 7);
  return {nbWeeks, nbDays, firstDay};
}

struct Skills_Parsed parse_skills(const std::string &s) {
  int nbSkills;
  std::map<string, int> skillToInt;
  std::vector<std::string> intToSkill;
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  int i = 0;
  while (line.length() > 0) {
    boost::trim(line);
    skillToInt[line] = i;
    intToSkill.push_back(line);
    std::getline(X, line);
    i++;
  }
  nbSkills = i;
  return {nbSkills, skillToInt, intToSkill};
}

struct Shifts_Parsed parse_shifts(std::string s) {
  int nbShifts;
  vector<double> hoursInShift;
  vector<string> intToShift;
  std::map<string, int> shiftToInt;
  vector<Tools::Time> start_time;
  // for the rest shift
  start_time.push_back(Tools::Time());
  vector<Tools::Time> end_time;
  end_time.push_back(Tools::Time());
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  vector<PAbstractShift> pShifts;

  hoursInShift.push_back(0.0);
  shiftToInt[REST_SHIFT] = 0;
  intToShift.push_back(REST_SHIFT);
  int i = 1;
  while (line.length() > 0) {
    // Parser la ligne de shift:
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    boost::trim(tokens[0]);
    shiftToInt[tokens[0]] = i;
    intToShift.push_back(tokens[0]);
    start_time.push_back(Tools::readHourFromStr(tokens[1]));
    end_time.push_back(Tools::readHourFromStr(tokens[2]));
    hoursInShift.push_back(end_time[i].diff(start_time[i]));

    std::getline(X, line);
    i++;
  }
  nbShifts = i;
  // Create the allowed successors
  vector2D<int> successors(nbShifts);
  for (int i = 0; i < nbShifts; i++) {
    for (int j = 0; j < nbShifts; j++) {
      if (j != i) {
        successors[i].push_back(j);
      }
    }
  }
  return {nbShifts, hoursInShift, intToShift, shiftToInt};
}

struct ShiftTypes_Parsed parse_shiftsType(const std::string &current_s,
                                          const std::map<std::string,
                                                         int> shiftToInt,
                                          const vector<string> &intToShift,
                                          const vector<double> &hoursInShift) {
  vector<PShift> pShifts;
  int nbShiftsType;
  vector<string> intToShiftType;
  std::map<string, int> shiftTypeToInt;
  vector2D<int> shiftTypeIDToShiftID;
  vector<int> shiftIDToShiftTypeID(shiftToInt.size(), -1);

  std::map<string, int> inputToShiftGenre;
  int nbShifts = shiftToInt.size();
  shiftTypeIDToShiftID = vector2D<int>();
  // Any shift can come after any other shift (possibly at the
  // cost of a penalty)
  vector<int> successorList;
  successorList.reserve(nbShifts);
  for (int i = 0; i < nbShifts; i++) {
    successorList.push_back(i);
  }

  // Add all shifts to the shift genre map:
  for (const auto &myPair : shiftToInt) {
    inputToShiftGenre[myPair.first] = 0;
  }

  std::stringstream X(current_s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);

  // IMPORTANT : INSERT REST SHIFT !!!!!!
  intToShiftType.push_back(REST_SHIFT);
  shiftTypeToInt[REST_SHIFT] = 0;
  // initialize the shift
  shiftIDToShiftTypeID[0] = 0;
  shiftTypeIDToShiftID.push_back(vector<int>(1, 0));

  int i = 1;
  // TODO(Flore): gerer le cas ou cette categorie est vide
  while (line.length() > 0) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    int nb_param = stoi(tokens[1]);
    boost::trim(tokens[0]);
    intToShiftType.push_back(tokens[0]);
    shiftTypeToInt[tokens[0]] = i;
    shiftTypeIDToShiftID.push_back(vector<int>());
    for (int j = 0; j < nb_param; j++) {
      boost::trim(tokens[2 + j]);
      shiftTypeIDToShiftID[i].push_back(shiftToInt.at(tokens[2 + j]));
      shiftIDToShiftTypeID[shiftToInt.at(tokens[2 + j])] = i;
    }

    std::getline(X, line);
    i++;
  }

  // TODO(Flore): double check qe ca fonctionne correctement
  int toAdd = i - 1;
  for (int s = 0; s < nbShifts; s++) {
    if (shiftIDToShiftTypeID[s] == -1) {
      toAdd++;
      intToShiftType.push_back(intToShift[s]);
      shiftTypeToInt[intToShift[s]] = toAdd;
      shiftTypeIDToShiftID.push_back(vector<int>(1, s));
      shiftIDToShiftTypeID[s] = toAdd;
    }
  }
  nbShiftsType = intToShiftType.size();
  for (int l = 0; l < nbShifts; l++) {
    pShifts.push_back(
        std::make_shared<Shift>(intToShift[l], l, shiftIDToShiftTypeID[l],
                                ceil(hoursInShift[l]), successorList));
  }

  // Add all shifts to the shift genre map:
  for (const auto &myPair : shiftTypeToInt) {
    inputToShiftGenre[myPair.first] = 1;
  }
  return {inputToShiftGenre, pShifts, nbShiftsType,
          intToShiftType, shiftTypeToInt, shiftTypeIDToShiftID,
          shiftIDToShiftTypeID};
}
std::map<std::string, vector<PBaseResource>> parse_contracts(
    std::string s, ShiftsFactory shiftFactory,
    std::map<string, int> inputToShiftGenre, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays) {
  // Removing the first 2 lines
  s.erase(0, s.find("\n") + 1);
  s.erase(0, s.find("\n") + 1);
  std::map<string, vector<PBaseResource>> ruleSets;
  std::regex reg("\\}\n\\{");
  std::sregex_token_iterator iter(s.begin(), s.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  int i = 0;
  for (const string &block : tokens) {
    struct Single_Contract contract = parse_single_contract(
        block, shiftFactory, inputToShiftGenre, shiftToInt, skillToInt,
        shiftTypeToInt, pWeights, nDays);
    ruleSets[contract.name] = contract.ressources;
    i++;
  }
  return ruleSets;
}
struct Single_Contract parse_single_contract(
    std::string contract, ShiftsFactory shiftFactory,
    std::map<string, int> inputToShiftGenre, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays) {
  std::vector<PBaseResource> cons;
  std::stringstream X(contract);
  bool parsing_constraints = false;
  int nb_constraints = 0;
  string name;
  for (std::string line; std::getline(X, line);) {
    if (line.find("contractName") != -1) {
      std::regex reg(",");
      std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
      std::sregex_token_iterator end;
      vector<string> tokens(iter, end);
      boost::trim(tokens[1]);
      name = tokens[1];
    } else if (line.find("constraints") != -1) {
      parsing_constraints = true;
    } else if (parsing_constraints) {
      vector<PBaseResource> new_cons = parse_contract_constraints(
          line, shiftFactory, inputToShiftGenre, shiftToInt, skillToInt,
          shiftTypeToInt, pWeights, nDays);
      cons.insert(cons.end(), std::begin(new_cons), std::end(new_cons));
      nb_constraints += new_cons.size();
    }
  }
  return {name, cons};
}
vector<PBaseResource> parse_contract_constraints(
    std::string contract, ShiftsFactory shiftFactory,
    std::map<string, int> inputShiftToInt, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays) {
  vector<PBaseResource> cons;
  std::regex reg(",");
  std::sregex_token_iterator iter(contract.begin(), contract.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  // TODO(Flore): max duration and average duration should not be hard coded
  int maxDuration = 8;
  int avgDuration = 6;
  for (string &s : tokens) {
    boost::trim(s);
  }
  if (tokens[0].find("NumberOfFreeDaysAfterShift") != -1) {
    int shiftTypeInt;
    PAbstractShift pAShift;

    double cost;
    if (tokens[3].find("Work") != -1) {
      pAShift = shiftFactory.pAnyWorkShift();
    } else if (tokens[3].find("Rest") != -1) {
      pAShift = shiftFactory.pAnyRestShift();
    } else if (inputShiftToInt.at(tokens[3]) == 0) {
      int shiftId = shiftToInt.at(tokens[3]);
      pAShift = shiftFactory.pAnyTypeShift(shiftFactory.pShift(shiftId)->type);
    } else if (inputShiftToInt.at(tokens[3]) == 1) {
      shiftTypeInt = shiftTypeToInt.at(tokens[3]);
      pAShift = shiftFactory.pAnyTypeShift(shiftTypeInt);
    }
    if (tokens[2] == "hard") {
      cost = XLARGE_SCORE;
    } else {
      cost = std::stod(tokens[2]);
    }
    cons.push_back(std::make_shared<SoftFreeDaysAfterShiftResource>(
        pAShift, std::stoi(tokens[1]), cost));

  } else if (tokens[0].find("CompleteWeekends") != -1) {
    double cost;
    if (tokens[1] == "hard") {
      cost = XLARGE_SCORE;
    } else {
      cost = std::stod(tokens[1]);
    }
    cons.push_back(std::make_shared<SoftIdentWeekendResource>(
        std::make_shared<ShiftWorkComparator>(), cost, SATURDAY, SUNDAY, "WE"));
  } else if (tokens[0].find("MinMaxConsecutiveShiftType") != -1) {
    int shiftTypeInt;
    PAbstractShift pAShift;
    if (tokens[5].find("Work") != -1) {
      pAShift = shiftFactory.pAnyWorkShift();
    } else if (tokens[5].find("Rest") != -1) {
      pAShift = shiftFactory.pAnyRestShift();
    } else if (inputShiftToInt.at(tokens[5]) == 0) {
      int shiftId = shiftToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftFactory.pShift(shiftId)->type);
    } else if (inputShiftToInt.at(tokens[5]) == 1) {
      shiftTypeInt = shiftTypeToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftTypeInt);
    }

    if (std::stoi(tokens[1]) > std::stoi(tokens[3])) {
      std::string errorMessage =
          std::string("Error: Min and Max are not consistent for Total number "
                      "of WE in 4 weeks ");
      throw std::runtime_error(errorMessage);
    }
    double costMin, costMax;
    if (tokens[2] == "hard") {
      costMin = XLARGE_SCORE;
    } else {
      costMin = std::stod(tokens[2]);
    }
    if (tokens[4] == "hard") {
      costMax = XLARGE_SCORE;
    } else {
      costMax = std::stod(tokens[4]);
    }
    cons.push_back(std::make_shared<SoftConsShiftResource>(
        std::stoi(tokens[1]), std::stoi(tokens[3]), costMin, costMax, pAShift,
        CONS_SHIFTS_COST, nDays, 0));
  } else if (tokens[0].find("MinMaxConsecutiveWeekends") != -1) {
    if (std::stoi(tokens[1]) > std::stoi(tokens[3])) {
      std::string errorMessage =
          std::string("Error: Min and Max are not consistent for Total number "
                      "of WE in 4 weeks ");
      throw std::runtime_error(errorMessage);
    }
    double costMin, costMax;
    if (tokens[2] == "hard") {
      costMin = XLARGE_SCORE;
    } else {
      costMin = std::stod(tokens[2]);
    }
    if (tokens[4] == "hard") {
      costMax = XLARGE_SCORE;
    } else {
      costMax = std::stod(tokens[4]);
    }
    cons.push_back(std::make_shared<SoftConsWeekendShiftResource>(
        std::stoi(tokens[1]), std::stoi(tokens[3]), costMin, costMax,
        shiftFactory.pAnyWorkShift(), nDays));

  } else if (tokens[0].find("MinMaxHoursInFourWeeks") != -1) {
    PAbstractShift pAShift;
    int shiftTypeInt;
    if (tokens[5].find("Work") != -1) {
      pAShift = shiftFactory.pAnyWorkShift();
    } else if (tokens[5].find("Rest") != -1) {
      pAShift = shiftFactory.pAnyRestShift();
    } else if (inputShiftToInt.at(tokens[5]) == 0) {
      int shiftId = shiftToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftFactory.pShift(shiftId)->type);
    } else if (inputShiftToInt.at(tokens[5]) == 1) {
      shiftTypeInt = shiftTypeToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftTypeInt);
    }
    if (std::stoi(tokens[1]) > std::stoi(tokens[3])) {
      std::string errorMessage =
          std::string("Error: Min and Max are not consistent for Total number "
                      "of WE in 4 weeks ");
      throw std::runtime_error(errorMessage);
    }
    double costMin, costMax;
    if (tokens[2] == "hard") {
      costMin = XLARGE_SCORE;
    } else {
      costMin = std::stod(tokens[2]);
    }
    if (tokens[4] == "hard") {
      costMax = XLARGE_SCORE;
    } else {
      costMax = std::stod(tokens[4]);
    }
    cons.push_back(std::make_shared<SoftTotalShiftDurationResource>(
        std::stoi(tokens[1]), std::stoi(tokens[3]), costMin, costMax, pAShift,
        nDays, maxDuration, avgDuration));
  } else if (tokens[0].find("MinMaxNumAssignmentsInFourWeeks") != -1) {
    PAbstractShift pAShift;
    int shiftTypeInt;
    if (tokens[5].find("Work") != -1) {
      pAShift = shiftFactory.pAnyWorkShift();
    } else if (tokens[5].find("Rest") != -1) {
      pAShift = shiftFactory.pAnyRestShift();
    } else if (inputShiftToInt.at(tokens[5]) == 0) {
      int shiftId = shiftToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftFactory.pShift(shiftId)->type);
    } else if (inputShiftToInt.at(tokens[5]) == 1) {
      shiftTypeInt = shiftTypeToInt.at(tokens[5]);
      pAShift = shiftFactory.pAnyTypeShift(shiftTypeInt);
    }
    if (std::stoi(tokens[1]) > std::stoi(tokens[3])) {
      std::string errorMessage =
          std::string("Error: Min and Max are not consistent for Total number "
                      "of WE in 4 weeks ");
      throw std::runtime_error(errorMessage);
    }
    double costMin, costMax;
    if (tokens[2] == "hard") {
      costMin = XLARGE_SCORE;
    } else {
      costMin = std::stod(tokens[2]);
    }
    if (tokens[4] == "hard") {
      costMax = XLARGE_SCORE;
    } else {
      costMax = std::stod(tokens[4]);
    }
    cons.push_back(std::make_shared<SoftTotalShiftDurationResource>(
        std::stoi(tokens[1]), std::stoi(tokens[3]), costMin, costMax, pAShift,
        nDays, 1, 1));
  } else if (tokens[0].find("IdentShiftTypesDuringWeekend") != -1) {
    if (tokens[1] == "hard") {
      cons.push_back(std::make_shared<SoftIdentWeekendResource>(
          std::make_shared<ShiftComparator>(), XLARGE_SCORE));
    } else {
      cons.push_back(std::make_shared<SoftIdentWeekendResource>(
          std::make_shared<ShiftComparator>(), std::stod(tokens[1])));
    }
  } else if (tokens[0].find("TotalWeekendsInFourWeeks") != -1) {
    if (std::stoi(tokens[1]) > std::stoi(tokens[3])) {
      std::string errorMessage =
          std::string("Error: Min and Max are not consistent for Total number "
                      "of WE in 4 weeks ");
      throw std::runtime_error(errorMessage);
    }
    if (tokens[2] == "hard") {
      cons.push_back(std::make_shared<HardTotalWeekendsResource>(
          std::stoi(tokens[1]), shiftFactory.pAnyWorkShift(), nDays));
    } else {
      cons.push_back(std::make_shared<SoftTotalWeekendsResource>(
          std::stoi(tokens[1]), std::stod(tokens[2]),
          shiftFactory.pAnyWorkShift(), nDays));
      if (tokens[4] == "hard") {
        cons.push_back(std::make_shared<HardTotalWeekendsResource>(
            std::stoi(tokens[3]), shiftFactory.pAnyWorkShift(), nDays));
      } else {
        cons.push_back(std::make_shared<SoftTotalWeekendsResource>(
            std::stoi(tokens[3]), std::stod(tokens[4]),
            shiftFactory.pAnyWorkShift(), nDays));
      }
    }
  } else if (tokens[0].find("unwantedPatterns") != -1) {
    double cost = XLARGE_SCORE;
    if (tokens[1].find("hard") == -1)
      cost = std::stod(tokens[1]);
    int nbPat = std::stoi(tokens[2]);
    for (int i = 0; i < nbPat; i++) {
      vector<PAbstractDay> days;
      vector<PAbstractShift> shifts;

      std::regex regPat(";");
      std::sregex_token_iterator iterPat(tokens[3 + i].begin(),
                                         tokens[3 + i].end(), regPat, -1);
      std::sregex_token_iterator endPat;
      vector<string> patItems(iterPat, endPat);

      std::regex regD("\\|");
      std::sregex_token_iterator iterD(patItems[0].begin(), patItems[0].end(),
                                       regD, -1);
      std::sregex_token_iterator endD;
      vector<string> daysStr(iterD, endD);

      std::sregex_token_iterator iterS(patItems[1].begin(), patItems[1].end(),
                                       regD, -1);
      std::sregex_token_iterator endS;
      vector<string> shiftsStr(iterS, endS);
      int patSiz = daysStr.size();
      for (int j = 0; j < patSiz; j++) {
        int shiftTypeInt;
        PAbstractShift pAShift;
        boost::trim(shiftsStr[j]);
        if (inputShiftToInt.at(shiftsStr[j]) == 0) {
          int shiftId = shiftToInt.at(shiftsStr[j]);
          pAShift =
              shiftFactory.pAnyTypeShift(shiftFactory.pShift(shiftId)->type);
        } else if (inputShiftToInt.at(shiftsStr[j]) == 1) {
          shiftTypeInt = shiftTypeToInt.at(shiftsStr[j]);
          pAShift = shiftFactory.pAnyTypeShift(shiftTypeInt);
        }
        days.push_back(std::make_shared<Day>(daysOfWeekByName.at(daysStr[j])));
        shifts.push_back(pAShift);
      }
      cons.push_back(std::make_shared<SoftForbiddenPatternResource>(
          Pattern(shifts, days), cost));
    }

  } else if (tokens[0].find("unwantedSkills") != -1) {
    boost::trim(tokens[2]);
    double cost = XLARGE_SCORE;
    if (tokens[1].find("hard") == -1)
      cost = std::stod(tokens[1]);
    int nbSkills = std::stoi(tokens[2]);
    for (int i = 0; i < nbSkills; i++) {
      int skill = skillToInt.at(tokens[3 + i]);
      cons.push_back(std::make_shared<AlternativeShiftResource>(skill, cost));
    }
  }
  return cons;
}
std::map<std::string, std::vector<PBaseResource>> parse_group_contracts(
    const std::string &s,
    const std::map<std::string, std::vector<PBaseResource>> &ruleSets) {
  std::map<std::string, std::vector<PBaseResource>> new_ruleSets = ruleSets;
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  int id = 0;
  while (line.size() > 0) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    boost::trim(tokens[0]);
    vector<PBaseResource> cons;
    int nbSet = std::stoi(tokens[1]);
    for (int j = 0; j < nbSet; j++) {
      boost::trim(tokens[j + 2]);
      cons.insert(cons.end(), std::begin(ruleSets.at(tokens[j + 2])),
                  std::end(ruleSets.at(tokens[j + 2])));
    }

    new_ruleSets[tokens[0]] = cons;
    id++;
    getline(X, line);
  }
  return new_ruleSets;
}

struct Nurses_Parsed parse_nurses(const std::string &s,
                                  const std::map<string,
                                                 PContract> &pContractsByName,
                                  const vector<PContract> &pContracts,
                                  const std::map<std::string,
                                                 std::vector<PBaseResource>>
                                  &ruleSets,
                                  const vector<PNurse> &pNurses,
                                  vector<int> skills,
                                  vector<int> shifts) {
  // Second, creer les objects contrats et les objects nurse, les ajouter au
  // scenario.
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  std::map<string, PContract> new_pContractsByName = pContractsByName;
  vector<PContract> new_pContracts = pContracts;
  vector<PNurse> new_pNurses = pNurses;
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    int id = stoi(tokens[0]);
    boost::trim(tokens[1]);
    boost::trim(tokens[3]);
    string name = tokens[1];
    std::cout << name << std::endl;
    int nbContracts = std::stoi(tokens[2]);
    int nbContractGroups = std::stoi(tokens[nbContracts + 3]);
    vector<PBaseResource> cons;
    string contract_name = "";
    for (int i = 3; i < nbContracts + 3; i++) {
      std::cout << tokens[i] << std::endl;
      contract_name.append("_");
      boost::trim(tokens[i]);
      auto consAdd = ruleSets.at(tokens[i]);
      cons.insert(cons.end(), std::begin(consAdd), std::end(consAdd));
      contract_name.append(tokens[i]);
    }
    for (int i = nbContracts + 4; i < nbContractGroups + nbContracts + 4; i++) {
      contract_name.append("_");
      boost::trim(tokens[i]);
      std::cout << tokens[i] << std::endl;
      auto consAdd = ruleSets.at(tokens[i]);
      cons.insert(cons.end(), std::begin(consAdd), std::end(consAdd));
      contract_name.append(tokens[i]);
    }
    PContract pc = std::make_shared<Contract>(id, contract_name, cons);
    new_pContractsByName[contract_name] = pc;
    new_pContracts.push_back(pc);

    new_pNurses.push_back(std::make_shared<Nurse>(
        id, name, shifts.size(), skills.size(), skills, shifts, pc));
    getline(X, line);
  }
  return {new_pContractsByName, new_pContracts, new_pNurses};
}

PDemand parse_demand(const std::string &s,
                     std::map<string, int> shiftToInt,
                     std::map<string, int> skillToInt, int nbDays,
                     const tm &startDay) {
  string name = "standard";
  vector<int> days_vec(nbDays);
  vector<vector<int>> days_skill_vec;
  vector3D<int> minDemand;
  vector3D<int> optDemand;
  Tools::initVector3D(&minDemand, nbDays, shiftToInt.size(), skillToInt.size(),
                      0);
  Tools::initVector3D(&optDemand, nbDays, shiftToInt.size(), skillToInt.size(),
                      0);
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);

    boost::trim(tokens[0]);
    boost::trim(tokens[1]);
    boost::trim(tokens[2]);
    const tm tmDem(*Tools::readDateFromStr(tokens[0]));
    int dayId(tmDem.tm_yday - startDay.tm_yday);
    int shiftId = shiftToInt.at(tokens[1]);
    int skillId = skillToInt.at(tokens[2]);
    int min = stoi(tokens[3]);
    int opt = stoi(tokens[5]);

    minDemand[dayId][shiftId][skillId] = min;
    optDemand[dayId][shiftId][skillId] = opt;
    std::getline(X, line);
  }
  PDemand pDemand = std::make_shared<Demand>(nbDays,
                                             0,
                                             shiftToInt.size(),
                                             skillToInt.size(),
                                             name,
                                             minDemand,
                                             optDemand);
  pDemand->preprocessMinDemand();
  pDemand->preprocessOptDemand();
  return pDemand;
}
void parse_preferences(const string &s, const PScenario &pScenario) {
  // For each nurse, maps the day to the set of shifts that he/she wants to have
  // off
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  PPreferences pPref = std::make_shared<Preferences>(
      pScenario->nNurses(), pScenario->nDays(), pScenario->nShifts());
  std::cout << "nbDays : " << pScenario->nDays() << std::endl;
  ShiftsFactory shiftFactory(pScenario->pShifts());
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);

    double pref;
    boost::trim(tokens[0]);
    boost::trim(tokens[3]);
    boost::trim(tokens[4]);
    const tm tmRequest(*Tools::readDateFromStr(tokens[0]));
    int dayId(tmRequest.tm_yday - pScenario->startDate().tm_yday);

    if (tokens[4] == "hard") {
      pref = LARGE_SCORE;
    } else if (std::stod(tokens[4]) < 5) {
      pref = WEAK;
    } else if (std::stod(tokens[4]) <= 8) {
      pref = MODERATE;
    } else {
      pref = STRONG;
    }

    if (tokens[2] == "OFF") {
      if (tokens[3] == "Any") {
        pPref->addShiftOff(std::stoi(tokens[1]), dayId,
                           pScenario->shiftsFactory().pAnyWorkShift(), pref);
      } else {
        pPref->addShiftOff(std::stoi(tokens[1]), dayId,
                           pScenario->pShift(tokens[3]), pref);
      }
    } else if (tokens[2] == "ON") {
      if (tokens[3] == "Any") {
        pPref->addShiftOn(std::stoi(tokens[1]), dayId, pScenario->pShift(1),
                          pref);
      } else {
        pPref->addShiftOn(std::stoi(tokens[1]), dayId,
                          pScenario->pShift(tokens[3]), pref);
      }
    }
    getline(X, line);
  }

  pScenario->linkWithPreferences(pPref);
}
struct HistoryPeriod parse_history_period(const string &s) {
  int nbDays_h;
  std::tm firstDay;
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  std::regex reg(",");
  std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  string start_time;

  firstDay = *Tools::readDateFromStr(tokens[0]);
  nbDays_h = Tools::readDateFromStr(tokens[1])->tm_yday - firstDay.tm_yday + 1;
  return {firstDay, nbDays_h};
}

std::vector<State> parse_history(const string &s,
                                 const PScenario &pScenario,
                                 struct HistoryPeriod history) {
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  ShiftsFactory shiftFactory(pScenario->pShifts());
  std::vector<State> initialState;
  std::vector<std::vector<PShift>> hist;
  const PShift &pNoneShift =
      pScenario->shiftsFactory().pNoneShift()->pIncludedShifts().front();
  for (int n = 0; n < pScenario->nNurses(); n++) {
    State nurseState(
        0,
        0,
        0,
        0,
        0,
        0,
        pNoneShift);
    initialState.push_back(nurseState);
    hist.push_back(std::vector<PShift>(history.nbDays_history, pNoneShift));
  }
  std::vector<std::vector<PShift>>();
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    boost::trim(tokens[0]);
    boost::trim(tokens[1]);
    boost::trim(tokens[2]);
    const tm tmRequest(*Tools::readDateFromStr(tokens[0]));
    int dayId = tmRequest.tm_yday - history.firstDay_history.tm_yday;
    int nurse_id = stoi(tokens[1]);
    int shift_id = pScenario->shiftType(tokens[2]);
    PShift shift = pScenario->shiftsFactory().pShift(shift_id);
    hist.at(nurse_id).at(dayId) = shift;
    getline(X, line);
  }
  for (int n = 0; n < pScenario->nNurses(); n++) {
    for (int j = 0; j < history.nbDays_history; j++) {
      initialState[n].addDayToState(initialState[n], hist[n][j]);
    }
  }
  return initialState;
}
