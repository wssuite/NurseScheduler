//
// Copyright 2023 Flore Caye
//

#include <chrono>// NOLINT [build/c++11]
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>// NOLINT [build/c++11]
#include <string>
#include <utility>

#include "solvers/mp/sp/rcspp/resources/UnwantedShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"
#include "ParseUI.h"
#include "tools/Tools.h"
#include "data/Nurse.h"
#include "data/Shift.h"

using std::string;
using std::vector;
using std::map;


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

#define B_CHECK_LINE(line, label) else if (line == label) { try {  // NOLINT
#define E_CHECK_LINE(label) } catch(...) { \
  std::cerr << "Exception thrown while processing ui input block " \
            << label << std::endl; \
  rethrow_exception(std::current_exception()); } }

PScenario readScenarioUI(const std::string &fileName, const string &name) {
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
  std::map<string, PAbstractShift> pShiftGroups;
  // "input" -> 0 if input is a shift, 1 if input is a shift type,
  // 2 if input if a shift group
  map<string, int> inputToShiftGenre;
  RulesSet ruleSets;
  PDemand pDemand;
  PPreferences pPref;
  vector<State> initialState;
  PScenario pScenario;
  ShiftsFactory shiftsFactory;
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
  vector<PDemand> pDemands;

  for (const string &s : tokens) {
    std::stringstream X(s);
    std::string line;
    std::getline(X, line);
    Tools::trim(&line);
    if (line == "HEADERS") {
      header = s;
    } B_CHECK_LINE(line, "SCHEDULING_PERIOD")
        const struct Sched_Period t = parse_scheduling_period(s);
        nWeeks = t.nWeeks;
        nDays = t.nDays;
        startDay = t.startDay;
        header.append(s);
        Day::setFirstDayOfWeek(
                (DayOfWeek) (startDay.tm_wday >= 0 ? startDay.tm_wday - 1 : 6));
    E_CHECK_LINE("SCHEDULING_PERIOD")
    B_CHECK_LINE(line, "SKILLS")
        const struct Skills_Parsed t = parse_skills(s);
        nSkills = t.nbSkills;
        skillToInt = t.skillToInt;
        intToSkill = t.intToSkill;
    E_CHECK_LINE("SKILLS")
    B_CHECK_LINE(line, "SHIFTS")
        const struct Shifts_Parsed t = parse_shifts(s);
        nShifts = t.nbShifts;
        shiftDurations = t.hoursInShift;
        intToShift = t.intToShift;
        shiftToInt = t.shiftToInt;
    E_CHECK_LINE("SHIFTS")
    B_CHECK_LINE(line, "SHIFT_TYPES")
        const struct ShiftTypes_Parsed
                t = parse_shiftsType(s, shiftToInt, intToShift, shiftDurations);
        inputToShiftGenre = t.inputToShiftGenre;
        pShifts = t.pShifts;
        nShiftsType = t.nbShiftsType;
        intToShiftType = t.intToShiftType;
        shiftTypeToInt = t.shiftTypeToInt;
        shiftTypeIDToShiftID = t.shiftTypeIDToShiftID;
        shiftIDToShiftTypeID = t.shiftIDToShiftTypeID;
        shiftsFactory = ShiftsFactory(pShifts, intToShiftType);
    E_CHECK_LINE("SHIFT_TYPES")
    B_CHECK_LINE(line, "SHIFT_GROUPS")
        pShiftGroups = parse_shiftGroups(
                s, shiftToInt, shiftTypeToInt, shiftsFactory);
    E_CHECK_LINE("SHIFT_GROUPS")
    B_CHECK_LINE(line, "CONTRACTS")
        ruleSets = parse_contracts(s,
                                   shiftsFactory,
                                   inputToShiftGenre,
                                   shiftToInt,
                                   skillToInt,
                                   shiftTypeToInt,
                                   pShiftGroups,
                                   nDays,
                                   startDay);
        contracts_processed = ruleSets.size();
    E_CHECK_LINE("CONTRACTS")
    B_CHECK_LINE(line, "CONTRACT_GROUPS")
        if (!contracts_processed)
          Tools::throwError("Cannot parse CONTRACT_GROUPS block, as CONTRACTS "
                            "block has not be defined before.");
        ruleSets = parse_group_contracts(s, ruleSets);
    E_CHECK_LINE("CONTRACT_GROUPS")
    B_CHECK_LINE(line, "EMPLOYEES")
        // As I go thought the employee,
        // I craft them contracts based on the rule Sets as needed
        vector<int> allSkills, allShifts;
        for (int i = 0; i < nSkills; i++) {
          allSkills.push_back(i);
        }
        for (int i = 0; i < nShifts; i++) {
          allShifts.push_back(i);
        }
        struct Nurses_Parsed nurses = parse_nurses(
                s,
                ruleSets,
                allSkills,
                allShifts);

        pScenario = std::make_shared<Scenario>(name,
                                               nWeeks,
                                               intToSkill,
                                               skillToInt,
                                               pShifts,
                                               intToShift,
                                               nurses.pContracts,
                                               nurses.pNurses,
                                               header);
    E_CHECK_LINE("EMPLOYEES")
    B_CHECK_LINE(line, "HOSPITAL_DEMAND")
        auto tm = pScenario->startDate(), now = std::tm();
        // check if tm is equals to either now if not already set
        // or the new date
        if ((tm.tm_year != now.tm_year || tm.tm_yday != now.tm_yday) &&
            (tm.tm_year != startDay.tm_year || tm.tm_yday != startDay.tm_yday))
          Tools::throwError("Demand starts at different days.");
        pScenario->setStartDate(startDay);
        auto pDs = parse_demands(s, shiftToInt, skillToInt, nDays, startDay);
        pDemands.insert(pDemands.end(), pDs.begin(), pDs.end());
        pScenario->linkWithDemand(pDemands);
    E_CHECK_LINE("HOSPITAL_DEMAND")
    B_CHECK_LINE(line, "PREFERENCES")
        // parse s - line
        parse_preferences(s, pScenario);
    E_CHECK_LINE("PREFERENCES")
    B_CHECK_LINE(line, "HISTORY")
        initialState = parse_history(s, pScenario, startDay);
    E_CHECK_LINE("HISTORY")
  }

  // Initialize the history of every nurse to empty state
  // Add a fictitious shift just for the initial state if no history given
  if (initialState.empty()) {
    for (int n = 0; n < nNurses; n++)
      initialState.push_back(State::noneState(pScenario->shiftsFactory()));
  }
  pScenario->setInitialState(initialState);
  return pScenario;
}

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
    Tools::trim(&line);
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
  start_time.emplace_back();
  vector<Tools::Time> end_time;
  end_time.emplace_back();
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  vector<PAbstractShift> pShifts;

  hoursInShift.push_back(0.0);
  shiftToInt[REST_SHIFT] = 0;
  intToShift.emplace_back(REST_SHIFT);
  int i = 1;
  while (line.length() > 0) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    Tools::trim(&tokens[0]);
    shiftToInt[tokens[0]] = i;
    intToShift.push_back(tokens[0]);
    start_time.push_back(Tools::readHourFromStr(tokens[1]));
    end_time.push_back(Tools::readHourFromStr(tokens[2]));
    double h = end_time[i].diff(start_time[i]);
    if (h < 1e-3) {
      std::cout << "Warning: shift " << tokens[0] << " has a 0 ("
                << tokens[1] << " -> " << tokens[2] << ") duration, "
                << "thus it is replaced by a default 1 hour duration."
                << std::endl;
      h = 1;
    }
    hoursInShift.push_back(h);

    std::getline(X, line);
    i++;
  }
  return {i, hoursInShift, intToShift, shiftToInt};
}

struct ShiftTypes_Parsed parse_shiftsType(
        const std::string &current_s,
        const std::map<std::string, int> shiftToInt,
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
    Tools::trim(&tokens[0]);
    intToShiftType.push_back(tokens[0]);
    shiftTypeToInt[tokens[0]] = i;
    shiftTypeIDToShiftID.emplace_back();
    for (int j = 0; j < nb_param; j++) {
      Tools::trim(&tokens[2 + j]);
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
      shiftTypeIDToShiftID.emplace_back(1, s);
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

std::map<string, PAbstractShift>
parse_shiftGroups(const std::string &current_s,
                  const std::map<std::string, int> &shiftToInt,
                  const std::map<string, int> &shiftTypeToInt,
                  const ShiftsFactory &shiftFactory) {
  std::map<string, PAbstractShift> shiftGroups;
  std::stringstream X(current_s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);

  while (line.length() > 0) {
    vector<string> tokens = Tools::tokenize<string>(line, ',');
    Tools::trim(&tokens[0]);
    if (tokens[0] == "Rest" || tokens[0] == "Work") {
      std::cout << "Shift group " << tokens[0]
                << " is a reserved name -> group ignored." << std::endl;
      std::getline(X, line);
      continue;
    }
    int nShifts = stoi(tokens[1]), nShiftTypes = stoi(tokens[2+nShifts]);
    std::vector<PAbstractShift> allPAShifts;
    for (int i=0; i < nShifts; i++)
      allPAShifts.push_back(shiftFactory.pShift(shiftToInt.at(tokens[2+i])));
    for (int i=0; i < nShiftTypes; i++)
      allPAShifts.push_back(shiftFactory.pAnyTypeShift(
              shiftTypeToInt.at(tokens[3+nShifts])));
    // add rest shift, if no shift present
    if (allPAShifts.empty()) allPAShifts = {shiftFactory.pAnyRestShift()};
    shiftGroups[tokens[0]] = std::make_shared<Shifts>(allPAShifts);
    std::getline(X, line);
  }

  return shiftGroups;
}

RulesSet parse_contracts(
    std::string s,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputToShiftGenre,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay) {
  // Removing the first 2 lines
  s.erase(0, s.find("\n") + 1);
  s.erase(0, s.find("\n") + 1);
  RulesSet ruleSets;
  std::regex reg("\\}\n\\{");
  std::sregex_token_iterator iter(s.begin(), s.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  int i = 0;
  for (const string &block : tokens) {
    struct Single_Contract contract = parse_single_contract(
        block, shiftFactory, inputToShiftGenre, shiftToInt, skillToInt,
        shiftTypeToInt, pShiftGroups, nDays, startDay);
    ruleSets[contract.name] = contract;
    i++;
  }
  return ruleSets;
}

struct Single_Contract parse_single_contract(
    std::string contract,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputToShiftGenre,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay) {
  std::vector<PBaseResource> cons;
  std::stringstream X(contract);
  std::regex reg(",");
  bool parsing_constraints = false;
  string name;
  std::map<int, double> unwantedSkills;
  for (std::string line; std::getline(X, line);) {
    if (line.find("contractName") != -1) {
      std::regex reg(",");
      std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
      std::sregex_token_iterator end;
      vector<string> tokens(iter, end);
      Tools::trim(&tokens[1]);
      name = tokens[1];
    } else if (line.find("constraints") != -1) {
      parsing_constraints = true;
    } else if (parsing_constraints) {
      // try to parse UnwantedSkills
      std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
      std::sregex_token_iterator end;
      vector<string> tokens(iter, end);
      if (tokens[0].find("UnwantedSkills") != -1) {
        Tools::trim(&tokens[2]);
        double cost = getCost(tokens[1]);
        int nSkills = std::stoi(tokens[2]);
        for (int i = 0; i < nSkills; i++) {
          int sk = skillToInt.at(tokens[3 + i]);
          if (unwantedSkills[sk] > 0) {
            std::cout << "WARNING: skills " << tokens[3 + i]
                      << " is marked as unwanted at least twice" << std::endl;
            if (unwantedSkills[sk] < cost)
              unwantedSkills[sk] = cost;
          } else {
            unwantedSkills[sk] = cost;
          }
        }
      } else {
        // otherwise try to parse the rest
        vector<PBaseResource> new_cons = parse_contract_constraints(
                line, shiftFactory, inputToShiftGenre, shiftToInt, skillToInt,
                shiftTypeToInt, pShiftGroups, nDays, startDay);
        cons.insert(cons.end(), std::begin(new_cons), std::end(new_cons));
      }
    }
  }
  return {name, cons, unwantedSkills};
}

vector<PBaseResource> parse_contract_constraints(
    std::string contract,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputShiftToInt,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay) {
  vector<PBaseResource> cons;
  std::regex reg(",");
  std::sregex_token_iterator iter(contract.begin(), contract.end(), reg, -1);
  std::sregex_token_iterator end;
  vector<string> tokens(iter, end);
  double ratio4W = nDays / (4.0 * 7);
  for (string &s : tokens)
    Tools::trim(&s);
  if (tokens[0].find("NumberOfFreeDaysAfterShift") != -1) {
    PAbstractShift pAShift = findAShift(
            tokens[3],
            shiftFactory,
            inputShiftToInt,
            shiftToInt,
            shiftTypeToInt,
            pShiftGroups,
            true);
    double cost = getCost(tokens[2], HARD_COST);
    cons.push_back(std::make_shared<SoftFreeDaysAfterShiftResource>(
            pAShift, std::stoi(tokens[1]), cost));
  } else if (tokens[0].find("CompleteWeekends") != -1) {
    double cost = getCost(tokens[1], INFEAS_COST);
    cons.push_back(std::make_shared<SoftIdentWeekendResource>(
            std::make_shared<ShiftWorkComparator>(), cost));
  } else if (tokens[0].find("MinMaxConsecutiveShiftType") != -1) {
    PAbstractShift pAShift = findAShift(
            tokens[5],
            shiftFactory,
            inputShiftToInt,
            shiftToInt,
            shiftTypeToInt,
            pShiftGroups,
            true);
    double costMin = getCost(tokens[2]),
            costMax = getCost(tokens[4]);
    int lb = std::stoi(tokens[1]), ub = std::stoi(tokens[3]);
    checkBounds(lb, ub, "MinMaxConsecutiveShiftType");
    CostType cT = pAShift->isWork() ?
            (pAShift->isAnyWork() ? CONS_WORK_COST : CONS_SHIFTS_COST) :
            CONS_REST_COST;
    if (isHardCost(costMax))
      cons.push_back(std::make_shared<HardConsShiftResource>(
              isSoftCost(costMin) ? 0 : lb, ub,
              pAShift,
              cT, nDays, 0));
    if (isSoftCost(costMin) || isSoftCost(costMax))
      cons.push_back(std::make_shared<SoftConsShiftResource>(
              lb, ub, costMin, costMax, pAShift,
              cT, nDays, 0));
  } else if (tokens[0].find("MinMaxConsecutiveWeekends") != -1) {
    double costMin = getCost(tokens[2], INFEAS_COST),
            costMax = getCost(tokens[4], INFEAS_COST);
    int lb = std::stoi(tokens[1]), ub = std::stoi(tokens[3]);
    checkBounds(lb, ub, "MinMaxConsecutiveWeekends");
    cons.push_back(std::make_shared<SoftConsWeekendShiftResource>(
            lb, ub, costMin, costMax,
            shiftFactory.pAnyWorkShift(), nDays));

  } else if (tokens[0].find("MinMaxHoursInFourWeeks") != -1) {
    PAbstractShift pAShift = findAShift(
            tokens[5],
            shiftFactory,
            inputShiftToInt,
            shiftToInt,
            shiftTypeToInt,
            pShiftGroups,
            false);
    int maxDuration = 0;
    for (const auto &pS : pAShift->pIncludedShifts())
      if (pS->duration > maxDuration) maxDuration = pS->duration;
    double costMin = getCost(tokens[2]),
            costMax = getCost(tokens[4]);
    int lb = static_cast<int>(round(std::stoi(tokens[1]) * ratio4W)),
            ub = static_cast<int>(round(std::stoi(tokens[3]) * ratio4W));
    checkBounds(lb, ub, "MinMaxHoursInFourWeeks");
    if (isHardCost(costMax))
      cons.push_back(std::make_shared<HardTotalShiftDurationResource>(
              isSoftCost(costMin) ? 0 : lb, ub,
              pAShift, nDays, false, maxDuration));
    if (isSoftCost(costMin) || isSoftCost(costMax))
      cons.push_back(std::make_shared<SoftTotalShiftDurationResource>(
              lb, ub, costMin, costMax, pAShift, nDays, false, maxDuration));
  } else if (tokens[0].find("MinMaxNumAssignmentsInFourWeeks") != -1) {
    PAbstractShift pAShift = findAShift(
            tokens[5],
            shiftFactory,
            inputShiftToInt,
            shiftToInt,
            shiftTypeToInt,
            pShiftGroups,
            false);
    double costMin = getCost(tokens[2]),
            costMax = getCost(tokens[4]);
    int lb = static_cast<int>(round(std::stoi(tokens[1]) * ratio4W)),
            ub = static_cast<int>(round(std::stoi(tokens[3]) * ratio4W));
    checkBounds(lb, ub, "MinMaxNumAssignmentsInFourWeeks");
    if (isHardCost(costMax))
      cons.push_back(std::make_shared<HardTotalShiftDurationResource>(
              isSoftCost(costMin) ? 0 : lb, ub, pAShift, nDays, true));
    if (isSoftCost(costMin) || isSoftCost(costMax))
      cons.push_back(std::make_shared<SoftTotalShiftDurationResource>(
              lb, ub, costMin, costMax, pAShift, nDays, true));
  } else if (tokens[0].find("IdentShiftTypesDuringWeekend") != -1) {
    double cost = getCost(tokens[1], INFEAS_COST);
    cons.push_back(std::make_shared<SoftIdentWeekendResource>(
            std::make_shared<ShiftComparator>(), cost));
  } else if (tokens[0].find("TotalWeekendsInFourWeeks") != -1) {
    PAbstractShift pAS = findAShift(
            tokens[5],
            shiftFactory,
            inputShiftToInt,
            shiftToInt,
            shiftTypeToInt,
            pShiftGroups,
            true);
    double costMin = getCost(tokens[2]),
            costMax = getCost(tokens[4]);
    int lb = static_cast<int>(round(std::stoi(tokens[1]) * ratio4W)),
            ub = static_cast<int>(round(std::stoi(tokens[3]) * ratio4W));
    checkBounds(lb, ub, "TotalWeekendsInFourWeeks");
    if (isHardCost(costMax))
      cons.push_back(std::make_shared<HardTotalWeekendsResource>(
              isSoftCost(costMin) ? 0 : lb, ub, pAS, nDays));
    if (isSoftCost(costMin) || isSoftCost(costMax))
      cons.push_back(std::make_shared<SoftTotalWeekendsResource>(
              lb, ub, costMin, costMax, pAS, nDays));
  } else if (tokens[0].find("UnwantedPatterns") != -1) {
    double cost = getCost(tokens[1]);
    int nbPat = std::stoi(tokens[2]);
    vector<PAbstractDay> patDays;
    vector<PAbstractShift> patShifts;
    for (int i = 0; i < nbPat; i++) {
      vector<PAbstractShift> shifts;

      std::regex regPat(";");
      std::sregex_token_iterator iterPat(tokens[3 + i].begin(),
                                         tokens[3 + i].end(), regPat, -1);
      std::sregex_token_iterator endPat;
      vector<string> patItems(iterPat, endPat);

      // create one abstract days for all of them
      std::regex regD("\\|");
      std::sregex_token_iterator iterD(
              patItems[0].begin(), patItems[0].end(), regD, -1);
      std::sregex_token_iterator endD;
      vector<string> daysStr(iterD, endD);
      // empty string or * fits all day
      if (daysStr.size() == 1 && (daysStr[0].empty() || daysStr[0] == "*")) {
        patDays.push_back(std::make_shared<AnyDay>());
      } else {
        vector<PAbstractDay> days;
        for (const auto &d : daysStr) {
          auto it = daysOfWeekByName.find(d);
          if (it != daysOfWeekByName.end()) {
            days.push_back(std::make_shared<WeekDay>(it->second));
          } else {
            const tm tmRequest(*Tools::readDateFromStr(d));
            int dayId(tmRequest.tm_yday - startDay.tm_yday);
            days.push_back(std::make_shared<Day>(dayId));
          }
        }
        auto pADay = std::make_shared<Days>(days);
        patDays.push_back(pADay);
      }

      // create one abstract shift for all shift types
      std::sregex_token_iterator iterS(
              patItems[1].begin(), patItems[1].end(), regD, -1);
      std::sregex_token_iterator endS;
      vector<string> shiftsStr(iterS, endS);
      for (auto s : shiftsStr) {
        Tools::trim(&s);
        PAbstractShift pAShift = findAShift(
                s,
                shiftFactory,
                inputShiftToInt,
                shiftToInt,
                shiftTypeToInt,
                pShiftGroups,
                true);
        shifts.push_back(pAShift);
      }
      patShifts.push_back(std::make_shared<Shifts>(shifts));
    }

    // create the pattern and add the constraint
    auto pat = Pattern(patShifts, patDays);
    if (isHardCost(cost)) {
      cons.push_back(std::make_shared<HardForbiddenPatternResource>(pat));
    } else {
      cons.push_back(std::make_shared<SoftForbiddenPatternResource>(pat, cost));
    }

  } else if (tokens[0].find("UnwantedShift") != -1) {
    Tools::trim(&tokens[2]);
    double cost = getCost(tokens[1]);
    vector<PShift> pAltShifts;
    if (inputShiftToInt.at(tokens[2]) == 1) {
      int t = shiftTypeToInt.at(tokens[2]);
      pAltShifts = shiftFactory.pAnyTypeShift(t)->pIncludedShifts();
    } else {
      int s = shiftToInt.at(tokens[2]);
      pAltShifts = {shiftFactory.pShift(s)};
    }
    cons.push_back(std::make_shared<UnwantedShiftResource>(
            shiftToInt.size(), pAltShifts, cost));
  } else if (tokens[0].find('}') == string::npos) {
    Tools::throwException("Constraint %s is nor recognized: %s",
                          tokens[0].c_str(), contract.c_str());
  }
  return cons;
}

PAbstractShift findAShift(
        std::string name,
        const ShiftsFactory &shiftFactory,
        const std::map<string, int> &inputShiftToInt,
        const std::map<string, int> &shiftToInt,
        const std::map<string, int> &shiftTypeToInt,
        const std::map<string, PAbstractShift> &pShiftGroups,
        bool getShiftTypeForShift) {
  Tools::trim(&name);
  if (name.find("Work") != -1) {
    return shiftFactory.pAnyWorkShift();
  } else if (name.find("Rest") != -1) {
    return shiftFactory.pAnyRestShift();
  } else if (pShiftGroups.find(name) != pShiftGroups.end()) {
    auto pAS = pShiftGroups.at(name);
    if (!getShiftTypeForShift) return pAS;
    // transform pShift into their associated shift types
    vector<PAbstractShift> newShiftsGroup;
    for (const auto &pAS2 : pAS->pIncludedShifts()) {
      auto pS = std::dynamic_pointer_cast<Shift>(pAS2);
      if (pS) {
        const auto &pAS3 = shiftFactory.pAnyTypeShift(pS->type);
        if (pAS3->pIncludedShifts().size() > 1) {
          std::cout << "WARNING: cannot use shift (" << name << ") in patterns,"
                    << " but only shift type. Will use the associated shift "
                    << "type instead that includes more than one shift: "
                    << pAS->name << std::endl;
        }
        newShiftsGroup.push_back(pAS3);
      } else {
        newShiftsGroup.push_back(pAS2);
      }
    }
    return std::make_shared<Shifts>(newShiftsGroup);
  } else if (inputShiftToInt.at(name) == 0) {
    int shiftId = shiftToInt.at(name);
    const auto &pS = shiftFactory.pShift(shiftId);
    if (!getShiftTypeForShift) return pS;
    auto pAS = shiftFactory.pAnyTypeShift(pS->type);
    if (pAS->pIncludedShifts().size() > 1) {
      std::cout << "WARNING: cannot use shift (" << name << ") in patterns, "
                << "but only shift type. Will use the associated shift "
                << "type instead that includes more than one shift: "
                << pAS->name << std::endl;
    }
    return pAS;
  } else {
    int shiftTypeInt = shiftTypeToInt.at(name);
    return shiftFactory.pAnyTypeShift(shiftTypeInt);
  }
}

double getCost(const std::string &c, double hardCost) {
  size_t i = Tools::toLowerCase(c).find("hard");
  return i != string::npos ? hardCost : std::stod(c);
}

void checkBounds(int lb, int ub, const std::string &consName) {
  if (lb > ub) {
    Tools::throwException("Min (%d) and Max (%d) are not consistent for %s",
                          lb, ub, consName.c_str());
  }
}

RulesSet parse_group_contracts(const std::string &s, const RulesSet &ruleSets) {
  RulesSet new_ruleSets = ruleSets;
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
    Tools::trim(&tokens[0]);
    struct Single_Contract gCont = {tokens[0], {}, {}};
    int nbSet = std::stoi(tokens[1]);
    for (int j = 0; j < nbSet; j++) {
      Tools::trim(&tokens[j + 2]);
      auto cont = ruleSets.at(tokens[j + 2]);
      gCont.resources.insert(gCont.resources.end(),
              std::begin(cont.resources), std::end(cont.resources));
      for (const auto &p : cont.unwantedSkills) {
        if (gCont.unwantedSkills[p.first] < p.second)
          gCont.unwantedSkills[p.first] = p.second;
      }
    }
    new_ruleSets[gCont.name] = gCont;
    id++;
    getline(X, line);
  }
  return new_ruleSets;
}

struct Nurses_Parsed parse_nurses(
        const std::string &s,
        const RulesSet &ruleSets,
        const vector<int> &skills,
        const vector<int> &shifts) {
  // Second, creer les objects contrats et les objects nurse, les ajouter au
  // scenario.
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  std::map<string, PContract> new_pContractsByName;
  vector<PContract> new_pContracts;
  vector<PNurse> new_pNurses;
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    int id = stoi(tokens[0]);
    Tools::trim(&tokens[1]);
    Tools::trim(&tokens[3]);
    string name = tokens[1];
    int nbContracts = std::stoi(tokens[2]);
    int nbContractGroups = std::stoi(tokens[nbContracts + 3]);
    vector<PBaseResource> cons;
    std::map<int, double> unwantedSkills;
    string contract_name;
    for (int i = 3; i < nbContracts + 3; i++) {
      contract_name.append("_");
      Tools::trim(&tokens[i]);
      auto consAdd = ruleSets.at(tokens[i]);
      cons.insert(cons.end(), std::begin(consAdd.resources),
                  std::end(consAdd.resources));
      for (const auto &p : consAdd.unwantedSkills)
        if (unwantedSkills[p.first] < p.second)
          unwantedSkills[p.first] = p.second;
      contract_name.append(tokens[i]);
    }
    for (int i = nbContracts + 4; i < nbContractGroups + nbContracts + 4; i++) {
      contract_name.append("_");
      Tools::trim(&tokens[i]);
      auto consAdd = ruleSets.at(tokens[i]);
      cons.insert(cons.end(), std::begin(consAdd.resources),
                  std::end(consAdd.resources));
      for (const auto &p : consAdd.unwantedSkills)
        if (unwantedSkills[p.first] < p.second)
          unwantedSkills[p.first] = p.second;
      contract_name.append(tokens[i]);
    }
    // create a new contract if needed
    PContract pc;
    auto it = new_pContractsByName.find(contract_name);
    if (it != new_pContractsByName.end()) {
      pc = it->second;
    } else {
      vector<int> altSkills;
      for (const auto &p : unwantedSkills) altSkills.push_back(p.first);
      std::stable_sort(altSkills.begin(), altSkills.end());
      vector<double> altSkillCosts;
      for (int sk : altSkills) altSkillCosts.push_back(unwantedSkills[sk]);
      pc = std::make_shared<Contract>(
              new_pContracts.size(), contract_name,
              cons, altSkills, altSkillCosts);
      new_pContractsByName[contract_name] = pc;
      new_pContracts.push_back(pc);
    }

    vector<int> nurseSkills;
    for (int sk : skills)
      if (unwantedSkills.find(sk) == unwantedSkills.end())
        nurseSkills.push_back(sk);

    new_pNurses.push_back(std::make_shared<Nurse>(
        id, name, skills.size(), nurseSkills, shifts, pc));
    getline(X, line);
  }
  return {new_pContracts, new_pNurses};
}

std::vector<PDemand> parse_demands(
        const std::string &s,
        std::map<string, int> shiftToInt,
        std::map<string, int> skillToInt,
        int nbDays,
        const tm &startDay) {
  vector<int> days_vec(nbDays);
  vector<vector<int>> days_skill_vec;
  int maxLb, maxUb;
  double maxLbCost, maxUbCost;
  vector3D<int> lbDemand, ubDemand;
  vector3D<double> lbCosts, ubCosts;
  Tools::initVector3D(
          &lbDemand, nbDays, shiftToInt.size(), skillToInt.size(), 0);
  Tools::initVector3D(
          &ubDemand, nbDays, shiftToInt.size(), skillToInt.size(), 0);
  Tools::initVector3D(
          &lbCosts, nbDays, shiftToInt.size(), skillToInt.size(), .0);
  Tools::initVector3D(
          &ubCosts, nbDays, shiftToInt.size(), skillToInt.size(), .0);
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);

    Tools::trim(&tokens[0]);
    Tools::trim(&tokens[1]);
    Tools::trim(&tokens[2]);
    const tm tmDem(*Tools::readDateFromStr(tokens[0]));
    int dayId(tmDem.tm_yday - startDay.tm_yday);
    int shiftId = shiftToInt.at(tokens[1]);
    int skillId = skillToInt.at(tokens[2]);
    int lb = stoi(tokens[3]), ub = stoi(tokens[5]);
    double lbCost = getCost(tokens[4]), ubCost = getCost(tokens[6]);
    lbDemand[dayId][shiftId][skillId] = lb;
    lbCosts[dayId][shiftId][skillId] = lbCost;
    ubDemand[dayId][shiftId][skillId] = ub;
    ubCosts[dayId][shiftId][skillId] = ubCost;
    if (lb > maxLb) maxLb = lb;
    if (ub > maxUb) maxUb = ub;
    if (lbCost > maxLbCost) maxLbCost = lbCost;
    if (ubCost > maxUbCost) maxUbCost = ubCost;
    std::getline(X, line);
  }

  vector<PDemand> pDemands;
  if (maxLb > 0 && maxLbCost >= 1e-3) {
    char buff[20];
    snprintf(buff, sizeof(buff), "LB (max: %d, %.0f)", maxLb, maxLbCost);
    auto pD = std::make_shared<Demand>(
            nbDays, 0, shiftToInt.size(), skillToInt.size(), buff,
            lbDemand, D_GE, maxLbCost);
    pD->setCosts(lbCosts);
    pD->preprocess();
    pDemands.push_back(pD);
  }

  if (maxUbCost >= 1e-3) {
    char buff[20];
    snprintf(buff, sizeof(buff), "UB (max: %d, %.0f)", maxUb, maxUbCost);
    auto pD = std::make_shared<Demand>(
            nbDays, 0, shiftToInt.size(), skillToInt.size(), buff,
            ubDemand, D_LE, maxUbCost);
    pD->setCosts(ubCosts);
    pD->preprocess();
    pDemands.push_back(pD);
  }

  // only return active demand
  return pDemands;
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
    Tools::trim(&tokens[0]);
    Tools::trim(&tokens[3]);
    Tools::trim(&tokens[4]);
    const tm tmRequest(*Tools::readDateFromStr(tokens[0]));
    int dayId(tmRequest.tm_yday - pScenario->startDate().tm_yday);
    double pref = getCost(tokens[4]);
    PAbstractShift pAS;
    if (tokens[3] == "Any" || tokens[3] == "Work") {
      pAS = pScenario->shiftsFactory().pAnyWorkShift();
    } else {
      try {
        int st = pScenario->shiftType(tokens[3]);
        pAS = pScenario->pAnyTypeShift(st);
      } catch (...) {
        pAS = pScenario->pShift(tokens[3]);
      }
    }
    if (tokens[2] == "OFF") {
      pPref->addShiftOff(std::stoi(tokens[1]), dayId, pAS, pref);
    } else if (tokens[2] == "ON") {
      pPref->addShiftOn(std::stoi(tokens[1]), dayId, pAS, pref);
    }
    getline(X, line);
  }

  pScenario->linkWithPreferences(pPref);
}

std::vector<State> parse_history(const string &s,
                                 const PScenario &pScenario,
                                 std::tm startDate) {
  std::stringstream X(s);
  std::string line;
  std::getline(X, line);
  std::getline(X, line);
  // store for each day the associated shift read
  int firstDay = startDate.tm_yday;
  std::vector<std::map<int, PShift>> hist(pScenario->nNurses());
  while (!line.empty()) {
    std::regex reg(",");
    std::sregex_token_iterator iter(line.begin(), line.end(), reg, -1);
    std::sregex_token_iterator end;
    vector<string> tokens(iter, end);
    Tools::trim(&tokens[0]);
    Tools::trim(&tokens[1]);
    Tools::trim(&tokens[2]);
    std::tm tmRequest(*Tools::readDateFromStr(tokens[0]));
    int dayId = tmRequest.tm_yday;
    if (tmRequest.tm_year < startDate.tm_year) dayId -= 365;
    int nurse_id = stoi(tokens[1]);
    int shift_id = pScenario->shiftType(tokens[2]);
    PShift shift = pScenario->shiftsFactory().pShift(shift_id);
    hist.at(nurse_id)[dayId] = shift;
    if (dayId <= firstDay) firstDay = dayId;
    getline(X, line);
  }

  // build the initial state of each nurse
  std::vector<State> initialState;
  initialState.reserve(pScenario->nNurses());
  const PShift &pRestShift = pScenario->pRestShift();
  int histLength = startDate.tm_yday - firstDay + 1;
  for (int n = 0; n < pScenario->nNurses(); n++) {
    // if no history is provided for the nurse, add None shift
    if (hist.at(n).empty()) {
      initialState.push_back(State::noneState(pScenario->shiftsFactory()));
    } else {
      // otherwise, built initial state
      State nurseState(pRestShift, -histLength);
      for (int d = firstDay; d < startDate.tm_yday; ++d) {
        // if no shift provided, suppose a rest shift
        auto it = hist.at(n).find(d);
        nurseState.addDayToState(
                nurseState, it != hist.at(n).end() ? it->second : pRestShift);
      }
      nurseState.resetTotal();
      initialState.push_back(nurseState);
    }
  }
  return initialState;
}
