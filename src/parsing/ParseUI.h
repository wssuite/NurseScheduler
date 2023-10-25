//
// Copyright 2023 Flore Caye
//

#ifndef SRC_PARSING_PARSEUI_H_
#define SRC_PARSING_PARSEUI_H_

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <streambuf>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/resources/UnwantedShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

#include "tools/Tools.h"
#include "data/Scenario.h"

struct Sched_Period {
  int nWeeks;
  int nDays;
  std::tm startDay;
};

struct Skills_Parsed {
  int nbSkills;
  std::map<string, int> skillToInt;
  std::vector<std::string> intToSkill;
};

struct Shifts_Parsed {
  int nbShifts;
  vector<double> hoursInShift;
  vector<string> intToShift;
  std::map<string, int> shiftToInt;
};
struct ShiftTypes_Parsed {
  std::map<string, int> inputToShiftGenre;
  vector<PShift> pShifts;
  int nbShiftsType;
  vector<string> intToShiftType;
  std::map<string, int> shiftTypeToInt;
  vector2D<int> shiftTypeIDToShiftID;
  vector<int> shiftIDToShiftTypeID;
};

struct Single_Contract {
  string name;
  vector<PBaseResource> resources;
  std::map<int, double> unwantedSkills;
};

typedef std::map<std::string, struct Single_Contract> RulesSet;

struct Nurses_Parsed {
  vector<PContract> pContracts;
  vector<PNurse> pNurses;
};

struct HistoryPeriod {
  tm firstDay_history;
  int nbDays_history;
};

// Read the scenario file and store the content in a Scenario instance
PScenario readScenarioUI(const std::string &fileName, const string &name);

// Write the output custom file from values in the scenario and the solution
// instances
void writeUI(const string &path, const PScenario &pScenario);

bool string_starts_with(std::string s, std::string pattern);

struct Sched_Period parse_scheduling_period(const std::string &sch_per);

struct Skills_Parsed parse_skills(const std::string &s);

struct Shifts_Parsed parse_shifts(std::string s);

struct ShiftTypes_Parsed
parse_shiftsType(const std::string &current_s,
                 std::map<std::string, int> shiftToInt,
                 const vector<string> &intToShift,
                 const vector<double> &hoursInShift);

std::map<string, PAbstractShift>
parse_shiftGroups(const std::string &current_s,
                 const std::map<std::string, int> &shiftToInt,
                 const std::map<string, int> &shiftTypeToInt,
                 const ShiftsFactory &shiftFactory);

RulesSet parse_contracts(
    std::string s,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputShiftToInt,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay);

struct Single_Contract parse_single_contract(
    std::string contract,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputShiftToInt,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay);

vector<PBaseResource> parse_contract_constraints(
    std::string contract,
    const ShiftsFactory &shiftFactory,
    const std::map<string, int> &inputShiftToInt,
    const std::map<string, int> &shiftToInt,
    const std::map<string, int> &skillToInt,
    const std::map<string, int> &shiftTypeToInt,
    const std::map<string, PAbstractShift> &pShiftGroups,
    int nDays, const tm &startDay);

PAbstractShift findAShift(
        std::string name,
        const ShiftsFactory &shiftFactory,
        const std::map<string, int> &inputShiftToInt,
        const std::map<string, int> &shiftToInt,
        const std::map<string, int> &shiftTypeToInt,
        const std::map<string, PAbstractShift> &pShiftGroups,
        bool getShiftTypeForShift);

static double getCost(const std::string &c, double hardCost = HARD_COST);

static void checkBounds(int lb, int ub, const std::string &consName);

// void make_contracts_from_groups(const std::string& s,
// std::map<std::string,std::vector<PBaseResource>>& ruleSets,
// vector<PContract>& pContracts, std::map<string, PContract>&
// pContractsByName); Contract craftContract(int id,const std::string&
// name,std::map<std::string,std::vector<PBaseResource>> ruleSet); void
// setGeneralRules(std::map<std::string,std::vector<PBaseResource>>
// generalRules, Scenario& pScenario);
struct Nurses_Parsed parse_nurses(
        const std::string &s,
        const RulesSet &ruleSets,
        const vector<int> &skills,
        const vector<int> &shifts);

RulesSet parse_group_contracts(const std::string &s, const RulesSet &ruleSets);

// PConstContract craftFromRules(int id,int nbApplicableRuleSets,
// std::vector<int> applicableRuleSets, std::vector<RuleSet> ruleSets, Scenario&
// scenario);
std::vector<PDemand> parse_demands(
        const std::string &s, std::map<string, int> shiftToInt,
        std::map<string, int> skillToInt, int nbDays,
        const tm &firstDay);

void parse_preferences(const std::string &s, const PScenario &pScenario);
std::vector<State> parse_history(const string &s, const PScenario &pScenario);
struct HistoryPeriod parse_history_period(const string &s);
std::vector<State> parse_history(const string &s,
                                 const PScenario &pScenario,
                                 std::tm startDate);
#endif  // SRC_PARSING_PARSEUI_H_
