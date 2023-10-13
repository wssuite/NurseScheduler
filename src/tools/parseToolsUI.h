//
// Copyright 2023 Flore Caye
//

#ifndef SRC_TOOLS_PARSETOOLSUI_H_
#define SRC_TOOLS_PARSETOOLSUI_H_

#include <map>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <streambuf>
#include <string>
#include <ctime>
#include <vector>

#include <boost/variant.hpp>
#include "solvers/mp/sp/rcspp/resources/AlternativeShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/ForbiddenPatternResource.h"
#include "solvers/mp/sp/rcspp/resources/FreeDaysAfterShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/IdentWeekendResource.h"
#include "solvers/mp/sp/rcspp/resources/PreferenceResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

#include "Tools.h"
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
  vector<PBaseResource> ressources;
};

struct Nurses_Parsed {
  std::map<string, PContract> pContractsByName;
  vector<PContract> pContracts;
  vector<PNurse> pNurses;
};

struct HistoryPeriod {
  tm firstDay_history;
  int nbDays_history;
};

bool string_starts_with(std::string s, std::string pattern);

struct Sched_Period parse_scheduling_period(const std::string &sch_per);

struct Skills_Parsed parse_skills(const std::string &s);

struct Shifts_Parsed parse_shifts(std::string s);

struct ShiftTypes_Parsed
parse_shiftsType(const std::string &current_s,
                 std::map<std::string, int> shiftToInt,
                 const vector<string> &intToShift,
                 const vector<double> &hoursInShift);

std::map<std::string, std::vector<PBaseResource>> parse_contracts(
    std::string s, ShiftsFactory shiftFactory,
    std::map<string, int> inputToShiftGenre, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays);

struct Single_Contract parse_single_contract(
    std::string contract, ShiftsFactory shiftFactory,
    std::map<string, int> inputToShiftGenre, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays);

std::vector<PBaseResource> parse_contract_constraints(
    std::string contract, ShiftsFactory shiftFactory,
    std::map<string, int> inputShiftToInt, std::map<string, int> shiftToInt,
    std::map<string, int> skillToInt, std::map<string, int> shiftTypeToInt,
    const PWeights &pWeights, int nDays);

// void make_contracts_from_groups(const std::string& s,
// std::map<std::string,std::vector<PBaseResource>>& ruleSets,
// vector<PContract>& pContracts, std::map<string, PContract>&
// pContractsByName); Contract craftContract(int id,const std::string&
// name,std::map<std::string,std::vector<PBaseResource>> ruleSet); void
// setGeneralRules(std::map<std::string,std::vector<PBaseResource>>
// generalRules, Scenario& pScenario);
struct Nurses_Parsed parse_nurses(const std::string &s,
                                  const std::map<string,
                                                 PContract> &pContractsByName,
                                  const vector<PContract> &pContracts,
                                  const std::map<std::string,
                                                 std::vector<PBaseResource>>
                                  &ruleSets,
                                  const vector<PNurse> &pNurses,
                                  vector<int> skills,
                                  vector<int> shifts);

std::map<std::string, std::vector<PBaseResource>> parse_group_contracts(
    const std::string &s,
    const std::map<std::string, std::vector<PBaseResource>> &ruleSets);

// PConstContract craftFromRules(int id,int nbApplicableRuleSets,
// std::vector<int> applicableRuleSets, std::vector<RuleSet> ruleSets, Scenario&
// scenario);
PDemand parse_demand(const std::string &s, std::map<string, int> shiftToInt,
                     std::map<string, int> skillToInt, int nbDays,
                     const tm &firstDay);

void parse_preferences(const std::string &s, const PScenario &pScenario);
std::vector<State> parse_history(const string &s, const PScenario &pScenario);
struct HistoryPeriod parse_history_period(const string &s);
std::vector<State> parse_history(const string &s,
                                 const PScenario &pScenario,
                                 struct HistoryPeriod history);
#endif  // SRC_TOOLS_PARSETOOLSUI_H_
