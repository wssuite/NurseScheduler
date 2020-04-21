//
// Created by antoine legrain on 2020-04-19.
//

#include "RosterSP.h"

using std::string;
using std::vector;
using std::map;

// Constructors and destructor
RosterSP::RosterSP() {}

RosterSP::RosterSP(Scenario* scenario, int nbDays, const Contract* contract, vector<State>* pInitState):
    SubProblem(scenario, nbDays, contract,  pInitState) {

  nLabels_ = 5;
}

RosterSP::~RosterSP(){}