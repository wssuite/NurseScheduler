//
// Created by Flore Caye on 2022-08-11.
//

#include "data.h"



#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>

#include "jlcxx/jlcxx.hpp"
#include "data/Scenario.h"


JLCXX_MODULE init_scenario_module(jlcxx::Module& mod_scenario) {
    mod_scenario.add_type<PScenario>("PScenario")
            .constructor<>(true);
}
