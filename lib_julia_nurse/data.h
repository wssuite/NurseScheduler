//
// Created by Flore Caye on 2022-08-11.
//

#ifndef JULIA_WRAP_DATA_H
#define JULIA_WRAP_DATA_H

#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>

#include "jlcxx/jlcxx.hpp"
JLCXX_MODULE init_scenario_module(jlcxx::Module& mod_scenario);

#endif //JULIA_WRAP_DATA_H
