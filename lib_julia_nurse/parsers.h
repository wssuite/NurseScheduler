//
// Created by Flore Caye on 2022-08-11.
//

#ifndef JULIA_WRAP_PARSERS_H
#define JULIA_WRAP_PARSERS_H
#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "ParseArguments.h"

JLCXX_MODULE init_parsearg_module(jlcxx::Module& mod_parsearg);

#endif //JULIA_WRAP_PARSERS_H
