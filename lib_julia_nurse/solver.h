//
// Created by Flore Caye on 2022-08-11.
//

#ifndef JULIA_WRAP_SOLVER_H
#define JULIA_WRAP_SOLVER_H
#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"

#include "tools/InputPaths.h"
#include "data/Scenario.h"
#include "solvers/InitializeSolver.h"

JLCXX_MODULE init_initsolver_module(jlcxx::Module& mod_initsolver);
JLCXX_MODULE init_detersolver_module(jlcxx::Module& mod_detersolver);
#endif //JULIA_WRAP_SOLVER_H
