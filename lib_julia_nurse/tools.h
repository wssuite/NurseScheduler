//
// Created by Flore Caye on 2022-08-11.
//

#ifndef JULIA_WRAP_TOOLS_H
#define JULIA_WRAP_TOOLS_H


#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "tools/InputPaths.h"
#include "tools/Tools.h"

JLCXX_MODULE init_Tools_module(jlcxx::Module& mod_tools);
JLCXX_MODULE init_inputpaths_module(jlcxx::Module& mod_inputpaths);
#endif //JULIA_WRAP_TOOLS_H
