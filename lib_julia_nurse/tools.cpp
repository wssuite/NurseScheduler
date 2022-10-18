//
// Created by Flore Caye on 2022-08-11.
//

#include "tools.h"

#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "tools/InputPaths.h"
#include "tools/Tools.h"
#include "tools/GlobalStats.h"

JLCXX_MODULE init_Tools_module(jlcxx::Module& mod_tools){
    using namespace Tools;
    mod_tools.method("initializeRandomGenerator", [](int i){return initializeRandomGenerator(i);});
    mod_tools.add_type<LogOutput>("LogOutput");
    mod_tools.add_type<struct GlobalStats>("GlobalStats")
            .method("statStream", &GlobalStats::lnsStatsToString)
            .method("statStream", &GlobalStats::toString);
}
JLCXX_MODULE init_inputpaths_module(jlcxx::Module& mod_inputpaths) {
   //std::cout << "Lib included!! Julia, here we come!!! " <<  std::endl;
    mod_inputpaths.add_type<InputPaths>("InputPaths")
            .constructor<const std::string ,
                             const std::string ,
                             int,
                             std::vector<int>,
                             const std::string ,
                             const std::string ,
                             const std::string,
                             int ,
                             int ,
                             int ,
                             const std::string,
                             int,
                             const std::string,
                             int ,
                             int>()
            .constructor<const std::string ,
                    const std::string ,
                    std::vector<int>,
                    const std::string ,
                    const std::string ,
                    const std::string,
                    int ,
                    int ,
                    int ,
                    const std::string,
                    int,
                    const std::string,
                    int ,
                    int>()
                    .method("solutionPaths", [](InputPaths &inputPaths)-> const std::string& {return inputPaths.solutionPath();});
}