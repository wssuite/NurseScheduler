//
// Created by Flore Caye on 2022-08-11.
//

#include "solver.h"

#include<julia.h>
#include <type_traits>
#include <string>
#include <memory>
#include <iostream>
#include <functional>
#include <vector>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/array.hpp"

#include "tools/InputPaths.h"
#include "data/Scenario.h"
#include "solvers/InitializeSolver.h"
#include "solvers/DeterministicSolver.h"
#include "solvers/Solver.h"


JLCXX_MODULE init_initsolver_module(jlcxx::Module& mod_initsolver) {
    mod_initsolver.method("initializeMultipleWeeksINRC2",[](const InputPaths &inputPaths, const std::string& logPath="") -> PScenario {return initializeMultipleWeeksINRC2(inputPaths, "");});
    mod_initsolver.method("initializeScenarioINRC2", &initializeScenarioINRC2);
}
JLCXX_MODULE init_detersolver_module(jlcxx::Module& mod_detersolver) {
    mod_detersolver.add_type<DeterministicSolver>("DeterministicSolver")
            .constructor<const PScenario&,
                const InputPaths>()
                .method("solve", &DeterministicSolver::solve)
                .method("status", [](DeterministicSolver pSolver) -> bool {return pSolver.status();})
                .method("displaySolutionMultipleWeeks", &Solver::displaySolutionMultipleWeeks)
                .method("getGlobalStats", &DeterministicSolver::getGlobalStat);
}

