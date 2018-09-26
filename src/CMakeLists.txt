project(nursescheduler)

set(SOURCE_FILES
        BcpModeler.cpp
#        CbcModeler.cpp
        Demand.cpp
        DemandGenerator.cpp
        DeterministicSolver.cpp
        GlobalStats.cpp
        Greedy.cpp
        InitializeSolver.cpp
        InputPaths.cpp
        MasterProblem.cpp
        MyTools.cpp
        Nurse.cpp
        ReadWrite.cpp
        Roster.cpp
        RotationPricer.cpp
        Scenario.cpp
        Solver.cpp
        StochasticSolver.cpp
        SubProblem.cpp
        TreeManager.cpp)

#ifdef USE_CBC
#include "CbcModeler.h"
#endif

# build the dynamic library
add_library(${NURSE_SCHEDULER_LIB} SHARED ${SOURCE_FILES})
target_link_libraries(${NURSE_SCHEDULER_LIB} ${LIBS})
