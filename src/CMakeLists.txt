project(nursescheduler)

set(SOURCE_FILES
        BcpModeler.cpp
        Demand.cpp
        DeterministicMain_test.cpp
        DeterministicSolver.cpp
        GlobalStats.cpp
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
#        StochasticSolver.cpp
        SubProblem.cpp
        TreeManager.cpp)

# build the dynamic library
add_library(${NURSE_SCHEDULER_LIB} SHARED ${SOURCE_FILES})
target_link_libraries(${NURSE_SCHEDULER_LIB} ${LIBS})