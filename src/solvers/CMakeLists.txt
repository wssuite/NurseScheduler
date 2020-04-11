add_subdirectory(mp)

set(LIB_NAME solvers)
set(SOURCE_FILES
        DeterministicSolver.cpp
        Greedy.cpp
        InitializeSolver.cpp
        MasterProblem.cpp
        Solver.cpp
        StochasticSolver.cpp
        )

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} mp ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
