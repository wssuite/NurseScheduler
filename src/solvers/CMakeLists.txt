add_subdirectory(mp)

set(LIB_NAME solvers)

set(SOURCE_FILES
        DeterministicSolver.cpp
        InitializeSolver.cpp
        Solver.cpp
        StochasticSolver.cpp
        DynamicWeights.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} mp tools data ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
