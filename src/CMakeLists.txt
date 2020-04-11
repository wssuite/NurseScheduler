add_subdirectory(data)
add_subdirectory(solvers)
add_subdirectory(tools)

set(SOURCE_FILES
        DeterministicMain_test.cpp
        )

# build the dynamic library
add_library(${NURSE_SCHEDULER_LIB} SHARED ${SOURCE_FILES})
target_link_libraries(${NURSE_SCHEDULER_LIB} data solvers tools ${LIBS})