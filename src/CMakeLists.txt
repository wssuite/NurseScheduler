# Add custom target and methods to always relink the library
# Custom target will always cause its dependencies to be evaluated and is
# run by default
add_custom_target(relink ALL DEPENDS always_relink)

SET(RM_LIB_NAME ${NURSE_SCHEDULER_LIB})

# RELINK will always be rebuilt because it depends on always_rebuild
if (WIN_32)
    add_custom_command(
            TARGET relink
            COMMAND bash -c "rd ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${RM_LIB_NAME}*"
            VERBATIM COMMENT)
else()
    add_custom_command(
            TARGET relink
            COMMAND bash -c "rm -rf ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${RM_LIB_NAME}*"
            VERBATIM COMMENT)
endif()

# Dummy output which is never actually produced. Anything that depends on.
# this will always be rebuilt.
add_custom_command(
        OUTPUT always_relink
        COMMAND cmake -E echo
)

# add the subdirectories
add_subdirectory(data)
add_subdirectory(solvers)
add_subdirectory(tools)

# build the dynamic library
if (DEFINED RCSPP_SHARED)
    MESSAGE("Linking rcspp as a shared library.")
    # shared library for tools and data
    add_library(nsbase SHARED Parameters.cpp ParseArguments.cpp)
    target_link_libraries(nsbase data tools ${LIBS})
    # shared library for rcspp
    add_library(${NURSE_SCHEDULER_LIB} SHARED ReadWrite.cpp InitializeInstance.cpp)
    target_link_libraries(${NURSE_SCHEDULER_LIB} solvers nsbase rcspp ${LIBS})
    # ensure that relink is always build before nsbase and NURSE_SCHEDULER_LIB
    add_dependencies(nsbase relink)
else()
    add_library(${NURSE_SCHEDULER_LIB} SHARED
            Parameters.cpp ParseArguments.cpp ReadWrite.cpp InitializeInstance.cpp solvers/DeterministicSolver.cpp)
    target_link_libraries(${NURSE_SCHEDULER_LIB} data solvers tools ${LIBS})
    # ensure that relink is always build before NURSE_SCHEDULER_LIB
    add_dependencies(${NURSE_SCHEDULER_LIB} relink)
endif()
MESSAGE("MAIN LIB: ${NURSE_SCHEDULER_LIB}")