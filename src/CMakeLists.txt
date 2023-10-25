# add the subdirectories
add_subdirectory(data)
add_subdirectory(parsing)
add_subdirectory(solvers)
add_subdirectory(tools)

# build the dynamic library
if (DEFINED RCSPP_SHARED)
    MESSAGE("Linking rcspp as a shared library.")
    # shared library for tools and data
    add_library(nsbase SHARED Parameters.cpp ParseArguments.cpp)
    target_link_libraries(nsbase data parsing tools ${LIBS})
    # shared library for rcspp
    add_library(${NURSE_SCHEDULER_LIB} SHARED InitializeInstance.cpp)
    target_link_libraries(${NURSE_SCHEDULER_LIB} solvers nsbase rcspp ${LIBS})
else()
    add_library(${NURSE_SCHEDULER_LIB} SHARED
            Parameters.cpp ParseArguments.cpp InitializeInstance.cpp solvers/DeterministicSolver.cpp)
    target_link_libraries(${NURSE_SCHEDULER_LIB} data parsing solvers tools ${LIBS})
endif()

MESSAGE("MAIN LIB: ${NURSE_SCHEDULER_LIB}")


# Add custom target and methods to always relink the library when changing build type

# Dummy output to rebuild anything that depends on when changing build type.
# this will always be rebuilt when changing config.
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    add_custom_command(
            OUTPUT always_relink
            COMMAND rm -f ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/.debug.cmake && touch ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/.release.cmake
    )
else ()
    add_custom_command(
            OUTPUT ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/debug.lock
            COMMAND touch ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/debug.lock
    )
endif ()

# Custom target will always cause its dependencies to be evaluated and is
# run by default
SET(REL_LOCK ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/.release.lock)
SET(DBG_LOCK ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/.debug.lock)

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(LOCK ${REL_LOCK})
    SET(UNLOCK ${DBG_LOCK})
else ()
    SET(LOCK ${DBG_LOCK})
    SET(UNLOCK ${REL_LOCK})
endif ()

# add lock. The lock will trigger relink when removed
add_custom_command(
        OUTPUT ${LOCK}
        COMMAND bash -c "touch ${LOCK}"
)

# RELINK will always be rebuilt once trigger by ${LOCK}
add_custom_target(relink ALL DEPENDS ${LOCK})
SET(RM_LIB_NAME ${NURSE_SCHEDULER_LIB})
if (WIN_32)
    add_custom_command(
            TARGET relink
            COMMAND bash -c "rd ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${RM_LIB_NAME}*"
            VERBATIM COMMENT)
    # remove other lock once nurse scheduler has been rebuild
    add_custom_command(
            TARGET ${NURSE_SCHEDULER_LIB} POST_BUILD
            COMMAND bash -c "rd -f ${UNLOCK}"
    )
else()
    add_custom_command(
            TARGET relink
            COMMAND bash -c "rm -rf ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${RM_LIB_NAME}*"
            VERBATIM COMMENT)
    # remove other lock once nurse scheduler has been rebuild
    add_custom_command(
            TARGET ${NURSE_SCHEDULER_LIB} POST_BUILD
            COMMAND bash -c "rm -f ${UNLOCK}"
    )
endif()

# ensure that relink is always build before NURSE_SCHEDULER_LIB
add_dependencies(${NURSE_SCHEDULER_LIB} relink)