cmake_minimum_required(VERSION 3.5)
project(nursescheduler)

# set my own definitions
SET(CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
SET(NURSE_SCHEDULER_LIB nursescheduler)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# read user definitions
include(CMakeDefinitionsLists.txt)

MESSAGE("Type: ${CMAKE_BUILD_TYPE}")
set(BCP_DIR ${BCPDIROPT})
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set(BCP_DIR ${BCPDIRDBG})
        add_definitions(-fstack-protector)
        SET(CMAKE_VERBOSE_MAKEFILE on)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
        #    add_definitions(-D_GLIBCXX_DEBUG)
        # add asan flags (help to detect memory runtime errors)
        if(${ASAN})
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
            set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
        endif()
        add_definitions(-DDBG)
endif()
add_definitions(-pedantic -Wall -Wno-sign-compare)

# add boost
if(DEFINED BOOST_DIR)
    set(Boost_INCLUDE_DIR ${BOOST_DIR})
else()
    FIND_PACKAGE( Boost 1.48 COMPONENTS graph REQUIRED )
endif()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
MESSAGE("Boost directory: ${Boost_INCLUDE_DIR}")

# add BCP
MESSAGE("BCP directory: ${BCP_DIR}")
INCLUDE_DIRECTORIES(${BCP_DIR}/include/coin)
LINK_DIRECTORIES(${BCP_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SET(LIBS Bcp Osi OsiClp ClpSolver Clp CoinUtils bz2 z lapack blas m pthread)

# build main library
INCLUDE_DIRECTORIES(src)
add_subdirectory(src)

# exec
add_executable(staticscheduler src/DeterministicMain.cpp)
target_link_libraries(staticscheduler ${NURSE_SCHEDULER_LIB})

add_executable(dynamicscheduler src/DynamicMain.cpp)
target_link_libraries(dynamicscheduler ${NURSE_SCHEDULER_LIB})

add_executable(pricer src/PricerMain.cpp)
target_link_libraries(pricer ${NURSE_SCHEDULER_LIB})