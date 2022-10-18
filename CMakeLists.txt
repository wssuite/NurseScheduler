cmake_minimum_required(VERSION 3.10)
project(nursescheduler)

# set my own definitions
SET(CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
SET(NURSE_SCHEDULER_LIB nursescheduler)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# read user definitions
include(CMakeDefinitionsLists.txt)

MESSAGE("Type: ${CMAKE_BUILD_TYPE} ${CMAKE_SYSTEM_PROCESSOR}")
set(BCP_DIR ${BCPDIROPT})
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set(BCP_DIR ${BCPDIRDBG})
    endif()
    #add_definitions(-fstack-protector)
    SET(CMAKE_VERBOSE_MAKEFILE on)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-3")
    #    add_definitions(-D_GLIBCXX_DEBUG)
    # add asan flags (help to detect memory runtime errors)
    if(${ASAN})
        MESSAGE("Add address sanitizer.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fsanitize=leak")
        set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fsanitize=leak")
    endif()
    add_definitions(-DDBG)
endif()
add_definitions(-pedantic -Wno-deprecated -Wno-sign-compare)

# add boost
if(DEFINED BOOST_DIR)
    message("Boost ${BOOST_DIR}")
    set(Boost_INCLUDE_DIRS ${BOOST_DIR})
else()
    FIND_PACKAGE(Boost 1.48 COMPONENTS graph REQUIRED )
endif()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
MESSAGE("Boost directory: ${Boost_INCLUDE_DIRS}")

# set default BLAS and LAPACK libs
if(NOT DEFINED BLAS_LIBS)
    set(BLAS_LIBS lapack blas)
endif()

# add BCP
LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SET(LIBS ${BLAS_LIBS} bz2 z m pthread)

# find coin or depedencies with pkg-config
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
    pkg_search_module(COINUTILS coinutils)
    pkg_search_module(OSI osi)
    pkg_search_module(OSI-CLP osi-clp)
    pkg_search_module(CLP clp)
    pkg_search_module(OSI-CBC osi-cbc)
    pkg_search_module(CGL cgl)
    pkg_search_module(CBC cbc)
    if (CBC_FOUND)
        add_definitions(-DUSE_CBC)
    endif()
    pkg_search_module(BCP bcp)
    if (BCP_FOUND)
        Message("BCP modified version from the submodules is used.")
        add_definitions(-DMBCP)
    endif()
    pkg_search_module(OSI-GUROBI osi-gurobi)
    if (OSI-GUROBI_FOUND)
        Message("Add gurobi library: -l${GUROBI_LIB}")
        include_directories(${GUROBI_INCLUDE_DIR})
        add_definitions(-DUSE_GUROBI)
        SET(LIBS ${LIBS} ${GUROBI_LIB})
    endif()
    include_directories(${COINUTILS_INCLUDE_DIRS} ${OSI_INCLUDE_DIRS}
            ${OSI-CLP_INCLUDE_DIRS} ${CLP_INCLUDE_DIRS} ${CBC_INCLUDE_DIRS}
            ${BCP_INCLUDE_DIRS} ${OSI-GUROBI_INCLUDE_DIRS})
    link_directories(${COINUTILS_LIBRARY_DIRS} ${OSI_LIBRARY_DIRS}
            ${OSI-CLP_LIBRARY_DIRS} ${CLP_LIBRARY_DIRS}
            ${OSI-CBC_LIBRARY_DIRS} ${CBC_LIBRARY_DIRS}
            ${BCP_LIBRARY_DIRS} ${OSI-GUROBI_LIBRARY_DIRS})
    SET(LIBS ${LIBS} ${COINUTILS_LIBRARIES} ${OSI_LIBRARIES}
            ${OSI-CLP_LIBRARIES} ${CLP_LIBRARIES} ${OSI-CBC_LIBRARIES}
            ${CBC_LIBRARIES} ${BCP_LIBRARIES} ${OSI-GUROBI_LIBRARIES})
endif()

# if pkg_config is not available or the pkg_config files are not available
if(NOT PKG_CONFIG_FOUND OR NOT BCP_FOUND)
    MESSAGE("BCP directory: ${BCP_DIR}")
    INCLUDE_DIRECTORIES(${BCP_DIR}/include/coin)
    LINK_DIRECTORIES(${BCP_DIR}/lib)
    SET(LIBS ${LIBS} CoinUtils Osi OsiClp Clp Bcp)
endif()

if(DEFINED CTR)
    add_definitions(-DCTR)
endif()



# build librar
INCLUDE_DIRECTORIES(src)
add_subdirectory(src)

# build executables
include_directories(main)
add_subdirectory(main)

MESSAGE("Linked directories: ${LINK_DIRECTORIES}")
MESSAGE("Cmake linker flags: ${CMAKE_LINKER_FLAGS}")

if(DEFINED BUILD_JULIA)
    include_directories(lib_julia_nurse)
    add_subdirectory(lib_julia_nurse)
endif()