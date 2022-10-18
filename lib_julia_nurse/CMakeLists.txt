#project(Julia_wrap)
#cmake_minimum_required(VERSION 3.23)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_CXX_STANDARD 17)
include(./CMakeDefinitionsLists.txt)
set(CMAKE_MACOSX_RPATH 1)
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${NURSE_DIR};${JlCxx_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}; ${JlCxx_DIR})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

find_package(Julia)
find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")
include_directories(JlCxx)
include_directories(${Julia_INCLUDE_DIRS})
include_directories(${JlCxx_DIR})
include_directories(${JlCxx_DIR}/../libcxxwrap-julia/include)

find_library(NURSE_LIB nursescheduler)
message(STATUS "Found NURSE_LIB at ${NURSE_LIB}")
include_directories(${NURSE_DIR}/src)

add_library(julia_wrap SHARED tools.cpp tools.h data.cpp data.h solver.cpp solver.h parsers.cpp parsers.h)
target_link_libraries(julia_wrap ${NURSE_LIB}/libnursescheduler.dylib JlCxx::cxxwrap_julia)
