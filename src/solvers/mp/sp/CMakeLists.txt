add_subdirectory(rcspp)

set(LIB_NAME sp)
set(SOURCE_FILES
        SubProblem.cpp
        RosterSP.cpp
        RotationSP.cpp
        RCSPPSubProblem.cpp
        CyclicRosterSP.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} rcspp ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
