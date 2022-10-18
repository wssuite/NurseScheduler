add_subdirectory(rcspp)

set(LIB_NAME sp)
set(SOURCE_FILES
        SubProblem.cpp
        RosterSP.cpp
        RotationSP.cpp
        RCSPPSubProblem.cpp
        CyclicRosterSP.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
if(DEFINED RCSPP_SHARED)
    target_link_libraries(${LIB_NAME} mp ${LIBS})
else()
    target_link_libraries(${LIB_NAME} rcspp mp ${LIBS})
endif()
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
