set(LIB_NAME boost)
set(SOURCE_FILES
        LongRotationSP.cpp
        RosterSP.cpp
        RotationSP.cpp
        RCGraph.cpp
        PrincipalGraph.cpp
        RCSPP.cpp
        SubProblem.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
if(DEFINED RCSPP_SHARED)
    target_link_libraries(${LIB_NAME} nsbase ${LIBS})
else()
    target_link_libraries(${LIB_NAME} rcspp ${LIBS})
endif()
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)