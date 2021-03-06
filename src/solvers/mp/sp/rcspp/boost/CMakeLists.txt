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
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)