set(LIB_NAME rcspp)
set(SOURCE_FILES
        RCGraph.cpp
        SubProblem.cpp
        PrincipalGraph.cpp
        PriceLabelGraph.cpp
        ShortSP.cpp
        #RosterSP.cpp
        )

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
