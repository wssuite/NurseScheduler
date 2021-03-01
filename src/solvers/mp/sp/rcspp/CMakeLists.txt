add_subdirectory(resources)

set(LIB_NAME rcspp)
set(SOURCE_FILES
        MyRCLabel.cpp
        RCGraph.cpp
        MyRCGraph.cpp
        PrincipalGraph.cpp
        BoostRCSPP.cpp
        MyRCSPP.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} resources ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)