add_subdirectory(boost)
add_subdirectory(resources)

set(LIB_NAME rcspp)
set(SOURCE_FILES
        RCLabel.cpp
        RCGraph.cpp
        RCSPP.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} resources boost ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)