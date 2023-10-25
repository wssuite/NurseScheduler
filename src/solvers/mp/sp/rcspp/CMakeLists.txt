set(SUB_LIBS resources)
add_subdirectory(resources)

if (${USE_BOOST})
    set(SUB_LIBS ${SUB_LIBS} boost)
    add_subdirectory(boost)
endif ()

set(LIB_NAME rcspp)
set(SOURCE_FILES
        RCLabel.cpp
        RCGraph.cpp
        RCSPP.cpp)

if(DEFINED RCSPP_SHARED)
    # build the dynamic library
    add_library(${LIB_NAME} SHARED ${SOURCE_FILES})
    target_link_libraries(${LIB_NAME} ${SUB_LIBS} nsbase ${LIBS})
else()
    add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
    target_link_libraries(${LIB_NAME} ${SUB_LIBS} sp ${LIBS})
    set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()