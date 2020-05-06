set(LIB_NAME tools)
set(SOURCE_FILES
        DemandGenerator.cpp
        GlobalStats.cpp
        InputPaths.cpp
        MyTools.cpp
        ReadWrite.cpp
        )

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
