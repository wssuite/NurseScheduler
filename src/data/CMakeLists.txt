set(LIB_NAME data)
set(SOURCE_FILES
        Demand.cpp
        Nurse.cpp
        Roster.cpp
        Scenario.cpp
        Day.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} tools ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)