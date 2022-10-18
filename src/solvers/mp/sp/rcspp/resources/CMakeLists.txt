set(LIB_NAME resources)
set(SOURCE_FILES
        TotalWeekendsResource.cpp
        TotalShiftDurationResource.cpp
        ConsShiftResource.cpp
        ConsWeekendShiftResource.cpp
        ForbiddenPatternResource.cpp
        IdentWeekendResource.cpp
        PreferenceResource.cpp
        FreeDaysAfterShiftResource.cpp
        AlternativeShiftResource.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
if(DEFINED RCSPP_SHARED)
    target_link_libraries(${LIB_NAME} nsbase ${LIBS})
else()
    target_link_libraries(${LIB_NAME} rcspp ${LIBS})
endif()
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)