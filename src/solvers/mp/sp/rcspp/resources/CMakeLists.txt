set(LIB_NAME resources)
set(SOURCE_FILES
        TotalWeekendsResource.cpp
        TotalShiftResource.cpp
        ConsShiftResource.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)