add_subdirectory(sp)
add_subdirectory(modeler)

set(LIB_NAME mp)
set(SOURCE_FILES
        MasterProblem.cpp
        RotationMP.cpp
        RCPricer.cpp
        TreeManager.cpp
        RosterMP.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} sp modeler ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

