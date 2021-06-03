set(LIB_NAME constraints)
set(SOURCE_FILES
        ConstraintsMP.cpp
        AssignmentConstraints.cpp
        ResourceConstraints.cpp
        DynamicConstraints.cpp)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

