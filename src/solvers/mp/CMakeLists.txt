add_subdirectory(rcspp)
add_subdirectory(modeler)

set(LIB_NAME mp)
set(SOURCE_FILES
        RotationPricer.cpp
        TreeManager.cpp
        )

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} rcspp modeler ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

