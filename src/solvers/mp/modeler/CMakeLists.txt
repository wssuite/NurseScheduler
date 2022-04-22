set(LIB_NAME modeler)
set(SOURCE_FILES
        BcpModeler.cpp
        Stabilization.h
        Stabilization.cpp
        BcpBranchingCandidates.cpp
        BcpHeuristics.cpp)

if (CBC_FOUND)
    set(SOURCE_FILES ${SOURCE_FILES} CbcModeler.cpp)
endif()

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
