set(LIB_NAME modeler)
set(SOURCE_FILES
        BcpModeler.cpp
#        CbcModeler.cpp
        Stabilization.h
        Stabilization.cpp BcpBranchingCandidates.cpp BcpBranchingCandidates.h)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
