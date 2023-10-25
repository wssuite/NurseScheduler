set(LIB_NAME parsing)
set(SOURCE_FILES
        ParseUI.cpp
        ParseINRC2.cpp ParseINRC2.h ParseINRC.cpp ParseINRC.h ParseNRP.cpp ParseNRP.h)

add_library(${LIB_NAME} STATIC ${SOURCE_FILES})
target_link_libraries(${LIB_NAME} data tools ${LIBS})
set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
