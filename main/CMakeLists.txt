# exec
add_executable(staticscheduler DeterministicMain.cpp)
target_link_libraries(staticscheduler ${NURSE_SCHEDULER_LIB})

add_executable(dynamicscheduler DynamicMain.cpp)
target_link_libraries(dynamicscheduler ${NURSE_SCHEDULER_LIB})

add_executable(pricer PricerMain.cpp)
target_link_libraries(pricer ${NURSE_SCHEDULER_LIB})

# add all subdirectories (you can create your own main)
FILE(GLOB directories *)
MESSAGE("Directories:" ${directories})
FOREACH(d ${directories})
    IF(IS_DIRECTORY ${d})
        MESSAGE("Add subdirectory " ${d})
        add_subdirectory(${d})
    ENDIF()
ENDFOREACH()
