
add_airliner_app_def(pe
    FILE PE
    SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_app.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_config_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_events.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_msg.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_tbldefs.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_version.h
   
    INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/

    USER_DOCS   ${CMAKE_CURRENT_SOURCE_DIR}/user_doxy.in
    DESIGN_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/detail_doxy.in
)

add_airliner_app_unit_test(pe-ut
    UTASSERT
    SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_app_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_app_test.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_cmds_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_cmds_test.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_config_tbl_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_config_tbl_test.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_testrunner.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_test_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/pe_test_utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_app.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_app.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_config_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_events.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_msg.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_tbldefs.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pe_version.h
        
    INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/
        
    VALGRIND_SUPPRESSION_FILE
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit_test/valgrind.supp
)
