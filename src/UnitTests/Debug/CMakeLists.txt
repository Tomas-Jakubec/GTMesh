
## below are targets for each test
add_executable( VariableExportTest VariableExportTest.cpp )
target_compile_options( VariableExportTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( VariableExportTest ${GTEST_BOTH_LIBRARIES} )

add_executable( DBGVAR_JSONTest DBGVAR_JSONTest.cpp )
target_compile_options( DBGVAR_JSONTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( DBGVAR_JSONTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
add_test( VariableExportTest ${EXECUTABLE_OUTPUT_PATH}/VariableExportTest${CMAKE_EXECUTABLE_SUFFIX} )
add_test( DBGVAR_JSONTest ${EXECUTABLE_OUTPUT_PATH}/DBGVAR_JSONTest${CMAKE_EXECUTABLE_SUFFIX} )
