
## below are targets for each test
add_executable( BinarySerializerTest BinarySerializerTest.cpp )
target_compile_options( BinarySerializerTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( BinarySerializerTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
add_test( BinarySerializerTest ${EXECUTABLE_OUTPUT_PATH}/BinarySerializerTest${CMAKE_EXECUTABLE_SUFFIX} )
