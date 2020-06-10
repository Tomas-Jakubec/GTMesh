
## below are targets for each test
add_executable( TraitsTest TraitsTest.cpp )
target_compile_options( TraitsTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( TraitsTest ${GTEST_BOTH_LIBRARIES} )

add_executable( ArithmeticTraitsTest ArithmeticTraitsTest.cpp )
target_compile_options( ArithmeticTraitsTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( ArithmeticTraitsTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
add_test( TraitsTest ${EXECUTABLE_OUTPUT_PATH}/TraitsTest${CMAKE_EXECUTABLE_SUFFIX} )
add_test( ArithmeticTraitsTest ${EXECUTABLE_OUTPUT_PATH}/ArithmeticTraitsTest${CMAKE_EXECUTABLE_SUFFIX} )
