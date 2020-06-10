
## below are targets for each test
add_executable( MeshDataContainerTest MeshDataContainerTest.cpp )
target_compile_options( MeshDataContainerTest PRIVATE ${CXX_TESTS_FLAGS} )
target_link_libraries( MeshDataContainerTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
add_test( MeshDataContainerTest ${EXECUTABLE_OUTPUT_PATH}/MeshDataContainerTest${CMAKE_EXECUTABLE_SUFFIX} )
