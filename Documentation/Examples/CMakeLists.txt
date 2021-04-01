cmake_minimum_required(VERSION 3.12.2)

project(GTMesh_Examples)

# set cmake's include path so that we can include modules from
# the cmake directory in the TNL repository
#set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )

# set C++ standard
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )


include_directories(${CMAKE_SOURCE_DIR}/../../src)
add_subdirectory( Debug )
add_subdirectory( BinarySerializer )
#add_subdirectory( Traits )
#add_subdirectory( UnstructuredMesh )

## below are targets for each test
# add_executable( MyTest MyTest.cpp )
# target_compile_options( MyTest PRIVATE ${CXX_TESTS_FLAGS} )
# target_link_libraries( MyTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
# add_test( MyTest ${EXECUTABLE_OUTPUT_PATH}/MyTest${CMAKE_EXECUTABLE_SUFFIX} )
