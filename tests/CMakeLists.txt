cmake_minimum_required(VERSION 3.12.2)

project(GTMesh)

# set cmake's include path so that we can include modules from
# the cmake directory in the TNL repository
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )

# set C++ standard
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

# force colorized output in continuous integration
if( DEFINED ENV{CI_JOB_NAME} OR ${CMAKE_GENERATOR} STREQUAL "Ninja" )
   message(STATUS "Continuous integration or Ninja detected -- forcing compilers to produce colorized output.")
   if( CXX_COMPILER_NAME MATCHES "clang" )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics" )
   else()
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color" )
   endif()
endif()

# enable tests
enable_testing()
include(BuildGtest)
include_directories(${CMAKE_SOURCE_DIR}/../src)
include_directories(${CMAKE_SOURCE_DIR})
add_subdirectory( Debug )
add_subdirectory( BinarySerializer )
add_subdirectory( Traits )
add_subdirectory( UnstructuredMesh )

## below are targets for each test
# add_executable( MyTest MyTest.cpp )
# target_compile_options( MyTest PRIVATE ${CXX_TESTS_FLAGS} )
# target_link_libraries( MyTest ${GTEST_BOTH_LIBRARIES} )


## add all targets to the CMake's "test" target (i.e. `make test` runs all tests)
# add_test( MyTest ${EXECUTABLE_OUTPUT_PATH}/MyTest${CMAKE_EXECUTABLE_SUFFIX} )
