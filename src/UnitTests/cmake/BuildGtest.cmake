# Gtest developers recommend to build the gtest libraries directly from
# the projects' build systems, see
# https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

#find_package( GTest )
#if( GTEST_FOUND )
#   set( CXX_TESTS_FLAGS "${CXX_TESTS_FLAGS} -DHAVE_GTEST" )
#endif( GTEST_FOUND )


# compatibility with the GTest package
set( GTEST_BOTH_LIBRARIES gtest gtest_main )
set( CXX_TESTS_FLAGS ${CXX_TESTS_FLAGS} -DHAVE_GTEST )


# Download and unpack googletest at configure time
configure_file(cmake/Gtest.cmake.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
 RESULT_VARIABLE result
 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
 message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
 RESULT_VARIABLE result
 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
 message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                ${CMAKE_BINARY_DIR}/googletest-build)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()
