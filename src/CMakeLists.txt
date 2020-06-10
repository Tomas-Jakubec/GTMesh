# ASK TNL_TARGET_INCLUDE_DIRECTORY
INSTALL( DIRECTORY GTMesh/ DESTINATION ${TNL_TARGET_INCLUDE_DIRECTORY}
         MESSAGE_NEVER
         FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" )



# Note that it is important to start building examples as soon as possible,
# because they take the longest time and other stuff can be pipelined before
# they are finished (at least with Ninja).
# if( ${WITH_EXAMPLES} )
#    add_subdirectory( Examples )
# endif()

# if( ${WITH_BENCHMARKS} )
#    ADD_SUBDIRECTORY( Benchmarks )
# endif()

# if( ${WITH_TOOLS} )
#    ADD_SUBDIRECTORY( Tools )
# endif()

if( ${WITH_TESTS} )
   ADD_SUBDIRECTORY( UnitTests )
endif()

# if( ${WITH_PYTHON} )
#    ADD_SUBDIRECTORY( Python )
# endif()
