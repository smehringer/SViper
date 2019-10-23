# ===========================================================================
#                  SViper - Polish Structural Variants
# ===========================================================================

project (sviper CXX)

cmake_minimum_required (VERSION 3.0.0)

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

message("source dir: ${CMAKE_CURRENT_SOURCE_DIR}")

set (CMAKE_INCLUDE_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan/include
     ${CMAKE_INCLUDE_PATH})
set (CMAKE_PREFIX_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan/util/cmake
     ${CMAKE_PREFIX_PATH})
set (CMAKE_MODULE_PATH
     ${CMAKE_CURRENT_SOURCE_DIR}/extern/seqan/util/cmake
     ${CMAKE_MODULE_PATH})

# Search SeqAn and select dependencies.
find_package(OpenMP QUIET)
find_package(ZLIB   QUIET)
find_package(SeqAn  QUIET REQUIRED CONFIG)

message(STATUS "These dependencies were found:")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     ZLIB       ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
	message (WARNING "WARNING: OpenMP not found. SViper will be built without multi-threading! "
    "This is probably not what you want! Use GCC >= 4.9.1, Clang >= 3.8.0 or ICC >= 16.0.2\n")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
	message (WARNING "WARNING: Zlib not found. Building SViper without support for gzipped input and output (this includes support for .bam).")
endif (NOT ZLIB_FOUND)

message(STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Enable global exception handler for all seqan apps.
set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic")

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic")

add_executable(sviper "src/sviper.cpp")
target_include_directories(sviper PUBLIC "include")
target_link_libraries (sviper ${SEQAN_LIBRARIES})


# Add utility targets
# currently: utilities_compare_vcf.cpp  utilities_get_supporting_reads.cpp  utilities_merge_split_alignments.cpp  utilities_stats.cpp

add_executable(merge_split_alignments "src/utilities_merge_split_alignments.cpp")
target_include_directories(merge_split_alignments PUBLIC "include")
target_link_libraries (merge_split_alignments ${SEQAN_LIBRARIES})

add_executable(compare_vcf "src/utilities_compare_vcf.cpp")
target_include_directories(compare_vcf PUBLIC "include")
target_link_libraries (compare_vcf ${SEQAN_LIBRARIES})

add_executable(get_supporting_reads "src/utilities_get_supporting_reads.cpp")
target_include_directories(get_supporting_reads PUBLIC "include")
target_link_libraries (get_supporting_reads ${SEQAN_LIBRARIES})

add_executable(stats "src/utilities_stats.cpp")
target_include_directories(stats PUBLIC "include")
target_link_libraries (stats ${SEQAN_LIBRARIES})

# ===========================================================================
# Include test subdirectory
# ===========================================================================

#message (STATUS "Configuring tests")
#add_subdirectory (tests)

