cmake_minimum_required( VERSION 3.11 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 17)
project(HiggsAnalysisCombinedLimit VERSION 0.0.1)

option( MODIFY_ROOTMAP "Modify generated Rootmap to take out classes already bundled in StatAnalysis" FALSE )
option( INSTALL_PYTHON "Install the Python library and scripts" TRUE )
option( USE_VDT "Use VDT (fast and vectorisable mathematical functions)" TRUE )

# Can build with CMake after e.g. setting up StatAnalysis release like this:
# export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
# source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh; asetup StatAnalysis,0.3.2
# mkdir build; cd build
# cmake path/to/source # change this path to where-ever you cloned Combine repo to
# make -j4

find_package(ROOT REQUIRED COMPONENTS MathMore RooFitCore RooFit RooStats HistFactory)
find_package(Eigen3 REQUIRED)
find_package(Vdt)
find_package(LCG QUIET) # only used for FindBoost in StatAnalysis
find_package(Boost CONFIG REQUIRED COMPONENTS program_options filesystem)

message(STATUS "Using ROOT From: ${ROOT_INCLUDE_DIRS}")
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

set(LIBNAME HiggsAnalysisCombinedLimit)

file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} interface/*.h*)
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.c*)

# includes require "HiggsAnalysis/CombinedLimit" prefix in many places
# so create a symlink in the build dir
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/HiggsAnalysis)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/HiggsAnalysis/CombinedLimit" )
include_directories(${CMAKE_CURRENT_BINARY_DIR})

ROOT_GENERATE_DICTIONARY(G__${LIBNAME} HiggsAnalysis/CombinedLimit/src/classes.h LINKDEF src/classes_def.xml
        MODULE ${LIBNAME}
        OPTIONS ${ROOTCLING_OPTIONS})
add_library(${LIBNAME} SHARED ${SOURCES} G__${LIBNAME}.cxx)
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries (${LIBNAME} Eigen3::Eigen ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

if(NOT USE_VDT)
  target_compile_definitions(${LIBNAME} PUBLIC COMBINE_NO_VDT)
else()
  target_link_libraries (${LIBNAME} VDT::VDT)
endif()

add_executable(combine bin/combine.cpp)
target_link_libraries(combine PUBLIC ${LIBNAME})

if(MODIFY_ROOTMAP)
        # edit the generated rootmap in-situ before installation
        
        # Define the path to the generated file
        set(GENERATED_FILE_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap")

        # Define the custom command to search and replace in-place
        add_custom_command(
          OUTPUT "${GENERATED_FILE_PATH}.bak"  # Declare the generated file as output
          DEPENDS "${GENERATED_FILE_PATH}"  # Ensure it depends on the original generated file
          COMMAND sed -i.bak -e "/class RooParamKeysPdf/d"
                            -e "/class RooStarMomentMorph/d"
                            -e "/class RooStats::HistFactory::RooBSpline/d"
                            -e "/class RooStats::HistFactory::RooBSplineBases/d" 
                            -e "/class ResponseFunction/d"
                            -e "/class RooTwoSidedCBShape/d"
                            "${GENERATED_FILE_PATH}"
          COMMENT "Removing conflicting classes from generated rootmap"
        )

        # Define a custom target that depends on the custom command output
        add_custom_target(
                ModifyRootmapFile  # Name of the custom target
                DEPENDS "${GENERATED_FILE_PATH}.bak"  # Ensure the custom command is executed
        )

        # Add the custom target as a dependency to another target to ensure it's built
        add_dependencies(combine ModifyRootmapFile)
endif()

#####################
# Installation part #
#####################


# Install the dictionaries.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}_rdict.pcm
        ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
        DESTINATION lib)

# Install the libraries and header files.
install(TARGETS ${LIBNAME}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/HiggsAnalysis/CombinedLimit/interface
        )
# dictionary requires classes.h to be in include path exactly as it was specified
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/HiggsAnalysis/CombinedLimit/src/classes.h
        DESTINATION include/HiggsAnalysis/CombinedLimit/src)

# Install the "combine" executable in the bin directory.
install(TARGETS combine DESTINATION bin)


if(INSTALL_PYTHON)
  # In conda-forge cross-compiling builds we hardcode Python_SITELIB because
  # find_package does not ever find the right version, unfortunately.
  # In any other environment, find_package should work fine
  if(NOT DEFINED Python_SITELIB)
          find_package(Python REQUIRED COMPONENTS Interpreter)
  endif()

  # The the Python library installation directory relative to the install prefix.
  file(RELATIVE_PATH Python_SITELIB_IN_PREFIX ${CMAKE_INSTALL_PREFIX} ${Python_SITELIB})
  message (STATUS "Using Python install location:" ${Python_SITELIB_IN_PREFIX})

  # Install the scripts like "text2workspace.py" to the bin directory.
  install(DIRECTORY scripts/ DESTINATION bin USE_SOURCE_PERMISSIONS)

  # The python package will be installed in such a way that the original
  # CMSSW-style directory structure is kept, for maximal compatibility.
  install(DIRECTORY python/ DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis/CombinedLimit)

  # Create empty __init__.py file in the build directory that will be installed
  # in the Python library directories.
  set(empty_init_py "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")
  file(TOUCH ${empty_init_py})
  INSTALL(FILES ${empty_init_py} DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis)
  INSTALL(FILES ${empty_init_py} DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis/CombinedLimit)
endif()
