cmake_minimum_required( VERSION 3.11 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 17)
project(HiggsAnalysisCombinedLimit VERSION 10.2.1)

option( MODIFY_ROOTMAP "Modify generated Rootmap to take out classes already bundled in StatAnalysis" FALSE )
option( USE_VDT "Use VDT (fast and vectorisable mathematical functions)" TRUE )
option( BUILD_TESTS "Build tests" FALSE )


# Use the same mechanism as ROOT for the install directories of components (see RootInstallDirs.cmake)
if(NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "user executables (bin)")
endif()
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "object code libraries (lib)")
endif()
if(NOT DEFINED CMAKE_INSTALL_PYTHONDIR)
  set(CMAKE_INSTALL_PYTHONDIR "python" CACHE PATH "python libraries and modules (python)")
endif()
if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
  set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "C header files (include)")
endif()

message (STATUS "Install prefix               : " ${CMAKE_INSTALL_PREFIX})
message (STATUS "Executables installed to     : " ${CMAKE_INSTALL_BINDIR})
message (STATUS "Libraries installed to       : " ${CMAKE_INSTALL_LIBDIR})
message (STATUS "Header files installed to    : " ${CMAKE_INSTALL_INCLUDEDIR})
message (STATUS "Python module installed to   : " ${CMAKE_INSTALL_PYTHONDIR})

# Can build with CMake after e.g. setting up StatAnalysis release like this:
# export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
# source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh; asetup StatAnalysis,0.3.2
# mkdir build; cd build
# cmake path/to/source # change this path to where-ever you cloned Combine repo to
# make -j4

find_package(ROOT REQUIRED COMPONENTS MathMore RooFitCore RooFit RooStats HistFactory)
find_package(Eigen3 REQUIRED)
if(USE_VDT)
  find_package(Vdt REQUIRED)
endif()
find_package(LCG QUIET) # only used for FindBoost in StatAnalysis
find_package(Boost CONFIG REQUIRED COMPONENTS program_options filesystem)

message(STATUS "Using ROOT From: ${ROOT_INCLUDE_DIRS}")
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

set(LIBNAME HiggsAnalysisCombinedLimit)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")

file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} interface/*.h*)
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.c*)

# includes require "HiggsAnalysis/CombinedLimit" prefix in many places
# so create a symlink in the build dir
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/HiggsAnalysis)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/HiggsAnalysis/CombinedLimit" )
include_directories(${CMAKE_CURRENT_BINARY_DIR})

ROOT_GENERATE_DICTIONARY(G__${LIBNAME} HiggsAnalysis/CombinedLimit/src/classes.h LINKDEF src/classes_def.xml
        MODULE ${LIBNAME}
        NOINSTALL
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
set_target_properties(combine PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(MODIFY_ROOTMAP)
        # edit the generated rootmap in-situ before installation

        # Define the path to the generated file
        set(GENERATED_FILE_PATH "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${LIBNAME}.rootmap")

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


# Copy the Python files over to the build directory and add the missing
# __init__.py files, so that text2workspace can be used from the build directory
# Assume we want to copy ${CMAKE_SOURCE_DIR}/assets to ${CMAKE_BINARY_DIR}/assets
set(python_build_dir "${CMAKE_BINARY_DIR}/python")
file(GLOB_RECURSE python_module_files "${CMAKE_SOURCE_DIR}/python/*")
add_custom_command(
    OUTPUT "${python_build_dir}/HiggsAnalysis/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/python" "${python_build_dir}/HiggsAnalysis/CombinedLimit"
    COMMAND ${CMAKE_COMMAND} -E touch "${python_build_dir}/HiggsAnalysis/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E touch "${python_build_dir}/HiggsAnalysis/CombinedLimit/__init__.py"
    DEPENDS ${python_module_files}
    COMMENT "Copying Python module files to build directory"
)
# Create a target that depends on the copy, so it actually happens
add_custom_target(copy_python_module DEPENDS "${python_build_dir}/HiggsAnalysis/__init__.py")
add_dependencies(combine copy_python_module)


# Copy the python scripts into the build tree
set(python_scripts
    ${CMAKE_SOURCE_DIR}/scripts/combineCards.py
    ${CMAKE_SOURCE_DIR}/scripts/plotGof.py
    ${CMAKE_SOURCE_DIR}/scripts/pruneUncerts.py
    ${CMAKE_SOURCE_DIR}/scripts/combineTool.py
    ${CMAKE_SOURCE_DIR}/scripts/plotImpacts.py
    ${CMAKE_SOURCE_DIR}/scripts/text2workspace.py
    ${CMAKE_SOURCE_DIR}/scripts/commentUncerts.py
    ${CMAKE_SOURCE_DIR}/scripts/plot1DScan.py
    ${CMAKE_SOURCE_DIR}/scripts/plotLimitGrid.py
    ${CMAKE_SOURCE_DIR}/scripts/plotBSMxsBRLimit.py
    ${CMAKE_SOURCE_DIR}/scripts/plotLimits.py
)
add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/.stamp_python_scripts"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/bin"
    COMMAND ${CMAKE_COMMAND} -E copy ${python_scripts} "${CMAKE_BINARY_DIR}/bin"
    COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_BINARY_DIR}/.stamp_python_scripts"
    DEPENDS ${python_scripts}
    COMMENT "Copying Python scripts to ${CMAKE_BINARY_DIR}/bin"
)
add_custom_target(copy_python_scripts DEPENDS "${CMAKE_BINARY_DIR}/.stamp_python_scripts")
add_dependencies(combine copy_python_scripts)



#####################
# Installation part #
#####################


# Install the dictionaries.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/lib${LIBNAME}_rdict.pcm
        ${CMAKE_CURRENT_BINARY_DIR}/lib/lib${LIBNAME}.rootmap
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

# Install the scripts like "text2workspace.py" to the bin directory.
install(DIRECTORY scripts/ DESTINATION ${CMAKE_INSTALL_BINDIR} USE_SOURCE_PERMISSIONS)

# The python package will be installed in such a way that the original
# CMSSW-style directory structure is kept, for maximal compatibility.
install(DIRECTORY python/ DESTINATION ${CMAKE_INSTALL_PYTHONDIR}/HiggsAnalysis/CombinedLimit)

# Create empty __init__.py file in the build directory that will be installed
# in the Python library directories.
set(empty_init_py "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")
file(TOUCH ${empty_init_py})
INSTALL(FILES ${empty_init_py} DESTINATION ${CMAKE_INSTALL_PYTHONDIR}/HiggsAnalysis)
INSTALL(FILES ${empty_init_py} DESTINATION ${CMAKE_INSTALL_PYTHONDIR}/HiggsAnalysis/CombinedLimit)


if(BUILD_TESTS)
  find_package(GTest 1.10 REQUIRED)

  enable_testing()
  add_subdirectory(test)
endif()
