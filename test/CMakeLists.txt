cmake_minimum_required(VERSION 3.11)

# Detect whether this is being built standalone
if(NOT PROJECT_NAME)
    project(HiggsAnalysisCombinedLimit LANGUAGES CXX)
    set(REPO "${CMAKE_SOURCE_DIR}/..")
    set(standalone_tests 1)
    # Check if the ROOTSYS environment variable is defined
    if(DEFINED ENV{ROOTSYS})
        message(STATUS "ROOTSYS found at $ENV{ROOTSYS}")
        include_directories($ENV{ROOTSYS}/include)
    else()
        message(STATUS "ROOTSYS not defined — skipping ROOT includes")
    endif()
else()
    set(REPO ${CMAKE_SOURCE_DIR})
endif()

enable_testing()

# Loot for gtest. If not found, download it from the web.
# Taken from ROOTs SearchInstalledSoftware.cmake
find_package(GTest 1.10)
if(NOT GTEST_FOUND)
  include(ExternalProject)

  # Add googletest
  # http://stackoverflow.com/questions/9689183/cmake-googletest

  set(_gtest_byproduct_binary_dir
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-prefix/src/googletest-build)
  set(_gtest_byproducts
    ${_gtest_byproduct_binary_dir}/lib/libgtest.a
    ${_gtest_byproduct_binary_dir}/lib/libgtest_main.a
    ${_gtest_byproduct_binary_dir}/lib/libgmock.a
    ${_gtest_byproduct_binary_dir}/lib/libgmock_main.a
    )

  set(GTEST_CXX_FLAGS "${ROOT_EXTERNAL_CXX_FLAGS}")
  if(MSVC)
     if(winrtdebug)
      set(GTEST_BUILD_TYPE Debug)
    else()
      set(GTEST_BUILD_TYPE Release)
    endif()
    set(_gtest_byproducts
      ${_gtest_byproduct_binary_dir}/lib/gtest.lib
      ${_gtest_byproduct_binary_dir}/lib/gtest_main.lib
      ${_gtest_byproduct_binary_dir}/lib/gmock.lib
      ${_gtest_byproduct_binary_dir}/lib/gmock_main.lib
    )
    if(CMAKE_GENERATOR MATCHES Ninja)
      set(GTEST_BUILD_COMMAND "BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>")
    else()
      set(GTEST_BUILD_COMMAND "BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config ${GTEST_BUILD_TYPE}")
    endif()
    if(asan)
      if(NOT winrtdebug)
        set(gtestbuild "RelWithDebInfo")
      endif()
      set(GTEST_CXX_FLAGS "${ROOT_EXTERNAL_CXX_FLAGS} ${ASAN_EXTRA_CXX_FLAGS}")
    endif()
    set(EXTRA_GTEST_OPTS
      -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${_gtest_byproduct_binary_dir}/lib/
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL:PATH=${_gtest_byproduct_binary_dir}/lib/
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${_gtest_byproduct_binary_dir}/lib/
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO:PATH=${_gtest_byproduct_binary_dir}/lib/
      -Dgtest_force_shared_crt=ON
      ${GTEST_BUILD_COMMAND})
  else()
    set(GTEST_BUILD_TYPE Release)
  endif()
  if(APPLE)
    set(EXTRA_GTEST_OPTS
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT})
  endif()

  ExternalProject_Add(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_SHALLOW 1
    GIT_TAG release-1.12.1
    UPDATE_COMMAND ""
    # # Force separate output paths for debug and release builds to allow easy
    # # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    # CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
    #            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
    #            -Dgtest_force_shared_crt=ON
    CMAKE_ARGS -G ${CMAKE_GENERATOR}
                  -DCMAKE_BUILD_TYPE=${GTEST_BUILD_TYPE}
                  # Don't inherit from main project, because it doesn't use C (only C++).
                  #-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DCMAKE_CXX_FLAGS=${GTEST_CXX_FLAGS}
                  -DCMAKE_AR=${CMAKE_AR}
                  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                  ${EXTRA_GTEST_OPTS}
    # Disable install step
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${_gtest_byproducts}
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_OUTPUT_ON_FAILURE ON
    TIMEOUT 600
  )

  # Specify include dirs for gtest and gmock
  ExternalProject_Get_Property(googletest source_dir)
  set(GTEST_INCLUDE_DIR ${source_dir}/googletest/include)
  set(GMOCK_INCLUDE_DIR ${source_dir}/googlemock/include)
  # Create the directories. Prevents bug https://gitlab.kitware.com/cmake/cmake/issues/15052
  file(MAKE_DIRECTORY ${GTEST_INCLUDE_DIR} ${GMOCK_INCLUDE_DIR})

  # Libraries
  ExternalProject_Get_Property(googletest binary_dir)
  set(_G_LIBRARY_PATH ${binary_dir}/lib/)

  # Use gmock_main instead of gtest_main because it initializes gtest as well.
  # Note: The libraries are listed in reverse order of their dependencies.
  foreach(lib gtest gtest_main gmock gmock_main)
    add_library(${lib} IMPORTED STATIC GLOBAL)
    set_target_properties(${lib} PROPERTIES
      IMPORTED_LOCATION "${_G_LIBRARY_PATH}${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
    add_dependencies(${lib} googletest)
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND
        ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER_EQUAL 9)
      target_compile_options(${lib} INTERFACE -Wno-deprecated-copy)
    endif()
  endforeach()
  target_include_directories(gtest INTERFACE ${GTEST_INCLUDE_DIR})
  target_include_directories(gmock INTERFACE ${GMOCK_INCLUDE_DIR})

  set_property(TARGET gtest PROPERTY IMPORTED_LOCATION ${_G_LIBRARY_PATH}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX})
  set_property(TARGET gtest_main PROPERTY IMPORTED_LOCATION ${_G_LIBRARY_PATH}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX})
  set_property(TARGET gmock PROPERTY IMPORTED_LOCATION ${_G_LIBRARY_PATH}/${CMAKE_STATIC_LIBRARY_PREFIX}gmock${CMAKE_STATIC_LIBRARY_SUFFIX})
  set_property(TARGET gmock_main PROPERTY IMPORTED_LOCATION ${_G_LIBRARY_PATH}/${CMAKE_STATIC_LIBRARY_PREFIX}gmock_main${CMAKE_STATIC_LIBRARY_SUFFIX})

  # Starting from cmake 3.23, the GTest targets will have stable names.
  # ROOT was updated to use those, but for older CMake versions, we have to declare the aliases:
  foreach(LIBNAME gtest_main gmock_main gtest gmock)
    if(NOT TARGET GTest::${LIBNAME} AND TARGET ${LIBNAME})
      add_library(GTest::${LIBNAME} ALIAS ${LIBNAME})
    endif()
  endforeach()
endif()

include(CombineTestMacros.cmake)

# Set the location of text2workspace.py, assuming installation to a virtual
# environment inside the source directory.
set(t2w ${CMAKE_BINARY_DIR}/bin/text2workspace.py)

# Counting datacard
ADD_COMBINE_TEST(counting_datacard
  COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/toy-hgg-125.txt
  T2W_COMMAND
    text2workspace.py toy-hgg-125.txt -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
  COMBINE_COMMANDS
    "combine -M MultiDimFit toy-hgg-125.root  --setParameterRanges r=-1,1"
)

# Counting datacard Fixed Point from csv - combineTool.py
COMBINE_ADD_TEST(counting_datacard_from_csv
    COMMAND combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/fixed.csv
    FIXTURES_REQUIRED counting_datacard
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    # The output is unfortunatley not predictable on all platforms: on
    # CMSSW_11_3_4 - ROOT 6.22.09 and conda with ROOT 6.26.04, it will not
    # print the combine commands that the tool is running.
    # TODO: investigate this.
    # CHECKOUT OUTPUT counting_datacard_from_csv.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard_from_csv.out
    ENVIRONMENT PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH} # So the combineTool.py finds the combine executable
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
)
# set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv >> references/counting_datacard_from_csv.out")

# Parametric analysis
ADD_COMBINE_TEST(parametric_analysis
  COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
                   ${REPO}/data/tutorials/CAT23001/parametric-analysis-datacard-input.root
  T2W_COMMAND
    text2workspace.py datacard-3-parametric-analysis.txt --mass 125
  COMBINE_COMMANDS
    "combine -M MultiDimFit datacard-3-parametric-analysis.root --algo singles --setParameterRanges r=-2,1"
)

# TODO: make these tests work too
#
# # Template analysis CMSHistFunc
# ADD_COMBINE_TEST(cmshistfunc
#   COPY_TO_BUILDDIR ${REPO}/data/ci/template-analysis_shapeInterp.txt
#                    ${REPO}/data/ci/htt_input.root
#   T2W_COMMAND
#     text2workspace.py template-analysis_shapeInterp.txt -o ws_template-analysis.root --mass 200
#   COMBINE_COMMANDS
#     "combine -M MultiDimFit ws_template-analysis.root --algo singles  --setParameterRanges r=-1,1"
#     "combine -M FitDiagnostics ws_template-analysis.root  -t -1 --setParameters r=1 --plots  --setParameterRanges r=-1,1"
# )
# 
# # Template analysis CMSHistFunc shapeN
# ADD_COMBINE_TEST(cmshistfunc_shapeN
#   COPY_TO_BUILDDIR ${REPO}/data/ci/template-analysis_shapeNInterp.txt
#                    ${REPO}/data/ci/htt_input.root
#   T2W_COMMAND
#     text2workspace.py template-analysis_shapeNInterp.txt -o ws_template-analysis-shapeN.root --mass 200
#   COMBINE_COMMANDS
#     "combine -M MultiDimFit ws_template-analysis-shapeN.root --algo singles  --setParameterRanges r=-1,1"
#     "combine -M FitDiagnostics ws_template-analysis-shapeN.root  -t -1 --setParameters r=1 --plots  --setParameterRanges r=-1,1"
# )

# Template analysis CMSHistSum
ADD_COMBINE_TEST(cmshistsum
  COPY_TO_BUILDDIR ${REPO}/data/ci/template-analysis_shapeInterp.txt
                   ${REPO}/data/ci/htt_input.root
  T2W_COMMAND
    text2workspace.py template-analysis_shapeInterp.txt --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit template-analysis_shapeInterp.root --algo singles  --setParameterRanges r=-1,1 --X-rtd FAST_VERTICAL_MORPH"
)

# Template analysis CMSHistSum with shapeN
ADD_COMBINE_TEST(cmshistsum_shapeN
  COPY_TO_BUILDDIR ${REPO}/data/ci/template-analysis_shapeNInterp.txt
                   ${REPO}/data/ci/htt_input.root
  T2W_COMMAND
    text2workspace.py template-analysis_shapeNInterp.txt --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit template-analysis_shapeNInterp.root --algo singles  --setParameterRanges r=-1,1 --X-rtd FAST_VERTICAL_MORPH"
)

# Template-analysis datacard -> text2workspace
COMBINE_ADD_TEST(template_analysis-text2workspace
    COMMAND text2workspace.py template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    COPY_TO_BUILDDIR
        ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt
        ${REPO}/data/ci/htt_input.root
    ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
                PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
    FIXTURES_SETUP template_analysis_workspace
)


# Build and run the testCreateNLL helper on the generated workspace
COMBINE_ADD_GTEST(template-analysis-testCreateNLL
    testCreateNLL.cxx
    LIBRARIES HiggsAnalysisCombinedLimit
)
set_property(TEST gtest-template-analysis-testCreateNLL
    PROPERTY FIXTURES_REQUIRED template_analysis_workspace
)


# test_interference.py
COMBINE_ADD_TEST(test_interference_py
    COMMAND python3 test_interference.py
    COPY_TO_BUILDDIR ${REPO}/test/test_interference.py
    ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
                PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    # The test output is not deterministic yet, so unfortunately we can't do
    # the comparison with the reference file for now.
    # CHECKOUT OUTPUT test_interference_py.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/test_interference_py.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "python test_interference.py >> references/test_interference_py.out")


# Write the script to produce the reference files. After building combine, you
# can produce the reference files as follows:
#
#   cd build
#   ctest -j8
#   cd test
#   sh create_reference_files.sh
#
# You have now a new reference files directory in build/test/references that
# you can use to replace test/references in the repository.

get_property(all_combine_commands GLOBAL PROPERTY ALL_COMBINE_COMMANDS)
message(STATUS "All combine commands:")
set(_outfile "${CMAKE_BINARY_DIR}/test/create_reference_files.sh")
file(WRITE "${_outfile}" "#!/usr/bin/env bash\n\n")
file(APPEND "${_outfile}" "mkdir references\n\n")
file(APPEND "${_outfile}" "PATH=${CMAKE_BINARY_DIR}/bin:$PATH\n")
file(APPEND "${_outfile}" "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$LD_LIBRARY_PATH\n")
file(APPEND "${_outfile}" "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$PYTHONPATH\n\n")
foreach(cmd IN LISTS all_combine_commands)
    file(APPEND "${_outfile}" "${cmd}\n")
endforeach()
message(STATUS "Wrote script to produce test reference files ${_outfile}")
