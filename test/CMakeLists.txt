cmake_minimum_required(VERSION 3.11)

# Detect whether this is being built standalone
if(NOT PROJECT_NAME)
    project(HiggsAnalysisCombinedLimit LANGUAGES CXX)
    set(REPO "${CMAKE_SOURCE_DIR}/..")
    set(standalone_tests 1)
else()
    set(REPO ${CMAKE_SOURCE_DIR})
endif()

enable_testing()

# Look for gtest. If not found, download it from the web.
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
  T2W_COMMAND
    text2workspace.py ${REPO}/data/tutorials/multiDim/toy-hgg-125.txt -o counting_datacard.root -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
  COMBINE_COMMANDS
    "combine -M MultiDimFit counting_datacard.root  --setParameterRanges r=-5,5"
)

# Counting datacard Fixed Point from csv - combineTool.py
COMBINE_ADD_TEST(counting_datacard_from_csv
    COMMAND combineTool.py -M MultiDimFit counting_datacard.root --fromfile fixed.csv
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/fixed.csv
    FIXTURES_REQUIRED counting_datacard
    # The output is unfortunatley not predictable on all platforms: on
    # CMSSW_11_3_4 - ROOT 6.22.09 and conda with ROOT 6.26.04, it will not
    # print the combine commands that the tool is running.
    # TODO: investigate this.
    # CHECKOUT OUTPUT counting_datacard_from_csv.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard_from_csv.out
)
# set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv >> references/counting_datacard_from_csv.out")

# Parametric analysis
ADD_COMBINE_TEST(parametric_analysis
  T2W_COMMAND
    text2workspace.py ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -o datacard-3-parametric-analysis.root --mass 125
  COMBINE_COMMANDS
    "combine -M MultiDimFit datacard-3-parametric-analysis.root --algo singles --setParameterRanges r=-2,1"
)

# RooHistPdf
set(CMSSW_VERSION $ENV{CMSSW_VERSION})
if(NOT CMSSW_VERSION STREQUAL "CMSSW_11_3_4")
  # For CMSSW_11_3_4 with ROOT 6.22, this test crashes because of some TList
  # memory problem in the end:
  #   Error in <TList::Clear>: A list is accessing an object (0x152333b0) already deleted (list name = TList)
  ADD_COMBINE_TEST(histpdf
    T2W_COMMAND
      text2workspace.py ${REPO}/data/ci/datacard_RooHistPdf.txt.gz -o ws_RooHistPdf.root
    COMBINE_COMMANDS
      "combine -M MultiDimFit ws_RooHistPdf.root --algo singles -v -2  --setParameterRanges r=-1,2"
  )
endif()

# Template analysis CMSHistFunc
ADD_COMBINE_TEST(cmshistfunc
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/template-analysis_shapeInterp.txt -o ws_template-analysis.root --mass 200
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis.root --algo singles  --setParameterRanges r=-5,5"
)

# Template analysis CMSHistFunc with channel masks
ADD_COMBINE_TEST(cmshistfunc_channel_masks
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/htt_multiple_regions.txt  -o ws_template-analysis-channel_masks.root --mass 125 --channel-masks
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis-channel_masks.root --algo singles  --setParameterRanges r=-5,5 --setParameters mask_htt_tt_2_8TeV=1"
)

# Template analysis CMSHistFunc shapeN
ADD_COMBINE_TEST(cmshistfunc_shapeN
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/template-analysis_shapeNInterp.txt -o ws_template-analysis-shapeN.root --mass 200
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis-shapeN.root --algo singles  --setParameterRanges r=-5,5"
)

# Template analysis CMSHistSum
ADD_COMBINE_TEST(cmshistsum
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/template-analysis_shapeInterp.txt -o cmshistsum.root --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit cmshistsum.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
)

# Template analysis CMSHistSum with shapeN
ADD_COMBINE_TEST(cmshistsum_shapeN
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/template-analysis_shapeNInterp.txt -o cmshistsum_shapeN.root --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit cmshistsum_shapeN.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
)

# Template analysis with large integrals
ADD_COMBINE_TEST(template_analysis_large_integrals
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/templ_datacard_largeYields.txt -o template_analysis_large_integrals.root
  COMBINE_COMMANDS
    "combine -M MultiDimFit template_analysis_large_integrals.root --algo singles  --setParameterRanges r=-5,5"
)

# Template analysis with large integrals using CMSHistSum
ADD_COMBINE_TEST(template_analysis_large_integrals_cmshistsum
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/templ_datacard_largeYields.txt -o template_analysis_large_integrals_cmshistsum.root  --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit template_analysis_large_integrals_cmshistsum.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
)

# RooParametricHist
# Unfortunately, the result of this test is not deterministic, so we can't enable it yet.
# TODO: understand why this is, enable the test here, and remove corresponding test from cvmfs-ci.yml
# ADD_COMBINE_TEST(parametric_hist
#   T2W_COMMAND
#   text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO map=.*/*hcc*:r[1,-500,500] --PO map=.*/zcc:z[1,-5,5] ${REPO}/data/ci/datacard_RooParametricHist.txt -o ws_RooParametricHist.root
#   COMBINE_COMMANDS
#     "combine -M MultiDimFit ws_RooParametricHist.root --algo singles"
# )

# Template analysis CMSHistFunc with channel masks
ADD_COMBINE_TEST(cmshistsum_channel_masks
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/htt_multiple_regions.txt  -o ws_template-analysis-channel_masks.root --mass 125 --channel-masks --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis-channel_masks.root --algo singles  --setParameterRanges r=-5,5 --setParameters mask_htt_tt_2_8TeV=1 --X-rtd FAST_VERTICAL_MORPH"
)

# Template-analysis datacard -> text2workspace
COMBINE_ADD_TEST(template_analysis-text2workspace
    COMMAND text2workspace.py ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    FIXTURES_SETUP template_analysis_workspace
)


if(BUILD_TESTS AND NOT DEFINED standalone_tests)
    # Build and run the testCreateNLL helper on the generated workspace
    COMBINE_ADD_GTEST(template-analysis-testCreateNLL
        testCreateNLL.cxx
        LIBRARIES HiggsAnalysisCombinedLimit
    )
    set_property(TEST gtest-template-analysis-testCreateNLL
        PROPERTY FIXTURES_REQUIRED template_analysis_workspace
    )
endif()


# test_interference.py
COMBINE_ADD_TEST(test_interference_py
    COMMAND python3 test_interference.py
    COPY_TO_BUILDDIR ${REPO}/test/test_interference.py
    # The test output is not deterministic yet, so unfortunately we can't do
    # the comparison with the reference file for now.
    # CHECKOUT OUTPUT test_interference_py.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/test_interference_py.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "python test_interference.py >> references/test_interference_py.out")


COMBINE_ADD_TEST(simple-counting-experiment-text2workspace
    COMMAND text2workspace.py ${REPO}/data/tutorials/counting/simple-counting-experiment.txt -o simple-counting-experiment.root
    # No output, so nothing to compare. We just check that it runs.
)

# Test that setting a parameter outside its range throws an error
COMBINE_ADD_TEST(setParameters-out-of-range
    COMMAND combine -M MultiDimFit ${REPO}/data/tutorials/counting/simple-counting-experiment.txt --setParameters r=-1.05
    PASSRC 1
    PASSREGEX "outside its range"
)

# Test that setting a parameter outside its range via regex throws an error
COMBINE_ADD_TEST(setParameters-out-of-range-regex
    COMMAND combine -M MultiDimFit ${REPO}/data/tutorials/counting/simple-counting-experiment.txt --setParameters "rgx{r}=-1.05"
    PASSRC 1
    PASSREGEX "outside its range"
)

COMBINE_ADD_TEST(simple-shapes-TH1-text2workspace
    COMMAND text2workspace.py ${REPO}/data/tutorials/shapes/simple-shapes-TH1.txt -o simple-shapes-TH1.root
    # No output, so nothing to compare. We just check that it runs.
)

COMBINE_ADD_TEST(LHC-Higgs-coupling-models-text2workspace
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/toy-hgg-125.txt
    COMMAND text2workspace.py toy-hgg-125.txt -P HiggsAnalysis.CombinedLimit.LHCHCGModels:K3 -m 125
    CHECKOUT OUTPUT LHC-Higgs-coupling-models-text2workspace.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-Higgs-coupling-models-text2workspace.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "text2workspace.py toy-hgg-125.txt -P HiggsAnalysis.CombinedLimit.LHCHCGModels:K3 -m 125 > references/LHC-Higgs-coupling-models-text2workspace.out"
)

COMBINE_ADD_TEST(hzz4l-significance
    COMMAND combine -M Significance ${REPO}/test/inputs/hzz4l_datacard.txt
    CHECKOUT OUTPUT hzz4l-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/hzz4l-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine -M Significance hzz4l_datacard.txt > references/hzz4l-significance.out"
)

COMBINE_ADD_TEST(LHC-limits
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01
    CHECKOUT OUTPUT LHC-limits.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-limits.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01 > references/LHC-limits.out"
)

COMBINE_ADD_TEST(LHC-significance
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125
    CHECKOUT OUTPUT LHC-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125 > references/LHC-significance.out"
)

COMBINE_ADD_TEST(parametric-analysis-significance
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M Significance --mass 125
    CHECKOUT OUTPUT parametric-analysis-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric-analysis-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-3-parametric-analysis.txt -M Significance --mass 125 > references/parametric-analysis-significance.out"
)

COMBINE_ADD_TEST(MarkovChainMC
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100
    CHECKOUT OUTPUT MarkovChainMC.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/MarkovChainMC.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100 > references/MarkovChainMC.out"
)

COMBINE_ADD_TEST(multi-signal-text2workspace
    COMMAND text2workspace.py ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125
    CHECKOUT OUTPUT multi-signal-text2workspace.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal-text2workspace.out
    FIXTURES_SETUP multi-signal-text2workspace
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "text2workspace.py datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125 > references/multi-signal-text2workspace.out"
)

COMBINE_ADD_TEST(multi-signal
    COMMAND combine datacard-5-multi-signal.root -M MultiDimFit --algo singles --mass 125
    CHECKOUT OUTPUT multi-signal.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal.out
    FIXTURES_REQUIRED multi-signal-text2workspace
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-5-multi-signal.root -M MultiDimFit --algo singles --mass 125 > references/multi-signal.out"
)

COMBINE_ADD_TEST(multi-signal-ChannelCompatibilityCheck
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt -M ChannelCompatibilityCheck --mass 125
    CHECKOUT OUTPUT multi-signal-ChannelCompatibilityCheck.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal-ChannelCompatibilityCheck.out
    FIXTURES_REQUIRED multi-signal
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-5-multi-signal.txt -M ChannelCompatibilityCheck --mass 125 > references/multi-signal-ChannelCompatibilityCheck.out"
)


# Write the script to produce the reference files. After building combine, you
# can produce the reference files as follows:
#
#   cd build
#   ctest -j$(nproc)
#   cd test
#   sh create_reference_files.sh
#
# You have now a new reference files directory in build/test/references that
# you can use to replace test/references in the repository.

get_property(all_combine_commands GLOBAL PROPERTY ALL_COMBINE_COMMANDS)
message(STATUS "All combine commands:")
set(_outfile "${CMAKE_BINARY_DIR}/test/create_reference_files.sh")
file(WRITE "${_outfile}" "#!/usr/bin/env bash\n\n")
file(APPEND "${_outfile}" "rm -rf references\n\n")
file(APPEND "${_outfile}" "mkdir references\n\n")
if(DEFINED standalone_tests)
    file(APPEND "${_outfile}" "export COMBINE_SRC=${CMAKE_SOURCE_DIR}/..\n\n")
else()
    file(APPEND "${_outfile}" "export COMBINE_SRC=${CMAKE_SOURCE_DIR}\n\n")
endif()
file(APPEND "${_outfile}" "PATH=${CMAKE_BINARY_DIR}/bin:$PATH\n")
file(APPEND "${_outfile}" "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$LD_LIBRARY_PATH\n")
file(APPEND "${_outfile}" "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$PYTHONPATH\n\n")
foreach(cmd IN LISTS all_combine_commands)
    file(APPEND "${_outfile}" "${cmd}\n")
endforeach()
message(STATUS "Wrote script to produce test reference files ${_outfile}")
