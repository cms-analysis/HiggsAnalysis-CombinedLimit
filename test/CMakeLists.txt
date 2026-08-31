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

# Asymptotic limits
ADD_COMBINE_TEST(asymptotic_limits
  T2W_COMMAND
    text2workspace.py ${REPO}/data/tutorials/counting/realistic-counting-experiment.txt -o asymptotic_limits.root -m 125
  COMBINE_COMMANDS
    "combine -M AsymptoticLimits asymptotic_limits.root -m 125 -n AsymptoticLimitsRegression"
  DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.001
)
set_tests_properties(asymptotic_limits-text2workspace asymptotic_limits
    PROPERTIES LABELS AsymptoticLimits TIMEOUT 120
)

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
    # The output is unfortunately not predictable on all platforms: on
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

# RooMultiPdf
if(NOT CMSSW_VERSION MATCHES "^CMSSW_11")
  # This test does not work with the ROOT 6.22 of CMSSW_11; the condition is
  # inherited from the former cvmfs-ci.yml step.
  ADD_COMBINE_TEST(multipdf
    T2W_COMMAND
      text2workspace.py ${REPO}/data/ci/datacard_RooMultiPdf.txt.gz -o ws_RooMultiPdf.root
    # Both commands run with the lowest verbosity: with ROOT >= 6.38, Minos
    # stepping to the boundary of r makes RooFit print range-clipping warnings
    # on every likelihood evaluation, and on some platforms Migrad probes a
    # point where the Chebychev transfer functions are negative, which
    # produces an evaluation error log. Neither affects the fit result.
    COMBINE_COMMANDS
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=2 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles ws_RooMultiPdf.root"
      # Also run Hesse on the result.
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=2 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles --saveFitResult ws_RooMultiPdf.root"
      # Start the discrete profiling from a different pdf index. The result
      # must not depend on the starting index: the parameters of the pdfs that
      # are inactive at the start have to be floated again when the scan
      # activates them (see CascadeMinimizer::minimize).
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=0 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles ws_RooMultiPdf.root"
    # The last digit of the Minos uncertainties differs between ROOT versions,
    # so compare floating point numbers with a relative tolerance of 1 %.
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
  )
endif()

# Template analysis CMSHistFunc
ADD_COMBINE_TEST(cmshistfunc
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/template-analysis_shapeInterp.txt -o ws_template-analysis.root --mass 200
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis.root --algo singles  --setParameterRanges r=-5,5"
)

# FastScan on the template-analysis workspace (combineTool.py)
#
# The output mixes Python prints with RooFit messages written from C++. When
# stdout is not a terminal, both sides are block-buffered independently, so
# the lines would come out in an arbitrary order (a RooFit message can even
# be spliced into the middle of a Python line). Force line buffering on both
# sides to get the lines in program order: PYTHONUNBUFFERED for Python and
# stdbuf for the C stdio buffer used by RooFit.
set(_fastscan_command "PYTHONUNBUFFERED=1 stdbuf -oL combineTool.py -M FastScan -w ws_template-analysis.root:w")
COMBINE_ADD_TEST(cmshistfunc-FastScan
    COMMAND bash -c "${_fastscan_command}"
    FIXTURES_REQUIRED cmshistfunc
    CHECKOUT OUTPUT cmshistfunc-FastScan.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistfunc-FastScan.out
    # Compare floating point numbers with a relative tolerance of 1 % plus an
    # absolute tolerance of 0.1. The absolute tolerance is needed because the
    # "Better minimum" lines contain NLL differences at the level of numerical
    # noise, and for nearly flat parameters the best scan point can move by
    # one step (0.07 for the parameters with a -7..7 range) between ROOT
    # versions. The lines that are not reproducible (the RooFit timing message
    # for the NLL creation, the printout of the NLL object and the RooFit
    # numerical derivator error) are removed by filter_test_output.sh.
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01 --atol 0.1
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "${_fastscan_command} > references/cmshistfunc-FastScan.out"
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
ADD_COMBINE_TEST(parametric_hist
  COPY_TO_BUILDDIR ${REPO}/data/ci/shapes_RooParametricHist.root
  T2W_COMMAND
  text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO map=.*/*hcc*:r[1,-500,500] --PO map=.*/zcc:z[1,-5,5] ${REPO}/data/ci/datacard_RooParametricHist.txt -o ws_RooParametricHist.root
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_RooParametricHist.root --algo singles"
  # The last digits of the Minos uncertainties differ between ROOT versions,
  # so compare floating point numbers with a relative tolerance of 1 %.
  DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
)

# Template analysis CMSHistFunc with channel masks
ADD_COMBINE_TEST(cmshistsum_channel_masks
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/htt_multiple_regions.txt  -o ws_template-analysis-channel_masks.root --mass 125 --channel-masks --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit ws_template-analysis-channel_masks.root --algo singles  --setParameterRanges r=-5,5 --setParameters mask_htt_tt_2_8TeV=1 --X-rtd FAST_VERTICAL_MORPH"
)

# FitDiagnostics shapes from a CMSHistErrorPropagator and a CMSHistSum workspace
# of the same datacard must agree: the content and uncertainty of every bin of
# every process, pre-fit and after the background-only and s+b fits.
COMBINE_ADD_TEST(cmshistsum_vs_cmshistfunc_fitdiagnostics-text2workspace
    COMMAND bash -c "text2workspace.py ${REPO}/data/ci/template-analysis_shapeInterp.txt -o fitdiag_cmshistfunc.root --mass 200 && text2workspace.py ${REPO}/data/ci/template-analysis_shapeInterp.txt -o fitdiag_cmshistsum.root --mass 200 --use-histsum"
    FIXTURES_SETUP cmshistsum_vs_cmshistfunc_fitdiagnostics
)
COMBINE_ADD_TEST(cmshistsum_vs_cmshistfunc_fitdiagnostics
    COMMAND bash -c "combine -M FitDiagnostics fitdiag_cmshistfunc.root -m 200 -n cmshistfunc --saveShapes --saveWithUncertainties --numToysForShapes 20 && combine -M FitDiagnostics fitdiag_cmshistsum.root -m 200 -n cmshistsum --saveShapes --saveWithUncertainties --numToysForShapes 20 && python3 ${REPO}/test/printFitDiagnosticsShapes.py fitDiagnosticscmshistfunc.root > shapes_cmshistfunc.txt && python3 ${REPO}/test/printFitDiagnosticsShapes.py fitDiagnosticscmshistsum.root > shapes_cmshistsum.txt && python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 1e-5 shapes_cmshistfunc.txt shapes_cmshistsum.txt"
    FIXTURES_REQUIRED cmshistsum_vs_cmshistfunc_fitdiagnostics
)

# Template-analysis datacard -> text2workspace
COMBINE_ADD_TEST(template_analysis-text2workspace
    COMMAND text2workspace.py ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    FIXTURES_SETUP template_analysis_workspace
)


if(BUILD_TESTS AND NOT DEFINED standalone_tests)
    COMBINE_ADD_GTEST(hzz4l-readfile
        testHZZ4LReadFile.cxx
        LIBRARIES HiggsAnalysisCombinedLimit
        LABELS HZZ4L
    )
    COMBINE_ADD_GTEST(roopower
        testRooPower.cxx
        LIBRARIES HiggsAnalysisCombinedLimit
        LABELS RooPower
    )
    # Build and run the testCreateNLL helper on the generated workspace
    COMBINE_ADD_GTEST(template-analysis-testCreateNLL
        testCreateNLL.cxx
        LIBRARIES HiggsAnalysisCombinedLimit
    )
    set_property(TEST gtest-template-analysis-testCreateNLL
        PROPERTY FIXTURES_REQUIRED template_analysis_workspace
    )
endif()

# Test NLL consistency between backends
if(NOT DEFINED standalone_tests)
    set(nll_consistency_backends combine)
    if(ROOT_VERSION VERSION_GREATER_EQUAL "6.30.00")
        list(APPEND nll_consistency_backends cpu codegen)
    endif()
    # Check that the NLLs computed with different backends are consistent for various cases
    foreach(backend IN LISTS nll_consistency_backends)
        set(nll_consistency_cases extra_constraints multipdf_penalty
            ignore_constraint_poi_selection ignore_constraint_selected_nuisance ignored_extra_constraint)
        # The generic constraint-only regression uses an empty-channel model;
        # the Combine backend requires at least one channel.
        if(NOT backend STREQUAL "combine")
            list(APPEND nll_consistency_cases constraint_only_prefit)
        endif()
        foreach(case IN LISTS nll_consistency_cases)
            COMBINE_ADD_TEST(nll_consistency_${case}_${backend}
                COMMAND python3 ${REPO}/test/test_nll_backend_consistency.py
                    --case ${case} --backend ${backend} --combine $<TARGET_FILE:combine>
                LABELS NLLConsistency
                TIMEOUT 300
            )
        endforeach()
    endforeach()
    # Check special cases for FitDiagnostics and ChannelCompatibilityCheck
    foreach(case fitdiagnostics_prefit channel_compatibility_constraints)
        COMBINE_ADD_TEST(nll_consistency_${case}_combine
            COMMAND python3 ${REPO}/test/test_nll_backend_consistency.py
                --case ${case} --backend combine --combine $<TARGET_FILE:combine>
            LABELS NLLConsistency
            TIMEOUT 240
        )
    endforeach()
endif()

# ModelTools.py
COMBINE_ADD_TEST(test_model_tools_py
    COMMAND python3 ${REPO}/test/test_model_tools.py
)

# Copies of the vertical-morphing pdfs must follow their own nuisance parameter
COMBINE_ADD_TEST(test_morphing_copy_py
    COMMAND python3 ${REPO}/test/test_morphing_copy.py
)

# RooBernsteinFast v1 -> v2 schema evolution
COMBINE_ADD_TEST(roo_bernstein_fast_schema_evolution
    COMMAND python3 test_roo_bernstein_fast_schema_evolution.py
    COPY_TO_BUILDDIR ${REPO}/test/test_roo_bernstein_fast_schema_evolution.py
                    ${REPO}/data/ci/test_bernstein_fast_v1.root
    PASSREGEX "All checks passed"
)

# HybridNew must reach the relative accuracy (issue #1232)
COMBINE_ADD_TEST(hybridnew_accuracy
    COMMAND python3 test_hybridnew_accuracy.py
    COPY_TO_BUILDDIR ${REPO}/test/test_hybridnew_accuracy.py
                    ${REPO}/data/ci/hybridnew_small_limit_counting.txt
    PASSREGEX "All checks passed"
)

# test_interference.py
COMBINE_ADD_TEST(test_interference_py
    COMMAND python3 test_interference.py
    COPY_TO_BUILDDIR ${REPO}/test/test_interference.py
    CHECKOUT OUTPUT test_interference_py.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/test_interference_py.out
    # The output contains the best-fit values of fits to toys, so compare
    # floating point numbers with a relative tolerance of 1 %.
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "python3 ${REPO}/test/test_interference.py > references/test_interference_py.out"
)


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

# With autoMCStats, the Barlow-Beeston parameters are constant during the
# analytic minimization and floating again when the fit result is saved. Check
# that printing that fit result (done for -v 2 and above) does not crash.
COMBINE_ADD_TEST(multidimfit-verbose-analytic-bb
    COMMAND combine -M MultiDimFit ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt --algo singles -v 2
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
    # No meaningful output. We just check that it runs.
)

COMBINE_ADD_TEST(hzz4l-significance
    COMMAND combine -M Significance ${REPO}/test/inputs/hzz4l_datacard.txt
    CHECKOUT OUTPUT hzz4l-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/hzz4l-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine -M Significance ${REPO}/test/inputs/hzz4l_datacard.txt > references/hzz4l-significance.out"
)

COMBINE_ADD_TEST(LHC-limits
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01
    CHECKOUT OUTPUT LHC-limits.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-limits.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine ${REPO}/data/tutorials/CAT23001/datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01 > references/LHC-limits.out"
)

COMBINE_ADD_TEST(LHC-significance
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125
    CHECKOUT OUTPUT LHC-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125 > references/LHC-significance.out"
)

COMBINE_ADD_TEST(parametric-analysis-significance
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M Significance --mass 125
    CHECKOUT OUTPUT parametric-analysis-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric-analysis-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt -M Significance --mass 125 > references/parametric-analysis-significance.out"
)

COMBINE_ADD_TEST(MarkovChainMC
    COMMAND combine ${REPO}/data/tutorials/CAT23001/datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100
    CHECKOUT OUTPUT MarkovChainMC.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/MarkovChainMC.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine ${REPO}/data/tutorials/CAT23001/datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100 > references/MarkovChainMC.out"
)

COMBINE_ADD_TEST(multi-signal-text2workspace
    COMMAND text2workspace.py ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125
    CHECKOUT OUTPUT multi-signal-text2workspace.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal-text2workspace.out
    FIXTURES_SETUP multi-signal-text2workspace
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "text2workspace.py ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125 > references/multi-signal-text2workspace.out"
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
    "combine ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt -M ChannelCompatibilityCheck --mass 125 > references/multi-signal-ChannelCompatibilityCheck.out"
)


# AsymptoticLimits directly on a text datacard. This is deliberately close to
# the asymptotic_limits test above: the point is to exercise the path where
# combine itself converts the text datacard to a workspace on the fly, which
# is the command that motivated the pandas dependency in pixi.toml. The -n is
# only there to give the output file a unique name in the shared build
# directory; it does not change the printed output.
COMBINE_ADD_TEST(realistic-counting-experiment-AsymptoticLimits
    COMMAND combine -M AsymptoticLimits -d ${REPO}/data/tutorials/counting/realistic-counting-experiment.txt -n RealisticCountingExperiment
    CHECKOUT OUTPUT realistic-counting-experiment-AsymptoticLimits.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/realistic-counting-experiment-AsymptoticLimits.out
    # The last digits of the limits can differ between ROOT versions, so
    # compare floating point numbers with a relative tolerance of 1 %.
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine -M AsymptoticLimits -d ${REPO}/data/tutorials/counting/realistic-counting-experiment.txt -n RealisticCountingExperiment > references/realistic-counting-experiment-AsymptoticLimits.out"
)

# Check that the combine executable and the executable scripts can at least
# be loaded (i.e. all their imports work) and print their help.
foreach(executable
    combine
    combineCards.py
    combineTool.py
    commentUncerts.py
    plot1DScan.py
    plotBSMxsBRLimit.py
    plotGof.py
    plotImpacts.py
    plotLimitGrid.py
    plotLimits.py
    pruneUncerts.py
    text2workspace.py
)
    COMBINE_ADD_TEST(help-${executable}
        COMMAND ${executable} --help
        PASSREGEX "[Uu]sage:"
    )
endforeach()

# Write the script to produce the reference files. After building combine, you
# can produce the reference files as follows:
#
#   cd build
#   ctest -j$(nproc)
#   sh create_reference_files.sh
#   cp references/*.out /path/to/source/test/references/
#
# The script is placed in CMAKE_CURRENT_BINARY_DIR, which is the same
# directory where CTest generates the workspace .root files. This is
# build/test/ for non-standalone builds and build/ for standalone builds.

get_property(all_combine_commands GLOBAL PROPERTY ALL_COMBINE_COMMANDS)
message(STATUS "All combine commands:")
set(_outfile "${CMAKE_CURRENT_BINARY_DIR}/create_reference_files.sh")
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
    # Apply the same output filter as COMBINE_ADD_TEST does before comparing
    # the test output to the reference file.
    string(REGEX REPLACE " (>>?) (references/[^ ]+)$"
           " | sh ${CMAKE_CURRENT_SOURCE_DIR}/filter_test_output.sh \\1 \\2" cmd "${cmd}")
    file(APPEND "${_outfile}" "${cmd}\n")
endforeach()
message(STATUS "Wrote script to produce test reference files ${_outfile}")
