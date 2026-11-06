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

# Pre-fit nuisance parameter values and uncertainties from FitDiagnostics, for
# a model with a flat (lnU) constraint (issue #600)
ADD_COMBINE_TEST(fitdiagnostics_prefit_nuisances
  T2W_COMMAND
    text2workspace.py ${REPO}/data/ci/counting_lnU.txt -o ws_counting_lnU.root
  COMBINE_COMMANDS
    "combine -M FitDiagnostics ws_counting_lnU.root -n _counting_lnU"
    "python3 ${REPO}/test/checkPrefitNuisances.py fitDiagnostics_counting_lnU.root lumi=1 bkg2N=4 bkgU=0"
  # The last digits of the fitted values differ between platforms.
  DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
)

# ModelTools.py
COMBINE_ADD_TEST(test_model_tools_py
    COMMAND python3 ${REPO}/test/test_model_tools.py
)

# CombineHarvester python interface and datacard validation
COMBINE_ADD_TEST(test_combineharvester_py
    COMMAND python3 ${REPO}/test/test_combineharvester.py ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt
    PASSREGEX "All checks passed"
)
COMBINE_ADD_TEST(test_combineharvester_shapes_py
    COMMAND python3 ${REPO}/test/test_combineharvester_shapes.py ${REPO}/data/ci/template-analysis_shapeInterp.txt
    PASSREGEX "All checks passed"
)
COMBINE_ADD_TEST(test_combineharvester_gmn_py
    COMMAND python3 ${REPO}/test/test_combineharvester_gmn.py ${REPO}/data/ci/counting_gmN.txt
    PASSREGEX "All checks passed"
)
COMBINE_ADD_TEST(test_combineharvester_directives_py
    COMMAND python3 ${REPO}/test/test_combineharvester_directives.py
    PASSREGEX "All checks passed"
)
COMBINE_ADD_TEST(combineharvester-validate-datacards
    COMMAND ValidateDatacards.py ${REPO}/data/ci/template-analysis_shape_autoMCStats.txt --jsonFile validation.json
    PASSREGEX "Validation results"
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

# KS and AD goodness of fit must match each category's data, model and output branch
COMBINE_ADD_TEST(gof_ks_ad_categories
    COMMAND python3 test_gof_ks_ad_categories.py
    COPY_TO_BUILDDIR ${REPO}/test/test_gof_ks_ad_categories.py
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

# CombineHarvester round trip
#
# The "-ch" tests repeat the tests above on the datacards after CH has parsed
# and written them, and compare the output to the same reference files. The
# ch-write-* tests write the datacards to ch_roundtrip/ch_cards/ with
# ch_roundtrip.py. All of them run in ch_roundtrip/ and have the CHRoundTrip
# label. LHC-limits and LHC-significance are not repeated: the parameters are
# in another order in the datacards written by CH, so the toys differ and the
# results would only agree within their uncertainties.
set(ch_dir ${CMAKE_CURRENT_BINARY_DIR}/ch_roundtrip)
file(MAKE_DIRECTORY ${ch_dir})
set(ch_roundtrip python3 ${CMAKE_CURRENT_SOURCE_DIR}/ch_roundtrip.py)

# Datacards without $MASS in their shapes lines
foreach(card
    data/tutorials/counting/realistic-counting-experiment.txt
    data/tutorials/counting/simple-counting-experiment.txt
    data/tutorials/multiDim/toy-hgg-125.txt
    data/tutorials/shapes/simple-shapes-TH1.txt
    data/tutorials/CAT23001/datacard-1-counting-experiment.txt
    data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
    data/tutorials/CAT23001/datacard-5-multi-signal.txt
    data/ci/templ_datacard_largeYields.txt
    data/ci/datacard_RooParametricHist.txt
    data/ci/template-analysis_shape_autoMCStats.txt
    test/inputs/hzz4l_datacard.txt
)
    get_filename_component(name ${card} NAME)
    string(REGEX REPLACE "\\.txt$" "" name ${name})
    COMBINE_ADD_TEST(ch-write-${name}
        COMMAND ${ch_roundtrip} ${REPO}/${card} ch_cards/${name}.txt
        WORKING_DIR ${ch_dir}
        FIXTURES_SETUP ch-card-${name}
        LABELS CHRoundTrip
    )
endforeach()
# CH can't read gzipped datacards, so these are unzipped to ch_roundtrip/
# first, with their shapes files next to them
foreach(name RooHistPdf RooMultiPdf)
    COMBINE_ADD_TEST(ch-write-datacard_${name}
        PRECMD sh -c "cp ${REPO}/data/ci/shapes_${name}.root . && gzip -dc ${REPO}/data/ci/datacard_${name}.txt.gz > datacard_${name}.txt"
        COMMAND ${ch_roundtrip} datacard_${name}.txt ch_cards/datacard_${name}.txt
        WORKING_DIR ${ch_dir}
        FIXTURES_SETUP ch-card-datacard_${name}
        LABELS CHRoundTrip
    )
endforeach()
# Datacards with $MASS in their shapes lines, written for the mass of the tests
COMBINE_ADD_TEST(ch-write-template-analysis_shapeInterp
    COMMAND ${ch_roundtrip} --mass=200 ${REPO}/data/ci/template-analysis_shapeInterp.txt ch_cards/template-analysis_shapeInterp.txt
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP ch-card-template-analysis_shapeInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(ch-write-template-analysis_shapeNInterp
    COMMAND ${ch_roundtrip} --mass=200 ${REPO}/data/ci/template-analysis_shapeNInterp.txt ch_cards/template-analysis_shapeNInterp.txt
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP ch-card-template-analysis_shapeNInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(ch-write-htt_multiple_regions
    COMMAND ${ch_roundtrip} --mass=125 ${REPO}/data/ci/htt_multiple_regions.txt ch_cards/htt_multiple_regions.txt
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP ch-card-htt_multiple_regions
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(asymptotic_limits-ch-text2workspace
    COMMAND text2workspace.py ch_cards/realistic-counting-experiment.txt -o asymptotic_limits.root -m 125
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP asymptotic_limits-ch
    FIXTURES_REQUIRED ch-card-realistic-counting-experiment
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(asymptotic_limits-ch
    COMMAND bash -c "combine -M AsymptoticLimits asymptotic_limits.root -m 125 -n AsymptoticLimitsRegression"
    FIXTURES_REQUIRED asymptotic_limits-ch
    CHECKOUT OUTPUT asymptotic_limits-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/asymptotic_limits.out
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.001
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(counting_datacard-ch-text2workspace
    COMMAND text2workspace.py ch_cards/toy-hgg-125.txt -o counting_datacard.root -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP counting_datacard-ch
    FIXTURES_REQUIRED ch-card-toy-hgg-125
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(counting_datacard-ch
    COMMAND bash -c "combine -M MultiDimFit counting_datacard.root  --setParameterRanges r=-5,5"
    FIXTURES_REQUIRED counting_datacard-ch
    CHECKOUT OUTPUT counting_datacard-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(counting_datacard_from_csv-ch
    COMMAND combineTool.py -M MultiDimFit counting_datacard.root --fromfile fixed.csv
    PRECMD ${CMAKE_COMMAND} -E copy ${REPO}/data/tutorials/multiDim/fixed.csv .
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED counting_datacard-ch
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(parametric_analysis-ch-text2workspace
    COMMAND text2workspace.py ch_cards/datacard-3-parametric-analysis.txt -o datacard-3-parametric-analysis.root --mass 125
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP parametric_analysis-ch
    FIXTURES_REQUIRED ch-card-datacard-3-parametric-analysis
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(parametric_analysis-ch
    COMMAND bash -c "combine -M MultiDimFit datacard-3-parametric-analysis.root --algo singles --setParameterRanges r=-2,1"
    FIXTURES_REQUIRED parametric_analysis-ch
    CHECKOUT OUTPUT parametric_analysis-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric_analysis.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

if(NOT CMSSW_VERSION STREQUAL "CMSSW_11_3_4")
  COMBINE_ADD_TEST(histpdf-ch-text2workspace
      COMMAND text2workspace.py ch_cards/datacard_RooHistPdf.txt -o ws_RooHistPdf.root
      WORKING_DIR ${ch_dir}
      FIXTURES_SETUP histpdf-ch
      FIXTURES_REQUIRED ch-card-datacard_RooHistPdf
      LABELS CHRoundTrip
  )
  COMBINE_ADD_TEST(histpdf-ch
      COMMAND bash -c "combine -M MultiDimFit ws_RooHistPdf.root --algo singles -v -2  --setParameterRanges r=-1,2"
      FIXTURES_REQUIRED histpdf-ch
      CHECKOUT OUTPUT histpdf-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/histpdf.out
      WORKING_DIR ${ch_dir}
      LABELS CHRoundTrip
  )
endif()

if(NOT CMSSW_VERSION MATCHES "^CMSSW_11")
  COMBINE_ADD_TEST(multipdf-ch-text2workspace
      COMMAND text2workspace.py ch_cards/datacard_RooMultiPdf.txt -o ws_RooMultiPdf.root
      WORKING_DIR ${ch_dir}
      FIXTURES_SETUP multipdf-ch
      FIXTURES_REQUIRED ch-card-datacard_RooMultiPdf
      LABELS CHRoundTrip
  )
  string(JOIN " && " _multipdf_ch_command
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=2 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles ws_RooMultiPdf.root"
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=2 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles --saveFitResult ws_RooMultiPdf.root"
      "combine -v -1 -M MultiDimFit -m 125.38 --setParameters pdf_index_ggh=0 --freezeParameters MH --cminDefaultMinimizerTolerance 0.001 --cminDefaultMinimizerStrategy 1 --X-rtd FAST_VERTICAL_MORPH --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_maskChannels=2 --algo singles ws_RooMultiPdf.root"
  )
  COMBINE_ADD_TEST(multipdf-ch
      COMMAND bash -c "${_multipdf_ch_command}"
      FIXTURES_REQUIRED multipdf-ch
      CHECKOUT OUTPUT multipdf-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multipdf.out
      DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
      WORKING_DIR ${ch_dir}
      LABELS CHRoundTrip
  )
endif()

COMBINE_ADD_TEST(cmshistfunc-ch-text2workspace
    COMMAND text2workspace.py ch_cards/template-analysis_shapeInterp.txt -o ws_template-analysis.root --mass 200
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistfunc-ch
    FIXTURES_REQUIRED ch-card-template-analysis_shapeInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistfunc-ch
    COMMAND bash -c "combine -M MultiDimFit ws_template-analysis.root --algo singles  --setParameterRanges r=-5,5"
    FIXTURES_REQUIRED cmshistfunc-ch
    CHECKOUT OUTPUT cmshistfunc-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistfunc.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistfunc-FastScan-ch
    COMMAND bash -c "${_fastscan_command}"
    FIXTURES_REQUIRED cmshistfunc-ch
    CHECKOUT OUTPUT cmshistfunc-FastScan-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistfunc-FastScan.out
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01 --atol 0.1
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistfunc_channel_masks-ch-text2workspace
    COMMAND text2workspace.py ch_cards/htt_multiple_regions.txt  -o ws_template-analysis-channel_masks.root --mass 125 --channel-masks
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistfunc_channel_masks-ch
    FIXTURES_REQUIRED ch-card-htt_multiple_regions
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistfunc_channel_masks-ch
    COMMAND bash -c "combine -M MultiDimFit ws_template-analysis-channel_masks.root --algo singles  --setParameterRanges r=-5,5 --setParameters mask_htt_tt_2_8TeV=1"
    FIXTURES_REQUIRED cmshistfunc_channel_masks-ch
    CHECKOUT OUTPUT cmshistfunc_channel_masks-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistfunc_channel_masks.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistfunc_shapeN-ch-text2workspace
    COMMAND text2workspace.py ch_cards/template-analysis_shapeNInterp.txt -o ws_template-analysis-shapeN.root --mass 200
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistfunc_shapeN-ch
    FIXTURES_REQUIRED ch-card-template-analysis_shapeNInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistfunc_shapeN-ch
    COMMAND bash -c "combine -M MultiDimFit ws_template-analysis-shapeN.root --algo singles  --setParameterRanges r=-5,5"
    FIXTURES_REQUIRED cmshistfunc_shapeN-ch
    CHECKOUT OUTPUT cmshistfunc_shapeN-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistfunc_shapeN.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistsum-ch-text2workspace
    COMMAND text2workspace.py ch_cards/template-analysis_shapeInterp.txt -o cmshistsum.root --mass 200 --for-fits --no-wrappers --use-histsum
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistsum-ch
    FIXTURES_REQUIRED ch-card-template-analysis_shapeInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistsum-ch
    COMMAND bash -c "combine -M MultiDimFit cmshistsum.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
    FIXTURES_REQUIRED cmshistsum-ch
    CHECKOUT OUTPUT cmshistsum-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistsum.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistsum_shapeN-ch-text2workspace
    COMMAND text2workspace.py ch_cards/template-analysis_shapeNInterp.txt -o cmshistsum_shapeN.root --mass 200 --for-fits --no-wrappers --use-histsum
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistsum_shapeN-ch
    FIXTURES_REQUIRED ch-card-template-analysis_shapeNInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistsum_shapeN-ch
    COMMAND bash -c "combine -M MultiDimFit cmshistsum_shapeN.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
    FIXTURES_REQUIRED cmshistsum_shapeN-ch
    CHECKOUT OUTPUT cmshistsum_shapeN-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistsum_shapeN.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(template_analysis_large_integrals-ch-text2workspace
    COMMAND text2workspace.py ch_cards/templ_datacard_largeYields.txt -o template_analysis_large_integrals.root
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP template_analysis_large_integrals-ch
    FIXTURES_REQUIRED ch-card-templ_datacard_largeYields
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(template_analysis_large_integrals-ch
    COMMAND bash -c "combine -M MultiDimFit template_analysis_large_integrals.root --algo singles  --setParameterRanges r=-5,5"
    FIXTURES_REQUIRED template_analysis_large_integrals-ch
    CHECKOUT OUTPUT template_analysis_large_integrals-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/template_analysis_large_integrals.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(template_analysis_large_integrals_cmshistsum-ch-text2workspace
    COMMAND text2workspace.py ch_cards/templ_datacard_largeYields.txt -o template_analysis_large_integrals_cmshistsum.root  --for-fits --no-wrappers --use-histsum
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP template_analysis_large_integrals_cmshistsum-ch
    FIXTURES_REQUIRED ch-card-templ_datacard_largeYields
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(template_analysis_large_integrals_cmshistsum-ch
    COMMAND bash -c "combine -M MultiDimFit template_analysis_large_integrals_cmshistsum.root --algo singles  --setParameterRanges r=-5,5 --X-rtd FAST_VERTICAL_MORPH"
    FIXTURES_REQUIRED template_analysis_large_integrals_cmshistsum-ch
    CHECKOUT OUTPUT template_analysis_large_integrals_cmshistsum-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/template_analysis_large_integrals_cmshistsum.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

# The extArg lines still read shapes_RooParametricHist.root in the working
# directory
COMBINE_ADD_TEST(parametric_hist-ch-text2workspace
    COMMAND text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO map=.*/*hcc*:r[1,-500,500] --PO map=.*/zcc:z[1,-5,5] ch_cards/datacard_RooParametricHist.txt -o ws_RooParametricHist.root
    PRECMD ${CMAKE_COMMAND} -E copy ${REPO}/data/ci/shapes_RooParametricHist.root .
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP parametric_hist-ch
    FIXTURES_REQUIRED ch-card-datacard_RooParametricHist
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(parametric_hist-ch
    COMMAND bash -c "combine -M MultiDimFit ws_RooParametricHist.root --algo singles"
    FIXTURES_REQUIRED parametric_hist-ch
    CHECKOUT OUTPUT parametric_hist-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric_hist.out
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistsum_channel_masks-ch-text2workspace
    COMMAND text2workspace.py ch_cards/htt_multiple_regions.txt  -o ws_template-analysis-channel_masks.root --mass 125 --channel-masks --for-fits --no-wrappers --use-histsum
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistsum_channel_masks-ch
    FIXTURES_REQUIRED ch-card-htt_multiple_regions
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistsum_channel_masks-ch
    COMMAND bash -c "combine -M MultiDimFit ws_template-analysis-channel_masks.root --algo singles  --setParameterRanges r=-5,5 --setParameters mask_htt_tt_2_8TeV=1 --X-rtd FAST_VERTICAL_MORPH"
    FIXTURES_REQUIRED cmshistsum_channel_masks-ch
    CHECKOUT OUTPUT cmshistsum_channel_masks-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/cmshistsum_channel_masks.out
    WORKING_DIR ${ch_dir}
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(cmshistsum_vs_cmshistfunc_fitdiagnostics-text2workspace-ch
    COMMAND bash -c "text2workspace.py ch_cards/template-analysis_shapeInterp.txt -o fitdiag_cmshistfunc.root --mass 200 && text2workspace.py ch_cards/template-analysis_shapeInterp.txt -o fitdiag_cmshistsum.root --mass 200 --use-histsum"
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP cmshistsum_vs_cmshistfunc_fitdiagnostics-ch
    FIXTURES_REQUIRED ch-card-template-analysis_shapeInterp
    LABELS CHRoundTrip
)
COMBINE_ADD_TEST(cmshistsum_vs_cmshistfunc_fitdiagnostics-ch
    COMMAND bash -c "combine -M FitDiagnostics fitdiag_cmshistfunc.root -m 200 -n cmshistfunc --saveShapes --saveWithUncertainties --numToysForShapes 20 && combine -M FitDiagnostics fitdiag_cmshistsum.root -m 200 -n cmshistsum --saveShapes --saveWithUncertainties --numToysForShapes 20 && python3 ${REPO}/test/printFitDiagnosticsShapes.py fitDiagnosticscmshistfunc.root > shapes_cmshistfunc.txt && python3 ${REPO}/test/printFitDiagnosticsShapes.py fitDiagnosticscmshistsum.root > shapes_cmshistsum.txt && python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 1e-5 shapes_cmshistfunc.txt shapes_cmshistsum.txt"
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED cmshistsum_vs_cmshistfunc_fitdiagnostics-ch
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(template_analysis-text2workspace-ch
    COMMAND text2workspace.py ch_cards/template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-template-analysis_shape_autoMCStats
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(simple-counting-experiment-text2workspace-ch
    COMMAND text2workspace.py ch_cards/simple-counting-experiment.txt -o simple-counting-experiment.root
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-simple-counting-experiment
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(setParameters-out-of-range-ch
    COMMAND combine -M MultiDimFit ch_cards/simple-counting-experiment.txt --setParameters r=-1.05
    PASSRC 1
    PASSREGEX "outside its range"
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-simple-counting-experiment
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(multidimfit-verbose-analytic-bb-ch
    COMMAND combine -M MultiDimFit ch_cards/template-analysis_shape_autoMCStats.txt --algo singles -v 2
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-template-analysis_shape_autoMCStats
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(setParameters-out-of-range-regex-ch
    COMMAND combine -M MultiDimFit ch_cards/simple-counting-experiment.txt --setParameters "rgx{r}=-1.05"
    PASSRC 1
    PASSREGEX "outside its range"
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-simple-counting-experiment
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(simple-shapes-TH1-text2workspace-ch
    COMMAND text2workspace.py ch_cards/simple-shapes-TH1.txt -o simple-shapes-TH1.root
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-simple-shapes-TH1
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(LHC-Higgs-coupling-models-text2workspace-ch
    COMMAND text2workspace.py ch_cards/toy-hgg-125.txt -o toy-hgg-125.root -P HiggsAnalysis.CombinedLimit.LHCHCGModels:K3 -m 125
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-toy-hgg-125
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(hzz4l-significance-ch
    COMMAND combine -M Significance ch_cards/hzz4l_datacard.txt
    CHECKOUT OUTPUT hzz4l-significance-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/hzz4l-significance.out
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-hzz4l_datacard
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(parametric-analysis-significance-ch
    COMMAND combine ch_cards/datacard-3-parametric-analysis.txt -M Significance --mass 125
    CHECKOUT OUTPUT parametric-analysis-significance-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric-analysis-significance.out
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-datacard-3-parametric-analysis
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(MarkovChainMC-ch
    COMMAND combine ch_cards/datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100
    CHECKOUT OUTPUT MarkovChainMC-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/MarkovChainMC.out
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-datacard-1-counting-experiment
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(multi-signal-text2workspace-ch
    COMMAND text2workspace.py ch_cards/datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125
    CHECKOUT OUTPUT multi-signal-text2workspace-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal-text2workspace.out
    WORKING_DIR ${ch_dir}
    FIXTURES_SETUP multi-signal-text2workspace-ch
    FIXTURES_REQUIRED ch-card-datacard-5-multi-signal
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(multi-signal-ch
    COMMAND combine datacard-5-multi-signal.root -M MultiDimFit --algo singles --mass 125
    CHECKOUT OUTPUT multi-signal-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal.out
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED multi-signal-text2workspace-ch
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(multi-signal-ChannelCompatibilityCheck-ch
    COMMAND combine ch_cards/datacard-5-multi-signal.txt -M ChannelCompatibilityCheck --mass 125
    CHECKOUT OUTPUT multi-signal-ChannelCompatibilityCheck-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/multi-signal-ChannelCompatibilityCheck.out
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED multi-signal-ch ch-card-datacard-5-multi-signal
    LABELS CHRoundTrip
)

COMBINE_ADD_TEST(realistic-counting-experiment-AsymptoticLimits-ch
    COMMAND combine -M AsymptoticLimits -d ch_cards/realistic-counting-experiment.txt -n RealisticCountingExperiment
    CHECKOUT OUTPUT realistic-counting-experiment-AsymptoticLimits-ch.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/realistic-counting-experiment-AsymptoticLimits.out
    DIFFCMD python3 ${CMAKE_CURRENT_SOURCE_DIR}/numdiff.py --rtol 0.01
    WORKING_DIR ${ch_dir}
    FIXTURES_REQUIRED ch-card-realistic-counting-experiment
    LABELS CHRoundTrip
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
