cmake_minimum_required(VERSION 3.11)

# Detect whether this is being built standalone
if(NOT PROJECT_NAME)
    project(HiggsAnalysisCombinedLimit LANGUAGES CXX)
    set(REPO "${CMAKE_SOURCE_DIR}/..")
    set(standalone_tests 1)
else()
    set(REPO ${CMAKE_SOURCE_DIR})
endif()

if(NOT DEFINED standalone_tests)
  find_package(GTest 1.10 REQUIRED)
endif()

enable_testing()

include(CombineTestMacros.cmake)

# Set the location of text2workspace.py, assuming installation to a virtual
# environment inside the source directory.
set(t2w ${CMAKE_BINARY_DIR}/bin/text2workspace.py)

# Counting datacard
ADD_COMBINE_TEST(counting_datacard
  COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/toy-hgg-125.txt
  T2W_COMMAND
    text2workspace.py toy-hgg-125.txt -o counting_datacard.root -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
  COMBINE_COMMANDS
    "combine -M MultiDimFit counting_datacard.root  --setParameterRanges r=-1,1"
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
    FIXTURES_SETUP template_analysis_workspace
)


if(NOT DEFINED standalone_tests)
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
    COMMAND text2workspace.py simple-counting-experiment.txt
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/counting/simple-counting-experiment.txt
    # No output, so nothing to compare. We just check that it runs.
)

COMBINE_ADD_TEST(simple-shapes-TH1-text2workspace
    COMMAND text2workspace.py simple-shapes-TH1.txt
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/shapes/simple-shapes-TH1.txt
                     ${REPO}/data/tutorials/shapes/simple-shapes-TH1_input.root
    # No output, so nothing to compare. We just check that it runs.
)

COMBINE_ADD_TEST(LHC-Higgs-coupling-models-text2workspace
    COMMAND text2workspace.py toy-hgg-125.txt -P HiggsAnalysis.CombinedLimit.LHCHCGModels:K3 -m 125
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/multiDim/toy-hgg-125.txt
    CHECKOUT OUTPUT LHC-Higgs-coupling-models-text2workspace.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-Higgs-coupling-models-text2workspace.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "text2workspace.py toy-hgg-125.txt -P HiggsAnalysis.CombinedLimit.LHCHCGModels:K3 -m 125 > references/LHC-Higgs-coupling-models-text2workspace.out"
)

COMBINE_ADD_TEST(hzz4l-significance
    COMMAND combine -M Significance hzz4l_datacard.txt
    COPY_TO_BUILDDIR ${REPO}/test/inputs/hzz4l_datacard.txt
                     ${REPO}/test/inputs/Workspace_m4l.root
    CHECKOUT OUTPUT hzz4l-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/hzz4l-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine -M Significance hzz4l_datacard.txt > references/hzz4l-significance.out"
)

COMBINE_ADD_TEST(LHC-limits
    COMMAND combine datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-2-template-analysis.txt
                     ${REPO}/data/tutorials/CAT23001/template-analysis-datacard-input.root
    # TODO: fix this test with ROOT master. Right now, it throws an exception because of variable clipping.
    # CHECKOUT OUTPUT LHC-limits.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-limits.out
)
# set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
#     "combine datacard-2-template-analysis.txt -M HybridNew --LHCmode LHC-limits --rMax 2.0 --clsAcc 0.01 > references/LHC-limits.out"
# )

COMBINE_ADD_TEST(LHC-significance
    COMMAND combine datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
    CHECKOUT OUTPUT LHC-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/LHC-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-3-parametric-analysis.txt -M HybridNew --LHCmode LHC-significance -T 5000 --mass 125 > references/LHC-significance.out"
)

COMBINE_ADD_TEST(parametric-analysis-significance
    COMMAND combine datacard-3-parametric-analysis.txt -M Significance --mass 125
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
    CHECKOUT OUTPUT parametric-analysis-significance.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/parametric-analysis-significance.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-3-parametric-analysis.txt -M Significance --mass 125 > references/parametric-analysis-significance.out"
)

COMBINE_ADD_TEST(MarkovChainMC
    COMMAND combine datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-1-counting-experiment.txt
    CHECKOUT OUTPUT MarkovChainMC.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/MarkovChainMC.out
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS
    "combine datacard-1-counting-experiment.txt -M MarkovChainMC --tries 100 > references/MarkovChainMC.out"
)

COMBINE_ADD_TEST(multi-signal-text2workspace
    COMMAND text2workspace.py datacard-5-multi-signal.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH -o datacard-5-multi-signal.root --mass 125
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt
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
    COMMAND combine datacard-5-multi-signal.txt -M ChannelCompatibilityCheck --mass 125
    COPY_TO_BUILDDIR ${REPO}/data/tutorials/CAT23001/datacard-5-multi-signal.txt
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
