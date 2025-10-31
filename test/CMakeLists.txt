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
    CHECKOUT OUTPUT counting_datacard_from_csv.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard_from_csv.out
    ENVIRONMENT PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH} # So the combineTool.py finds the combine executable
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv >> references/counting_datacard_from_csv.out")

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
