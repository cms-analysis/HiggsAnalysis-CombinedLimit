function(COMBINE_ADD_GTEST test_suite)
  cmake_parse_arguments(ARG
    ""
    "TIMEOUT;FAILREGEX"
    "COPY_TO_BUILDDIR;LIBRARIES;LABELS;INCLUDE_DIRS;ENVIRONMENT" ${ARGN})

  ROOT_GET_SOURCES(source_files . ${ARG_UNPARSED_ARGUMENTS})
  ROOT_EXECUTABLE(${test_suite} ${source_files} LIBRARIES ${ARG_LIBRARIES} GTest::gtest GTest::gmock GTest::gtest_main GTest::gmock_main)
  target_include_directories(${test_suite} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  if (ARG_INCLUDE_DIRS)
    target_include_directories(${test_suite} PRIVATE ${ARG_INCLUDE_DIRS})
  endif(ARG_INCLUDE_DIRS)

  ROOT_PATH_TO_STRING(name_with_path ${test_suite} PATH_SEPARATOR_REPLACEMENT "-")
  string(REPLACE "-test-" "-" clean_name_with_path ${name_with_path})
  ROOT_ADD_TEST(
    gtest${clean_name_with_path}
    COMMAND ${test_suite} ${extra_command}
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    COPY_TO_BUILDDIR "${ARG_COPY_TO_BUILDDIR}"
    ${willfail}
    TIMEOUT "${ARG_TIMEOUT}"
    LABELS "${ARG_LABELS}"
    FAILREGEX "${ARG_FAILREGEX}"
    ENVIRONMENT "${ARG_ENVIRONMENT}"
  )
endfunction()


# Arguments:
#   TEST_BASENAME           Base name of the test (e.g. "parametric_analysis" or "cmshistsum")
#   COPY_TO_BUILDDIR        List of required files to copy (e.g. datacard and ROOT inputs)
#   T2W_COMMAND
#   COMBINE_COMMANDS
function(ADD_COMBINE_TEST TEST_BASENAME)
    cmake_parse_arguments(ARG "" "" "COPY_TO_BUILDDIR;T2W_COMMAND;COMBINE_COMMANDS" ${ARGN})

  # --- text2workspace test ---
  ROOT_ADD_TEST(${TEST_BASENAME}-text2workspace
      COMMAND ${CMAKE_BINARY_DIR}/bin/${ARG_T2W_COMMAND}
      COPY_TO_BUILDDIR ${ARG_COPY_TO_BUILDDIR}
      FIXTURES_SETUP ${TEST_BASENAME}
      ENVIRONMENT
          LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
  )

  # Combine multiple commands into a single shell command chain
  set(_combined_command "")
  foreach(cmd IN LISTS ARG_COMBINE_COMMANDS)
      string(APPEND _combined_command "${cmd} && ")
      set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "${cmd} >> references/${TEST_BASENAME}.out")
  endforeach()
  string(REGEX REPLACE " && $" "" _combined_command "${_combined_command}")  # remove trailing &&

  # --- combine test ---
  ROOT_ADD_TEST(${TEST_BASENAME}
      COMMAND bash -c "${_combined_command}"
      FIXTURES_REQUIRED ${TEST_BASENAME} # requires corresponding text2workspace run
      # We compare the output to reference files to validate the best-fit parameter values
      WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
      CHECKOUT OUTPUT ${TEST_BASENAME}.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/${TEST_BASENAME}.out
      # For this test, we are not interested in the standard error, but you could
      # compare this too:
      # CHECKERR
      # ERROR ${TEST_BASENAME}.err
      # ERRREF ${CMAKE_CURRENT_SOURCE_DIR}/references/${TEST_BASENAME}.err
      ENVIRONMENT PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
  )
endfunction()

# Set the location of text2workspace.py, assuming installation to a virtual
# environment inside the source directory.
set(t2w ${CMAKE_BINARY_DIR}/bin/text2workspace.py)

# Counting datacard
ADD_COMBINE_TEST(counting_datacard
  COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/tutorials/multiDim/toy-hgg-125.txt
  T2W_COMMAND
    text2workspace.py toy-hgg-125.txt -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
  COMBINE_COMMANDS
    "combine -M MultiDimFit toy-hgg-125.root  --setParameterRanges r=-1,1"
)

# Counting datacard Fixed Point from csv - combineTool.py
ROOT_ADD_TEST(counting_datacard_from_csv
    COMMAND ${CMAKE_BINARY_DIR}/bin/combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv
    COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/tutorials/multiDim/fixed.csv
    FIXTURES_REQUIRED counting_datacard
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    CHECKOUT OUTPUT counting_datacard_from_csv.out OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard_from_csv.out
    ENVIRONMENT PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH} # So the combineTool.py finds the combine executable
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
)
set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "combineTool.py -M MultiDimFit toy-hgg-125.root --fromfile fixed.csv >> references/counting_datacard_from_csv.out")

# Parametric analysis
ADD_COMBINE_TEST(parametric_analysis
  COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
                   ${CMAKE_SOURCE_DIR}/data/tutorials/CAT23001/parametric-analysis-datacard-input.root
  T2W_COMMAND
    text2workspace.py datacard-3-parametric-analysis.txt --mass 125
  COMBINE_COMMANDS
    "combine -M MultiDimFit datacard-3-parametric-analysis.root --algo singles --setParameterRanges r=-2,1"
)

# TODO: make these tests work too
#
# # Template analysis CMSHistFunc
# ADD_COMBINE_TEST(cmshistfunc
#   COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shapeInterp.txt
#                    ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
#   T2W_COMMAND
#     text2workspace.py template-analysis_shapeInterp.txt -o ws_template-analysis.root --mass 200
#   COMBINE_COMMANDS
#     "combine -M MultiDimFit ws_template-analysis.root --algo singles  --setParameterRanges r=-1,1"
#     "combine -M FitDiagnostics ws_template-analysis.root  -t -1 --setParameters r=1 --plots  --setParameterRanges r=-1,1"
# )
# 
# # Template analysis CMSHistFunc shapeN
# ADD_COMBINE_TEST(cmshistfunc_shapeN
#   COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shapeNInterp.txt
#                    ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
#   T2W_COMMAND
#     text2workspace.py template-analysis_shapeNInterp.txt -o ws_template-analysis-shapeN.root --mass 200
#   COMBINE_COMMANDS
#     "combine -M MultiDimFit ws_template-analysis-shapeN.root --algo singles  --setParameterRanges r=-1,1"
#     "combine -M FitDiagnostics ws_template-analysis-shapeN.root  -t -1 --setParameters r=1 --plots  --setParameterRanges r=-1,1"
# )

# Template analysis CMSHistSum
ADD_COMBINE_TEST(cmshistsum
  COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shapeInterp.txt
                   ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
  T2W_COMMAND
    text2workspace.py template-analysis_shapeInterp.txt --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit template-analysis_shapeInterp.root --algo singles  --setParameterRanges r=-1,1 --X-rtd FAST_VERTICAL_MORPH"
)

# Template analysis CMSHistSum with shapeN
ADD_COMBINE_TEST(cmshistsum_shapeN
  COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shapeNInterp.txt
                   ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
  T2W_COMMAND
    text2workspace.py template-analysis_shapeNInterp.txt --mass 200 --for-fits --no-wrappers --use-histsum
  COMBINE_COMMANDS
    "combine -M MultiDimFit template-analysis_shapeNInterp.root --algo singles  --setParameterRanges r=-1,1 --X-rtd FAST_VERTICAL_MORPH"
)

# Template-analysis datacard -> text2workspace
ROOT_ADD_TEST(template_analysis-text2workspace
    COMMAND ${t2w} template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    COPY_TO_BUILDDIR
        ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shape_autoMCStats.txt
        ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
    ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
                PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
    FIXTURES_SETUP template_analysis_workspace
)

# Build and run the testCreateNLL helper on the generated workspace
COMBINE_ADD_GTEST(template_analysis-testCreateNLL
    testCreateNLL.cxx
    LIBRARIES HiggsAnalysisCombinedLimit
)
set_property(TEST gtest-template-analysis-testCreateNLL
    PROPERTY FIXTURES_REQUIRED template_analysis_workspace
)


# Write the script to produce the reference files. After building combine, you
# can produce the reference files as follows:
#
#   cd build
#   ctest -j8 .
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
