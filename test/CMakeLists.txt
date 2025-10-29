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

# Set the location of text2workspace.py, assuming installation to a virtual
# environment inside the source directory.
set(t2w ${CMAKE_SOURCE_DIR}/virtualenv/bin/text2workspace.py)

# Counting datacard - text2workspace
ROOT_ADD_TEST(counting_datacard-t2w
    COMMAND ${t2w} toy-hgg-125.txt -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingXSHiggs --PO modes=ggH,qqH
    COPY_TO_BUILDDIR ${CMAKE_SOURCE_DIR}/data/tutorials/multiDim/toy-hgg-125.txt
    FIXTURES_SETUP counting_datacard
    # Append the build directory to the library Path to we can use Combine classes
    ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}
)

# Counting datacard - combine
ROOT_ADD_TEST(counting_datacard
    COMMAND combine -M MultiDimFit toy-hgg-125.root  --setParameterRanges r=-1,1
    FIXTURES_REQUIRED counting_datacard # requires corresponding text2workspace run
    # We compare the output to reference files to validate the best-fit parameter values
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    CHECKOUT
    OUTPUT counting_datacard.out
    OUTREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard.out
    # For this test, we are not interested in the standard error, but you could
    # compare this too:
    # CHECKERR
    # ERROR counting_datacard.err
    # ERRREF ${CMAKE_CURRENT_SOURCE_DIR}/references/counting_datacard.err
)

# Template-analysis datacard -> text2workspace
ROOT_ADD_TEST(template_analysis-t2w
    COMMAND ${t2w} template-analysis_shape_autoMCStats.txt -o template-analysis_shape_autoMCStats.root
    COPY_TO_BUILDDIR
        ${CMAKE_SOURCE_DIR}/data/ci/template-analysis_shape_autoMCStats.txt
        ${CMAKE_SOURCE_DIR}/data/ci/htt_input.root
    WORKING_DIR ${CMAKE_CURRENT_BINARY_DIR}
    ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}
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
