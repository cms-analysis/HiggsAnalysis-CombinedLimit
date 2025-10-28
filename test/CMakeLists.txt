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
