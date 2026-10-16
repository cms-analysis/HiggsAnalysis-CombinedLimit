#!/bin/sh
# Filter applied to the standard output of every test that is compared against
# a reference file in test/references/, both when running the test and when
# (re)generating the references with create_reference_files.sh. It removes lines
# that are expected to change between Combine versions but that carry no
# information about the correctness of the result, so the reference files don't
# have to be touched on every release.
#
# Removed lines:
#  - the Combine version banner
#  - the RooFit timing message for the NLL creation
#  - the printout of a RooAbsReal (e.g. the NLL object in FastScan): its value
#    is the offset NLL, i.e. numerical noise around zero, and the parameter
#    order depends on RooFit internals
#  - the RooFit numerical derivator error, whose presence and wording depend
#    on the ROOT version
sed -e '/^ <<< v[^ ]* >>>$/d' \
    -e '/Creation of NLL object took/d' \
    -e '/^RooAbsReal::/d' \
    -e '/error in numerical derivator/d'
