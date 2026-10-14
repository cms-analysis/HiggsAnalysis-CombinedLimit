#!/bin/sh
# Filter applied to the standard output of every test that is compared against
# a reference file in test/references/, both when running the test and when
# (re)generating the references with create_reference_files.sh. It removes lines
# that are expected to change between Combine versions but that carry no
# information about the correctness of the result, so the reference files don't
# have to be touched on every release.
sed -e '/^ <<< v[^ ]* >>>$/d'
