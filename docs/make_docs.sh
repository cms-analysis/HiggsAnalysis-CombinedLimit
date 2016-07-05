#!/bin/bash -e

doxygen Doxyfile
# We need a custom page resize function
# because of the modifications to the
# page layout
cp docs/resize.js docs/html/
