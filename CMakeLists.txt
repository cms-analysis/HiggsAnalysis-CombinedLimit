# CMakeLists.txt file for cms-analysis/HiggsAnalysis-CombinedLimit
# Authors: Jonas Rembser, CERN 2022

cmake_minimum_required(VERSION 3.12)
project(HiggsAnalysis-CombinedLimit CXX)

##########################################
# Define configuration and build targets #
##########################################

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -m64 -O2 -rdynamic -g -Wall")

find_package(ROOT REQUIRED COMPONENTS Core MathMore RooFitCore RooFit RooStats Minuit HistFactory)
include(${ROOT_USE_FILE}) # inherit ROOT's compile options, including C++ standard

find_package(Boost 1.45.0 COMPONENTS program_options filesystem system)
add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS)

find_package(Eigen3 3.0 REQUIRED)

# To get the Python library install directory into Python_SITELIB
find_package(Python COMPONENTS Interpreter Development)

set(LIBNAME HiggsAnalysisCombinedLimit)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(interface)

file(GLOB SOURCE_FILES src/*.c*)
file(GLOB HEADER_FILES interface/*.h*)

ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ../src/classes.h LINKDEF src/classes_def.xml
                         MODULE ${LIBNAME}
                         OPTIONS --deep)

add_library(${LIBNAME} SHARED ${SOURCE_FILES} G__${LIBNAME}.cxx)
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")

target_link_libraries(${LIBNAME} PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(${LIBNAME} PUBLIC ${Boost_LIBRARIES})

add_executable(combine bin/combine.cpp)
target_link_libraries(combine PUBLIC ${LIBNAME})

# Create empty __init__.py file in the build directory that will be installed
# in the Python library directories.
set(empty_init_py "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")


#####################
# Installation part #
#####################


# Install the dictionaries.
install(FILES ${PROJECT_BINARY_DIR}/lib${LIBNAME}_rdict.pcm
              ${PROJECT_BINARY_DIR}/lib${LIBNAME}.rootmap
        DESTINATION lib)

# Install the libraries and header files.
install(TARGETS ${LIBNAME}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/HiggsAnalysis/CombinedLimit/interface
)

# Install the "combine" executable in the bin directory.
install(TARGETS combine DESTINATION bin)

# Install the scripts like "text2workspace" to the bin directory.
install(DIRECTORY scripts/ DESTINATION bin)

# Check if the Python library installation directory is outside the install
# prefix. If it is, we error out because CMake should not install files outside
# the prefix. In the future, one can imagine to let the user choose where the
# Python libraries get installed in the prefix with a CMake configuration flag.
cmake_path(IS_PREFIX CMAKE_INSTALL_PREFIX "${Python_SITELIB}" sitelib_in_prefix)
if(NOT ${sitelib_in_prefix})
    message( FATAL_ERROR "Your Python library installation directory ${Python_SITELIB} "
                         "is outside the install prefix ${CMAKE_INSTALL_PREFIX}! "
                         "This is not supported for now. Consider changing the install prefix "
                         "with the -DCMAKE_INSTALL_PREFIX:PATH=<path> cmake configuration option.")
endif()

# The the Python library installation directory relative to the install prefix.
file(RELATIVE_PATH Python_SITELIB_IN_PREFIX ${CMAKE_INSTALL_PREFIX} ${Python_SITELIB})

# The python package will be installed in such a way that the original
# CMSSW-style directory structure is kept, for maximal compatibility.
install(DIRECTORY python/ DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis/CombinedLimit)

# Create empty __init__.py files in the Python package subdirectories such that
# the Python imports work.
file(TOUCH ${empty_init_py})
INSTALL(FILES ${empty_init_py} DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis/__init__.py)
INSTALL(FILES ${empty_init_py} DESTINATION ${Python_SITELIB_IN_PREFIX}/HiggsAnalysis/CombinedLimit/__init__.py)
