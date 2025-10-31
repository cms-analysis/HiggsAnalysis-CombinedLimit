#----------------------------------------------------------------------------
# Based on ROOT_ADD_TEST
# function COMBINE_ADD_TEST( <name> COMMAND cmd [arg1... ]
#                        [PRECMD cmd [arg1...]] [POSTCMD cmd [arg1...]]
#                        [OUTPUT outfile] [ERROR errfile] [INPUT infile]
#                        [ENVIRONMENT var1=val1 var2=val2 ...
#                        [DEPENDS test1 ...]
#                        [RUN_SERIAL]
#                        [TIMEOUT seconds]
#                        [DEBUG]
#                        [SOURCE_DIR dir] [BINARY_DIR dir]
#                        [WORKING_DIR dir] [COPY_TO_BUILDDIR files]
#                        [BUILD target] [PROJECT project]
#                        [PASSREGEX exp] [FAILREGEX epx]
#                        [PASSRC code]
#                        [RESOURCE_LOCK lock]
#                        [FIXTURES_SETUP ...] [FIXTURES_CLEANUP ...] [FIXTURES_REQUIRED ...]
#                        [LABELS label1 label2]
#                        [PROPERTIES prop1 value1 prop2 value2...]
#                       )
#
function(COMBINE_ADD_TEST test)
  CMAKE_PARSE_ARGUMENTS(ARG "DEBUG;WILLFAIL;CHECKOUT;CHECKERR;RUN_SERIAL"
                            "TIMEOUT;BUILD;INPUT;OUTPUT;ERROR;SOURCE_DIR;BINARY_DIR;WORKING_DIR;PROJECT;PASSRC;RESOURCE_LOCK"
                            "COMMAND;COPY_TO_BUILDDIR;DIFFCMD;OUTCNV;OUTCNVCMD;PRECMD;POSTCMD;ENVIRONMENT;DEPENDS;PASSREGEX;OUTREF;ERRREF;FAILREGEX;LABELS;PYTHON_DEPS;FIXTURES_SETUP;FIXTURES_CLEANUP;FIXTURES_REQUIRED;PROPERTIES"
                            ${ARGN})

  #- Handle COMMAND argument
  list(LENGTH ARG_COMMAND _len)
  if(_len LESS 1)
    if(NOT ARG_BUILD)
      message(FATAL_ERROR "COMBINE_ADD_TEST: command is mandatory (without build)")
    endif()
  else()
    list(GET ARG_COMMAND 0 _prg)
    list(REMOVE_AT ARG_COMMAND 0)

    if(TARGET ${_prg})                                 # if command is a target, get the actual executable
      set(_prg "$<TARGET_FILE:${_prg}>")
      set(_cmd ${_prg} ${ARG_COMMAND})
    else()
      find_program(_exe ${_prg})
      if(_exe)                                         # if the command is found in the system, use it
        set(_cmd ${_exe} ${ARG_COMMAND})
      elseif(NOT IS_ABSOLUTE ${_prg})                  # if not absolute, assume is found in current binary dir
        set(_prg ${CMAKE_CURRENT_BINARY_DIR}/${_prg})
        set(_cmd ${_prg} ${ARG_COMMAND})
      else()                                           # take as it is
        set(_cmd ${_prg} ${ARG_COMMAND})
      endif()
      unset(_exe CACHE)
    endif()

    string(REPLACE ";" "^" _cmd "${_cmd}")
  endif()

  set(_command ${CMAKE_COMMAND} -DCMD=${_cmd})

  #- Handle PRE and POST commands
  if(ARG_PRECMD)
    string(REPLACE ";" "^" _pre "${ARG_PRECMD}")
    set(_command ${_command} -DPRE=${_pre})
  endif()

  if(ARG_POSTCMD)
    string(REPLACE ";" "^" _post "${ARG_POSTCMD}")
    set(_command ${_command} -DPOST=${_post})
  endif()

  #- Handle INPUT, OUTPUT, ERROR, DEBUG arguments
  if(ARG_INPUT)
    set(_command ${_command} -DIN=${ARG_INPUT})
  endif()

  if(ARG_OUTPUT)
    set(_command ${_command} -DOUT=${ARG_OUTPUT})
  endif()

  if(ARG_OUTREF)
    set(_command ${_command} -DOUTREF=${ARG_OUTREF})
  endif()

  if(ARG_ERRREF)
    set(_command ${_command} -DERRREF=${ARG_ERRREF})
  endif()

  if(ARG_ERROR)
    set(_command ${_command} -DERR=${ARG_ERROR})
  endif()

  if(ARG_WORKING_DIR)
    set(_command ${_command} -DCWD=${ARG_WORKING_DIR})
  endif()

  if(ARG_DEBUG)
    set(_command ${_command} -DDBG=ON)
  endif()

  if(ARG_PASSRC)
    set(_command ${_command} -DRC=${ARG_PASSRC})
  endif()

  if(ARG_OUTCNVCMD)
    string(REPLACE ";" "^" _outcnvcmd "${ARG_OUTCNVCMD}")
    string(REPLACE "=" "@" _outcnvcmd "${_outcnvcmd}")
    set(_command ${_command} -DCNVCMD=${_outcnvcmd})
  endif()

  if(ARG_OUTCNV)
    string(REPLACE ";" "^" _outcnv "${ARG_OUTCNV}")
    set(_command ${_command} -DCNV=${_outcnv})
  endif()

  if(ARG_DIFFCMD)
    string(REPLACE ";" "^" _diff_cmd "${ARG_DIFFCMD}")
    set(_command ${_command} -DDIFFCMD=${_diff_cmd})
  endif()

  if(ARG_CHECKOUT)
    set(_command ${_command} -DCHECKOUT=true)
  endif()

  if(ARG_CHECKERR)
    set(_command ${_command} -DCHECKERR=true)
  endif()

  set(_command ${_command} -DSYS=${ROOTSYS})

  #- Handle ENVIRONMENT argument
  if(ARG_ENVIRONMENT)
    string(REPLACE ";" "#" _env "${ARG_ENVIRONMENT}")
    set(_command ${_command} -DENV=${_env})
  endif()

  #- Copy files to the build directory.
  if(ARG_COPY_TO_BUILDDIR)
    string(REPLACE ";" "^" _copy_files "${ARG_COPY_TO_BUILDDIR}")
    set(_command ${_command} -DCOPY=${_copy_files})
  endif()

  set(_command ${_command} -P ${CMAKE_CURRENT_SOURCE_DIR}/CombineTestDriver.cmake)

  if(ARG_WILLFAIL)
    set(test ${test}_WILL_FAIL)
  endif()

  #- Now we can actually add the test
  if(ARG_BUILD)
    if(NOT ARG_SOURCE_DIR)
      set(ARG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    if(NOT ARG_BINARY_DIR)
      set(ARG_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
    endif()
    if(NOT ARG_PROJECT)
       if(NOT PROJECT_NAME STREQUAL "ROOT")
         set(ARG_PROJECT ${PROJECT_NAME})
       else()
         set(ARG_PROJECT ${ARG_BUILD})
       endif()
    endif()
    add_test(NAME ${test} COMMAND ${CMAKE_CTEST_COMMAND}
      --build-and-test  ${ARG_SOURCE_DIR} ${ARG_BINARY_DIR}
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${ARG_BUILD}
      --build-project ${ARG_PROJECT}
      --build-config $<CONFIGURATION>
      --build-noclean
      --test-command ${_command} )
  else()
    add_test(NAME ${test} COMMAND ${_command})
  endif()

  #- provided fixtures and resource lock are set here
  if (ARG_FIXTURES_SETUP)
    set_property(TEST ${test} PROPERTY
      FIXTURES_SETUP ${ARG_FIXTURES_SETUP})
  endif()

  if (ARG_FIXTURES_CLEANUP)
    set_property(TEST ${test} PROPERTY
      FIXTURES_CLEANUP ${ARG_FIXTURES_CLEANUP})
  endif()

  if (ARG_FIXTURES_REQUIRED)
    set_property(TEST ${test} PROPERTY
      FIXTURES_REQUIRED ${ARG_FIXTURES_REQUIRED})
  endif()

  if (ARG_RESOURCE_LOCK)
    set_property(TEST ${test} PROPERTY
      RESOURCE_LOCK ${ARG_RESOURCE_LOCK})
  endif()

  #- Handle TIMEOUT and DEPENDS arguments
  if(ARG_TIMEOUT)
    set_property(TEST ${test} PROPERTY TIMEOUT ${ARG_TIMEOUT})
  endif()

  if(ARG_DEPENDS)
    set_property(TEST ${test} PROPERTY DEPENDS ${ARG_DEPENDS})
  endif()

  if(ARG_PASSREGEX)
    set_property(TEST ${test} PROPERTY PASS_REGULAR_EXPRESSION ${ARG_PASSREGEX})
  endif()

  if(ARG_FAILREGEX)
    set_property(TEST ${test} PROPERTY FAIL_REGULAR_EXPRESSION ${ARG_FAILREGEX})
  endif()

  if(ARG_WILLFAIL)
    set_property(TEST ${test} PROPERTY WILL_FAIL true)
  endif()

  if(ARG_LABELS)
    set_tests_properties(${test} PROPERTIES LABELS "${ARG_LABELS}")
  endif()

  if(ARG_RUN_SERIAL)
    set_property(TEST ${test} PROPERTY RUN_SERIAL true)
  endif()

  # Pass PROPERTIES argument to the set_tests_properties as-is
  if(ARG_PROPERTIES)
    set_tests_properties(${test} PROPERTIES ${ARG_PROPERTIES})
  endif()

endfunction()

# Based on ROOT_ADD_GTEST
function(COMBINE_ADD_GTEST test_suite)
  cmake_parse_arguments(ARG
    ""
    "TIMEOUT;FAILREGEX"
    "COPY_TO_BUILDDIR;LIBRARIES;LABELS;ENVIRONMENT" ${ARGN})

  set(source_files ${ARG_UNPARSED_ARGUMENTS})
  add_executable(${test_suite} ${source_files})
  target_link_libraries(${test_suite} PRIVATE ${ARG_LIBRARIES} GTest::gtest GTest::gmock GTest::gtest_main GTest::gmock_main)
  target_include_directories(${test_suite} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

  COMBINE_ADD_TEST(
    gtest-${test_suite}
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
  COMBINE_ADD_TEST(${TEST_BASENAME}-text2workspace
      COMMAND ${ARG_T2W_COMMAND}
      COPY_TO_BUILDDIR ${ARG_COPY_TO_BUILDDIR}
      FIXTURES_SETUP ${TEST_BASENAME}
      ENVIRONMENT
          LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
          PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}
          PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}
  )

  # Combine multiple commands into a single shell command chain
  set(_combined_command "")
  foreach(cmd IN LISTS ARG_COMBINE_COMMANDS)
      string(APPEND _combined_command "${cmd} && ")
      set_property(GLOBAL APPEND PROPERTY ALL_COMBINE_COMMANDS "${cmd} >> references/${TEST_BASENAME}.out")
  endforeach()
  string(REGEX REPLACE " && $" "" _combined_command "${_combined_command}")  # remove trailing &&

  # --- combine test ---
  COMBINE_ADD_TEST(${TEST_BASENAME}
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
