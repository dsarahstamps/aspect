# CMake generated Testfile for 
# Source directory: /Users/q/Documents/aspect/tests
# Build directory: /Users/q/Documents/aspect/build_deal.ii.app/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(quick_mpi "/Applications/CMake.app/Contents/bin/cmake" "-DBINARY_DIR=/Users/q/Documents/aspect/build_deal.ii.app" "-DTESTNAME=tests.quick_mpi" "-DERROR=\"Test quick_mpi failed\"" "-P" "/Users/q/Documents/aspect/tests/run_test.cmake")
set_tests_properties(quick_mpi PROPERTIES  TIMEOUT "600" WORKING_DIRECTORY "/Users/q/Documents/aspect/build_deal.ii.app/tests")
