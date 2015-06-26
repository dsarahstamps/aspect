# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/q/Documents/aspect

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/q/Documents/aspect/build_deal.ii.app

# Utility rule file for tests.quick_mpi.

# Include the progress variables for this target.
include tests/CMakeFiles/tests.quick_mpi.dir/progress.make

tests/CMakeFiles/tests.quick_mpi: tests/output-quick_mpi/screen-output

tests/output-quick_mpi/screen-output: aspect
tests/output-quick_mpi/screen-output: tests/quick_mpi.x.prm
tests/output-quick_mpi/screen-output: tests/libquick_mpi.so
tests/output-quick_mpi/screen-output: ../tests/quick_mpi.sh
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/q/Documents/aspect/build_deal.ii.app/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating output-quick_mpi/screen-output"
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && for i in /Users/q/Documents/aspect/build_deal.ii.app/tests/output-quick_mpi/* ; do if echo $i | grep -q -v .cmp.notime ; then rm -f $i ; fi ; done
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && mpirun -np 2 /Users/q/Documents/aspect/build_deal.ii.app/tests/../aspect /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm > /Users/q/Documents/aspect/build_deal.ii.app/tests/output-quick_mpi/screen-output.tmp
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && ../../tests/quick_mpi.sh screen-output </Users/q/Documents/aspect/build_deal.ii.app/tests/output-quick_mpi/screen-output.tmp >/Users/q/Documents/aspect/build_deal.ii.app/tests/output-quick_mpi/screen-output
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && /usr/bin/perl -pi /Users/q/Documents/aspect/tests/normalize.pl /Users/q/Documents/aspect/build_deal.ii.app/tests/output-quick_mpi/*

tests/quick_mpi.x.prm: ../tests/quick_mpi.prm
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/q/Documents/aspect/build_deal.ii.app/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating quick_mpi.x.prm"
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && cp /Users/q/Documents/aspect/tests/quick_mpi.prm /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && echo '' >> /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && echo 'set Output directory = output-quick_mpi' >> /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && echo ''set Additional shared libraries = ./libquick_mpi.so'' >> /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && /usr/bin/perl -pi -e 's!\@SOURCE_DIR\@!/Users/q/Documents/aspect/tests!g ' /Users/q/Documents/aspect/build_deal.ii.app/tests/quick_mpi.x.prm

tests.quick_mpi: tests/CMakeFiles/tests.quick_mpi
tests.quick_mpi: tests/output-quick_mpi/screen-output
tests.quick_mpi: tests/quick_mpi.x.prm
tests.quick_mpi: tests/CMakeFiles/tests.quick_mpi.dir/build.make
.PHONY : tests.quick_mpi

# Rule to build all files generated by this target.
tests/CMakeFiles/tests.quick_mpi.dir/build: tests.quick_mpi
.PHONY : tests/CMakeFiles/tests.quick_mpi.dir/build

tests/CMakeFiles/tests.quick_mpi.dir/clean:
	cd /Users/q/Documents/aspect/build_deal.ii.app/tests && $(CMAKE_COMMAND) -P CMakeFiles/tests.quick_mpi.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/tests.quick_mpi.dir/clean

tests/CMakeFiles/tests.quick_mpi.dir/depend:
	cd /Users/q/Documents/aspect/build_deal.ii.app && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/q/Documents/aspect /Users/q/Documents/aspect/tests /Users/q/Documents/aspect/build_deal.ii.app /Users/q/Documents/aspect/build_deal.ii.app/tests /Users/q/Documents/aspect/build_deal.ii.app/tests/CMakeFiles/tests.quick_mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/tests.quick_mpi.dir/depend

