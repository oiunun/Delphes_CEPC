# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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
CMAKE_COMMAND = /cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-centos7-gcc49/external/CMake/2.8.5/install/bin/cmake

# The command to remove a file.
RM = /cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-centos7-gcc49/external/CMake/2.8.5/install/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /workfs2/bes/lig/higgs/FSClasser

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workfs2/bes/lig/higgs/FSClasser/build

# Utility rule file for ExperimentalSubmit.

CMakeFiles/ExperimentalSubmit:
	/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-centos7-gcc49/external/CMake/2.8.5/install/bin/ctest -D ExperimentalSubmit

ExperimentalSubmit: CMakeFiles/ExperimentalSubmit
ExperimentalSubmit: CMakeFiles/ExperimentalSubmit.dir/build.make
.PHONY : ExperimentalSubmit

# Rule to build all files generated by this target.
CMakeFiles/ExperimentalSubmit.dir/build: ExperimentalSubmit
.PHONY : CMakeFiles/ExperimentalSubmit.dir/build

CMakeFiles/ExperimentalSubmit.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ExperimentalSubmit.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ExperimentalSubmit.dir/clean

CMakeFiles/ExperimentalSubmit.dir/depend:
	cd /workfs2/bes/lig/higgs/FSClasser/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workfs2/bes/lig/higgs/FSClasser /workfs2/bes/lig/higgs/FSClasser /workfs2/bes/lig/higgs/FSClasser/build /workfs2/bes/lig/higgs/FSClasser/build /workfs2/bes/lig/higgs/FSClasser/build/CMakeFiles/ExperimentalSubmit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ExperimentalSubmit.dir/depend

