# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build

# Include any dependencies generated for this target.
include CMakeFiles/mm_simd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mm_simd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mm_simd.dir/flags.make

CMakeFiles/mm_simd.dir/common/helper.cpp.o: CMakeFiles/mm_simd.dir/flags.make
CMakeFiles/mm_simd.dir/common/helper.cpp.o: ../common/helper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mm_simd.dir/common/helper.cpp.o"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mm_simd.dir/common/helper.cpp.o -c /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/helper.cpp

CMakeFiles/mm_simd.dir/common/helper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mm_simd.dir/common/helper.cpp.i"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/helper.cpp > CMakeFiles/mm_simd.dir/common/helper.cpp.i

CMakeFiles/mm_simd.dir/common/helper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mm_simd.dir/common/helper.cpp.s"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/helper.cpp -o CMakeFiles/mm_simd.dir/common/helper.cpp.s

CMakeFiles/mm_simd.dir/common/helper.cpp.o.requires:

.PHONY : CMakeFiles/mm_simd.dir/common/helper.cpp.o.requires

CMakeFiles/mm_simd.dir/common/helper.cpp.o.provides: CMakeFiles/mm_simd.dir/common/helper.cpp.o.requires
	$(MAKE) -f CMakeFiles/mm_simd.dir/build.make CMakeFiles/mm_simd.dir/common/helper.cpp.o.provides.build
.PHONY : CMakeFiles/mm_simd.dir/common/helper.cpp.o.provides

CMakeFiles/mm_simd.dir/common/helper.cpp.o.provides.build: CMakeFiles/mm_simd.dir/common/helper.cpp.o


CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o: CMakeFiles/mm_simd.dir/flags.make
CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o: ../common/compute_kernels_default.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o -c /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/compute_kernels_default.cpp

CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.i"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/compute_kernels_default.cpp > CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.i

CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.s"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/common/compute_kernels_default.cpp -o CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.s

CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.requires:

.PHONY : CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.requires

CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.provides: CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.requires
	$(MAKE) -f CMakeFiles/mm_simd.dir/build.make CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.provides.build
.PHONY : CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.provides

CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.provides.build: CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o


CMakeFiles/mm_simd.dir/src/main.cpp.o: CMakeFiles/mm_simd.dir/flags.make
CMakeFiles/mm_simd.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mm_simd.dir/src/main.cpp.o"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mm_simd.dir/src/main.cpp.o -c /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/src/main.cpp

CMakeFiles/mm_simd.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mm_simd.dir/src/main.cpp.i"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/src/main.cpp > CMakeFiles/mm_simd.dir/src/main.cpp.i

CMakeFiles/mm_simd.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mm_simd.dir/src/main.cpp.s"
	/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel/19.0.5-gcc-uglchea/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/src/main.cpp -o CMakeFiles/mm_simd.dir/src/main.cpp.s

CMakeFiles/mm_simd.dir/src/main.cpp.o.requires:

.PHONY : CMakeFiles/mm_simd.dir/src/main.cpp.o.requires

CMakeFiles/mm_simd.dir/src/main.cpp.o.provides: CMakeFiles/mm_simd.dir/src/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/mm_simd.dir/build.make CMakeFiles/mm_simd.dir/src/main.cpp.o.provides.build
.PHONY : CMakeFiles/mm_simd.dir/src/main.cpp.o.provides

CMakeFiles/mm_simd.dir/src/main.cpp.o.provides.build: CMakeFiles/mm_simd.dir/src/main.cpp.o


# Object files for target mm_simd
mm_simd_OBJECTS = \
"CMakeFiles/mm_simd.dir/common/helper.cpp.o" \
"CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o" \
"CMakeFiles/mm_simd.dir/src/main.cpp.o"

# External object files for target mm_simd
mm_simd_EXTERNAL_OBJECTS =

mm_simd: CMakeFiles/mm_simd.dir/common/helper.cpp.o
mm_simd: CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o
mm_simd: CMakeFiles/mm_simd.dir/src/main.cpp.o
mm_simd: CMakeFiles/mm_simd.dir/build.make
mm_simd: /dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/qt/5.14.2-gcc-34o5ww4/lib/libQt5Widgets.so.5.14.2
mm_simd: simd/libconv_simd.a
mm_simd: /dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/qt/5.14.2-gcc-34o5ww4/lib/libQt5Gui.so.5.14.2
mm_simd: /dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/qt/5.14.2-gcc-34o5ww4/lib/libQt5Core.so.5.14.2
mm_simd: CMakeFiles/mm_simd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable mm_simd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mm_simd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mm_simd.dir/build: mm_simd

.PHONY : CMakeFiles/mm_simd.dir/build

CMakeFiles/mm_simd.dir/requires: CMakeFiles/mm_simd.dir/common/helper.cpp.o.requires
CMakeFiles/mm_simd.dir/requires: CMakeFiles/mm_simd.dir/common/compute_kernels_default.cpp.o.requires
CMakeFiles/mm_simd.dir/requires: CMakeFiles/mm_simd.dir/src/main.cpp.o.requires

.PHONY : CMakeFiles/mm_simd.dir/requires

CMakeFiles/mm_simd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mm_simd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mm_simd.dir/clean

CMakeFiles/mm_simd.dir/depend:
	cd /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build /dss/dsshome1/lxc04/di29wab/hpclab/A1/src/2dconv/build/CMakeFiles/mm_simd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mm_simd.dir/depend

