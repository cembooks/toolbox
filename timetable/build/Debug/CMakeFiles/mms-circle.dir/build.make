# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/siarhei/Documents/fem/cem/dealii/static/mms

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug

# Include any dependencies generated for this target.
include CMakeFiles/mms-circle.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mms-circle.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mms-circle.dir/flags.make

CMakeFiles/mms-circle.dir/src/main.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/src/main.cpp.o: ../../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mms-circle.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/src/main.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/mms/src/main.cpp

CMakeFiles/mms-circle.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/mms/src/main.cpp > CMakeFiles/mms-circle.dir/src/main.cpp.i

CMakeFiles/mms-circle.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/mms/src/main.cpp -o CMakeFiles/mms-circle.dir/src/main.cpp.s

CMakeFiles/mms-circle.dir/src/solver.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/src/solver.cpp.o: ../../src/solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mms-circle.dir/src/solver.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/src/solver.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/mms/src/solver.cpp

CMakeFiles/mms-circle.dir/src/solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/src/solver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/mms/src/solver.cpp > CMakeFiles/mms-circle.dir/src/solver.cpp.i

CMakeFiles/mms-circle.dir/src/solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/src/solver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/mms/src/solver.cpp -o CMakeFiles/mms-circle.dir/src/solver.cpp.s

CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o: ../../src/exact_solution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/mms/src/exact_solution.cpp

CMakeFiles/mms-circle.dir/src/exact_solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/src/exact_solution.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/mms/src/exact_solution.cpp > CMakeFiles/mms-circle.dir/src/exact_solution.cpp.i

CMakeFiles/mms-circle.dir/src/exact_solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/src/exact_solution.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/mms/src/exact_solution.cpp -o CMakeFiles/mms-circle.dir/src/exact_solution.cpp.s

CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o: ../../src/static_scalar_input.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/mms/src/static_scalar_input.cpp

CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/mms/src/static_scalar_input.cpp > CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.i

CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/mms/src/static_scalar_input.cpp -o CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.s

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o: /home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp > CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.i

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp -o CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.s

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o: CMakeFiles/mms-circle.dir/flags.make
CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o: /home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o -c /home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp > CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.i

CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp -o CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.s

# Object files for target mms-circle
mms__circle_OBJECTS = \
"CMakeFiles/mms-circle.dir/src/main.cpp.o" \
"CMakeFiles/mms-circle.dir/src/solver.cpp.o" \
"CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o" \
"CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o" \
"CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o" \
"CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o"

# External object files for target mms-circle
mms__circle_EXTERNAL_OBJECTS =

../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/src/main.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/src/solver.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/src/exact_solution.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/src/static_scalar_input.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/misc.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/home/siarhei/Documents/fem/cem/dealii/static/shared/src/constants.cpp.o
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/build.make
../../bin/Debug/mms-circle: /home/siarhei/dealii-9.6.0/lib/libdeal_II.g.so.9.7.0-pre
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libtbb.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libz.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_system.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_thread.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_regex.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libumfpack.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libcholmod.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libccolamd.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libcolamd.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libcamd.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libamd.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libmetis.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libarpack.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/liblapack.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libblas.so
../../bin/Debug/mms-circle: /usr/lib/x86_64-linux-gnu/libgmsh.so
../../bin/Debug/mms-circle: CMakeFiles/mms-circle.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable ../../bin/Debug/mms-circle"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mms-circle.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mms-circle.dir/build: ../../bin/Debug/mms-circle

.PHONY : CMakeFiles/mms-circle.dir/build

CMakeFiles/mms-circle.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mms-circle.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mms-circle.dir/clean

CMakeFiles/mms-circle.dir/depend:
	cd /home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/siarhei/Documents/fem/cem/dealii/static/mms /home/siarhei/Documents/fem/cem/dealii/static/mms /home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug /home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug /home/siarhei/Documents/fem/cem/dealii/static/mms/build/Debug/CMakeFiles/mms-circle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mms-circle.dir/depend

