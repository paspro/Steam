# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_SOURCE_DIR = /home/panos/Development/Steam/SteamTesting

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/panos/Development/Steam/build/SteamTesting

# Include any dependencies generated for this target.
include CMakeFiles/steam-testing.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/steam-testing.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/steam-testing.dir/flags.make

CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o: CMakeFiles/steam-testing.dir/flags.make
CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o: /home/panos/Development/Steam/SteamTesting/src/SteamTablesTesting.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/panos/Development/Steam/build/SteamTesting/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/panos/Development/Steam/SteamTesting/src/SteamTablesTesting.f -o CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o

CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.i"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/panos/Development/Steam/SteamTesting/src/SteamTablesTesting.f > CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.i

CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.s"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/panos/Development/Steam/SteamTesting/src/SteamTablesTesting.f -o CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.s

CMakeFiles/steam-testing.dir/src/SteamTesting.f.o: CMakeFiles/steam-testing.dir/flags.make
CMakeFiles/steam-testing.dir/src/SteamTesting.f.o: /home/panos/Development/Steam/SteamTesting/src/SteamTesting.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/panos/Development/Steam/build/SteamTesting/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/steam-testing.dir/src/SteamTesting.f.o"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/panos/Development/Steam/SteamTesting/src/SteamTesting.f -o CMakeFiles/steam-testing.dir/src/SteamTesting.f.o

CMakeFiles/steam-testing.dir/src/SteamTesting.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/steam-testing.dir/src/SteamTesting.f.i"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/panos/Development/Steam/SteamTesting/src/SteamTesting.f > CMakeFiles/steam-testing.dir/src/SteamTesting.f.i

CMakeFiles/steam-testing.dir/src/SteamTesting.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/steam-testing.dir/src/SteamTesting.f.s"
	/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/panos/Development/Steam/SteamTesting/src/SteamTesting.f -o CMakeFiles/steam-testing.dir/src/SteamTesting.f.s

# Object files for target steam-testing
steam__testing_OBJECTS = \
"CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o" \
"CMakeFiles/steam-testing.dir/src/SteamTesting.f.o"

# External object files for target steam-testing
steam__testing_EXTERNAL_OBJECTS =

steam-testing: CMakeFiles/steam-testing.dir/src/SteamTablesTesting.f.o
steam-testing: CMakeFiles/steam-testing.dir/src/SteamTesting.f.o
steam-testing: CMakeFiles/steam-testing.dir/build.make
steam-testing: ../SteamPolynomials/libSteamPolynomials.so
steam-testing: CMakeFiles/steam-testing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/panos/Development/Steam/build/SteamTesting/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran executable steam-testing"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/steam-testing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/steam-testing.dir/build: steam-testing

.PHONY : CMakeFiles/steam-testing.dir/build

CMakeFiles/steam-testing.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/steam-testing.dir/cmake_clean.cmake
.PHONY : CMakeFiles/steam-testing.dir/clean

CMakeFiles/steam-testing.dir/depend:
	cd /home/panos/Development/Steam/build/SteamTesting && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/panos/Development/Steam/SteamTesting /home/panos/Development/Steam/SteamTesting /home/panos/Development/Steam/build/SteamTesting /home/panos/Development/Steam/build/SteamTesting /home/panos/Development/Steam/build/SteamTesting/CMakeFiles/steam-testing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/steam-testing.dir/depend

