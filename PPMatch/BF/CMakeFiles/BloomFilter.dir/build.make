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
CMAKE_SOURCE_DIR = /home/lyu/Desktop/graph-homomorphism

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lyu/Desktop/graph-homomorphism

# Include any dependencies generated for this target.
include BF/CMakeFiles/BloomFilter.dir/depend.make

# Include the progress variables for this target.
include BF/CMakeFiles/BloomFilter.dir/progress.make

# Include the compile flags for this target's objects.
include BF/CMakeFiles/BloomFilter.dir/flags.make

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o: BF/CMakeFiles/BloomFilter.dir/flags.make
BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o: BF/BloomFilter.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lyu/Desktop/graph-homomorphism/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o -c /home/lyu/Desktop/graph-homomorphism/BF/BloomFilter.cpp

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BloomFilter.dir/BloomFilter.cpp.i"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lyu/Desktop/graph-homomorphism/BF/BloomFilter.cpp > CMakeFiles/BloomFilter.dir/BloomFilter.cpp.i

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BloomFilter.dir/BloomFilter.cpp.s"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lyu/Desktop/graph-homomorphism/BF/BloomFilter.cpp -o CMakeFiles/BloomFilter.dir/BloomFilter.cpp.s

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.requires:

.PHONY : BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.requires

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.provides: BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.requires
	$(MAKE) -f BF/CMakeFiles/BloomFilter.dir/build.make BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.provides.build
.PHONY : BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.provides

BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.provides.build: BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o


BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o: BF/CMakeFiles/BloomFilter.dir/flags.make
BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o: BF/MurmurHash3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lyu/Desktop/graph-homomorphism/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o -c /home/lyu/Desktop/graph-homomorphism/BF/MurmurHash3.cpp

BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.i"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lyu/Desktop/graph-homomorphism/BF/MurmurHash3.cpp > CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.i

BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.s"
	cd /home/lyu/Desktop/graph-homomorphism/BF && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lyu/Desktop/graph-homomorphism/BF/MurmurHash3.cpp -o CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.s

BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.requires:

.PHONY : BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.requires

BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.provides: BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.requires
	$(MAKE) -f BF/CMakeFiles/BloomFilter.dir/build.make BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.provides.build
.PHONY : BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.provides

BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.provides.build: BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o


# Object files for target BloomFilter
BloomFilter_OBJECTS = \
"CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o" \
"CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o"

# External object files for target BloomFilter
BloomFilter_EXTERNAL_OBJECTS =

BF/libBloomFilter.a: BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o
BF/libBloomFilter.a: BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o
BF/libBloomFilter.a: BF/CMakeFiles/BloomFilter.dir/build.make
BF/libBloomFilter.a: BF/CMakeFiles/BloomFilter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lyu/Desktop/graph-homomorphism/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libBloomFilter.a"
	cd /home/lyu/Desktop/graph-homomorphism/BF && $(CMAKE_COMMAND) -P CMakeFiles/BloomFilter.dir/cmake_clean_target.cmake
	cd /home/lyu/Desktop/graph-homomorphism/BF && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BloomFilter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
BF/CMakeFiles/BloomFilter.dir/build: BF/libBloomFilter.a

.PHONY : BF/CMakeFiles/BloomFilter.dir/build

BF/CMakeFiles/BloomFilter.dir/requires: BF/CMakeFiles/BloomFilter.dir/BloomFilter.cpp.o.requires
BF/CMakeFiles/BloomFilter.dir/requires: BF/CMakeFiles/BloomFilter.dir/MurmurHash3.cpp.o.requires

.PHONY : BF/CMakeFiles/BloomFilter.dir/requires

BF/CMakeFiles/BloomFilter.dir/clean:
	cd /home/lyu/Desktop/graph-homomorphism/BF && $(CMAKE_COMMAND) -P CMakeFiles/BloomFilter.dir/cmake_clean.cmake
.PHONY : BF/CMakeFiles/BloomFilter.dir/clean

BF/CMakeFiles/BloomFilter.dir/depend:
	cd /home/lyu/Desktop/graph-homomorphism && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lyu/Desktop/graph-homomorphism /home/lyu/Desktop/graph-homomorphism/BF /home/lyu/Desktop/graph-homomorphism /home/lyu/Desktop/graph-homomorphism/BF /home/lyu/Desktop/graph-homomorphism/BF/CMakeFiles/BloomFilter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : BF/CMakeFiles/BloomFilter.dir/depend
