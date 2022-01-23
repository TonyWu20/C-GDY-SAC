# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.22.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.22.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC"

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/opt/homebrew/Cellar/cmake/3.21.3_1/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/opt/homebrew/Cellar/cmake/3.22.1/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC/CMakeFiles" "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC//CMakeFiles/progress.marks"
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC/CMakeFiles" 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named my_maths

# Build rule for target.
my_maths: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 my_maths
.PHONY : my_maths

# fast build rule for target.
my_maths/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/build
.PHONY : my_maths/fast

#=============================================================================
# Target rules for targets named datastruct

# Build rule for target.
datastruct: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 datastruct
.PHONY : datastruct

# fast build rule for target.
datastruct/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/build
.PHONY : datastruct/fast

#=============================================================================
# Target rules for targets named parser

# Build rule for target.
parser: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 parser
.PHONY : parser

# fast build rule for target.
parser/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/parser.dir/build.make CMakeFiles/parser.dir/build
.PHONY : parser/fast

#=============================================================================
# Target rules for targets named assemble

# Build rule for target.
assemble: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 assemble
.PHONY : assemble

# fast build rule for target.
assemble/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/assemble.dir/build.make CMakeFiles/assemble.dir/build
.PHONY : assemble/fast

#=============================================================================
# Target rules for targets named test.o

# Build rule for target.
test.o: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 test.o
.PHONY : test.o

# fast build rule for target.
test.o/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/build
.PHONY : test.o/fast

# target to build an object file
src/assemble.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/assemble.dir/build.make CMakeFiles/assemble.dir/src/assemble.o
.PHONY : src/assemble.o

# target to preprocess a source file
src/assemble.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/assemble.dir/build.make CMakeFiles/assemble.dir/src/assemble.i
.PHONY : src/assemble.i

# target to generate assembly for a file
src/assemble.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/assemble.dir/build.make CMakeFiles/assemble.dir/src/assemble.s
.PHONY : src/assemble.s

# target to build an object file
src/atom.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/atom.o
.PHONY : src/atom.o

# target to preprocess a source file
src/atom.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/atom.i
.PHONY : src/atom.i

# target to generate assembly for a file
src/atom.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/atom.s
.PHONY : src/atom.s

# target to build an object file
src/lattice.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/lattice.o
.PHONY : src/lattice.o

# target to preprocess a source file
src/lattice.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/lattice.i
.PHONY : src/lattice.i

# target to generate assembly for a file
src/lattice.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/lattice.s
.PHONY : src/lattice.s

# target to build an object file
src/main.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.o
.PHONY : src/main.o

# target to preprocess a source file
src/main.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.i
.PHONY : src/main.i

# target to generate assembly for a file
src/main.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.s
.PHONY : src/main.s

# target to build an object file
src/molecule.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/molecule.o
.PHONY : src/molecule.o

# target to preprocess a source file
src/molecule.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/molecule.i
.PHONY : src/molecule.i

# target to generate assembly for a file
src/molecule.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/datastruct.dir/build.make CMakeFiles/datastruct.dir/src/molecule.s
.PHONY : src/molecule.s

# target to build an object file
src/my_maths.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.o
.PHONY : src/my_maths.o

# target to preprocess a source file
src/my_maths.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.i
.PHONY : src/my_maths.i

# target to generate assembly for a file
src/my_maths.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.s
.PHONY : src/my_maths.s

# target to build an object file
src/parser.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/parser.dir/build.make CMakeFiles/parser.dir/src/parser.o
.PHONY : src/parser.o

# target to preprocess a source file
src/parser.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/parser.dir/build.make CMakeFiles/parser.dir/src/parser.i
.PHONY : src/parser.i

# target to generate assembly for a file
src/parser.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/parser.dir/build.make CMakeFiles/parser.dir/src/parser.s
.PHONY : src/parser.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... assemble"
	@echo "... datastruct"
	@echo "... my_maths"
	@echo "... parser"
	@echo "... test.o"
	@echo "... src/assemble.o"
	@echo "... src/assemble.i"
	@echo "... src/assemble.s"
	@echo "... src/atom.o"
	@echo "... src/atom.i"
	@echo "... src/atom.s"
	@echo "... src/lattice.o"
	@echo "... src/lattice.i"
	@echo "... src/lattice.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/molecule.o"
	@echo "... src/molecule.i"
	@echo "... src/molecule.s"
	@echo "... src/my_maths.o"
	@echo "... src/my_maths.i"
	@echo "... src/my_maths.s"
	@echo "... src/parser.o"
	@echo "... src/parser.i"
	@echo "... src/parser.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

