# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.21.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.21.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/C-GDY-SAC"

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/opt/homebrew/Cellar/cmake/3.21.2/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/opt/homebrew/Cellar/cmake/3.21.2/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

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
# Target rules for targets named main.o

# Build rule for target.
main.o: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 main.o
.PHONY : main.o

# fast build rule for target.
main.o/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.o.dir/build.make CMakeFiles/main.o.dir/build
.PHONY : main.o/fast

#=============================================================================
# Target rules for targets named MyMaths

# Build rule for target.
MyMaths: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 MyMaths
.PHONY : MyMaths

# fast build rule for target.
MyMaths/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyMaths.dir/build.make CMakeFiles/MyMaths.dir/build
.PHONY : MyMaths/fast

#=============================================================================
# Target rules for targets named msiParser

# Build rule for target.
msiParser: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 msiParser
.PHONY : msiParser

# fast build rule for target.
msiParser/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/build
.PHONY : msiParser/fast

# target to build an object file
src/main.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.o.dir/build.make CMakeFiles/main.o.dir/src/main.o
.PHONY : src/main.o

# target to preprocess a source file
src/main.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.o.dir/build.make CMakeFiles/main.o.dir/src/main.i
.PHONY : src/main.i

# target to generate assembly for a file
src/main.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/main.o.dir/build.make CMakeFiles/main.o.dir/src/main.s
.PHONY : src/main.s

# target to build an object file
src/maths/MyMaths.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyMaths.dir/build.make CMakeFiles/MyMaths.dir/src/maths/MyMaths.o
.PHONY : src/maths/MyMaths.o

# target to preprocess a source file
src/maths/MyMaths.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyMaths.dir/build.make CMakeFiles/MyMaths.dir/src/maths/MyMaths.i
.PHONY : src/maths/MyMaths.i

# target to generate assembly for a file
src/maths/MyMaths.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyMaths.dir/build.make CMakeFiles/MyMaths.dir/src/maths/MyMaths.s
.PHONY : src/maths/MyMaths.s

# target to build an object file
src/msiParser/mod_msi.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/mod_msi.o
.PHONY : src/msiParser/mod_msi.o

# target to preprocess a source file
src/msiParser/mod_msi.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/mod_msi.i
.PHONY : src/msiParser/mod_msi.i

# target to generate assembly for a file
src/msiParser/mod_msi.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/mod_msi.s
.PHONY : src/msiParser/mod_msi.s

# target to build an object file
src/msiParser/parse_base.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_base.o
.PHONY : src/msiParser/parse_base.o

# target to preprocess a source file
src/msiParser/parse_base.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_base.i
.PHONY : src/msiParser/parse_base.i

# target to generate assembly for a file
src/msiParser/parse_base.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_base.s
.PHONY : src/msiParser/parse_base.s

# target to build an object file
src/msiParser/parse_mol.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_mol.o
.PHONY : src/msiParser/parse_mol.o

# target to preprocess a source file
src/msiParser/parse_mol.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_mol.i
.PHONY : src/msiParser/parse_mol.i

# target to generate assembly for a file
src/msiParser/parse_mol.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_mol.s
.PHONY : src/msiParser/parse_mol.s

# target to build an object file
src/msiParser/parse_msi.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_msi.o
.PHONY : src/msiParser/parse_msi.o

# target to preprocess a source file
src/msiParser/parse_msi.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_msi.i
.PHONY : src/msiParser/parse_msi.i

# target to generate assembly for a file
src/msiParser/parse_msi.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msiParser.dir/build.make CMakeFiles/msiParser.dir/src/msiParser/parse_msi.s
.PHONY : src/msiParser/parse_msi.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... MyMaths"
	@echo "... main.o"
	@echo "... msiParser"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/maths/MyMaths.o"
	@echo "... src/maths/MyMaths.i"
	@echo "... src/maths/MyMaths.s"
	@echo "... src/msiParser/mod_msi.o"
	@echo "... src/msiParser/mod_msi.i"
	@echo "... src/msiParser/mod_msi.s"
	@echo "... src/msiParser/parse_base.o"
	@echo "... src/msiParser/parse_base.i"
	@echo "... src/msiParser/parse_base.s"
	@echo "... src/msiParser/parse_mol.o"
	@echo "... src/msiParser/parse_mol.i"
	@echo "... src/msiParser/parse_mol.s"
	@echo "... src/msiParser/parse_msi.o"
	@echo "... src/msiParser/parse_msi.i"
	@echo "... src/msiParser/parse_msi.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

