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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.22.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.22.2/bin/cmake -E rm -f

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
	/opt/homebrew/Cellar/cmake/3.22.1/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/opt/homebrew/Cellar/cmake/3.22.2/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
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
# Target rules for targets named database

# Build rule for target.
database: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 database
.PHONY : database

# fast build rule for target.
database/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/build
.PHONY : database/fast

#=============================================================================
# Target rules for targets named msi_modeling

# Build rule for target.
msi_modeling: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 msi_modeling
.PHONY : msi_modeling

# fast build rule for target.
msi_modeling/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/build
.PHONY : msi_modeling/fast

#=============================================================================
# Target rules for targets named castepSeedGen

# Build rule for target.
castepSeedGen: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 castepSeedGen
.PHONY : castepSeedGen

# fast build rule for target.
castepSeedGen/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/build
.PHONY : castepSeedGen/fast

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

src/atom.o: src/atom.c.o
.PHONY : src/atom.o

# target to build an object file
src/atom.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/atom.c.o
.PHONY : src/atom.c.o

src/atom.i: src/atom.c.i
.PHONY : src/atom.i

# target to preprocess a source file
src/atom.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/atom.c.i
.PHONY : src/atom.c.i

src/atom.s: src/atom.c.s
.PHONY : src/atom.s

# target to generate assembly for a file
src/atom.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/atom.c.s
.PHONY : src/atom.c.s

src/cell.o: src/cell.c.o
.PHONY : src/cell.o

# target to build an object file
src/cell.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/cell.c.o
.PHONY : src/cell.c.o

src/cell.i: src/cell.c.i
.PHONY : src/cell.i

# target to preprocess a source file
src/cell.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/cell.c.i
.PHONY : src/cell.c.i

src/cell.s: src/cell.c.s
.PHONY : src/cell.s

# target to generate assembly for a file
src/cell.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/cell.c.s
.PHONY : src/cell.c.s

src/database/ads_database.o: src/database/ads_database.c.o
.PHONY : src/database/ads_database.o

# target to build an object file
src/database/ads_database.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/ads_database.c.o
.PHONY : src/database/ads_database.c.o

src/database/ads_database.i: src/database/ads_database.c.i
.PHONY : src/database/ads_database.i

# target to preprocess a source file
src/database/ads_database.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/ads_database.c.i
.PHONY : src/database/ads_database.c.i

src/database/ads_database.s: src/database/ads_database.c.s
.PHONY : src/database/ads_database.s

# target to generate assembly for a file
src/database/ads_database.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/ads_database.c.s
.PHONY : src/database/ads_database.c.s

src/database/database.o: src/database/database.c.o
.PHONY : src/database/database.o

# target to build an object file
src/database/database.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/database.c.o
.PHONY : src/database/database.c.o

src/database/database.i: src/database/database.c.i
.PHONY : src/database/database.i

# target to preprocess a source file
src/database/database.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/database.c.i
.PHONY : src/database/database.c.i

src/database/database.s: src/database/database.c.s
.PHONY : src/database/database.s

# target to generate assembly for a file
src/database/database.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/database.c.s
.PHONY : src/database/database.c.s

src/database/lattice_database.o: src/database/lattice_database.c.o
.PHONY : src/database/lattice_database.o

# target to build an object file
src/database/lattice_database.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/lattice_database.c.o
.PHONY : src/database/lattice_database.c.o

src/database/lattice_database.i: src/database/lattice_database.c.i
.PHONY : src/database/lattice_database.i

# target to preprocess a source file
src/database/lattice_database.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/lattice_database.c.i
.PHONY : src/database/lattice_database.c.i

src/database/lattice_database.s: src/database/lattice_database.c.s
.PHONY : src/database/lattice_database.s

# target to generate assembly for a file
src/database/lattice_database.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/database.dir/build.make CMakeFiles/database.dir/src/database/lattice_database.c.s
.PHONY : src/database/lattice_database.c.s

src/lattice.o: src/lattice.c.o
.PHONY : src/lattice.o

# target to build an object file
src/lattice.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/lattice.c.o
.PHONY : src/lattice.c.o

src/lattice.i: src/lattice.c.i
.PHONY : src/lattice.i

# target to preprocess a source file
src/lattice.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/lattice.c.i
.PHONY : src/lattice.c.i

src/lattice.s: src/lattice.c.s
.PHONY : src/lattice.s

# target to generate assembly for a file
src/lattice.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/lattice.c.s
.PHONY : src/lattice.c.s

src/main.o: src/main.c.o
.PHONY : src/main.o

# target to build an object file
src/main.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.c.o
.PHONY : src/main.c.o

src/main.i: src/main.c.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.c.i
.PHONY : src/main.c.i

src/main.s: src/main.c.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test.o.dir/build.make CMakeFiles/test.o.dir/src/main.c.s
.PHONY : src/main.c.s

src/misc.o: src/misc.c.o
.PHONY : src/misc.o

# target to build an object file
src/misc.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/misc.c.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/misc.c.o
.PHONY : src/misc.c.o

src/misc.i: src/misc.c.i
.PHONY : src/misc.i

# target to preprocess a source file
src/misc.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/misc.c.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/misc.c.i
.PHONY : src/misc.c.i

src/misc.s: src/misc.c.s
.PHONY : src/misc.s

# target to generate assembly for a file
src/misc.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/misc.c.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/castepSeedGen.dir/build.make CMakeFiles/castepSeedGen.dir/src/misc.c.s
.PHONY : src/misc.c.s

src/molecule.o: src/molecule.c.o
.PHONY : src/molecule.o

# target to build an object file
src/molecule.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/molecule.c.o
.PHONY : src/molecule.c.o

src/molecule.i: src/molecule.c.i
.PHONY : src/molecule.i

# target to preprocess a source file
src/molecule.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/molecule.c.i
.PHONY : src/molecule.c.i

src/molecule.s: src/molecule.c.s
.PHONY : src/molecule.s

# target to generate assembly for a file
src/molecule.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/molecule.c.s
.PHONY : src/molecule.c.s

src/my_maths.o: src/my_maths.c.o
.PHONY : src/my_maths.o

# target to build an object file
src/my_maths.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.c.o
.PHONY : src/my_maths.c.o

src/my_maths.i: src/my_maths.c.i
.PHONY : src/my_maths.i

# target to preprocess a source file
src/my_maths.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.c.i
.PHONY : src/my_maths.c.i

src/my_maths.s: src/my_maths.c.s
.PHONY : src/my_maths.s

# target to generate assembly for a file
src/my_maths.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/my_maths.dir/build.make CMakeFiles/my_maths.dir/src/my_maths.c.s
.PHONY : src/my_maths.c.s

src/parser.o: src/parser.c.o
.PHONY : src/parser.o

# target to build an object file
src/parser.c.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/parser.c.o
.PHONY : src/parser.c.o

src/parser.i: src/parser.c.i
.PHONY : src/parser.i

# target to preprocess a source file
src/parser.c.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/parser.c.i
.PHONY : src/parser.c.i

src/parser.s: src/parser.c.s
.PHONY : src/parser.s

# target to generate assembly for a file
src/parser.c.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/msi_modeling.dir/build.make CMakeFiles/msi_modeling.dir/src/parser.c.s
.PHONY : src/parser.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... castepSeedGen"
	@echo "... database"
	@echo "... msi_modeling"
	@echo "... my_maths"
	@echo "... test.o"
	@echo "... src/atom.o"
	@echo "... src/atom.i"
	@echo "... src/atom.s"
	@echo "... src/cell.o"
	@echo "... src/cell.i"
	@echo "... src/cell.s"
	@echo "... src/database/ads_database.o"
	@echo "... src/database/ads_database.i"
	@echo "... src/database/ads_database.s"
	@echo "... src/database/database.o"
	@echo "... src/database/database.i"
	@echo "... src/database/database.s"
	@echo "... src/database/lattice_database.o"
	@echo "... src/database/lattice_database.i"
	@echo "... src/database/lattice_database.s"
	@echo "... src/lattice.o"
	@echo "... src/lattice.i"
	@echo "... src/lattice.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/misc.o"
	@echo "... src/misc.i"
	@echo "... src/misc.s"
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

