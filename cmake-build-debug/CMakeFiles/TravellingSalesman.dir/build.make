# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.26

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2023.2.2\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2023.2.2\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\teres\Downloads\TravellingSalesman

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/TravellingSalesman.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/TravellingSalesman.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/TravellingSalesman.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TravellingSalesman.dir/flags.make

CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj: CMakeFiles/TravellingSalesman.dir/flags.make
CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj: C:/Users/teres/Downloads/TravellingSalesman/src/main.cpp
CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj: CMakeFiles/TravellingSalesman.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj -MF CMakeFiles\TravellingSalesman.dir\src\main.cpp.obj.d -o CMakeFiles\TravellingSalesman.dir\src\main.cpp.obj -c C:\Users\teres\Downloads\TravellingSalesman\src\main.cpp

CMakeFiles/TravellingSalesman.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TravellingSalesman.dir/src/main.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\teres\Downloads\TravellingSalesman\src\main.cpp > CMakeFiles\TravellingSalesman.dir\src\main.cpp.i

CMakeFiles/TravellingSalesman.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TravellingSalesman.dir/src/main.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\teres\Downloads\TravellingSalesman\src\main.cpp -o CMakeFiles\TravellingSalesman.dir\src\main.cpp.s

CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj: CMakeFiles/TravellingSalesman.dir/flags.make
CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj: C:/Users/teres/Downloads/TravellingSalesman/src/classes/Script.cpp
CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj: CMakeFiles/TravellingSalesman.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj -MF CMakeFiles\TravellingSalesman.dir\src\classes\Script.cpp.obj.d -o CMakeFiles\TravellingSalesman.dir\src\classes\Script.cpp.obj -c C:\Users\teres\Downloads\TravellingSalesman\src\classes\Script.cpp

CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\teres\Downloads\TravellingSalesman\src\classes\Script.cpp > CMakeFiles\TravellingSalesman.dir\src\classes\Script.cpp.i

CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\teres\Downloads\TravellingSalesman\src\classes\Script.cpp -o CMakeFiles\TravellingSalesman.dir\src\classes\Script.cpp.s

CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj: CMakeFiles/TravellingSalesman.dir/flags.make
CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj: C:/Users/teres/Downloads/TravellingSalesman/src/classes/TSPSolver.cpp
CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj: CMakeFiles/TravellingSalesman.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj -MF CMakeFiles\TravellingSalesman.dir\src\classes\TSPSolver.cpp.obj.d -o CMakeFiles\TravellingSalesman.dir\src\classes\TSPSolver.cpp.obj -c C:\Users\teres\Downloads\TravellingSalesman\src\classes\TSPSolver.cpp

CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\teres\Downloads\TravellingSalesman\src\classes\TSPSolver.cpp > CMakeFiles\TravellingSalesman.dir\src\classes\TSPSolver.cpp.i

CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\teres\Downloads\TravellingSalesman\src\classes\TSPSolver.cpp -o CMakeFiles\TravellingSalesman.dir\src\classes\TSPSolver.cpp.s

CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj: CMakeFiles/TravellingSalesman.dir/flags.make
CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj: C:/Users/teres/Downloads/TravellingSalesman/src/classes/Menu.cpp
CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj: CMakeFiles/TravellingSalesman.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj -MF CMakeFiles\TravellingSalesman.dir\src\classes\Menu.cpp.obj.d -o CMakeFiles\TravellingSalesman.dir\src\classes\Menu.cpp.obj -c C:\Users\teres\Downloads\TravellingSalesman\src\classes\Menu.cpp

CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\teres\Downloads\TravellingSalesman\src\classes\Menu.cpp > CMakeFiles\TravellingSalesman.dir\src\classes\Menu.cpp.i

CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\teres\Downloads\TravellingSalesman\src\classes\Menu.cpp -o CMakeFiles\TravellingSalesman.dir\src\classes\Menu.cpp.s

# Object files for target TravellingSalesman
TravellingSalesman_OBJECTS = \
"CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj" \
"CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj" \
"CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj" \
"CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj"

# External object files for target TravellingSalesman
TravellingSalesman_EXTERNAL_OBJECTS =

TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/src/main.cpp.obj
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/src/classes/Script.cpp.obj
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/src/classes/TSPSolver.cpp.obj
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/src/classes/Menu.cpp.obj
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/build.make
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/linkLibs.rsp
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/objects1.rsp
TravellingSalesman.exe: CMakeFiles/TravellingSalesman.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable TravellingSalesman.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\TravellingSalesman.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TravellingSalesman.dir/build: TravellingSalesman.exe
.PHONY : CMakeFiles/TravellingSalesman.dir/build

CMakeFiles/TravellingSalesman.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\TravellingSalesman.dir\cmake_clean.cmake
.PHONY : CMakeFiles/TravellingSalesman.dir/clean

CMakeFiles/TravellingSalesman.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\teres\Downloads\TravellingSalesman C:\Users\teres\Downloads\TravellingSalesman C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug C:\Users\teres\Downloads\TravellingSalesman\cmake-build-debug\CMakeFiles\TravellingSalesman.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TravellingSalesman.dir/depend
