# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /nas2/ywu/anaconda2/bin/cmake

# The command to remove a file.
RM = /nas2/ywu/anaconda2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /nas2/ywu/Translocator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /nas2/ywu/Translocator/cmake-build-debug

# Include any dependencies generated for this target.
include lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/depend.make

# Include the progress variables for this target.
include lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/progress.make

# Include the compile flags for this target's objects.
include lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.o: ../lib/zlib-1.2.7/adler32.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/adler32.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/adler32.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/adler32.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/adler32.c > CMakeFiles/zlibstatic.dir/adler32.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/adler32.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/adler32.c -o CMakeFiles/zlibstatic.dir/adler32.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.o: ../lib/zlib-1.2.7/compress.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/compress.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/compress.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/compress.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/compress.c > CMakeFiles/zlibstatic.dir/compress.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/compress.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/compress.c -o CMakeFiles/zlibstatic.dir/compress.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.o: ../lib/zlib-1.2.7/crc32.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/crc32.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/crc32.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/crc32.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/crc32.c > CMakeFiles/zlibstatic.dir/crc32.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/crc32.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/crc32.c -o CMakeFiles/zlibstatic.dir/crc32.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.o: ../lib/zlib-1.2.7/deflate.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/deflate.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/deflate.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/deflate.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/deflate.c > CMakeFiles/zlibstatic.dir/deflate.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/deflate.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/deflate.c -o CMakeFiles/zlibstatic.dir/deflate.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.o: ../lib/zlib-1.2.7/gzclose.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/gzclose.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/gzclose.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/gzclose.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/gzclose.c > CMakeFiles/zlibstatic.dir/gzclose.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/gzclose.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/gzclose.c -o CMakeFiles/zlibstatic.dir/gzclose.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.o: ../lib/zlib-1.2.7/gzlib.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/gzlib.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/gzlib.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/gzlib.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/gzlib.c > CMakeFiles/zlibstatic.dir/gzlib.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/gzlib.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/gzlib.c -o CMakeFiles/zlibstatic.dir/gzlib.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.o: ../lib/zlib-1.2.7/gzread.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/gzread.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/gzread.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/gzread.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/gzread.c > CMakeFiles/zlibstatic.dir/gzread.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/gzread.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/gzread.c -o CMakeFiles/zlibstatic.dir/gzread.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.o: ../lib/zlib-1.2.7/gzwrite.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/gzwrite.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/gzwrite.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/gzwrite.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/gzwrite.c > CMakeFiles/zlibstatic.dir/gzwrite.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/gzwrite.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/gzwrite.c -o CMakeFiles/zlibstatic.dir/gzwrite.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.o: ../lib/zlib-1.2.7/inflate.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/inflate.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/inflate.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/inflate.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/inflate.c > CMakeFiles/zlibstatic.dir/inflate.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/inflate.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/inflate.c -o CMakeFiles/zlibstatic.dir/inflate.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.o: ../lib/zlib-1.2.7/infback.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/infback.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/infback.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/infback.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/infback.c > CMakeFiles/zlibstatic.dir/infback.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/infback.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/infback.c -o CMakeFiles/zlibstatic.dir/infback.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.o: ../lib/zlib-1.2.7/inftrees.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/inftrees.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/inftrees.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/inftrees.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/inftrees.c > CMakeFiles/zlibstatic.dir/inftrees.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/inftrees.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/inftrees.c -o CMakeFiles/zlibstatic.dir/inftrees.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.o: ../lib/zlib-1.2.7/inffast.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/inffast.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/inffast.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/inffast.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/inffast.c > CMakeFiles/zlibstatic.dir/inffast.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/inffast.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/inffast.c -o CMakeFiles/zlibstatic.dir/inffast.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.o: ../lib/zlib-1.2.7/trees.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/trees.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/trees.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/trees.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/trees.c > CMakeFiles/zlibstatic.dir/trees.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/trees.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/trees.c -o CMakeFiles/zlibstatic.dir/trees.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.o: ../lib/zlib-1.2.7/uncompr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/uncompr.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/uncompr.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/uncompr.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/uncompr.c > CMakeFiles/zlibstatic.dir/uncompr.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/uncompr.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/uncompr.c -o CMakeFiles/zlibstatic.dir/uncompr.s

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.o: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/flags.make
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.o: ../lib/zlib-1.2.7/zutil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.o"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/zlibstatic.dir/zutil.o   -c /nas2/ywu/Translocator/lib/zlib-1.2.7/zutil.c

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlibstatic.dir/zutil.i"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /nas2/ywu/Translocator/lib/zlib-1.2.7/zutil.c > CMakeFiles/zlibstatic.dir/zutil.i

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlibstatic.dir/zutil.s"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /nas2/ywu/Translocator/lib/zlib-1.2.7/zutil.c -o CMakeFiles/zlibstatic.dir/zutil.s

# Object files for target zlibstatic
zlibstatic_OBJECTS = \
"CMakeFiles/zlibstatic.dir/adler32.o" \
"CMakeFiles/zlibstatic.dir/compress.o" \
"CMakeFiles/zlibstatic.dir/crc32.o" \
"CMakeFiles/zlibstatic.dir/deflate.o" \
"CMakeFiles/zlibstatic.dir/gzclose.o" \
"CMakeFiles/zlibstatic.dir/gzlib.o" \
"CMakeFiles/zlibstatic.dir/gzread.o" \
"CMakeFiles/zlibstatic.dir/gzwrite.o" \
"CMakeFiles/zlibstatic.dir/inflate.o" \
"CMakeFiles/zlibstatic.dir/infback.o" \
"CMakeFiles/zlibstatic.dir/inftrees.o" \
"CMakeFiles/zlibstatic.dir/inffast.o" \
"CMakeFiles/zlibstatic.dir/trees.o" \
"CMakeFiles/zlibstatic.dir/uncompr.o" \
"CMakeFiles/zlibstatic.dir/zutil.o"

# External object files for target zlibstatic
zlibstatic_EXTERNAL_OBJECTS =

lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/adler32.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/compress.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/crc32.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/deflate.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzclose.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzlib.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzread.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/gzwrite.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inflate.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/infback.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inftrees.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/inffast.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/trees.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/uncompr.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/zutil.o
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/build.make
lib/zlib-1.2.7/libz.a: lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/nas2/ywu/Translocator/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking C static library libz.a"
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && $(CMAKE_COMMAND) -P CMakeFiles/zlibstatic.dir/cmake_clean_target.cmake
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zlibstatic.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/build: lib/zlib-1.2.7/libz.a

.PHONY : lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/build

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/clean:
	cd /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 && $(CMAKE_COMMAND) -P CMakeFiles/zlibstatic.dir/cmake_clean.cmake
.PHONY : lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/clean

lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/depend:
	cd /nas2/ywu/Translocator/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /nas2/ywu/Translocator /nas2/ywu/Translocator/lib/zlib-1.2.7 /nas2/ywu/Translocator/cmake-build-debug /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7 /nas2/ywu/Translocator/cmake-build-debug/lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/zlib-1.2.7/CMakeFiles/zlibstatic.dir/depend

