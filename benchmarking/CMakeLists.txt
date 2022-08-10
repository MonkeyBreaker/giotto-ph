cmake_minimum_required(VERSION 3.1.0)
project(ripser)

set(TURF_PREFER_CPP11 TRUE CACHE BOOL "Prefer C++11 atomics, threads, etc. over native implementations")
set(JUNCTION_WITH_SAMPLES FALSE CACHE BOOL "Include all Junction samples in generated build system")
# Add Junction (and Turf) build targets to the CMake project.
# This defines JUNCTION_ALL_INCLUDE_DIRS and JUNCTION_ALL_LIBRARIES.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../gph/ext/junction junction)

# Add executable.
add_executable(ripser ripser.cpp)

target_compile_definitions(ripser PUBLIC NDEBUG)
target_compile_options(ripser PRIVATE -Ofast)
target_link_libraries(ripser -lbenchmark)

# Set include dirs and libraries for the executable.
include_directories(${JUNCTION_ALL_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../gph/inc/)

target_link_libraries(ripser ${JUNCTION_ALL_LIBRARIES} -lpthread)
