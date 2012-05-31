# cmake -DCMAKE_BUILD_TYPE=Release -DAEON_INCLUDE_DIR=../../src/ -DAEON_BUILD_DIR=../../build/ ../src

# CMake Required Version
cmake_minimum_required (VERSION 2.6)

# Project Name [EDIT]
# The project name must start with AeonExperiment
project (AeonExperimentDF)

# Add customized CMAKE modules
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${AEON_INCLUDE_DIR}/cmake/")

# Compiler options
add_definitions(-ansi -pedantic -Wno-long-long -Wall -Wextra -Wformat=2
                -Wredundant-decls -Wwrite-strings -Wmissing-include-dirs 
                -Wswitch-enum -Wuninitialized 
                -Wswitch-default -Winit-self -Wshadow -Wfloat-equal -fPIC)

# Debug Flags 
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb3 -DDEBUG ")

# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -finline-functions -fomit-frame-pointer -DLOGGING_DEBUG_DISABLED -DNDEBUG")

# Include TCMalloc
if (WITH_PROFILER)
    set (WITH_TCMALLOC 1)
endif()
if (WITH_TCMALLOC)
    find_package(Tcmalloc REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
    link_directories(${TCMALLOC_LIBRARIES})
    set(GOOGLE_LIBRARIES "tcmalloc")
    if (WITH_PROFILER)
        link_directories(${PROFILER_LIBRARIES})
        set(GOOGLE_LIBRARIES "tcmalloc profiler")
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
       	include(OptimizeForArchitecture)
       	OptimizeForArchitecture("auto")
endif ()


# Include dependencies
find_package(Boost 1.47.0)
include_directories("${AEON_INCLUDE_DIR}/framework/shared/" "${AEON_INCLUDE_DIR}/libraries/json-cpp/" "${AEON_INCLUDE_DIR}/libraries/utils/" ${Boost_INCLUDE_DIRS})
link_directories("${AEON_BUILD_DIR}/framework" ${Boost_LIBRARY_DIRS})

# Get the list of all the source code
file (GLOB_RECURSE
     EXPERIMENT_SRC
     *.h
     *.cpp)

# Create the library 
# [EDIT] edit the experiment library name "my-experiment"
# The library name must start with aeon-experiment
add_library(aeon-experiment-df SHARED ${EXPERIMENT_SRC})
target_link_libraries(aeon-experiment-df aeon-shared)
