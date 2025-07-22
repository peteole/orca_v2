cmake_minimum_required(VERSION 3.10)

# Project name
project(orca VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(src)

# Source files
set(SOURCES
    src/orca.cpp
)

# Create executable
add_executable(orca ${SOURCES})

# Set output directory
set_target_properties(orca PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Optional: Copy test data to build directory (if it exists)
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/test_data)
    file(COPY tests/test_data DESTINATION ${CMAKE_BINARY_DIR})
endif()

# Installation rules
install(TARGETS orca DESTINATION bin)
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/test_data)
    install(DIRECTORY tests/test_data/ DESTINATION share/orca/examples)
endif()