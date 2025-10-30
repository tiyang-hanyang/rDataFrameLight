# Minimum required CMake version
cmake_minimum_required(VERSION 3.14)

# Project name and C++ standard
project(rDataFrameLight)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE Debug)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(ROOT_DIR $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})
# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree Gpad ROOTDataFrame)

include(${ROOT_USE_FILE})

add_compile_options(-std=c++17)

include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/rDataFrameLight)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)

file(GLOB PROJECT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Root/*.cxx)
add_library(projectlib STATIC ${PROJECT_SOURCES})
target_link_libraries(projectlib stdc++fs)

# Create executable target
add_executable(skimSamples utils/skimSamples.cxx)
add_executable(collectHists utils/collectHists.cxx)
add_executable(plotHists utils/plotHists.cxx)
add_executable(appendArray utils/appendArray.cxx)

# Link against ROOT libraries
target_link_libraries(skimSamples projectlib ${ROOT_LIBRARIES})
target_link_libraries(collectHists projectlib ${ROOT_LIBRARIES})
target_link_libraries(plotHists projectlib ${ROOT_LIBRARIES})
target_link_libraries(appendArray projectlib ${ROOT_LIBRARIES})

# for simple install
file(WRITE ${CMAKE_BINARY_DIR}/setup.sh
"#!/bin/bash\n"
"export PATH=\"${CMAKE_BINARY_DIR}/bin:\$PATH\"\n"
"echo \"Environment setup: ${CMAKE_BINARY_DIR}/bin added to PATH\"\n"
)
