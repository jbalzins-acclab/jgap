cmake_minimum_required(VERSION 3.11)
project(standard_fit CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Linking jgap::jgap is all that is needed: it brings in jgap's headers and its public dependencies.
# Point CMAKE_PREFIX_PATH at jgap (build or install dir) and at its dependencies — see ../ReadMe.md.
find_package(jgap CONFIG REQUIRED)

add_executable(standard_fit StandardFit.cpp)
target_link_libraries(standard_fit PRIVATE jgap::jgap)
