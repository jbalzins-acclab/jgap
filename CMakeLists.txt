cmake_minimum_required(VERSION 3.11)
project(jgap LANGUAGES CXX OBJCXX)

option(JGAP_USE_MPS "Use Apple Metal Performance Shaders for GPU acceleration"
        ${APPLE})

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")

set(${PROJECT_SOURCE_DIR}/vcpkg_installed)

find_package(Eigen3 CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)
find_package(OpenBLAS CONFIG REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/src)

file(READ "${CMAKE_SOURCE_DIR}/vcpkg.json" VCPKG_JSON_CONTENT)
string(REGEX MATCH "\"version\"[ \t]*:[ \t]*\"([^\"]+)\"" _ ${VCPKG_JSON_CONTENT})
set(PACKAGE_VERSION ${CMAKE_MATCH_1})

configure_file(
        ${CMAKE_SOURCE_DIR}/src/Version.hpp.in
        ${CMAKE_BINARY_DIR}/generated/Version.hpp
        @ONLY
)
include_directories(${CMAKE_BINARY_DIR}/generated)

# Function to generate parser registry header
function(generate_parser_registry_header)
    file(GLOB_RECURSE HPP_FILES "${PROJECT_SOURCE_DIR}/src/*.hpp")
    set(PARSER_FILES)
    foreach (HPP_FILE ${HPP_FILES})
        file(READ "${HPP_FILE}" FILE_CONTENT)
        if (FILE_CONTENT MATCHES "REGISTER_PARSER")
            list(APPEND PARSER_FILES "${HPP_FILE}")
        endif ()
    endforeach ()
    set(HEADER_CONTENT "// Auto-generated parser registry\n")
    string(APPEND HEADER_CONTENT "#ifndef PARSER_REGISTRY_AUTO_HPP\n")
    string(APPEND HEADER_CONTENT "#define PARSER_REGISTRY_AUTO_HPP\n\n")
    string(APPEND HEADER_CONTENT "// Auto-include all parser headers\n")
    foreach (PARSER_FILE ${PARSER_FILES})
        file(RELATIVE_PATH REL_PATH "${PROJECT_SOURCE_DIR}/src" "${PARSER_FILE}")
        string(APPEND HEADER_CONTENT "#include \"${REL_PATH}\"\n")
    endforeach ()
    string(APPEND HEADER_CONTENT "\n#endif // PARSER_REGISTRY_AUTO_HPP\n")
    set(GENERATED_HEADER "${PROJECT_SOURCE_DIR}/src/ParserRegistryAuto.hpp")
    file(WRITE "${GENERATED_HEADER}" "${HEADER_CONTENT}")
    add_custom_command(
            OUTPUT "${GENERATED_HEADER}"
            DEPENDS ${PARSER_FILES}
            COMMAND ${CMAKE_COMMAND} -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" -DCMAKE_CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}" -P "${CMAKE_CURRENT_SOURCE_DIR}/generate_parsers.cmake"
            COMMENT "Regenerating parser registry header"
    )
    add_custom_target(generate_parser_header DEPENDS "${GENERATED_HEADER}")
    add_dependencies(jgap_lib generate_parser_header)
    list(LENGTH PARSER_FILES PARSER_COUNT)
    message(STATUS "Found ${PARSER_COUNT} parser files for auto-registration")
endfunction()

file(GLOB_RECURSE ALL_SRC_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/app/.*")
list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/Main\\.cpp")
set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP})

# ── MPS: include .mm source files alongside the .cpp ones ─────────────────────
if(JGAP_USE_MPS)
    file(GLOB_RECURSE ALL_SRC_MM "${PROJECT_SOURCE_DIR}/src/*.mm")
    list(APPEND LIBRARY_SOURCE_FILES ${ALL_SRC_MM})
endif()

add_library(jgap_lib STATIC ${LIBRARY_SOURCE_FILES})

generate_parser_registry_header()

file(GLOB APP_SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/app/*.cpp")
add_executable(jgap ${PROJECT_SOURCE_DIR}/src/Main.cpp)

target_link_libraries(jgap_lib PRIVATE
        Eigen3::Eigen
        OpenBLAS::OpenBLAS
)

target_link_libraries(jgap_lib PUBLIC
        pugixml::static pugixml::pugixml
        yaml-cpp::yaml-cpp
        TBB::tbb TBB::tbbmalloc
        HighFive
)

target_compile_definitions(jgap_lib PRIVATE
        EIGEN_FAST_MATH
        EIGEN_USE_BLAS
)

target_link_libraries(jgap PRIVATE jgap_lib)

add_custom_command(TARGET jgap_lib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/resources
        ${CMAKE_CURRENT_BINARY_DIR}/resources
)
add_custom_command(TARGET jgap_lib POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/test/resources
        ${CMAKE_CURRENT_BINARY_DIR}/test/resources
)

# ── MPS block ─────────────────────────────────────────────────────────────────
# All fixes vs. the previous version:
#   1. Targets jgap_lib (where the .mm source lives), not jgap.
#   2. Uses plain <Metal/Metal.h> ObjC headers — no metal-cpp download needed.
#   3. OBJCXX_STANDARD set on the library so clang picks up C++23 for .mm files.
#   4. The .mm source files are already added to LIBRARY_SOURCE_FILES above;
#      set_source_files_properties marks them OBJCXX so clang uses the right
#      front-end without renaming every surrounding .cpp.
if(JGAP_USE_MPS)
    # Tell CMake to compile *.mm files as Objective-C++
    file(GLOB_RECURSE MPS_MM_FILES "${PROJECT_SOURCE_DIR}/src/*.mm")
    set_source_files_properties(${MPS_MM_FILES} PROPERTIES LANGUAGE OBJCXX)

    # Propagate the define so the header switches to the GPU path
    target_compile_definitions(jgap_lib PUBLIC JGAP_USE_MPS)

    # ARC + C++23 for Obj-C++ TUs
    target_compile_options(jgap_lib PRIVATE
            $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc -std=c++23>
    )

    # Apple frameworks — linked into the library, transitively available to jgap
    target_link_libraries(jgap_lib PUBLIC
            "-framework Metal"
            "-framework MetalPerformanceShaders"
            "-framework Foundation"
    )
endif()

## TESTS
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(FetchContent)
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
            DOWNLOAD_EXTRACT_TIMESTAMP OLD
    )
    FetchContent_MakeAvailable(googletest)
    enable_testing()

    set(TEST_SOURCES
            ${PROJECT_SOURCE_DIR}/test/core/fit/TestQRGapFit.cpp
            ${PROJECT_SOURCE_DIR}/test/core/atomic/neighbours/TestNeighbourList.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/2b/TestTwoBodyTransformation.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/3b/TestAngle3bTransformation.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/eam/TestFSGenPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/eam/TestCoscutoffPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/eam/TestPolycutoffPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/core/transform/TestTransformationAggregatorImpl.cpp
            ${PROJECT_SOURCE_DIR}/test/core/kernels/TestSquaredExpKernel.cpp
            ${PROJECT_SOURCE_DIR}/test/core/potentials/gap/component/TestNBodyGapComponent.cpp
            ${PROJECT_SOURCE_DIR}/test/core/potentials/gap/component/TestManyBodyGapComponent.cpp
            ${PROJECT_SOURCE_DIR}/test/core/cutoff/TestCosCutoff.cpp
            ${PROJECT_SOURCE_DIR}/test/core/cutoff/TestPerriotPolynomialCutoff.cpp
    )

    add_executable(jgap_tests ${TEST_SOURCES})
    target_link_libraries(jgap_tests gtest_main jgap_lib)

    include(GoogleTest)
    gtest_discover_tests(jgap_tests)
endif ()