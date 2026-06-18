cmake_minimum_required(VERSION 3.11)
project(jgap)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED OFF)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")

find_package(Eigen3 CONFIG REQUIRED)
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
    # Find all .hpp files in include directory
    file(GLOB_RECURSE HPP_FILES "${PROJECT_SOURCE_DIR}/src/*.hpp")

    # Filter for files that likely contain parsers
    set(PARSER_FILES)
    foreach (HPP_FILE ${HPP_FILES})
        # Read file content
        file(READ "${HPP_FILE}" FILE_CONTENT)

        # Check if file contains REGISTER_PARSER macro
        if (FILE_CONTENT MATCHES "REGISTER_PARSER")
            list(APPEND PARSER_FILES "${HPP_FILE}")
        endif ()
    endforeach ()

    # Generate the header content
    set(HEADER_CONTENT "// Auto-generated parser registry\n")
    string(APPEND HEADER_CONTENT "#ifndef PARSER_REGISTRY_AUTO_HPP\n")
    string(APPEND HEADER_CONTENT "#define PARSER_REGISTRY_AUTO_HPP\n\n")
    string(APPEND HEADER_CONTENT "// Auto-include all parser headers\n")

    foreach (PARSER_FILE ${PARSER_FILES})
        file(RELATIVE_PATH REL_PATH "${PROJECT_SOURCE_DIR}/src" "${PARSER_FILE}")
        string(APPEND HEADER_CONTENT "#include \"${REL_PATH}\"\n")
    endforeach ()

    string(APPEND HEADER_CONTENT "\n#endif // PARSER_REGISTRY_AUTO_HPP\n")

    # Write the generated header
    set(GENERATED_HEADER "${PROJECT_SOURCE_DIR}/src/ParserRegistryAuto.hpp")
    file(WRITE "${GENERATED_HEADER}" "${HEADER_CONTENT}")

    # Make sure the header is regenerated when source files change
    add_custom_command(
            OUTPUT "${GENERATED_HEADER}"
            DEPENDS ${PARSER_FILES}
            COMMAND ${CMAKE_COMMAND} -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" -DCMAKE_CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}" -P "${CMAKE_CURRENT_SOURCE_DIR}/generate_parsers.cmake"
            COMMENT "Regenerating parser registry header"
    )

    # Add the generated header as a dependency
    add_custom_target(generate_parser_header DEPENDS "${GENERATED_HEADER}")
    add_dependencies(jgap_lib generate_parser_header)

    # Report how many parsers were found
    list(LENGTH PARSER_FILES PARSER_COUNT)
    message(STATUS "Found ${PARSER_COUNT} parser files for auto-registration")
endfunction()

file(GLOB_RECURSE ALL_SRC_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE ALL_SRC_C "${PROJECT_SOURCE_DIR}/src/*.c")

list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/Main\\.cpp")

set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP} ${ALL_SRC_C})
# Shared so that every serializer's REGISTER_SERIALIZATION static initializer is linked into the library
# and runs on load (a static archive would drop the unreferenced translation units), and so type_info is
# exported once and dynamic_cast (ValuePtr::as<>) works across the library boundary.
add_library(jgap_lib SHARED ${LIBRARY_SOURCE_FILES})
# Emit the library as libjgap.dylib (clients link -ljgap); the CMake target stays jgap_lib to avoid
# clashing with the jgap executable target.
set_target_properties(jgap_lib PROPERTIES OUTPUT_NAME jgap)

# Generate parser registry header
generate_parser_registry_header()

file(GLOB APP_SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/app/*.cpp")

add_executable(jgap
        ${PROJECT_SOURCE_DIR}/src/Main.cpp
        #${APP_SOURCE_FILES}
)

target_link_libraries(jgap_lib PRIVATE
        Eigen3::Eigen
        OpenBLAS::OpenBLAS
)

target_link_libraries(jgap_lib PUBLIC
        pugixml::static pugixml::pugixml
        TBB::tbb TBB::tbbmalloc
        HighFive
)

target_compile_definitions(jgap_lib PRIVATE
        EIGEN_FAST_MATH
        EIGEN_USE_BLAS
)

target_link_libraries(jgap PRIVATE
        jgap_lib
)

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
            ${PROJECT_SOURCE_DIR}/test/unit/core/fit/TestQRGapFit.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/atomic/neighbours/TestNeighbourList.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/2b/TestTwoBodyTransformation.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/3b/TestAngle3bTransformation.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/eam/TestFSGenPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/eam/TestCoscutoffPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/eam/TestPolycutoffPairFunction.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/transform/TestTransformationAggregatorImpl.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/kernels/TestSquaredExpKernel.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/potentials/gap/component/TestNBodyGapComponent.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/potentials/gap/component/TestManyBodyGapComponent.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/cutoff/TestCosCutoff.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/cutoff/TestPerriotPolynomialCutoff.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/potentials/zbl/TestZblPotential.cpp
            ${PROJECT_SOURCE_DIR}/test/unit/core/splines/TestSplines.cpp
    )

    add_executable(jgap_tests ${TEST_SOURCES})
    target_link_libraries(jgap_tests
            gtest_main
            jgap_lib
    )

    include(GoogleTest)
    gtest_discover_tests(jgap_tests)
endif ()