cmake_minimum_required(VERSION 3.30)
project(jgap)

set(CMAKE_TOOLCHAIN_FILE "$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
# cmake -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -DNDEBUG")
#set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
#set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

set(${PROJECT_SOURCE_DIR}/vcpkg_installed)

find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(ParallelSTL CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
#find_package(Boost REQUIRED COMPONENTS filesystem regex property_tree)

include_directories(${PROJECT_SOURCE_DIR}/include)

#[===[
set(LIBRARY_SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/src/core/JgapFit.cpp
        ${PROJECT_SOURCE_DIR}/src/io/log/SpdlogWrapper.cpp
        ${PROJECT_SOURCE_DIR}/src/io/params/JsonParamReader.cpp
        ${PROJECT_SOURCE_DIR}/src/core/potentials/QuipPotential.cpp
        ${PROJECT_SOURCE_DIR}/src/utils/utils.cpp
        ${PROJECT_SOURCE_DIR}/src/io/log/StdoutLogger.cpp
        ${PROJECT_SOURCE_DIR}/src/core/neighbours/NeighbourFinder.cpp
        ${PROJECT_SOURCE_DIR}/src/core/descriptors/TwoBodyDescriptor.cpp
        ${PROJECT_SOURCE_DIR}/src/core/descriptors/ThreeBodyDescriptor.cpp
        ${PROJECT_SOURCE_DIR}/src/core/descriptors/EamDescriptor.cpp
        ${PROJECT_SOURCE_DIR}/src/core/matrices/SimpleMatrixFactory.cpp
        ${PROJECT_SOURCE_DIR}/src/core/CutoffFunction.cpp
        ${PROJECT_SOURCE_DIR}/src/core/kernel/TwoBodySE.cpp
        ${PROJECT_SOURCE_DIR}/src/core/kernel/ThreeBodySE.cpp
        ${PROJECT_SOURCE_DIR}/src/core/kernel/EamKernelSE.cpp
)
]===]

file(GLOB_RECURSE ALL_SRC_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/main\\.cpp$")
set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP})
add_library(jgap STATIC ${LIBRARY_SOURCE_FILES})

set(APP_SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/src/main.cpp
)

add_executable(jgap_app ${APP_SOURCE_FILES})

target_link_libraries(jgap PRIVATE
        TBB::tbb TBB::tbbmalloc
        pstl::ParallelSTL
        nlohmann_json::nlohmann_json
        spdlog::spdlog
        Eigen3::Eigen
        #Boost::filesystem
        #Boost::regex
)

target_link_libraries(jgap_app PRIVATE
        jgap
        nlohmann_json::nlohmann_json
)

# TODO doesn't work in CLion's cmake-build-debug ?
add_custom_command(TARGET jgap_app POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/resources
        ${CMAKE_CURRENT_BINARY_DIR}/resources
)

## TESTS

if(TRUE)
    include(FetchContent)
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    FetchContent_MakeAvailable(googletest)
    enable_testing()

    set(TEST_SOURCES
            ${PROJECT_SOURCE_DIR}/test/core/TestInRamJgap.cpp
            ${PROJECT_SOURCE_DIR}/test/core/neighbours/TestNeighbourList.cpp
            ${PROJECT_SOURCE_DIR}/test/core/descriptors/TestTwoBodyDescriptor.cpp
            ${PROJECT_SOURCE_DIR}/test/core/descriptors/TestThreeBodyDescriptor.cpp
            ${PROJECT_SOURCE_DIR}/test/core/descriptors/TestEamDescriptor.cpp
            ${PROJECT_SOURCE_DIR}/test/core/descriptors/eam/TestEamDensityCalculator.cpp
            ${PROJECT_SOURCE_DIR}/test/core/potentials/TestZblPotential.cpp
            ${PROJECT_SOURCE_DIR}/test/core/potentials/TestIsolatedAtomPotential.cpp
    )

    add_executable(jgap_tests ${TEST_SOURCES})
    target_link_libraries(jgap_tests
            gtest_main
            jgap  # link your library into tests
            nlohmann_json::nlohmann_json
    )

    include(GoogleTest)
    gtest_discover_tests(jgap_tests)

    add_custom_command(TARGET jgap_tests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/test/resources
            ${CMAKE_CURRENT_BINARY_DIR}/test/resources
    )
endif ()