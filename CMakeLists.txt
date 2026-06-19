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

# BLAS only accelerates Eigen (EIGEN_USE_BLAS) and is OPTIONAL. Prefer an OpenBLAS CMake config (vcpkg /
# conda / spack ship one); otherwise fall back to CMake's generic BLAS finder (Homebrew / system /
# Accelerate / MKL); if nothing is found, Eigen uses its own (slower) routines.
find_package(OpenBLAS CONFIG QUIET)
if (OpenBLAS_FOUND)
    set(JGAP_BLAS_LIBS OpenBLAS::OpenBLAS)
else ()
    find_package(BLAS QUIET)
    if (BLAS_FOUND)
        set(JGAP_BLAS_LIBS ${BLAS_LIBRARIES})
    endif ()
endif ()
if (JGAP_BLAS_LIBS)
    message(STATUS "BLAS acceleration enabled (${JGAP_BLAS_LIBS})")
else ()
    # No BLAS: Eigen uses its own routines. Enable Eigen's OpenMP multi-threading so they are at least
    # parallel (a BLAS such as OpenBLAS is still strongly recommended for best performance).
    find_package(OpenMP QUIET)
    if (OpenMP_CXX_FOUND)
        message(STATUS "No BLAS found: using Eigen's built-in routines with OpenMP multi-threading "
                       "(installing OpenBLAS is strongly recommended for better performance)")
    else ()
        message(WARNING "No BLAS and no OpenMP found: Eigen will run single-threaded (slow). "
                        "Installing OpenBLAS (or enabling OpenMP) is strongly recommended.")
    endif ()
endif ()

# XML support (QUIP .xml conversion) is optional: it pulls in pugixml. When pugixml is not found, the
# converter sources and the convert app are simply left out of the build.
find_package(pugixml CONFIG)
if (pugixml_FOUND)
    message(STATUS "pugixml found: building with QUIP XML conversion support")
else ()
    message(STATUS "pugixml not found: building WITHOUT QUIP XML conversion support")
endif ()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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

file(GLOB_RECURSE ALL_SRC_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE ALL_SRC_C "${PROJECT_SOURCE_DIR}/src/*.c")

# The src/app/* executables are not part of the library.
list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/src/app/.*")

# The XML converter depends on pugixml; drop it from the library when pugixml is unavailable.
if (NOT pugixml_FOUND)
    list(FILTER ALL_SRC_CPP EXCLUDE REGEX ".*/io/convert/.*")
endif ()

set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP} ${ALL_SRC_C})
# Shared so that every serializer's REGISTER_SERIALIZATION static initializer is linked into the library
# and runs on load (a static archive would drop the unreferenced translation units), and so type_info is
# exported once and dynamic_cast (ValuePtr::as<>) works across the library boundary.
add_library(jgap_lib SHARED ${LIBRARY_SOURCE_FILES})
# Emit the library as libjgap (libjgap.so / libjgap.dylib depending on platform; clients link -ljgap);
# the CMake target stays jgap_lib to avoid clashing with the jgap executable target.
# OUTPUT_NAME -> libjgap.{so,dylib}; EXPORT_NAME -> the consumable target is jgap::jgap.
set_target_properties(jgap_lib PROPERTIES OUTPUT_NAME jgap EXPORT_NAME jgap)

# Hard guard: jgap MUST be a shared library. The serialization registry (REGISTER_SERIALIZATION static
# initializers + the per-type registry) and cross-boundary dynamic_cast (ValuePtr::as<>) rely on weak
# symbols coalescing to a single definition across the library/client boundary, which a static archive
# does not provide — registrations would silently vanish and casts return nullptr. Fail loudly if some
# build setting (e.g. BUILD_SHARED_LIBS=OFF + an edit to the add_library above) makes it static.
get_target_property(_jgap_lib_type jgap_lib TYPE)
if (NOT _jgap_lib_type STREQUAL "SHARED_LIBRARY")
    message(FATAL_ERROR
            "jgap_lib must be built as a SHARED library (it is currently ${_jgap_lib_type}). The "
            "serialization registry and cross-boundary dynamic_cast depend on shared-library symbol "
            "semantics; a static build breaks them silently.")
endif ()

add_executable(jgap
        ${PROJECT_SOURCE_DIR}/src/app/MainApp.cpp
)

# jgap's public headers include Eigen/TBB/HighFive (and pugixml for the converter), so those are PUBLIC:
# linking the jgap::jgap / jgap_lib target alone brings in their headers and libraries. BLAS is only used
# internally (Eigen's BLAS backend in .cpp), so it stays PRIVATE and is omitted when unavailable.
if (JGAP_BLAS_LIBS)
    target_link_libraries(jgap_lib PRIVATE ${JGAP_BLAS_LIBS})
elseif (OpenMP_CXX_FOUND)
    # Compiles jgap's sources with -fopenmp so Eigen's built-in routines run multi-threaded.
    target_link_libraries(jgap_lib PRIVATE OpenMP::OpenMP_CXX)
endif ()

target_link_libraries(jgap_lib PUBLIC
        Eigen3::Eigen
        TBB::tbb TBB::tbbmalloc
        HighFive
)

if (pugixml_FOUND)
    target_link_libraries(jgap_lib PUBLIC pugixml::static pugixml::pugixml)
endif ()

# Public include dirs so consumers of the (installed or in-tree) target get jgap's headers for free.
target_include_directories(jgap_lib PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>"
        "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_compile_definitions(jgap_lib PRIVATE EIGEN_FAST_MATH)
if (JGAP_BLAS_LIBS)
    target_compile_definitions(jgap_lib PRIVATE EIGEN_USE_BLAS)
endif ()

target_link_libraries(jgap PRIVATE
        jgap_lib
)

# QUIP XML -> .h5 converter app, built only when XML support is available (the converter it relies on is
# compiled into jgap_lib only then).
if (pugixml_FOUND)
    add_executable(jgap_convert ${PROJECT_SOURCE_DIR}/src/app/ConvertApp.cpp)
    target_link_libraries(jgap_convert PRIVATE jgap_lib)
endif ()

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

## INSTALL / EXPORT
# Makes jgap consumable as a package: after `cmake --install`, a downstream project only needs
#   find_package(jgap CONFIG REQUIRED)
#   target_link_libraries(my_app PRIVATE jgap::jgap)
# and gets jgap's headers and (transitively) Eigen/TBB/HighFive/pugixml automatically.
install(TARGETS jgap_lib EXPORT jgapTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
install(FILES ${CMAKE_BINARY_DIR}/generated/Version.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT jgapTargets
        NAMESPACE jgap::
        FILE jgapTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jgap
)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/jgapConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/jgapConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jgap
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/jgapConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/jgap
)

# Also export the build tree, so a consumer can `find_package(jgap)` against the BUILD directory
# (CMAKE_PREFIX_PATH=<jgap_build>) without installing. jgapConfig.cmake (generated above) sits in the
# build dir next to this jgapTargets.cmake and includes it.
export(EXPORT jgapTargets
        NAMESPACE jgap::
        FILE ${CMAKE_CURRENT_BINARY_DIR}/jgapTargets.cmake
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

    # All test sources are auto-discovered: any test/unit/**/Test*.cpp is picked up. CONFIGURE_DEPENDS
    # makes the build re-glob when files are added/removed, so no edit here is needed for a new test.
    file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS
            ${PROJECT_SOURCE_DIR}/test/unit/*.cpp)

    add_executable(jgap_tests ${TEST_SOURCES})
    target_link_libraries(jgap_tests
            gtest_main
            jgap_lib
    )

    include(GoogleTest)
    gtest_discover_tests(jgap_tests)
endif ()