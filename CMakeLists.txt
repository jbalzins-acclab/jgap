cmake_minimum_required(VERSION 3.20)
project(jgap)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED OFF)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native")

# Auto-detect pre-installed vcpkg dependencies in build / IDE trees if CMAKE_PREFIX_PATH is not specified
if (NOT CMAKE_PREFIX_PATH)
    file(GLOB _vcpkg_trees "${CMAKE_BINARY_DIR}/vcpkg_installed/*" "${CMAKE_SOURCE_DIR}/cmake-build-debug/vcpkg_installed/*")
    if (_vcpkg_trees)
        list(APPEND CMAKE_PREFIX_PATH ${_vcpkg_trees})
    endif ()
endif ()

find_package(Eigen3 CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)

find_package(TBB CONFIG QUIET)
if (TBB_FOUND)
    message(STATUS "Found TBB: multi-core parallelization will be enabled")
else ()
    message(STATUS "TBB not found: multi-core parallelization won't be enabled")
endif ()

# BLAS accelerates Eigen (EIGEN_USE_BLAS) when available.
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
    # No BLAS: enable Eigen's OpenMP multi-threading if available
    if (APPLE AND NOT OpenMP_CXX_FOUND)
        if (EXISTS "/opt/homebrew/opt/libomp")
            set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
            list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp")
        elseif (EXISTS "/usr/local/opt/libomp")
            set(OpenMP_ROOT "/usr/local/opt/libomp")
            list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/libomp")
        endif ()
    endif ()

    find_package(OpenMP QUIET)

    # Fallback hints for Apple Clang with Homebrew libomp
    if (APPLE AND NOT OpenMP_CXX_FOUND)
        if (EXISTS "/opt/homebrew/opt/libomp")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
            find_package(OpenMP QUIET)
        elseif (EXISTS "/usr/local/opt/libomp")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
            find_package(OpenMP QUIET)
        endif ()
    endif ()

    if (OpenMP_CXX_FOUND)
        message(STATUS "No BLAS found: using Eigen's built-in routines with OpenMP multi-threading")
    else ()
        message(STATUS "No BLAS and no OpenMP found: Eigen will run single-threaded")
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

# The src/app/* executables and src/pybind/* modules are not part of the library.
list(FILTER ALL_SRC_CPP EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/app/.*")
list(FILTER ALL_SRC_CPP EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/pybind/.*")

# The XML converter depends on pugixml; drop it from the library when pugixml is unavailable.
if (NOT pugixml_FOUND)
    list(FILTER ALL_SRC_CPP EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/jgap/io/convert/Quip.*")
endif ()

set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP} ${ALL_SRC_C})
# Shared so that every serializer's REGISTER_SERIALIZATION static initializer is linked into the library
# and runs on load (a static archive would drop the unreferenced translation units), and so type_info is
# exported once and dynamic_cast (ValuePtr::as<>) works across the library boundary.
add_library(jgap_lib SHARED ${LIBRARY_SOURCE_FILES})
add_library(jgap::jgap ALIAS jgap_lib)
# Emit the library as libjgap (libjgap.so / libjgap.dylib depending on platform; clients link -ljgap);
# the CMake target stays jgap_lib to avoid clashing with the jgap executable target.
# OUTPUT_NAME -> libjgap.{so,dylib}; EXPORT_NAME -> the consumable target is jgap::jgap.
set_target_properties(jgap_lib PROPERTIES
        OUTPUT_NAME jgap
        EXPORT_NAME jgap
        CXX_VISIBILITY_PRESET default
        VISIBILITY_INLINES_HIDDEN OFF
)

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

if (JGAP_BLAS_LIBS)
    target_link_libraries(jgap_lib PRIVATE ${JGAP_BLAS_LIBS})
endif ()

# Eigen, HighFive, and pugixml are internal implementation details used only in .cpp files;
# they stay PRIVATE so downstream consumers link against libjgap without needing third-party headers.
target_link_libraries(jgap_lib PRIVATE
        Eigen3::Eigen
        HighFive
)

if (pugixml_FOUND)
    target_link_libraries(jgap_lib PRIVATE pugixml::static pugixml::pugixml)
endif ()

if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Resolve libc++ symbols required by pre-compiled static dependencies when building with GCC on macOS
    target_link_libraries(jgap_lib PRIVATE c++)
endif ()

if (TBB_FOUND)
    target_compile_definitions(jgap_lib PRIVATE HAS_TBB)
    target_link_libraries(jgap_lib PRIVATE TBB::tbb TBB::tbbmalloc)
endif ()

find_package(ZLIB QUIET)
if (ZLIB_FOUND)
    target_link_libraries(jgap_lib PRIVATE ZLIB::ZLIB)
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
    target_link_libraries(jgap_convert PRIVATE jgap_lib pugixml::static pugixml::pugixml)
endif ()

## PYTHON / PYBIND11 MODULE (jgap_ase)
execute_process(
    COMMAND python3 -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if (PYBIND11_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
endif ()

find_package(pybind11 CONFIG QUIET)
if (pybind11_FOUND)
    message(STATUS "pybind11 found: building jgap_ase Python module")
    pybind11_add_module(jgap_ase "${PROJECT_SOURCE_DIR}/src/pybind/ASE/PyJGAPCalculator.cpp")
    target_link_libraries(jgap_ase PRIVATE jgap_lib)

    execute_process(
        COMMAND python3 -c "import sys; print(f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages')"
        OUTPUT_VARIABLE _PY_SITE_PACKAGES_REL
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if (NOT _PY_SITE_PACKAGES_REL)
        set(_PY_SITE_PACKAGES_REL "lib/python3/site-packages")
    endif ()

    install(TARGETS jgap_ase DESTINATION "${_PY_SITE_PACKAGES_REL}")
    if (EXISTS "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py")
        install(FILES "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py" DESTINATION "${_PY_SITE_PACKAGES_REL}")
    endif ()
else ()
    message(STATUS "pybind11 not found: skipping jgap_ase Python module build")
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
# and gets jgap's headers automatically.
install(TARGETS jgap_lib EXPORT jgapTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/jgap DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
install(FILES ${CMAKE_BINARY_DIR}/generated/Version.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jgap)

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