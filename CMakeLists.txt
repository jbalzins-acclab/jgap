cmake_minimum_required(VERSION 3.25)
project(jgap VERSION 7.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED OFF)

set(PACKAGE_VERSION ${PROJECT_VERSION})

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# 1. External & Fetched Dependencies
include(Dependencies)

# 2. Version Header Generation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

include_directories(${PROJECT_SOURCE_DIR}/src)

configure_file(
        ${CMAKE_SOURCE_DIR}/src/Version.hpp.in
        ${CMAKE_BINARY_DIR}/generated/Version.hpp
        @ONLY
)
include_directories(${CMAKE_BINARY_DIR}/generated)

# 3. Source Collection & Target Definitions
file(GLOB_RECURSE ALL_SRC_CPP "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE ALL_SRC_C "${PROJECT_SOURCE_DIR}/src/*.c")

# Exclude standalone executables and python modules from the main shared library
list(FILTER ALL_SRC_CPP EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/app/.*")
list(FILTER ALL_SRC_CPP EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/pybind/.*")

set(LIBRARY_SOURCE_FILES ${ALL_SRC_CPP} ${ALL_SRC_C})

add_library(jgap_lib SHARED ${LIBRARY_SOURCE_FILES})
add_library(jgap::jgap ALIAS jgap_lib)

set_target_properties(jgap_lib PROPERTIES
        OUTPUT_NAME jgap
        EXPORT_NAME jgap
        CXX_VISIBILITY_PRESET default
        VISIBILITY_INLINES_HIDDEN OFF
)

# Hard guard: jgap MUST be a shared library.
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

target_link_libraries(jgap_lib PRIVATE
        Eigen3::Eigen
        HighFive
        pugixml::pugixml
)

if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Resolve libc++ symbols required by pre-compiled static dependencies when building with GCC on macOS
    target_link_libraries(jgap_lib PRIVATE c++)
endif ()

if (OpenMP_CXX_FOUND)
    target_compile_definitions(jgap_lib PUBLIC HAS_OPENMP)
    target_link_libraries(jgap_lib PUBLIC OpenMP::OpenMP_CXX)
endif ()

# Public include dirs so consumers of the target get jgap's headers for free.
target_include_directories(jgap_lib PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>"
        "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_compile_definitions(jgap_lib PRIVATE EIGEN_FAST_MATH EIGEN_NO_DEBUG)
if (JGAP_BLAS_LIBS)
    target_compile_definitions(jgap_lib PRIVATE EIGEN_USE_BLAS)
endif ()

target_link_libraries(jgap PRIVATE jgap_lib)

# QUIP XML -> .h5 converter app
add_executable(jgap_convert ${PROJECT_SOURCE_DIR}/src/app/ConvertApp.cpp)
target_link_libraries(jgap_convert PRIVATE jgap_lib pugixml::pugixml)

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

# 4. Python Bindings
include(PythonBindings)

# 5. Installation Rules
include(Installation)

# 6. Unit Tests (Debug Mode)
include(Testing)