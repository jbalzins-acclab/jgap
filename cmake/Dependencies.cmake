include(FetchContent)

set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "Don't re-download fetched dependencies")

# ------------------------------------------------------------------------------
# 1. Eigen3 (Header-only Linear Algebra)
# ------------------------------------------------------------------------------
if (APPLE)
    if (EXISTS "/opt/homebrew/include/eigen3")
        list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
    elseif (EXISTS "/usr/local/include/eigen3")
        list(APPEND CMAKE_PREFIX_PATH "/usr/local")
    endif ()
endif ()

find_package(Eigen3 3.4 QUIET)

if (NOT TARGET Eigen3::Eigen)
    set(EIGEN_BUILD_DOC OFF CACHE INTERNAL "")
    set(BUILD_TESTING OFF CACHE INTERNAL "")
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE INTERNAL "")
    FetchContent_Declare(
        Eigen3
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        SYSTEM
    )
    FetchContent_MakeAvailable(Eigen3)
endif ()

# ------------------------------------------------------------------------------
# 2. pugixml (Fast C++ XML library) via FetchContent
# ------------------------------------------------------------------------------
FetchContent_Declare(
    pugixml
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    GIT_TAG v1.15
    GIT_SHALLOW TRUE
    SYSTEM
)
FetchContent_MakeAvailable(pugixml)
set(pugixml_FOUND TRUE)
if (TARGET pugixml-static)
    target_compile_definitions(pugixml-static PUBLIC "DBL_DIG=15")
endif ()

# ------------------------------------------------------------------------------
# 3. HDF5 (Host C library) & HighFive (C++ Header-only wrapper via FetchContent)
# ------------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C)
set(HIGHFIVE_USE_BOOST OFF CACHE INTERNAL "")
set(HIGHFIVE_BUILD_DOCS OFF CACHE INTERNAL "")
set(HIGHFIVE_UNIT_TESTS OFF CACHE INTERNAL "")
FetchContent_Declare(
    HighFive
    GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
    GIT_TAG v3.0.0-beta2
    GIT_SHALLOW TRUE
    SYSTEM
)
FetchContent_MakeAvailable(HighFive)

# ------------------------------------------------------------------------------
# 4. Multi-core Parallelization (GCD on Apple, OpenMP on Linux/others)
# ------------------------------------------------------------------------------
if (APPLE)
    message(STATUS "Apple platform detected: Grand Central Dispatch (GCD) enabled for multi-core parallelization")
else ()
    find_package(OpenMP QUIET)
    if (OpenMP_CXX_FOUND)
        message(STATUS "OpenMP enabled: multi-core parallelization will be active (${OpenMP_CXX_VERSION})")
    else ()
        message(STATUS "OpenMP not found: multi-core parallelization won't be enabled")
    endif ()
endif ()

# ------------------------------------------------------------------------------
# 5. BLAS / Acceleration Libraries from Host / System
# ------------------------------------------------------------------------------
# On macOS, check for OpenBLAS (e.g. Homebrew) first since it provides both BLAS
# and LAPACKE for multi-threaded QR acceleration, unlike Apple Accelerate.
if (APPLE AND NOT OpenBLAS_FOUND)
    if (EXISTS "/opt/homebrew/opt/openblas")
        set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas")
        list(PREPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/openblas")
    elseif (EXISTS "/usr/local/opt/openblas")
        set(OpenBLAS_ROOT "/usr/local/opt/openblas")
        list(PREPEND CMAKE_PREFIX_PATH "/usr/local/opt/openblas")
    else ()
        find_program(BREW_EXECUTABLE brew)
        if (BREW_EXECUTABLE)
            execute_process(
                COMMAND ${BREW_EXECUTABLE} --prefix openblas
                OUTPUT_VARIABLE _BREW_OPENBLAS_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            if (_BREW_OPENBLAS_PREFIX AND EXISTS "${_BREW_OPENBLAS_PREFIX}")
                set(OpenBLAS_ROOT "${_BREW_OPENBLAS_PREFIX}")
                list(PREPEND CMAKE_PREFIX_PATH "${_BREW_OPENBLAS_PREFIX}")
            endif ()
        endif ()
    endif ()
endif ()

find_package(OpenBLAS CONFIG QUIET)
if (NOT OpenBLAS_FOUND)
    find_package(OpenBLAS QUIET)
endif ()

if (OpenBLAS_FOUND)
    set(JGAP_BLAS_LIBS ${OpenBLAS_LIBRARIES})
    set(JGAP_BLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
    include_directories(${OpenBLAS_INCLUDE_DIRS})
elseif (APPLE)
    set(BLA_VENDOR Apple)
    find_package(BLAS QUIET)
    if (BLAS_FOUND)
        set(JGAP_BLAS_LIBS ${BLAS_LIBRARIES})
    endif ()
    unset(BLA_VENDOR)
else ()
    find_package(BLAS QUIET)
    if (BLAS_FOUND)
        set(JGAP_BLAS_LIBS ${BLAS_LIBRARIES})
    endif ()
endif ()

if (JGAP_BLAS_LIBS)
    message(STATUS "BLAS acceleration enabled (${JGAP_BLAS_LIBS})")
else ()
    if (APPLE)
        message(STATUS "No BLAS found: using Eigen's built-in routines (GCD multi-core for loops)")
    elseif (OpenMP_CXX_FOUND)
        message(STATUS "No BLAS found: using Eigen's built-in routines with OpenMP multi-threading")
    else ()
        message(STATUS "No BLAS and no OpenMP found: Eigen will run single-threaded")
    endif ()
endif ()

# ------------------------------------------------------------------------------
# 6. LAPACKE Acceleration for Eigen (e.g. Multi-threaded QR / Cholesky)
# ------------------------------------------------------------------------------
include(CheckCXXSourceCompiles)

find_path(LAPACKE_INCLUDE_DIR
    NAMES lapacke.h
    PATH_SUFFIXES openblas lapack lapacke
    HINTS
        ${OpenBLAS_INCLUDE_DIRS}
        ${OpenBLAS_ROOT}/include
        /opt/homebrew/opt/openblas/include
        /usr/local/opt/openblas/include
    PATHS
        /usr/include
        /usr/include/openblas
        /usr/local/include
        /usr/local/include/openblas
        /opt/homebrew/include
)

set(JGAP_LAPACK_LIBS "")
set(LAPACKE_FOUND FALSE)
if (LAPACKE_INCLUDE_DIR)
    set(_BACKUP_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(_BACKUP_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${LAPACKE_INCLUDE_DIR})

    # Stage 1: Check if BLAS library already provides LAPACKE (e.g. OpenBLAS with bundled LAPACKE)
    set(CMAKE_REQUIRED_LIBRARIES ${JGAP_BLAS_LIBS})
    check_cxx_source_compiles("
    #include <lapacke.h>
    int main() {
        double a[1] = {0.0};
        double tau[1] = {0.0};
        return LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 1, 1, a, 1, tau);
    }
    " LAPACKE_WORKS_BLAS_ALONE)

    if (LAPACKE_WORKS_BLAS_ALONE)
        set(LAPACKE_FOUND TRUE)
        message(STATUS "LAPACKE acceleration enabled via BLAS library (${LAPACKE_INCLUDE_DIR})")
    else ()
        # Stage 2: BLAS library does not bundle LAPACKE; search for standalone liblapacke
        find_library(LAPACKE_LIBRARY
            NAMES lapacke
            PATH_SUFFIXES lapack
        )
        if (LAPACKE_LIBRARY)
            set(CMAKE_REQUIRED_LIBRARIES ${JGAP_BLAS_LIBS} ${LAPACKE_LIBRARY})
            check_cxx_source_compiles("
            #include <lapacke.h>
            int main() {
                double a[1] = {0.0};
                double tau[1] = {0.0};
                return LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 1, 1, a, 1, tau);
            }
            " LAPACKE_WORKS_WITH_LIB)
            if (LAPACKE_WORKS_WITH_LIB)
                set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY})
                set(LAPACKE_FOUND TRUE)
                message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
            endif ()
        endif ()
    endif ()

    set(CMAKE_REQUIRED_INCLUDES ${_BACKUP_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${_BACKUP_CMAKE_REQUIRED_LIBRARIES})

    if (NOT LAPACKE_FOUND)
        message(STATUS "LAPACKE header found at ${LAPACKE_INCLUDE_DIR} but linking failed: Eigen QR will use built-in routines")
    endif ()
else ()
    message(STATUS "LAPACKE not found: Eigen QR and Cholesky will use built-in routines")
endif ()
