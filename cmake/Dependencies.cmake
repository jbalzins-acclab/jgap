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
if (APPLE)
    # On macOS, native Apple Accelerate framework provides peak AMX hardware-accelerated BLAS
    set(BLA_VENDOR Apple)
    find_package(BLAS QUIET)
    if (BLAS_FOUND)
        set(JGAP_BLAS_LIBS ${BLAS_LIBRARIES})
    endif ()
    unset(BLA_VENDOR)
else ()
    # Prioritize OpenMP-threaded OpenBLAS (openblaso) or pthreads (openblasp) over serial (openblas)
    find_library(JGAP_OPENBLAS_LIB
        NAMES openblaso openblasp openblas
        HINTS
            ${OpenBLAS_ROOT}/lib
            /usr/lib64
            /usr/lib
            /usr/local/lib
    )
    if (JGAP_OPENBLAS_LIB)
        set(JGAP_BLAS_LIBS ${JGAP_OPENBLAS_LIB})
        set(OpenBLAS_FOUND TRUE)
    else ()
        find_package(OpenBLAS CONFIG QUIET)
        if (NOT OpenBLAS_FOUND)
            find_package(OpenBLAS QUIET)
        endif ()

        if (OpenBLAS_FOUND)
            set(JGAP_BLAS_LIBS ${OpenBLAS_LIBRARIES})
            set(JGAP_BLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
            include_directories(${OpenBLAS_INCLUDE_DIRS})
        else ()
            find_package(BLAS QUIET)
            if (BLAS_FOUND)
                set(JGAP_BLAS_LIBS ${BLAS_LIBRARIES})
            endif ()
        endif ()
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
# 6. LAPACKE Acceleration for Eigen (Linux only)
# ------------------------------------------------------------------------------
# On macOS, native Apple Accelerate framework combined with Eigen's inlined blocked
# QR is up to 73% faster than external LAPACK without bundling external reference BLAS.
set(JGAP_LAPACK_LIBS "")
set(LAPACKE_FOUND FALSE)

if (APPLE)
    message(STATUS "macOS detected: native Apple Accelerate enabled for AMX hardware-accelerated BLAS with Eigen blocked QR")
else ()
    function(check_lapacke_candidate OUT_VAR INCLUDE_DIR LIBS)
        set(TEST_SRC "
#include <complex>
#ifndef lapack_complex_float
#define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
#define lapack_complex_double std::complex<double>
#endif
#include <lapacke.h>
int main() {
    double a[1] = {0.0};
    double tau[1] = {0.0};
    return LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 1, 1, a, 1, tau);
}
")
        file(WRITE "${CMAKE_BINARY_DIR}/lapacke_test.cpp" "${TEST_SRC}")
        try_compile(TEST_OK
            SOURCES "${CMAKE_BINARY_DIR}/lapacke_test.cpp"
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${INCLUDE_DIR}"
            LINK_LIBRARIES ${LIBS}
            OUTPUT_VARIABLE TEST_OUTPUT
        )
        if (TEST_OK)
            set(${OUT_VAR} TRUE PARENT_SCOPE)
        else ()
            set(${OUT_VAR} FALSE PARENT_SCOPE)
            if (JGAP_DEBUG_LAPACKE)
                message(STATUS "LAPACKE candidate test (${LIBS}) failed:\n${TEST_OUTPUT}")
            endif ()
        endif ()
    endfunction()

    find_path(LAPACKE_INCLUDE_DIR
        NAMES lapacke.h
        PATH_SUFFIXES lapacke lapack openblas
        HINTS
            ${OpenBLAS_INCLUDE_DIRS}
            ${OpenBLAS_ROOT}/include
        PATHS
            /usr/include/lapacke
            /usr/include/lapack
            /usr/include
            /usr/include/openblas
            /usr/local/include
            /usr/local/include/openblas
    )

    if (LAPACKE_INCLUDE_DIR)
        # Stage 1: Check if BLAS library already provides LAPACKE (e.g. Linux OpenBLAS built with bundled LAPACKE)
        if (JGAP_BLAS_LIBS)
            check_lapacke_candidate(LAPACKE_WORKS_BLAS_ALONE "${LAPACKE_INCLUDE_DIR}" "${JGAP_BLAS_LIBS}")
            if (LAPACKE_WORKS_BLAS_ALONE)
                set(LAPACKE_FOUND TRUE)
                message(STATUS "LAPACKE acceleration enabled via BLAS library (${LAPACKE_INCLUDE_DIR})")
            endif ()
        endif ()

        # Stage 2: BLAS library does not bundle LAPACKE; search for standalone liblapacke and liblapack
        if (NOT LAPACKE_FOUND)
            find_library(LAPACKE_LIBRARY
                NAMES lapacke
                PATH_SUFFIXES lapack
                HINTS
                    /usr/lib64
                    /usr/lib
                    /usr/local/lib
            )
            find_library(LAPACK_LIBRARY
                NAMES lapack
                PATH_SUFFIXES lapack
                HINTS
                    /usr/lib64
                    /usr/lib
                    /usr/local/lib
            )
            find_library(GFORTRAN_LIBRARY
                NAMES gfortran
                HINTS
                    /usr/lib64
                    /usr/lib
                    /usr/local/lib
            )

            if (LAPACKE_LIBRARY)
                # Try 1: LAPACKE + LAPACK + BLAS
                if (NOT LAPACKE_FOUND AND LAPACK_LIBRARY)
                    check_lapacke_candidate(TRY_1 "${LAPACKE_INCLUDE_DIR}" "${LAPACKE_LIBRARY};${LAPACK_LIBRARY};${JGAP_BLAS_LIBS}")
                    if (TRY_1)
                        set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY})
                        set(LAPACKE_FOUND TRUE)
                        message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY} and ${LAPACK_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
                    endif ()
                endif ()

                # Try 2: LAPACKE + LAPACK + BLAS + gfortran
                if (NOT LAPACKE_FOUND AND LAPACK_LIBRARY AND GFORTRAN_LIBRARY)
                    check_lapacke_candidate(TRY_2 "${LAPACKE_INCLUDE_DIR}" "${LAPACKE_LIBRARY};${LAPACK_LIBRARY};${JGAP_BLAS_LIBS};${GFORTRAN_LIBRARY}")
                    if (TRY_2)
                        set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} ${GFORTRAN_LIBRARY})
                        set(LAPACKE_FOUND TRUE)
                        message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY}, ${LAPACK_LIBRARY}, and ${GFORTRAN_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
                    endif ()
                endif ()

                # Try 3: LAPACKE + BLAS
                if (NOT LAPACKE_FOUND)
                    check_lapacke_candidate(TRY_3 "${LAPACKE_INCLUDE_DIR}" "${LAPACKE_LIBRARY};${JGAP_BLAS_LIBS}")
                    if (TRY_3)
                        set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY})
                        set(LAPACKE_FOUND TRUE)
                        message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
                    endif ()
                endif ()

                # Try 4: Standalone LAPACKE + LAPACK without BLAS
                if (NOT LAPACKE_FOUND AND LAPACK_LIBRARY)
                    check_lapacke_candidate(TRY_4 "${LAPACKE_INCLUDE_DIR}" "${LAPACKE_LIBRARY};${LAPACK_LIBRARY}")
                    if (TRY_4)
                        set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY})
                        set(LAPACKE_FOUND TRUE)
                        message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY} and ${LAPACK_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
                    endif ()
                endif ()

                # Try 5: Standalone LAPACKE alone
                if (NOT LAPACKE_FOUND)
                    check_lapacke_candidate(TRY_5 "${LAPACKE_INCLUDE_DIR}" "${LAPACKE_LIBRARY}")
                    if (TRY_5)
                        set(JGAP_LAPACK_LIBS ${LAPACKE_LIBRARY})
                        set(LAPACKE_FOUND TRUE)
                        message(STATUS "LAPACKE acceleration enabled via ${LAPACKE_LIBRARY} (${LAPACKE_INCLUDE_DIR})")
                    endif ()
                endif ()
            endif ()
        endif ()

        if (NOT LAPACKE_FOUND)
            message(STATUS "LAPACKE header found at ${LAPACKE_INCLUDE_DIR} but linking failed: Eigen QR will use built-in routines")
        endif ()
    else ()
        message(STATUS "LAPACKE not found: Eigen QR and Cholesky will use built-in routines")
    endif ()
endif ()
