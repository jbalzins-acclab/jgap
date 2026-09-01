include(FetchContent)

# ------------------------------------------------------------------------------
# 1. Eigen3 (Header-only Linear Algebra) via FetchContent
# ------------------------------------------------------------------------------
set(EIGEN_BUILD_DOC OFF CACHE INTERNAL "")
set(BUILD_TESTING OFF CACHE INTERNAL "")
set(EIGEN_BUILD_PKGCONFIG OFF CACHE INTERNAL "")
FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
    SYSTEM
)
FetchContent_MakeAvailable(Eigen3)

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
# 4. oneTBB (Multi-core Tasking) from Host / System or via FetchContent
# ------------------------------------------------------------------------------
find_package(TBB CONFIG QUIET)
if (NOT TBB_FOUND)
    find_package(TBB QUIET)
endif ()

if (NOT TBB_FOUND)
    message(STATUS "Host TBB not found, fetching oneTBB via FetchContent...")
    set(TBB_TEST OFF CACHE INTERNAL "" FORCE)
    set(TBB_STRICT OFF CACHE INTERNAL "" FORCE)
    set(TBB_EXAMPLES OFF CACHE INTERNAL "" FORCE)
    FetchContent_Declare(
        oneTBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG v2022.0.0
        GIT_SHALLOW TRUE
        SYSTEM
    )
    FetchContent_MakeAvailable(oneTBB)
    set(TBB_FOUND TRUE)
endif ()

if (TBB_FOUND)
    message(STATUS "TBB enabled: multi-core parallelization will be active")
endif ()

# ------------------------------------------------------------------------------
# 5. BLAS / OpenBLAS (Linear Algebra Acceleration) from Host / System
# ------------------------------------------------------------------------------
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
    # Fallback to OpenMP multi-threading in Eigen if available
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
