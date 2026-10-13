## PYTHON / PYBIND11 MODULE (_jgap & jgap package)

option(JGAP_BUILD_PYTHON "Build Python bindings (_jgap module and jgap package)" ON)
if (NOT JGAP_BUILD_PYTHON)
    message(STATUS "Python bindings disabled (JGAP_BUILD_PYTHON=OFF): skipping Python module build")
    return()
endif ()

# Prioritize active virtual environment (venv / conda)
set(Python3_FIND_VIRTUALENV FIRST)

find_package(Python3 COMPONENTS Interpreter Development.Module OPTIONAL_COMPONENTS Development QUIET)
if (NOT Python3_FOUND)
    find_package(Python3 COMPONENTS Interpreter QUIET)
endif ()

if (Python3_FOUND)
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if (PYBIND11_CMAKE_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
    endif ()

    find_package(pybind11 CONFIG QUIET)
    if (pybind11_FOUND)
        message(STATUS "pybind11 found for Python ${Python3_VERSION} (${Python3_EXECUTABLE}): building _jgap Python module")
        if (SKBUILD)
            if (APPLE)
                if (NOT JGAP_BLAS_LIBS)
                    message(FATAL_ERROR "Apple Accelerate BLAS acceleration is REQUIRED for macOS Python wheels (SKBUILD=ON), but BLAS was not found!")
                endif ()
            else ()
                if (NOT LAPACKE_FOUND)
                    message(FATAL_ERROR "LAPACKE acceleration is REQUIRED for Linux Python wheels (SKBUILD=ON), but LAPACKE was not found!")
                endif ()
            endif ()
            # Monolithic Python wheel module: compile library sources directly into _jgap to ensure a single,
            # self-contained binary with a single unified C++ runtime, preventing cross-DSO allocator collisions.
            pybind11_add_module(_jgap "${PROJECT_SOURCE_DIR}/src/pybind/PyJGAP.cpp" ${LIBRARY_SOURCE_FILES})
            target_link_libraries(_jgap PRIVATE
                Eigen3::Eigen
                HighFive
                pugixml::pugixml
                HDF5::HDF5
            )
            if (JGAP_LAPACK_LIBS)
                target_link_libraries(_jgap PRIVATE ${JGAP_LAPACK_LIBS})
            endif ()
            if (JGAP_BLAS_LIBS)
                target_link_libraries(_jgap PRIVATE ${JGAP_BLAS_LIBS})
            endif ()
            if (APPLE)
                target_compile_definitions(_jgap PUBLIC HAS_GCD)
            elseif (OpenMP_CXX_FOUND)
                target_compile_definitions(_jgap PUBLIC HAS_OPENMP)
                target_link_libraries(_jgap PUBLIC OpenMP::OpenMP_CXX)
            endif ()
            target_include_directories(_jgap PRIVATE
                "${PROJECT_SOURCE_DIR}/src"
                "${CMAKE_BINARY_DIR}/generated"
                ${HDF5_INCLUDE_DIRS}
            )
            target_compile_definitions(_jgap PRIVATE EIGEN_FAST_MATH EIGEN_NO_DEBUG "JGAP_VERSION=\"${PROJECT_VERSION}\"")
            if (JGAP_BLAS_LIBS)
                target_compile_definitions(_jgap PRIVATE EIGEN_USE_BLAS)
            endif ()
            if (LAPACKE_FOUND)
                target_include_directories(_jgap PRIVATE ${LAPACKE_INCLUDE_DIR})
                target_compile_definitions(_jgap PRIVATE EIGEN_USE_LAPACKE "lapack_complex_float=std::complex<float>" "lapack_complex_double=std::complex<double>")
            endif ()

            if (UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
                target_link_options(_jgap PRIVATE "-static-libstdc++" "-Wl,--exclude-libs,ALL")
            endif ()

            # scikit-build-core packaging for PyPI wheels: keep RPATH strictly local within package
            if (APPLE)
                set(_WHEEL_RPATH "@loader_path;@loader_path/.dylibs")
            else ()
                set(_WHEEL_RPATH "$ORIGIN;$ORIGIN/.libs")
            endif ()

            set_target_properties(_jgap PROPERTIES
                OUTPUT_NAME _jgap
                BUILD_WITH_INSTALL_RPATH TRUE
                INSTALL_RPATH "${_WHEEL_RPATH}"
            )
            install(TARGETS _jgap DESTINATION "jgap")
            install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/jgap" DESTINATION "."
                    FILES_MATCHING PATTERN "*.py")
        else ()
            pybind11_add_module(_jgap "${PROJECT_SOURCE_DIR}/src/pybind/PyJGAP.cpp")
            target_link_libraries(_jgap PRIVATE jgap_lib)
            target_compile_definitions(_jgap PRIVATE "JGAP_VERSION=\"${PROJECT_VERSION}\"")
            if (LAPACKE_FOUND)
                target_include_directories(_jgap PRIVATE ${LAPACKE_INCLUDE_DIR})
                target_compile_definitions(_jgap PRIVATE EIGEN_USE_LAPACKE "lapack_complex_float=std::complex<float>" "lapack_complex_double=std::complex<double>")
            endif ()

            # Legacy jgap_ase module for backward compatibility
            pybind11_add_module(jgap_ase "${PROJECT_SOURCE_DIR}/src/pybind/PyJGAP.cpp")
            target_link_libraries(jgap_ase PRIVATE jgap_lib)
            target_compile_definitions(jgap_ase PRIVATE "JGAP_VERSION=\"${PROJECT_VERSION}\"")
            if (LAPACKE_FOUND)
                target_include_directories(jgap_ase PRIVATE ${LAPACKE_INCLUDE_DIR})
                target_compile_definitions(jgap_ase PRIVATE EIGEN_USE_LAPACKE "lapack_complex_float=std::complex<float>" "lapack_complex_double=std::complex<double>")
            endif ()

            # Determine destination site-packages directory
            if (Python3_SITEARCH)
                file(RELATIVE_PATH _REL_CHECK "${CMAKE_INSTALL_PREFIX}" "${Python3_SITEARCH}")
                if (NOT _REL_CHECK MATCHES "^\\.\\.")
                    set(_PY_INSTALL_DIR "${_REL_CHECK}")
                else ()
                    set(_PY_INSTALL_DIR "${Python3_SITEARCH}")
                endif ()
            else ()
                execute_process(
                    COMMAND "${Python3_EXECUTABLE}" -c "import sys; print(f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages')"
                    OUTPUT_VARIABLE _PY_SITE_PACKAGES_REL
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                )
                if (NOT _PY_SITE_PACKAGES_REL)
                    set(_PY_SITE_PACKAGES_REL "lib/python3/site-packages")
                endif ()
                set(_PY_INSTALL_DIR "${_PY_SITE_PACKAGES_REL}")
            endif ()

            if (APPLE)
                set(_PY_INSTALL_RPATH "@loader_path/../../..;@loader_path/../../lib;@loader_path/..;@loader_path;${CMAKE_INSTALL_FULL_LIBDIR}")
            else ()
                set(_PY_INSTALL_RPATH "\$ORIGIN/../../..;\$ORIGIN/../../../${CMAKE_INSTALL_LIBDIR};\$ORIGIN/..;\$ORIGIN;${CMAKE_INSTALL_FULL_LIBDIR}")
            endif ()

            set_target_properties(_jgap PROPERTIES
                OUTPUT_NAME _jgap
                LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python/jgap"
                BUILD_WITH_INSTALL_RPATH FALSE
                BUILD_RPATH "${CMAKE_BINARY_DIR}"
                INSTALL_RPATH "${_PY_INSTALL_RPATH}"
            )

            set_target_properties(jgap_ase PROPERTIES
                OUTPUT_NAME jgap_ase
                BUILD_WITH_INSTALL_RPATH FALSE
                BUILD_RPATH "${CMAKE_BINARY_DIR}"
                INSTALL_RPATH "${_PY_INSTALL_RPATH}"
            )

            install(TARGETS _jgap DESTINATION "${_PY_INSTALL_DIR}/jgap")
            install(TARGETS jgap_ase DESTINATION "${_PY_INSTALL_DIR}")
            install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/jgap" DESTINATION "${_PY_INSTALL_DIR}"
                    FILES_MATCHING PATTERN "*.py")
            if (EXISTS "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py")
                install(FILES "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py" DESTINATION "${_PY_INSTALL_DIR}")
            endif ()
        endif ()
    else ()
        message(STATUS "pybind11 not found: skipping _jgap Python module build")
    endif ()
else ()
    message(STATUS "Python3 not found: skipping _jgap Python module build")
endif ()
