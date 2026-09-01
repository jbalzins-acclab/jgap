## PYTHON / PYBIND11 MODULE (_jgap & jgap package)

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
        pybind11_add_module(_jgap "${PROJECT_SOURCE_DIR}/src/pybind/PyJGAP.cpp")
        target_link_libraries(_jgap PRIVATE jgap_lib)

        # Place the built extension directly into python/jgap/ for in-tree use
        set_target_properties(_jgap PROPERTIES
            OUTPUT_NAME _jgap
            LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python/jgap"
            BUILD_WITH_INSTALL_RPATH FALSE
            BUILD_RPATH "${CMAKE_BINARY_DIR}"
            INSTALL_RPATH "@loader_path/../../lib;@loader_path/.."
        )

        # Legacy jgap_ase module for backward compatibility
        pybind11_add_module(jgap_ase "${PROJECT_SOURCE_DIR}/src/pybind/PyJGAP.cpp")
        target_link_libraries(jgap_ase PRIVATE jgap_lib)

        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import sys; print(f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages')"
            OUTPUT_VARIABLE _PY_SITE_PACKAGES_REL
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if (NOT _PY_SITE_PACKAGES_REL)
            set(_PY_SITE_PACKAGES_REL "lib/python3/site-packages")
        endif ()

        install(TARGETS _jgap DESTINATION "${_PY_SITE_PACKAGES_REL}/jgap")
        install(TARGETS jgap_ase DESTINATION "${_PY_SITE_PACKAGES_REL}")
        install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/jgap" DESTINATION "${_PY_SITE_PACKAGES_REL}"
                FILES_MATCHING PATTERN "*.py")
        if (EXISTS "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py")
            install(FILES "${PROJECT_SOURCE_DIR}/scripts/ASE/jgap_calculator.py" DESTINATION "${_PY_SITE_PACKAGES_REL}")
        endif ()
    else ()
        message(STATUS "pybind11 not found: skipping _jgap Python module build")
    endif ()
else ()
    message(STATUS "Python3 not found: skipping _jgap Python module build")
endif ()
