cmake_minimum_required(VERSION 3.25)

set(VALIDATION_COMMON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/common/ValidationUtils.cpp
)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/generated
)

# 1. QR Solver Variants Benchmark (flexible parameters, --max-m3b, etc.)
add_executable(ValidateQrVariants
    ValidateQrVariants.cpp
    ${VALIDATION_COMMON_SOURCES}
)
target_link_libraries(ValidateQrVariants PRIVATE jgap_lib)

# 2. Coefficient Verification against Reference QUIP Potentials
add_executable(ValidateQRGapFit
    ValidateQRGapFit.cpp
    ${VALIDATION_COMMON_SOURCES}
)
target_link_libraries(ValidateQRGapFit PRIVATE jgap_lib)

# 3. Test Set Prediction Verification against Reference QUIP XYZ
add_executable(ValidateGapEnergyEval
    ValidateGapEnergyEval.cpp
    ${VALIDATION_COMMON_SOURCES}
)
target_link_libraries(ValidateGapEnergyEval PRIVATE jgap_lib)

# 4. Tabulation Spline Grid Verification against Reference tabGAP Tables
add_executable(ValidateTabulation
    ValidateTabulation.cpp
    ${VALIDATION_COMMON_SOURCES}
)
target_link_libraries(ValidateTabulation PRIVATE jgap_lib HighFive)

# 5. TabGap Engine Verification against Reference LAMMPS Predictions
add_executable(ValidateTabGapEnergyEval
    ValidateTabGapEnergyEval.cpp
    ${VALIDATION_COMMON_SOURCES}
)
target_link_libraries(ValidateTabGapEnergyEval PRIVATE jgap_lib)

# 6. Physical Material Properties (Python + ASE + Elastic + JGAPCalculator)
# Compares a0, B, elastic constants C_ij, Ecoh, Evac_f against reference CSVs
if (Python3_EXECUTABLE)
    add_test(
        NAME ValidationElementalProperties
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_elemental_properties.py"
    )
    set_tests_properties(ValidationElementalProperties PROPERTIES
        LABELS "validation"
        ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/python:$ENV{PYTHONPATH}"
    )
endif ()

# CTest Registrations (Label: validation)
enable_testing()

add_test(NAME ValidationQrQuick COMMAND ValidateQrVariants --quick)
set_tests_properties(ValidationQrQuick PROPERTIES LABELS "validation")

add_test(NAME ValidationQRGapFit COMMAND ValidateQRGapFit)
set_tests_properties(ValidationQRGapFit PROPERTIES LABELS "validation")

add_test(NAME ValidationGapEnergyEval COMMAND ValidateGapEnergyEval)
set_tests_properties(ValidationGapEnergyEval PROPERTIES LABELS "validation")

add_test(NAME ValidationTabulation COMMAND ValidateTabulation)
set_tests_properties(ValidationTabulation PROPERTIES LABELS "validation")

add_test(NAME ValidationTabGapEnergyEval COMMAND ValidateTabGapEnergyEval)
set_tests_properties(ValidationTabGapEnergyEval PROPERTIES LABELS "validation")
