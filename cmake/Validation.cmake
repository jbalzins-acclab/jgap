## VALIDATION BENCHMARK SUITE
## Enabled via -DJGAP_BUILD_VALIDATION=ON

option(JGAP_BUILD_VALIDATION "Build jgap automated validation benchmark suite" OFF)

if (JGAP_BUILD_VALIDATION)
    message(STATUS "Configuring automated validation suite (JGAP_BUILD_VALIDATION=ON)")
    add_subdirectory(test/validation)
endif ()
