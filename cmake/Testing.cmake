## UNIT TESTS (Enabled by default in Debug, or via -DJGAP_BUILD_TESTS=ON / -DBUILD_TESTING=ON)

option(JGAP_BUILD_TESTS "Build jgap unit test suite" OFF)

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR JGAP_BUILD_TESTS OR BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP OLD
        SYSTEM
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
