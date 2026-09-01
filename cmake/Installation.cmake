## INSTALL / EXPORT
# Makes jgap consumable as a package: after `cmake --install`, a downstream project only needs:
#   (1) with CMake:
#       find_package(jgap CONFIG REQUIRED)
#       target_link_libraries(my_app PRIVATE jgap::jgap)
#   (2) without CMake:
#       c++ ... -ljgap ...
# and gets jgap's headers automatically.

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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
