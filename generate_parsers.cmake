# Generate parser registry header
set(HPP_FILES)
file(GLOB_RECURSE HPP_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")

# Filter for files that contain REGISTER_PARSER
set(PARSER_FILES)
foreach(HPP_FILE ${HPP_FILES})
    file(READ "${HPP_FILE}" FILE_CONTENT)

    if(FILE_CONTENT MATCHES "REGISTER_PARSER")
        list(APPEND PARSER_FILES "${HPP_FILE}")
    endif()
endforeach()

# Generate the header content
set(HEADER_CONTENT "// Auto-generated parser registry - DO NOT EDIT\n")
string(APPEND HEADER_CONTENT "#ifndef PARSER_REGISTRY_AUTO_HPP\n")
string(APPEND HEADER_CONTENT "#define PARSER_REGISTRY_AUTO_HPP\n\n")
string(APPEND HEADER_CONTENT "#include \"ParserRegistry.hpp\"\n\n")
string(APPEND HEADER_CONTENT "// Auto-include all parser headers\n")

foreach(PARSER_FILE ${PARSER_FILES})
    file(RELATIVE_PATH REL_PATH "${PROJECT_SOURCE_DIR}/include" "${PARSER_FILE}")
    string(APPEND HEADER_CONTENT "#include \"${REL_PATH}\"\n")
endforeach()

string(APPEND HEADER_CONTENT "\n#endif // PARSER_REGISTRY_AUTO_HPP\n")

# Write the generated header
set(GENERATED_HEADER "${CMAKE_CURRENT_BINARY_DIR}/ParserRegistryAuto.hpp")
file(WRITE "${GENERATED_HEADER}" "${HEADER_CONTENT}")

list(LENGTH PARSER_FILES PARSER_COUNT)
message(STATUS "Generated parser registry with ${PARSER_COUNT} parser headers")