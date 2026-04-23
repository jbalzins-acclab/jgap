#ifndef JGAP_LOGCONFIG_HPP
#define JGAP_LOGCONFIG_HPP

namespace jgap {
    enum class OutputRouting {
        NONE,                      // Do not print anywhere
        FILES_ONLY,                 // Print to log files only
        BOTH_STDOUT_AND_FILES,        // Print to both stdout and files for all levels
        MIXED_NON_DEBUG_STDOUT      // Default: non-debug to stdout; all levels to files
    };

    enum class MetadataVisibility {
        NONE,        // Do not include file:line/func anywhere
        FILES_ONLY,   // Default: include metadata in files only
        BOTH         // Include metadata in both stdout and files
    };

    struct LogConfig {
        OutputRouting routing = OutputRouting::MIXED_NON_DEBUG_STDOUT;
        MetadataVisibility metadata = MetadataVisibility::FILES_ONLY;
        bool stdoutLogDebug = false; // if routing sends to stdout, should DEBUG appear there?
    };
}

#endif
