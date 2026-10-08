#ifndef JGAP_LOGCONFIG_HPP
#define JGAP_LOGCONFIG_HPP

namespace jgap {
    enum class OutputRouting {
        None,                      // Do not print anywhere
        StdoutOnly,                // Default: print to stdout only (no file logging)
        FilesOnly,                 // Print to log files only
        BothStdoutAndFiles,        // Print to both stdout and files for all levels
        MixedNonDebugStdout        // Non-debug to stdout; all levels to files
    };

    enum class MetadataVisibility {
        None,        // Do not include file:line/func anywhere
        FilesOnly,   // Include metadata in files only
        Both         // Include metadata in both stdout and files
    };

    struct LogConfig {
        OutputRouting routing = OutputRouting::StdoutOnly;
        MetadataVisibility metadata = MetadataVisibility::None;
        bool stdout_log_debug = false; // if routing sends to stdout, should DEBUG appear there?
    };
}

#endif
