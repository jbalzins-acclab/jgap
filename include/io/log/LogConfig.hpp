#ifndef JGAP_LOGCONFIG_HPP
#define JGAP_LOGCONFIG_HPP

namespace jgap {
    enum class OutputRouting {
        None,                      // Do not print anywhere
        FilesOnly,                 // Print to log files only
        BothStdoutAndFiles,        // Print to both stdout and files for all levels
        MixedNonDebugToStdout      // Default: non-debug to stdout; all levels to files
    };

    enum class MetadataVisibility {
        None,        // Do not include file:line/func anywhere
        FilesOnly,   // Default: include metadata in files only
        Both         // Include metadata in both stdout and files
    };

    struct LogConfig {
        OutputRouting routing = OutputRouting::MixedNonDebugToStdout;
        MetadataVisibility metadata = MetadataVisibility::FilesOnly;
        bool stdoutLogDebug = false; // if routing sends to stdout, should DEBUG appear there?
    };
}

#endif
