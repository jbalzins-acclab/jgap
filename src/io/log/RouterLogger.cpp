#include "io/log/RouterLogger.hpp"


namespace jgap {

    RouterLogger::RouterLogger(LogConfig cfg, std::string file_path)
        : cfg(cfg) {
        // Setup stdout logger depending on routing
        if (cfg.routing == OutputRouting::BothStdoutAndFiles || cfg.routing == OutputRouting::MixedNonDebugStdout) {
            // For stdout, metadata visibility is FilesOnly -> only include metadata in files by default
            const auto stdoutMeta = (cfg.metadata == MetadataVisibility::Both) ? MetadataVisibility::Both : MetadataVisibility::None;
            stdout_logger = std::make_shared<StdoutLogger>(cfg.stdout_log_debug, stdoutMeta);
        }

        // Setup file logger depending on routing
        if (cfg.routing == OutputRouting::FilesOnly || cfg.routing == OutputRouting::BothStdoutAndFiles || cfg.routing == OutputRouting::MixedNonDebugStdout) {
            file_logger = std::make_shared<FileLogger>(std::move(file_path), cfg.metadata == MetadataVisibility::None ? MetadataVisibility::None : MetadataVisibility::FilesOnly);
        }
    }

    bool RouterLogger::sendToStdout(LogLevel level) const {
        if (!stdout_logger) return false;
        switch (cfg.routing) {
            case OutputRouting::None: return false;
            case OutputRouting::FilesOnly: return false;
            case OutputRouting::BothStdoutAndFiles: return true;
            case OutputRouting::MixedNonDebugStdout:
                return level != LogLevel::Debug || cfg.stdout_log_debug;
        }
        return false;
    }

    bool RouterLogger::sendToFile(LogLevel /*level*/) const {
        return file_logger != nullptr; // if constructed, always send
    }

    void RouterLogger::log(LogLevel level, std::string_view msg) {
        if (cfg.routing == OutputRouting::None) return;
        if (sendToStdout(level)) {
            stdout_logger->log(level, msg);
        }
        if (sendToFile(level)) {
            file_logger->log(level, msg);
        }
    }

    void RouterLogger::logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) {
        if (cfg.routing == OutputRouting::None) return;
        if (sendToStdout(level)) {
            stdout_logger->logWithSrc(level, msg, file, line, func);
        }
        if (sendToFile(level)) {
            file_logger->logWithSrc(level, msg, file, line, func);
        }
    }
}
