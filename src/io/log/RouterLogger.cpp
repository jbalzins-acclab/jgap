#include "io/log/RouterLogger.hpp"


namespace jgap {

    RouterLogger::RouterLogger(LogConfig cfg, std::string filePath)
        : _cfg(cfg) {
        // Setup stdout logger depending on routing
        if (_cfg.routing == OutputRouting::BOTH_STDOUT_AND_FILES || _cfg.routing == OutputRouting::MIXED_NON_DEBUG_STDOUT) {
            // For stdout, metadata visibility is FilesOnly -> only include metadata in files by default
            const auto stdoutMeta = (_cfg.metadata == MetadataVisibility::BOTH) ? MetadataVisibility::BOTH : MetadataVisibility::NONE;
            _stdoutLogger = std::make_shared<StdoutLogger>(_cfg.stdoutLogDebug, stdoutMeta);
        }

        // Setup file logger depending on routing
        if (_cfg.routing == OutputRouting::FILES_ONLY || _cfg.routing == OutputRouting::BOTH_STDOUT_AND_FILES || _cfg.routing == OutputRouting::MIXED_NON_DEBUG_STDOUT) {
            _fileLogger = std::make_shared<FileLogger>(std::move(filePath), _cfg.metadata == MetadataVisibility::NONE ? MetadataVisibility::NONE : MetadataVisibility::FILES_ONLY);
        }
    }

    bool RouterLogger::sendToStdout(LogLevel level) const {
        if (!_stdoutLogger) return false;
        switch (_cfg.routing) {
            case OutputRouting::NONE: return false;
            case OutputRouting::FILES_ONLY: return false;
            case OutputRouting::BOTH_STDOUT_AND_FILES: return true;
            case OutputRouting::MIXED_NON_DEBUG_STDOUT:
                return level != LogLevel::Debug || _cfg.stdoutLogDebug;
        }
        return false;
    }

    bool RouterLogger::sendToFile(LogLevel /*level*/) const {
        return _fileLogger != nullptr; // if constructed, always send
    }

    void RouterLogger::log(LogLevel level, std::string_view msg) {
        if (_cfg.routing == OutputRouting::NONE) return;
        if (sendToStdout(level)) {
            _stdoutLogger->log(level, msg);
        }
        if (sendToFile(level)) {
            _fileLogger->log(level, msg);
        }
    }

    void RouterLogger::logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) {
        if (_cfg.routing == OutputRouting::NONE) return;
        if (sendToStdout(level)) {
            _stdoutLogger->logWithSrc(level, msg, file, line, func);
        }
        if (sendToFile(level)) {
            _fileLogger->logWithSrc(level, msg, file, line, func);
        }
    }
}
