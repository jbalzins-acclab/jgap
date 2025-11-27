#include "io/log/RouterLogger.hpp"


namespace jgap {

    RouterLogger::RouterLogger(LogConfig cfg, std::string filePath)
        : _cfg(cfg) {
        // Setup stdout logger depending on routing
        if (_cfg.routing == OutputRouting::BothStdoutAndFiles || _cfg.routing == OutputRouting::MixedNonDebugToStdout) {
            // For stdout, metadata visibility is FilesOnly -> only include metadata in files by default
            const auto stdoutMeta = (_cfg.metadata == MetadataVisibility::Both) ? MetadataVisibility::Both : MetadataVisibility::None;
            _stdoutLogger = make_shared<StdoutLogger>(_cfg.stdoutLogDebug, stdoutMeta);
        }

        // Setup file logger depending on routing
        if (_cfg.routing == OutputRouting::FilesOnly || _cfg.routing == OutputRouting::BothStdoutAndFiles || _cfg.routing == OutputRouting::MixedNonDebugToStdout) {
            _fileLogger = make_shared<FileLogger>(std::move(filePath), _cfg.metadata == MetadataVisibility::None ? MetadataVisibility::None : MetadataVisibility::FilesOnly);
        }
    }

    bool RouterLogger::sendToStdout(LogLevel level) const {
        if (!_stdoutLogger) return false;
        switch (_cfg.routing) {
            case OutputRouting::None: return false;
            case OutputRouting::FilesOnly: return false;
            case OutputRouting::BothStdoutAndFiles: return true;
            case OutputRouting::MixedNonDebugToStdout:
                return level != LogLevel::Debug || _cfg.stdoutLogDebug;
        }
        return false;
    }

    bool RouterLogger::sendToFile(LogLevel /*level*/) const {
        return _fileLogger != nullptr; // if constructed, always send
    }

    void RouterLogger::log(LogLevel level, std::string_view msg) {
        if (_cfg.routing == OutputRouting::None) return;
        if (sendToStdout(level)) {
            _stdoutLogger->log(level, msg);
        }
        if (sendToFile(level)) {
            _fileLogger->log(level, msg);
        }
    }

    void RouterLogger::logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) {
        if (_cfg.routing == OutputRouting::None) return;
        if (sendToStdout(level)) {
            _stdoutLogger->logWithSrc(level, msg, file, line, func);
        }
        if (sendToFile(level)) {
            _fileLogger->logWithSrc(level, msg, file, line, func);
        }
    }
}
