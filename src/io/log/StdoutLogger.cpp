#include "io/log/StdoutLogger.hpp"

#include <iostream>
#include <chrono>
#include <format>


namespace jgap {

    StdoutLogger::StdoutLogger(const bool logDebug, MetadataVisibility meta)
        : _logDebug(logDebug), _metaVis(meta) {
        std::ios_base::sync_with_stdio(false);
        std::cin.tie(nullptr);
    }

    void StdoutLogger::log(LogLevel level, std::string_view msg) {
        std::string levelStr;

        std::ostream& out = level == LogLevel::Error ? std::cerr : std::cout;
        switch (level) {
            case LogLevel::Debug:
                if (!_logDebug) return;
                levelStr = "[DEBUG]"; break;
            case LogLevel::Info: levelStr = "[INFO ]"; break;
            case LogLevel::Warn:  levelStr = "[WARN ]"; break;
            case LogLevel::Error: levelStr = "[ERROR]"; break;
        }

        auto timestamp = std::chrono::system_clock::now();
        out << format("{:%Y-%m-%d %H:%M:%S} {} {}", timestamp, levelStr, msg) << std::endl;
    }

    void StdoutLogger::logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) {
        if (_metaVis == MetadataVisibility::BOTH) {
            const std::string withSrc = format("{} ({}:{}:{})", msg, file ? file : "?", line, func ? func : "?");
            log(level, withSrc);
        } else {
            log(level, msg);
        }
    }
}
