#include "io/log/FileLogger.hpp"

#include <chrono>
#include <format>
#include <filesystem>
#include <cstdio>

#if defined(_WIN32)
#  include <process.h>
static int jgap_getpid() { return _getpid(); }
#else
#  include <unistd.h>
static int jgap_getpid() { return static_cast<int>(getpid()); }
#endif


namespace fs = std::filesystem;

namespace jgap {

    static std::string timestamp_for_filename() {
        // Use local time in human-readable safe format: YYYYMMDD-HHMMSS
        auto now = chrono::system_clock::now();
        return format("{:%Y%m%d-%H%M%S}", now);
    }

    std::string FileLogger::autoFilePath() {
        const std::string baseDir = "logs/jgap";
        fs::create_directories(baseDir);
        const std::string name = format("jgap-{}-{}.log", timestamp_for_filename(), jgap_getpid());
        return (fs::path(baseDir) / name).std::string();
    }

    std::string FileLogger::levelTag(LogLevel level) {
        switch (level) {
            case LogLevel::Debug: return "[DEBUG]";
            case LogLevel::Info:  return "[INFO ]";
            case LogLevel::Warn:  return "[WARN ]";
            case LogLevel::Error: return "[ERROR]";
        }
        return "[INFO ]";
    }

    FileLogger::FileLogger(std::string filePath, MetadataVisibility meta)
        : _filePath(std::move(filePath)), _metaVis(meta) {
        if (_filePath.empty()) _filePath = autoFilePath();
        _out.open(_filePath, ios::out | ios::app);
    }

    FileLogger::~FileLogger() {
        if (_out.is_open()) _out.flush();
    }

    void FileLogger::log(LogLevel level, string_view msg) {
        lock_guard lock(_mtx);
        if (!_out.is_open()) return;
        const auto ts = chrono::system_clock::now();
        _out << format("{:%Y-%m-%d %H:%M:%S} {} {}\n", ts, levelTag(level), msg);
        _out.flush();
    }

    void FileLogger::logWithSrc(LogLevel level, string_view msg, const char* file, int line, const char* func) {
        if (_metaVis == MetadataVisibility::None) {
            log(level, msg);
            return;
        }
        const std::string withSrc = format("{} ({}:{}:{})", msg, file ? file : "?", line, func ? func : "?");
        log(level, withSrc);
    }
}
