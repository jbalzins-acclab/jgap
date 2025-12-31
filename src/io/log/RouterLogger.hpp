#ifndef JGAP_ROUTERLOGGER_HPP
#define JGAP_ROUTERLOGGER_HPP

#include <memory>
#include "io/log/Logger.hpp"
#include "io/log/LogConfig.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/log/FileLogger.hpp"

namespace jgap {
    class RouterLogger : public Logger {
    public:
        explicit RouterLogger(LogConfig cfg, std::string filePath = {});

        void log(LogLevel level, std::string_view msg) override;
        void logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) override;

        const LogConfig& config() const { return _cfg; }
        const std::string& filePath() const { return _fileLogger ? _fileLogger->filePath() : _filePathCached; }

    private:
        LogConfig _cfg{};
        std::shared_ptr<StdoutLogger> _stdoutLogger; // may be null if routing excludes stdout
        std::shared_ptr<FileLogger> _fileLogger;     // may be null if routing excludes files
        std::string _filePathCached;                 // in case file logger not constructed

        bool sendToStdout(LogLevel level) const;
        bool sendToFile(LogLevel level) const;
    };
}

#endif
