#ifndef JGAP_ROUTERLOGGER_HPP
#define JGAP_ROUTERLOGGER_HPP

#include <memory>
#include "jgap/io/log/Logger.hpp"
#include "jgap/io/log/LogConfig.hpp"
#include "jgap/io/log/StdoutLogger.hpp"
#include "jgap/io/log/FileLogger.hpp"

namespace jgap {
    class RouterLogger : public Logger {
    public:
        explicit RouterLogger(LogConfig cfg, std::string file_path = {});

        void log(LogLevel level, std::string_view msg) override;
        void logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) override;

        const LogConfig& config() const { return cfg; }
        const std::string& filePath() const { return file_logger ? file_logger->filePath() : file_path_cached; }

    private:
        LogConfig cfg{};
        std::shared_ptr<StdoutLogger> stdout_logger; // may be null if routing excludes stdout
        std::shared_ptr<FileLogger> file_logger;     // may be null if routing excludes files
        std::string file_path_cached;                 // in case file logger not constructed

        bool sendToStdout(LogLevel level) const;
        bool sendToFile(LogLevel level) const;
    };
}

#endif
