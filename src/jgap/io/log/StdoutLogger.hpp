#ifndef CONSOLELOGGER_HPP
#define CONSOLELOGGER_HPP

#include <iostream>
#include <memory>

#include "jgap/io/log/Logger.hpp"
#include "jgap/io/log/LogConfig.hpp"

namespace jgap {
    class StdoutLogger : public Logger {
    public:
        StdoutLogger(bool log_debug = true, MetadataVisibility meta = MetadataVisibility::FilesOnly);
        void log(LogLevel level, std::string_view msg) override;
        void logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) override;
    private:
        bool log_debug;
        MetadataVisibility meta_vis;
    };
}

#endif
