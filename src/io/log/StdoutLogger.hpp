#ifndef CONSOLELOGGER_HPP
#define CONSOLELOGGER_HPP

#include <iostream>
#include <memory>

#include "io/log/Logger.hpp"
#include "io/log/LogConfig.hpp"

namespace jgap {
    class StdoutLogger : public Logger {
    public:
        StdoutLogger(bool logDebug = true, MetadataVisibility meta = MetadataVisibility::FILES_ONLY);
        void log(LogLevel level, std::string_view msg) override;
        void logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) override;
    private:
        bool _logDebug;
        MetadataVisibility _metaVis;
    };
}

#endif
