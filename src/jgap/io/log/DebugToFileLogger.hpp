#ifndef JGAP_DEBUGTOFILELOGGER_HPP
#define JGAP_DEBUGTOFILELOGGER_HPP
#include "StdoutLogger.hpp"

namespace jgap {
    class DebugToFileLogger : public StdoutLogger {
    public:
        DebugToFileLogger(std::optional<std::string> debug_file_name);
        void log(LogLevel level, std::string_view msg) override;
    };
}

#endif