#ifndef CONSOLELOGGER_HPP
#define CONSOLELOGGER_HPP

#include <iostream>
#include <memory>

#include "io/log/Logger.hpp"

namespace jgap {
    class StdoutLogger : public Logger {
    public:
        explicit StdoutLogger(bool logDebug = true);
        void log(LogLevel level, std::string_view msg) override;
    private:
        bool _logDebug;
    };
}

#endif
