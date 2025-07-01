#ifndef CONSOLELOGGER_HPP
#define CONSOLELOGGER_HPP

#include <memory>

#include "io/log/Logger.hpp"

namespace jgap {
    class StdoutLogger : public Logger {
    public:
        explicit StdoutLogger(bool debug = true);

        void info(const string_view msg) override { print(msg, "INFO", false); }
        void debug(const string_view msg) override { if (_debug) print(msg, "DEBUG", false); }
        void warn(const string_view msg) override { print(msg, "WARN", true); }
        void error(const string_view msg) override { print(msg, "ERROR", true); }

        static void init(const bool debug = true) {
            logger = new StdoutLogger(debug);
        }

        static void initIfNotInitialized() {
            if (logger == nullptr) logger = new StdoutLogger();
        }

    private:
        bool _debug;
        void print(string_view msg, string_view type, bool isError);
    };
}

#endif
