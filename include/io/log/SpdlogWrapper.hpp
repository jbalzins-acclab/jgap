#ifndef SPDLOGWRAPPER_HPP
#define SPDLOGWRAPPER_HPP

#include "io/log/Logger.hpp"

#include <spdlog/spdlog.h>
#include <string_view>

using namespace std;

namespace jgap {
    class SpdlogWrapper : public Logger {
    public:
        explicit SpdlogWrapper(
            string_view logFilePrefix = "log-jgap-",
            string_view logDirectory = "logs/jgap",
            bool debug = false
            );
        ~SpdlogWrapper() override;

        void info(string_view msg) override;
        void debug(string_view msg) override;
        void warn(string_view msg) override;
        void error(string_view msg) override;

    private:
        shared_ptr<spdlog::logger> _logger;
    };
}

#endif