#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string_view>
#include <string>

using namespace std;

namespace jgap {
    class Logger {
    public:
        Logger() = default;
        virtual ~Logger() = default;

        virtual void info(string_view msg) = 0;
        virtual void debug(string_view msg) = 0;
        virtual void warn(string_view msg) = 0;
        virtual void error(string_view msg) = 0;
        virtual void error(const string_view msg, const bool throwRuntimeErr) {
            error(msg);
            if (throwRuntimeErr) {
                throw runtime_error(string(msg));
            }
        }

        inline static Logger* logger = nullptr;
    };
}

#endif
