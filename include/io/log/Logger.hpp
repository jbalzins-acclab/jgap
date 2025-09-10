#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string_view>
#include <string>
#include <format>   // C++20

using namespace std;

namespace jgap {
    enum class LogLevel { Debug, Info, Warn, Error };

    class Logger {
    public:
        virtual ~Logger() = default;

        // Core virtual
        virtual void log(LogLevel level, string_view msg) = 0;

        // Convenience overloads with formatting
        template <typename... Args>
        void debug(const string_view fmt, Args&&... args) {
            log(LogLevel::Debug, vformat(fmt, make_format_args(args...)));
        }

        template <typename... Args>
        void info(const string_view fmt, Args&&... args) {
            log(LogLevel::Info, vformat(fmt, make_format_args(args...)));
        }

        template <typename... Args>
        void warn(const string_view fmt, Args&&... args) {
            log(LogLevel::Warn, vformat(fmt, make_format_args(args...)));
        }

        template <typename... Args>
        void error(const string_view fmt, Args&&... args) {
            log(LogLevel::Error, vformat(fmt, make_format_args(args...)));
        }

        template <typename... Args>
        [[noreturn]] void logAndThrow(const string_view fmt, Args&&... args) {
            const string msg = vformat(fmt, make_format_args(args...));
            log(LogLevel::Error, msg);
            throw runtime_error(msg);
        }
    };
}

#endif
