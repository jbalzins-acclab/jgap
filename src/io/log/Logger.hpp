#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string_view>
#include <string>
#include <format>   // C++20
#include <stdexcept>


namespace jgap {
    enum class LogLevel { Debug, Info, Warn, Error };

    class Logger {
    public:
        virtual ~Logger() = default;

        // Core virtual
        virtual void log(LogLevel level, std::string_view msg) = 0;

        // Extended API with source location metadata (file/line/func)
        // Default implementation falls back to legacy log() if derived class doesn't override it
        virtual void logWithSrc(LogLevel level,
                                std::string_view msg,
                                const char* file,
                                int line,
                                const char* func) {
            (void)file; (void)line; (void)func; // suppress unused warnings if not used
            log(level, msg);
        }

        // Convenience overloads with formatting
        template <typename... Args>
        void debug(const std::string_view fmt, Args&&... args) {
            log(LogLevel::Debug, std::vformat(fmt, std::make_format_args(args...)));
        }

        template <typename... Args>
        void info(const std::string_view fmt, Args&&... args) {
            log(LogLevel::Info, std::vformat(fmt, std::make_format_args(args...)));
        }

        template <typename... Args>
        void warn(const std::string_view fmt, Args&&... args) {
            log(LogLevel::Warn, std::vformat(fmt, std::make_format_args(args...)));
        }

        template <typename... Args>
        void error(const std::string_view fmt, Args&&... args) {
            log(LogLevel::Error, std::vformat(fmt, std::make_format_args(args...)));
        }

        // Convenience helpers that also pass source location metadata
        template <typename... Args>
        void debugSrc(const char* file, int line, const char* func, const std::string_view fmt, Args&&... args) {
            logWithSrc(LogLevel::Debug, std::vformat(fmt, std::make_format_args(args...)), file, line, func);
        }

        template <typename... Args>
        void infoSrc(const char* file, int line, const char* func, const std::string_view fmt, Args&&... args) {
            logWithSrc(LogLevel::Info, std::vformat(fmt, std::make_format_args(args...)), file, line, func);
        }

        template <typename... Args>
        void warnSrc(const char* file, int line, const char* func, const std::string_view fmt, Args&&... args) {
            logWithSrc(LogLevel::Warn, std::vformat(fmt, std::make_format_args(args...)), file, line, func);
        }

        template <typename... Args>
        void errorSrc(const char* file, int line, const char* func, const std::string_view fmt, Args&&... args) {
            logWithSrc(LogLevel::Error, std::vformat(fmt, std::make_format_args(args...)), file, line, func);
        }

        template <typename... Args>
        [[noreturn]] void logAndThrow(const std::string_view fmt, Args&&... args) {
            const std::string msg = std::vformat(fmt, std::make_format_args(args...));
            log(LogLevel::Error, msg);
            throw std::runtime_error(msg);
        }

        // Same as logAndThrow but also passes source location metadata
        template <typename... Args>
        [[noreturn]] void logAndThrowSrc(const char* file, int line, const char* func, const std::string_view fmt, Args&&... args) {
            const std::string msg = std::vformat(fmt, std::make_format_args(args...));
            logWithSrc(LogLevel::Error, msg, file, line, func);
            throw std::runtime_error(msg);
        }
    };
}

#endif
