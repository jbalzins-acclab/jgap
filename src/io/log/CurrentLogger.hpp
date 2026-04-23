#ifndef JGAP_CURRENTLOGGER_HPP
#define JGAP_CURRENTLOGGER_HPP

#include "Logger.hpp"
#include "RouterLogger.hpp"
#include "LogConfig.hpp"
#include <memory>
#include <string>

namespace jgap {
    class CurrentLogger {
    public:
        static std::shared_ptr<Logger> get() {
            return _logger;
        }
        static void setLogger(const std::shared_ptr<Logger> &logger) {
            _logger = logger;
        }
        static void initDefault(const LogConfig& cfg = LogConfig{}, const std::string& filePath = {}) {
            _logger = std::make_shared<RouterLogger>(cfg, filePath);
        }
    private:
        inline static std::shared_ptr<Logger> _logger = std::make_shared<RouterLogger>(LogConfig{});
    };
}

#define JGAP_LOG CurrentLogger::get()

// Simple macros that automatically pass source location
#define JGAP_LOG_DEBUG(fmt, ...) jgap::CurrentLogger::get()->debugSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_INFO(fmt, ...)  jgap::CurrentLogger::get()->infoSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_WARN(fmt, ...)  jgap::CurrentLogger::get()->warnSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_ERROR(fmt, ...) jgap::CurrentLogger::get()->errorSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_AND_THROW(fmt, ...) jgap::CurrentLogger::get()->logAndThrowSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)

#endif
