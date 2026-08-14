#ifndef JGAP_CURRENTLOGGER_HPP
#define JGAP_CURRENTLOGGER_HPP

#include "jgap/core/io/log/Logger.hpp"
#include "jgap/core/io/log/RouterLogger.hpp"
#include "LogConfig.hpp"
#include <memory>
#include <string>
#include <exception>

namespace jgap {
    class CurrentLogger {
    public:
        static std::shared_ptr<Logger> get() {
            return logger;
        }
        static void setLogger(const std::shared_ptr<Logger> &new_logger) {
            logger = new_logger;
        }
        static void initDefault(const LogConfig& cfg = LogConfig{}, const std::string& file_path = {}) {
            logger = std::make_shared<RouterLogger>(cfg, file_path);
        }
    private:
        inline static std::shared_ptr<Logger> logger = std::make_shared<RouterLogger>(LogConfig{});
    };
}

#define JGAP_LOG jgap::CurrentLogger::get()

// Simple macros that automatically pass source location
#define JGAP_LOG_DEBUG(fmt, ...) jgap::CurrentLogger::get()->debugSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_INFO(fmt, ...)  jgap::CurrentLogger::get()->infoSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_WARN(fmt, ...)  jgap::CurrentLogger::get()->warnSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_ERROR(fmt, ...) jgap::CurrentLogger::get()->errorSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
#define JGAP_LOG_AND_THROW(fmt, ...) \
    {jgap::CurrentLogger::get()->logAndThrowSrc(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__); \
    throw std::runtime_error("Log and throw did't throw - CHECK YOUR LOGGER!");}

#endif
