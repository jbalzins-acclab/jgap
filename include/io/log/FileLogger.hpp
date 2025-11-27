#ifndef JGAP_FILELOGGER_HPP
#define JGAP_FILELOGGER_HPP

#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <mutex>

#include "io/log/Logger.hpp"
#include "io/log/LogConfig.hpp"

namespace jgap {
    class FileLogger : public Logger {
    public:
        // If filePath is empty, auto-generate in logs/jgap
        explicit FileLogger(std::string filePath = {}, MetadataVisibility meta = MetadataVisibility::FilesOnly);
        ~FileLogger() override;

        void log(LogLevel level, std::string_view msg) override;
        void logWithSrc(LogLevel level, std::string_view msg, const char* file, int line, const char* func) override;

        const std::string& filePath() const { return _filePath; }

    private:
        std::string _filePath;
        std::ofstream _out;
        std::mutex _mtx;
        MetadataVisibility _metaVis;

        static std::string autoFilePath();
        static std::string levelTag(LogLevel level);
    };
}

#endif
