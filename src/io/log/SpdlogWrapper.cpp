#include "io/log/SpdlogWrapper.hpp"

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace jgap {

    SpdlogWrapper::SpdlogWrapper(const string_view logFilePrefix, const string_view logDirectory, const bool debug) {

        string logFile = string(logDirectory);
        if (!logDirectory.empty()) {
            filesystem::create_directories(logDirectory);
            if (logFile.back() != '/') {
                logFile += "/";
            }
        }

        // Get the current time
        auto now = chrono::system_clock::now();
        auto time_t_now = chrono::system_clock::to_time_t(now);

        // Convert time to string in the format YYYY-MM-DD_HH-MM-SS
        tm tm_now = *localtime(&time_t_now);
        ostringstream oss;
        oss << put_time(&tm_now, "%Y-%m-%d_%H-%M-%S") << ".log";

        logFile += string(logFilePrefix) + "-" + oss.str();

        try {
            // Create a basic file sink for logging
            auto fileSink = make_shared<spdlog::sinks::basic_file_sink_mt>(logFile, true);

            // Create a console sink (optional, for logging to console)
            auto consoleSink = make_shared<spdlog::sinks::stdout_color_sink_mt>();

            // Create a logger and add both sinks
            _logger = make_shared<spdlog::logger>("log", spdlog::sinks_init_list{fileSink, consoleSink});

            cout << "Logging to " << logFile << " at level=";
            // Set the log level
            if (debug) {
                _logger->set_level(spdlog::level::debug);
                cout << "debug" << endl;
            } else {
                _logger->set_level(spdlog::level::info);
                cout << "info" << endl;
            }
        } catch (const spdlog::spdlog_ex& ex) {
            cerr << "Log initialization failed: " << ex.what() << endl;
            throw;
        }
    }

    SpdlogWrapper::~SpdlogWrapper() {
        spdlog::shutdown();
    }

    void SpdlogWrapper::info(const string_view msg) {
        _logger->info(msg);
    }

    void SpdlogWrapper::debug(const string_view msg) {
        _logger->debug(msg);
    }

    void SpdlogWrapper::warn(const string_view msg) {
        _logger->warn(msg);
    }

    void SpdlogWrapper::error(const string_view msg) {
        _logger->error(msg);
    }
}