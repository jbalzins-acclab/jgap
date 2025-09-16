#include "io/log/StdoutLogger.hpp"

#include <iostream>
#include <chrono>
#include <format>

using namespace std;

namespace jgap {

    StdoutLogger::StdoutLogger(const bool logDebug) : _logDebug(logDebug) {
        ios_base::sync_with_stdio(false);
        cin.tie(nullptr);
    }

    void StdoutLogger::log(LogLevel level, string_view msg) {
        string levelStr;

        ostream& out = level == LogLevel::Error ? cerr : cout;
        switch (level) {
            case LogLevel::Debug:
                if (!_logDebug) return;
                levelStr = "[DEBUG]"; break;
            case LogLevel::Info: levelStr = "[INFO ]"; break;
            case LogLevel::Warn:  levelStr = "[WARN ]"; break;
            case LogLevel::Error: levelStr = "[ERROR]"; break;
        }

        auto timestamp = chrono::system_clock::now();
        out << format("{:%Y-%m-%d %H:%M:%S} {} {}", timestamp, levelStr, msg) << endl;
    }
}
