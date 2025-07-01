#include "io/log/StdoutLogger.hpp"

#include <iostream>
#include <chrono>
#include <format>

using namespace std;

namespace jgap {

    StdoutLogger::StdoutLogger(const bool debug) : _debug(debug) {
        ios_base::sync_with_stdio(false);
        cin.tie(nullptr);
    }

    void StdoutLogger::print(string_view msg, string_view type, bool isError) {
        ostream& out = isError ? cerr : cout;

        auto timestamp = chrono::system_clock::now();

        out << format("{:%Y-%m-%d %H:%M:%S} [{}] {}", timestamp, type, msg) << endl;
    }
}
