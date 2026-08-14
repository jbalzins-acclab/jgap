#include "jgap/utils/Utils.hpp"

#include <array>
#include <chrono>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include "../core/io/xyz/XYZData.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "../core/io/log/CurrentLogger.hpp"

namespace jgap {
    bool getLine(std::istream &file, std::string &line) {
        if (!getline(file, line)) return false;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back(); // remove Windows carriage return
        }
        return true;
    }

    std::map<std::string, std::string> parseHeaderLine(const std::string &line) {
        std::map<std::string, std::string> header;

        try {
            size_t pos = 0;
            while (pos < line.size()) {
                if (isspace(line[pos])) {
                    pos++;
                    continue;
                }

                std::string property;
                while (line[pos] != '=' && !isspace(line[pos])) {
                    property += line[pos];
                    pos++;

                    if (pos >= line.size()) {
                        throw std::runtime_error("'=' not found after " + property);
                    }
                }

                if (header.contains(property)) {
                    throw std::runtime_error("Duplicate property: " + property);
                }

                std::string value;

                if (!isspace(line[pos])) {
                    pos++;

                    if (line[pos] == '"') {
                        pos++;
                        while (pos < line.size() && line[pos] != '"') {
                            value += line[pos];
                            pos++;
                        }
                        pos++;
                    } else {
                        while (pos < line.size() && !isspace(line[pos])) {
                            value += line[pos];
                            pos++;
                        }
                    }
                }

                header[property] = value;
            }
        } catch (std::exception &e) {
            JGAP_LOG_AND_THROW("Formatting error {} in : {}", e.what(), line);
        } catch (...) {
            JGAP_LOG_AND_THROW("Formatting error in: {}", line);
        }

        return header;
    }

    std::string uniqueStamp() {
        using namespace std::chrono;

        auto now = system_clock::now();
        auto tt = system_clock::to_time_t(now);

        std::tm local_tm{};
#ifdef _WIN32
        localtime_s(&local_tm, &tt);
#else
        localtime_r(&tt, &local_tm);
#endif

        auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

        std::ostringstream oss;
        // Human-readable and filename-friendly: YYYY-MM-DD_HH-MM-SS_mmm-pPID
        oss << std::put_time(&local_tm, "%Y-%m-%d_%H-%M-%S") << '_' << std::setw(3) << std::setfill('0') << ms.count();

#if _WIN32
        DWORD pid = GetCurrentProcessId();
        oss << "-p" << pid;
#else
        oss << "-p" << getpid();
#endif

        return oss.str();
    }

    double rms(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double res = 0.0;
        for (double i: x) {
            res += i * i;
        }
        res = sqrt(res / static_cast<double>(x.size()));
        return res;
    }

    double deviation(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double mean = 0.0;
        for (double i: x) {
            mean += i;
        }
        mean /= static_cast<double>(x.size());
        double variance = 0.0;
        for (double i: x) {
            variance += (i - mean) * (i - mean);
        }
        variance = variance / static_cast<double>(x.size());
        return sqrt(variance);
    }

    std::vector<std::string> split(const std::string &s, char delimiter) {
        std::vector<std::string> result;
        std::stringstream ss(s);
        std::string token;

        while (getline(ss, token, delimiter)) {
            result.push_back(token);
        }

        return result;
    }

    std::string join(const std::vector<std::string> &strs, char delimiter) {
        if (strs.empty()) return "";
        std::string res = strs[0];
        for (size_t i = 1; i < strs.size(); i++) {
            res += std::string{delimiter} + strs[i];
        }
        return res;
    }

    std::string withoutExtension(const std::string &s) {
        if (s.find('.') == std::string::npos) {
            return s;
        }

        auto after_split = split(s, '.');
        return join(std::vector(after_split.begin(), after_split.end() - 1), '.');
    }

    std::string vectorToString(const std::vector<double> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    std::string vectorToString(const std::vector<size_t> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }

    std::string vectorToString(const std::vector<std::string> &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
    }
}
