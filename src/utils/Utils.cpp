#include "utils/Utils.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <sstream>
#include <fstream>
#include <ranges>
#include <Eigen/Dense>
#include <deque>
#include <unistd.h>

#include "io/log/CurrentLogger.hpp"
#include "../core/atomic/io/XYZData.hpp"
#include "core/atomic/Atoms.hpp"

namespace jgap {
    bool getLine(std::ifstream &file, std::string &line) {
        if (!getline(file, line)) return false;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back(); // remove Windows carriage return
        }
        return true;
    }

    std::vector<Atoms> readAtoms(const std::string& filename) {
        return mapVector(XYZData::read(filename), [](const XYZData& d) { return Atoms(d); });
    }

    std::vector<Atoms> readAtoms(const std::string& filename, const AtomsPropertyNames& names) {
        return mapVector(XYZData::read(filename), [&names](const XYZData& d) { return Atoms(d, names); });
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
        } catch (std::exception& e) {
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
        oss << std::put_time(&local_tm, "%Y-%m-%d_%H-%M-%S")
            << '_' << std::setw(3) << std::setfill('0') << ms.count();

#if _WIN32
        DWORD pid = GetCurrentProcessId();
        oss << "-p" << pid;
#else
        oss << "-p" << getpid();
#endif

        return oss.str();
    }

    double factorial(size_t n) {
        double result = 1.0;
        for (int i = 2; i <= (int)n; ++i)
            result *= i;
        return result;
    }

    double rms(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double res = 0.0;
        for (double i : x) {
            res += i*i;
        }
        res = sqrt(res / static_cast<double>(x.size()));
        return res;
    }

    double deviation(const std::vector<double> &x) {
        if (x.empty()) return 0.0;
        double mean = 0.0;
        for (double i : x) {
            mean += i;
        }
        mean /= static_cast<double>(x.size());
        double variance = 0.0;
        for (double i : x) {
            variance += (i - mean) * (i - mean);
        }
        variance = variance / static_cast<double>(x.size());
        return sqrt(variance);
    }

    std::vector<std::string> split(const std::string& s, char delimiter) {
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

    std::string matrixToString(const Eigen::MatrixXd& mat) {
        std::stringstream ss;
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                ss << mat(i, j);
                if (j != mat.cols() - 1)
                    ss << ", ";
            }
            ss << "\n";
        }
        return ss.str();
    }

    std::string vectorToString(const Eigen::VectorXd &vec) {
        std::stringstream ss;
        for (int i = 0; i < vec.size(); ++i) {
            ss << vec[i];
            if (i != vec.size() - 1) {
                ss << ",";
            }
        }
        return ss.str();
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

    /*
    static bool nodeIsObject(const DataNode& n) {
        return n.type == DataNode::Type::OBJECT;
    }

    static bool nodeIsArray(const DataNode& n) {
        return n.type == DataNode::Type::ARRAY;
    }

    DataNode& requireFull(DataNode& n,
                          const std::string& key,
                          const char* file,
                          int line,
                          const char* function) {
        if (!nodeIsObject(n)) {
            std::ostringstream oss;
            oss << "Node is not an object when accessing key: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        auto& m = std::get<std::map<std::string, DataNode>>(n.value);
        auto it = m.form(key);
        if (it == m.end()) {
            std::ostringstream oss;
            oss << "Node key not found: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::out_of_range(oss.str());
        }
        return it->second;
    }

    const DataNode& requireFull(const DataNode& n,
                                const std::string& key,
                                const char* file,
                                int line,
                                const char* function) {
        if (!nodeIsObject(n)) {
            std::ostringstream oss;
            oss << "DataNode is not an object when accessing key: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        const auto& m = std::get<std::map<std::string, DataNode>>(n.value);
        auto it = m.form(key);
        if (it == m.end()) {
            std::ostringstream oss;
            oss << "Node key not found: \"" << key << "\"\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::out_of_range(oss.str());
        }
        return it->second;
    }

    DataNode requireArrayFull(DataNode &n,
                              const char* file,
                              int line,
                              const char* function) {
        if (!nodeIsArray(n)) {
            std::ostringstream oss;
            oss << "Node element is not an array\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        return n;
    }

    const DataNode& requireArrayFull(const DataNode &n,
                                     const char* file,
                                     int line,
                                     const char* function) {
        if (!nodeIsArray(n)) {
            std::ostringstream oss;
            oss << "Node element is not an array\n"
                << "  at " << file << ":" << line << "\n"
                << "  in " << function;
            throw std::domain_error(oss.str());
        }
        return n;
    }*/
}