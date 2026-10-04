#include "ValidationUtils.hpp"
#include <algorithm>

namespace jgap::validation {

    double computeRmse(const std::vector<double>& actual, const std::vector<double>& ref) {
        if (actual.empty() || actual.size() != ref.size()) return 0.0;
        double s = 0.0;
        for (size_t i = 0; i < actual.size(); ++i) {
            double d = actual[i] - ref[i];
            s += d * d;
        }
        return std::sqrt(s / static_cast<double>(actual.size()));
    }

    double computeNrmsePct(const std::vector<double>& actual, const std::vector<double>& ref) {
        if (actual.empty() || actual.size() != ref.size()) return 0.0;
        double num = 0.0;
        double denom = 0.0;
        for (size_t i = 0; i < actual.size(); ++i) {
            double d = actual[i] - ref[i];
            num += d * d;
            denom += ref[i] * ref[i];
        }
        if (denom <= 1e-30) return 0.0;
        return (std::sqrt(num) / std::sqrt(denom)) * 100.0;
    }

    double computeRmsnePct(const std::vector<double>& actual, const std::vector<double>& ref, double min_ref) {
        if (actual.empty() || actual.size() != ref.size()) return 0.0;
        double s = 0.0;
        size_t count = 0;
        for (size_t i = 0; i < actual.size(); ++i) {
            double r = std::abs(ref[i]);
            if (r >= min_ref) {
                double rel = (actual[i] - ref[i]) / r;
                s += rel * rel;
                count++;
            }
        }
        if (count == 0) return 0.0;
        return std::sqrt(s / static_cast<double>(count)) * 100.0;
    }

    double computeMaxRelPct(const std::vector<double>& actual, const std::vector<double>& ref, double min_ref) {
        if (actual.empty() || actual.size() != ref.size()) return 0.0;
        double max_rel = 0.0;
        for (size_t i = 0; i < actual.size(); ++i) {
            double r = std::abs(ref[i]);
            if (r >= min_ref) {
                double rel = std::abs(actual[i] - ref[i]) / r;
                if (rel > max_rel) max_rel = rel;
            }
        }
        return max_rel * 100.0;
    }

    double computeSigRelPct(const std::vector<double>& actual, const std::vector<double>& ref, double sig_fraction) {
        if (actual.empty() || actual.size() != ref.size()) return 0.0;
        double max_abs = 0.0;
        for (double v : ref) {
            if (std::abs(v) > max_abs) max_abs = std::abs(v);
        }
        double threshold = sig_fraction * max_abs;
        return computeMaxRelPct(actual, ref, threshold);
    }

    double computeCosineSim(const std::vector<double>& actual, const std::vector<double>& ref) {
        if (actual.empty() || actual.size() != ref.size()) return 1.0;
        double dot = 0.0, na = 0.0, nr = 0.0;
        for (size_t i = 0; i < actual.size(); ++i) {
            dot += actual[i] * ref[i];
            na += actual[i] * actual[i];
            nr += ref[i] * ref[i];
        }
        double denom = std::sqrt(na) * std::sqrt(nr);
        if (denom <= 1e-30) return 1.0;
        return dot / denom;
    }

    fs::path findRepoRoot() {
        fs::path p = fs::current_path();
        while (p.has_parent_path() && p != p.parent_path()) {
            if (fs::exists(p / "CMakeLists.txt") && fs::exists(p / "src" / "jgap")) {
                return p;
            }
            p = p.parent_path();
        }
        return fs::current_path();
    }

    fs::path getValidationResourceDir() {
        fs::path root = findRepoRoot();
        fs::path p = root / "test" / "resources" / "validation";
        if (fs::exists(p)) return p;
        // Check relative to cwd
        if (fs::exists("test/resources/validation")) return fs::path("test/resources/validation");
        return p;
    }

    bool ValidationReporter::assertLessOrEqual(
        const std::string& test_name,
        const std::string& metric_name,
        double actual,
        double threshold,
        const std::string& unit
    ) {
        bool pass = (actual <= threshold);
        if (pass) {
            pass_count_++;
            std::cout << "[PASS] " << test_name << " :: " << metric_name
                      << " = " << actual << unit << " (<= " << threshold << unit << ")\n";
        } else {
            fail_count_++;
            std::string msg = "[FAIL] " + test_name + " :: " + metric_name +
                              " = " + std::to_string(actual) + unit + " EXCEEDED threshold " +
                              std::to_string(threshold) + unit;
            std::cerr << msg << "\n";
            failure_messages_.push_back(msg);
        }
        return pass;
    }

    bool ValidationReporter::assertGreaterOrEqual(
        const std::string& test_name,
        const std::string& metric_name,
        double actual,
        double threshold,
        const std::string& unit
    ) {
        bool pass = (actual >= threshold);
        if (pass) {
            pass_count_++;
            std::cout << "[PASS] " << test_name << " :: " << metric_name
                      << " = " << actual << unit << " (>= " << threshold << unit << ")\n";
        } else {
            fail_count_++;
            std::string msg = "[FAIL] " + test_name + " :: " + metric_name +
                              " = " + std::to_string(actual) + unit + " BELOW threshold " +
                              std::to_string(threshold) + unit;
            std::cerr << msg << "\n";
            failure_messages_.push_back(msg);
        }
        return pass;
    }

    void ValidationReporter::recordFailure(const std::string& msg) {
        fail_count_++;
        std::cerr << "[FAIL] " << msg << "\n";
        failure_messages_.push_back(msg);
    }

    int ValidationReporter::summarizeAndExit() const {
        std::cout << "\n=======================================================\n";
        std::cout << "Validation Summary: " << suite_name_ << "\n";
        std::cout << "  Passed: " << pass_count_ << "\n";
        std::cout << "  Failed: " << fail_count_ << "\n";
        if (fail_count_ > 0) {
            std::cout << "  Failures:\n";
            for (const auto& f : failure_messages_) {
                std::cout << "    - " << f << "\n";
            }
        }
        std::cout << "=======================================================\n";
        return (fail_count_ == 0) ? 0 : 1;
    }

} // namespace jgap::validation
