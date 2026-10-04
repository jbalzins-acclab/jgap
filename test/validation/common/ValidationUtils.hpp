#ifndef JGAP_VALIDATION_UTILS_HPP
#define JGAP_VALIDATION_UTILS_HPP

#include <cmath>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <iomanip>

namespace jgap::validation {

    namespace fs = std::filesystem;

    // Numerical error metrics
    double computeRmse(const std::vector<double>& actual, const std::vector<double>& ref);
    double computeNrmsePct(const std::vector<double>& actual, const std::vector<double>& ref);
    double computeRmsnePct(const std::vector<double>& actual, const std::vector<double>& ref, double min_ref = 1e-12);
    double computeMaxRelPct(const std::vector<double>& actual, const std::vector<double>& ref, double min_ref = 1e-12);
    double computeSigRelPct(const std::vector<double>& actual, const std::vector<double>& ref, double sig_fraction = 0.05);
    double computeCosineSim(const std::vector<double>& actual, const std::vector<double>& ref);

    // Repository & resource location helpers
    fs::path findRepoRoot();
    fs::path getValidationResourceDir();

    // Assertion & threshold verification
    class ValidationReporter {
    public:
        ValidationReporter(std::string suite_name) : suite_name_(std::move(suite_name)) {}

        bool assertLessOrEqual(
            const std::string& test_name,
            const std::string& metric_name,
            double actual,
            double threshold,
            const std::string& unit = ""
        );

        bool assertGreaterOrEqual(
            const std::string& test_name,
            const std::string& metric_name,
            double actual,
            double threshold,
            const std::string& unit = ""
        );

        void recordFailure(const std::string& msg);

        int summarizeAndExit() const;

    private:
        std::string suite_name_;
        size_t pass_count_{0};
        size_t fail_count_{0};
        std::vector<std::string> failure_messages_;
    };

} // namespace jgap::validation

#endif // JGAP_VALIDATION_UTILS_HPP
