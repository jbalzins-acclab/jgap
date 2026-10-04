#ifndef JGAP_VALIDATION_CONFIG_HPP
#define JGAP_VALIDATION_CONFIG_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace jgap::validation {

    struct QrValidationOptions {
        size_t min_m3b{50};
        size_t max_m3b{100};
        size_t step_m3b{50};
        std::vector<size_t> sparse_seeds{42};
        bool quick{false};
        std::string output_csv{""};

        // Failure thresholds
        double max_allowed_nrmse_pct{1e-4};      // 0.0001%
        double max_allowed_sigrel_pct{0.05};     // 0.05%
        double min_allowed_cosine_sim{0.999999}; // 0.999999

        static QrValidationOptions parse(int argc, char* argv[]) {
            QrValidationOptions opt;
            for (int i = 1; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "--quick") {
                    opt.quick = true;
                    opt.min_m3b = 100;
                    opt.max_m3b = 100;
                    opt.sparse_seeds = {42};
                } else if (arg == "--max-m3b" && i + 1 < argc) {
                    opt.max_m3b = std::stoull(argv[++i]);
                } else if (arg == "--min-m3b" && i + 1 < argc) {
                    opt.min_m3b = std::stoull(argv[++i]);
                } else if (arg == "--step-m3b" && i + 1 < argc) {
                    opt.step_m3b = std::stoull(argv[++i]);
                } else if (arg == "--seeds" && i + 1 < argc) {
                    opt.sparse_seeds.clear();
                    std::stringstream ss(argv[++i]);
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        if (!token.empty()) {
                            opt.sparse_seeds.push_back(std::stoull(token));
                        }
                    }
                } else if (arg == "--csv" && i + 1 < argc) {
                    opt.output_csv = argv[++i];
                } else if (arg == "--threshold-nrmse" && i + 1 < argc) {
                    opt.max_allowed_nrmse_pct = std::stod(argv[++i]);
                } else if (arg == "-h" || arg == "--help") {
                    std::cout << "Usage: ValidateQrVariants [options]\n"
                              << "Options:\n"
                              << "  --quick               Run quick smoke test (N_3b=100, seed=42)\n"
                              << "  --min-m3b <N>         Minimum 3-body sparse points (default: 50)\n"
                              << "  --max-m3b <N>         Maximum 3-body sparse points (default: 100)\n"
                              << "  --step-m3b <N>        Step size for 3-body sparse points (default: 50)\n"
                              << "  --seeds <s1,s2,...>   Comma-separated list of random seeds\n"
                              << "  --csv <filepath>      Export detailed run metrics to CSV\n"
                              << "  --threshold-nrmse <P> Set maximum allowed NRMSE % (default: 1e-4)\n"
                              << "  -h, --help            Print help\n";
                    std::exit(0);
                }
            }
            return opt;
        }
    };

} // namespace jgap::validation

#endif // JGAP_VALIDATION_CONFIG_HPP
