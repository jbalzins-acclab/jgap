#include "PotentialLoader.hpp"

#include <filesystem>
#include <stdexcept>
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    ValuePtr<Potential> loadPotential(const std::vector<std::string>& pot_paths) {
        if (pot_paths.empty()) {
            throw std::invalid_argument("Empty paths list for loading potential");
        }
        if (pot_paths.size() == 1) {
            return loadPotential(pot_paths[0]);
        }
        return TabGapIO::read(pot_paths);
    }

    ValuePtr<Potential> loadPotential(const std::string& pot_path) {
        namespace fs = std::filesystem;
        fs::path p(pot_path);
        const bool has_extension = p.has_extension();

        if (!has_extension) {
            std::vector<std::string> tabgap_files;
            const std::string tabgap_h5 = pot_path + ".tabgap.h5";
            const std::string eam_fs = pot_path + ".eam.fs";

            if (fs::exists(tabgap_h5)) {
                tabgap_files.push_back(tabgap_h5);
            }
            if (fs::exists(eam_fs)) {
                tabgap_files.push_back(eam_fs);
            }

            if (!fs::exists(eam_fs)) {
                const fs::path dir = p.parent_path().empty() ? "." : p.parent_path();
                if (fs::exists(dir) && fs::is_directory(dir)) {
                    const std::string stem = p.filename().string();
                    for (const auto& entry: fs::directory_iterator(dir)) {
                        const std::string fname = entry.path().filename().string();
                        if (fname.starts_with(stem) && fname.ends_with(".eam.fs") && fname != eam_fs) {
                            tabgap_files.push_back(entry.path().string());
                        }
                    }
                }
            }

            if (!tabgap_files.empty()) {
                JGAP_LOG_INFO("Loading tabGAP potential from files: {}", utils::vectorToString(tabgap_files));
                return TabGapIO::read(tabgap_files);
            }

            if (fs::exists(pot_path + ".jgap.h5")) {
                return SerializationRegistry<Potential>::deserialize(pot_path + ".jgap.h5");
            }
            if (fs::exists(pot_path + ".h5")) {
                return SerializationRegistry<Potential>::deserialize(pot_path + ".h5");
            }
        }

        if (pot_path.ends_with(".tabgap.h5") || pot_path.ends_with(".eam.fs")) {
            return TabGapIO::read(std::vector<std::string>{pot_path});
        }

        return SerializationRegistry<Potential>::deserialize(pot_path);
    }
}
