#ifndef JGAP_TABGAPIO_HPP
#define JGAP_TABGAPIO_HPP

#include "core/potentials/TabGapPotential.hpp"

namespace jgap {

    using FileNames = std::vector<std::string>;

    class TabGapIO {
    public:
        static TabulationData read(const FileNames& fileNames);
        static FileNames write(const TabulationData& valuesTables,
                               const TabulationData& splineTables,
                               std::optional<std::string> outputFileNamePrefix);

    private:
        static std::string generateFileNamePrefix(const TabulationData& coeffs);

        static std::string writeH5(const TabulationData &valuesTables,
                              const TabulationData &splineTables,
                              const std::string &outputFileNamePrefix);
        static std::string writeEamFs(const TabulationData& valueTables, size_t index, const std::string &outputFileNamePrefix);

        static void readH5(const std::string& fileName, TabulationData& splineCoefficients);
        static void readEamFs(const std::string& fileName, TabulationData& splineCoefficients);
    };
}

#endif