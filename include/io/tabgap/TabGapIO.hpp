#ifndef JGAP_TABGAPIO_HPP
#define JGAP_TABGAPIO_HPP

#include "core/potentials/TabGapPotential.hpp"

namespace jgap {

    using FileNames = vector<string>;

    class TabGapIO {
    public:
        static TabulationData read(const FileNames& fileNames);
        static FileNames write(const TabulationData& valuesTables,
                               const TabulationData& splineTables,
                               optional<string> outputFileNamePrefix);

    private:
        static string generateFileNamePrefix(const TabulationData& coeffs);

        static string writeH5(const TabulationData &valuesTables,
                              const TabulationData &splineTables,
                              const string &outputFileNamePrefix);
        static string writeEamFs(const TabulationData& valueTables, size_t index, const string &outputFileNamePrefix);

        static void readH5(const string& fileName, TabulationData& splineCoefficients);
        static void readEamFs(const string& fileName, TabulationData& splineCoefficients);
    };
}

#endif