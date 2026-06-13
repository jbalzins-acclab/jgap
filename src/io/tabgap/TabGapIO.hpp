#ifndef JGAP_TABGAPIO_HPP
#define JGAP_TABGAPIO_HPP

#include <memory>
#include <optional>
#include <vector>
#include <highfive/H5File.hpp>

#include "core/potentials/tabgap/TabGapPotential.hpp"
#include "core/potentials/tabgap/components/TabGapComponent.hpp"
#include "core/potentials/tabgap/components/ThreeBodyTGComponent.hpp"
#include "core/potentials/tabgap/components/TwoBodyTGComponent.hpp"
#include "utils/ValuePtr.hpp"

namespace jgap {

    using Filenames = std::vector<std::string>;

    class TabGapIO {
    public:
        static TabGapPotential read(const Filenames& filenames);
        static Filenames write(const TabGapPotential& potential, const std::string &output_filename_prefix);

    private:
        static void write2b(HighFive::File& h5_file, const TwoBodyTGComponent& component);
        static void write3b(HighFive::File& h5_file, const ThreeBodyTGComponent& component);

        static std::string useSomeComponentsAndGenerateEamFs(const std::vector<Species>& all_species,
                                      std::map<SpeciesSet<2, Symmetric>, const TwoBodyTGComponent*>& pair_pots,
                                      std::multimap<Species, const EamTGComponent*>& eam_components);

        static void readH5(const std::string& filename, TabGapPotential& pot);
        static void readEamFs(const std::string& filename, TabGapPotential& pot);
    };
}

#endif