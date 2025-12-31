#ifndef JGAP_ATOMICSTRUCTURE_HPP
#define JGAP_ATOMICSTRUCTURE_HPP

#include <cassert>
#include <utility>

#include "NeighbourData.hpp"
#include "../Vector3.hpp"
#include "PredictionData.hpp"
#include "SpeciesData.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    class AtomicStructure {
    public:
        std::array<Vector3, 3> lattice_vectors;
        std::vector<Vector3> positions;
        std::vector<Species> species;

        std::optional<std::vector<NeighboursData>> neighbours_ascending_separation;

        std::map<std::string, std::string> properties;

        std::optional<double> energy;
        std::optional<Virials> virials;
        std::optional<std::vector<Vector3>> forces;

        std::optional<double> energy_sigma_inverse;
        std::optional<Virials> virial_sigmas_inverse;
        std::optional<std::vector<Vector3>> force_sigmas_inverse;

        AtomicStructure(const std::array<Vector3, 3> &lattice_vectors,
                        std::vector<Vector3> positions,
                        std::vector<Species> species,
                        std::map<std::string, std::string> properties = {},
                        const std::optional<double> energy = {},
                        const std::optional<Virials> &virials = {},
                        std::optional<std::vector<Vector3>> forces = {},
                        const std::optional<double> energy_sigma = {},
                        const std::optional<Virials> &virials_sigmas = {},
                        std::optional<std::vector<Vector3>> force_sigmas = {})
            : lattice_vectors(lattice_vectors),
              positions(std::move(positions)),
              species(std::move(species)),
              properties(std::move(properties)),
              energy(energy),
              virials(virials),
              forces(std::move(forces))
        {
            assert(positions.size() == species.size() && "Positions and species array size mismatch");
            assert(!forces.has_value() || forces.value().size() == positions.size()
                && "Forces and positions array size mismatch");

            encodeSpecies();
            setRegularizationParams(energy_sigma, virials_sigmas, std::move(force_sigmas));
        }

        void setRegularizationParams(const std::optional<double> energy_sigma = {},
                                     const std::optional<Virials> &virials_sigmas = {},
                                     std::optional<std::vector<Vector3>> force_sigmas = {}) {

            energy_sigma_inverse = energy_sigma.transform([](double val) -> double { return 1.0 / val; });
            virial_sigmas_inverse = virials_sigmas.transform([](Virials val) -> Virials { return {
                1.0 / val.xx, 1.0 / val.xy, 1.0 / val.xz, 1.0 / val.yy, 1.0 / val.yz, 1.0 / val.zz
            }; });
            force_sigmas_inverse = force_sigmas.transform([](std::vector<Vector3> vals) -> std::vector<Vector3> {
                std::vector<Vector3> result(vals.size());
                for (size_t i = 0; i < vals.size(); i++) {
                    result[i] = Vector3{
                        1.0 / vals[i].x,
                        1.0 / vals[i].y,
                        1.0 / vals[i].z
                    };
                }
                return result;
            });
        }

        struct AtomProxy {
            size_t index;
            AtomicStructure* structure;

            Vector3& position() const { return structure->positions[index]; }
            Species& species() const { return structure->species[index]; }
            EncodedSpecies& speciesEncoded() const { return structure->species_encoded[index]; }

            Vector3& force() const {
                if (!structure->forces) {
                    JGAP_LOG_AND_THROW("Forces not set");
                }
                return (*structure->forces)[index];
            }
            Vector3& forceSigmasInverse() const {
                if (!structure->force_sigmas_inverse) {
                    JGAP_LOG_AND_THROW("Force regularization not set");
                }
                return (*structure->force_sigmas_inverse)[index];
            }
            NeighboursData& neighbours() const {
                if (!structure->neighbours_ascending_separation) {
                    JGAP_LOG_AND_THROW("Neighbours not set");
                }
                return (*structure->neighbours_ascending_separation)[index];
            }
        };

        struct ConstAtomProxy {
            size_t index;
            const AtomicStructure* structure;

            Vector3 position() const { return structure->positions[index]; }
            Species species() const { return structure->species[index]; }
            EncodedSpecies speciesEncoded() const { return structure->species_encoded[index]; }

            Vector3 force() const {
                if (!structure->forces) {
                    JGAP_LOG_AND_THROW("Forces not set");
                }
                return (*structure->forces)[index];
            }
            Vector3 forceSigmasInverse() const {
                if (!structure->force_sigmas_inverse) {
                    JGAP_LOG_AND_THROW("Force regularization not set");
                }
                return (*structure->force_sigmas_inverse)[index];
            }
            NeighboursData neighboursAscendingSeparation() const {
                if (!structure->neighbours_ascending_separation) {
                    JGAP_LOG_AND_THROW("Neighbours not set");
                }
                return (*structure->neighbours_ascending_separation)[index];
            }
        };

        struct Iterator {
            AtomicStructure* structure;
            size_t index;

            Iterator& operator++() { ++index; return *this; }
            bool operator!=(const Iterator& other) const { return index != other.index; }
            AtomProxy operator*() const { return AtomProxy{index, structure}; }
        };

        struct ConstIterator {
            const AtomicStructure* structure;
            size_t index;

            ConstIterator& operator++() { ++index; return *this; }
            bool operator!=(const ConstIterator& other) const { return index != other.index; }
            ConstAtomProxy operator*() const { return ConstAtomProxy{index, structure}; }
        };

        Iterator begin() { return Iterator{this, 0}; }
        Iterator end() { return Iterator{this, positions.size()}; }

        ConstIterator begin() const { return ConstIterator{this, 0}; }
        ConstIterator end() const { return ConstIterator{this, positions.size()}; }

        AtomProxy operator[](const size_t i) { return AtomProxy{i, this}; }
        ConstAtomProxy operator[](const size_t i) const { return ConstAtomProxy{i, this}; }

        void setEnergyData(const Predictions& prediction);
        void adjust(const Predictions& prediction, bool subtract, bool setEmpty);
        AtomicStructure repeat(size_t a, size_t b, size_t c);

        double volume() const;
        size_t size() const { return positions.size(); }

    private:
        std::vector<EncodedSpecies> species_encoded;
        void encodeSpecies();
    };
}

#endif