#ifndef JGAP_ATOMICSTRUCTURE_HPP
#define JGAP_ATOMICSTRUCTURE_HPP

#include "NeighbourData.hpp"
#include "Vector3.hpp"
#include "PredictionData.hpp"
#include "SpeciesData.hpp"

namespace jgap {

    struct AtomicStructure {
        map<string, string> properties;
        array<Vector3, 3> lattice;
        vector<Vector3> positions;
        vector<Species> species;

        optional<double> neighbourListCutoff;
        optional<vector<NeighboursData>> neighbours;

        optional<double> energy;
        optional<vector<Vector3>> forces;
        optional<array<Vector3, 3>> virials;

        optional<double> energySigmaInverse;
        optional<vector<Vector3>> forceSigmasInverse;
        optional<array<Vector3, 3>> virialSigmasInverse;

        struct AtomProxy {
            size_t index;
            AtomicStructure* structure;

            Vector3& position() const { return structure->positions[index]; }
            Species& species() const { return structure->species[index]; }

            Vector3& force() const {
                if (!structure->forces)
                    throw std::runtime_error("Forces not set");
                return (*structure->forces)[index];
            }
            Vector3& forceSigmasInverse() const {
                if (!structure->forceSigmasInverse)
                    throw std::runtime_error("Forces not set");
                return (*structure->forceSigmasInverse)[index];
            }
            NeighboursData& neighbours() const {
                if (!structure->neighbours)
                    throw std::runtime_error("Neighbours not set");
                return (*structure->neighbours)[index];
            }
        };

        struct ConstAtomProxy {
            size_t index;
            const AtomicStructure* structure;

            Vector3 position() const { return structure->positions[index]; }
            Species species() const { return structure->species[index]; }

            Vector3 force() const {
                if (!structure->forces)
                    throw std::runtime_error("Forces not set");
                return (*structure->forces)[index];
            }
            Vector3 forceSigmasInverse() const {
                if (!structure->forceSigmasInverse)
                    throw std::runtime_error("Forces not set");
                return (*structure->forceSigmasInverse)[index];
            }
            NeighboursData neighbours() const {
                if (!structure->neighbours)
                    throw std::runtime_error("Neighbours not set");
                return (*structure->neighbours)[index];
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
    };
}

#endif