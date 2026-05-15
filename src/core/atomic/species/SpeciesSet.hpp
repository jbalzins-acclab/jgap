#ifndef JGAP_SPECIESSET_HPP
#define JGAP_SPECIESSET_HPP

#include "Species.hpp"
#include <array>
#include <algorithm>
#include <string>
#include <sstream>

namespace jgap {
    class SpeciesSet {
    public:
        static constexpr int16_t None = -1;
#ifdef MAX_SPECIES_IN_SET
        static constexpr size_t MAX_SPECIES = MAX_SPECIES_IN_SET;
#else
        static constexpr size_t MAX_SPECIES = 4;
#endif

        SpeciesSet(const SpeciesSet& other) = default;

        SpeciesSet(const Species& root) {
            representation_.fill(None);
            representation_[0] = root.id();
        }

        SpeciesSet& node(const Species& node_species) {
            // insertion sort, expecting MAX_SPECIES to be no more than 8
            int16_t id_to_insert = node_species.id();
            for (size_t i = 1; i < MAX_SPECIES; i++) {
                if (id_to_insert > representation_[i]) {
                    std::swap(representation_[i], id_to_insert);
                }
            }

            return *this;
        }

        bool operator==(const SpeciesSet& other) const {
            return representation_ == other.representation_;
        }

        bool operator<(const SpeciesSet &other) const {
            return representation_ < other.representation_;
        }

        int16_t operator[](size_t index) const {
            return representation_[index];
        }

        std::string toString() const {
            std::stringstream ss;
            bool first = true;
            for (size_t i = 0; i < MAX_SPECIES; i++) {
                if (representation_[i] != None) {
                    if (!first) {
                        ss << ",";
                    }
                    ss << representation_[i];
                    first = false;
                }
            }
            return ss.str();
        }

    private:
        std::array<int16_t, MAX_SPECIES> representation_{};
    };
}

#endif
