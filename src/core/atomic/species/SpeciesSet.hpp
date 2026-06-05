#ifndef JGAP_SPECIESSET_HPP
#define JGAP_SPECIESSET_HPP

#include "Species.hpp"
#include <array>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <initializer_list> // Required for std::initializer_list

namespace jgap {
    class SpeciesSet {
    public:
        static constexpr int16_t None = -1;
#ifdef MAX_SPECIES_IN_SET
        static constexpr size_t MaxSpecies = MAX_SPECIES_IN_SET;
#else
        static constexpr size_t MaxSpecies = 4;
#endif

        SpeciesSet(const SpeciesSet& other) = default;

        SpeciesSet(const Species& root) {
            representation.fill(None);
            representation[0] = root.getId();
        }

        SpeciesSet& node(const Species& node_species) {
            // insertion sort, expecting MaxSpecies to be no more than 8
            int16_t id_to_insert = node_species.getId();
            for (size_t i = 1; i < MaxSpecies; i++) {
                if (id_to_insert > representation[i]) {
                    std::swap(representation[i], id_to_insert);
                }
            }

            return *this;
        }

        SpeciesSet(const std::vector<std::string>& symbols) {
            representation.fill(None);
            if (symbols.size() > MaxSpecies) {
                JGAP_LOG_AND_THROW("More than MaxSpecies={} species in {}, consider re-compiling",
                                   MaxSpecies, symbols.size());
            }
            if (symbols.empty()) return;

            representation[0] = Species(symbols[0]).getId();
            for (size_t i = 1; i < symbols.size(); i++) {
                 node(Species(symbols[i]));
            }
        }

        // New constructor for initializer_list
        SpeciesSet(const std::initializer_list<std::string>& symbols) {
            representation.fill(None);
            if (symbols.size() > MaxSpecies) {
                JGAP_LOG_AND_THROW("More than MaxSpecies={} species in {}, consider re-compiling",
                                   MaxSpecies, symbols.size());
            }
            if (symbols.size() == 0) return;

            auto it = symbols.begin();
            representation[0] = Species(*it).getId();
            ++it;
            for (; it != symbols.end(); ++it) {
                 node(Species(*it));
            }
        }

        bool operator==(const SpeciesSet& other) const {
            return representation == other.representation;
        }

        bool operator<(const SpeciesSet &other) const {
            return representation < other.representation;
        }

        int16_t operator[](size_t index) const {
            return representation[index];
        }

        Species getRoot() const {
            return Species(representation[0]);
        }

        SpeciesSet rootOnly() const {
            return SpeciesSet(getRoot());
        }

        bool isRedundantTwoBody() const {
            if (representation[0] == None || representation[1] == None) return false;
            if (representation[0] <= representation[1]) return false;
            for (size_t i = 2; i < MaxSpecies; i++) {
                if (representation[i] != None) return false;
            }
            return true;
        }

        size_t clusterSize() const {
            size_t size = 0;
            for (size_t i = 0; i < MaxSpecies; i++) {
                if (representation[i] != None) {
                    size++;
                } else {
                    break;
                }
            }
            return size;
        }

        std::string toString() const {
            std::stringstream ss;
            bool first = true;
            for (size_t i = 0; i < MaxSpecies; i++) {
                if (representation[i] != None) {
                    if (!first) {
                        ss << ",";
                    }
                    ss << representation[i];
                    first = false;
                }
            }
            return ss.str();
        }

    private:
        std::array<int16_t, MaxSpecies> representation{};
    };
}

#endif