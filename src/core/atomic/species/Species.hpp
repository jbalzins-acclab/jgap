#ifndef JGAP_SPECIES_HPP
#define JGAP_SPECIES_HPP

#include <array>
#include <atomic>
#include <cassert>
#include <map>
#include <mutex>
#include <optional>
#include <string>

#include "../../Real.hpp"
#include "io/log/CurrentLogger.hpp"


namespace jgap {
    class Species {
    public:
        static constexpr size_t NumberOfElements = 118;

        static std::map<std::string, uint16_t> SymbolToAtomicNumber;
        static std::array<std::string, NumberOfElements+1> AtomicNumberToSymbol;
        static std::array<Real, NumberOfElements+1> Masses;

        static Species Anon() {
            static const Species anon("AnonymousSpecies");
            return anon;
        }

        static Species fromAtomicNumber(size_t Z) {
            return Species(AtomicNumberToSymbol[Z]);
        }

        Species(const Species& other) = default;
        Species& operator=(const Species& other) = default;

        Species(const std::string &symbol) {
            std::lock_guard lock(Mtx);

            if (SymbolIds.contains(symbol)) {
                id = SymbolIds[symbol];
            } else {
                id = static_cast<uint16_t>(SymbolIds.size());
                SymbolIds[symbol] = id;
                IdSymbols[id] = symbol;

                JGAP_LOG_DEBUG("Internally mapped {} with id = {}", symbol, id);
            }
        }

        Species(const char* symbol) : Species(std::string(symbol)) {}

        Species(uint16_t id) : id(id) {
            assert(IdSymbols.contains(id) && "Unknown species ID");
        }

        std::string symbol() const {
            return IdSymbols.at(id);
        }

        uint16_t getId() const {
            return id;
        }

        std::optional<size_t> atomicNumber() const {
            if (!SymbolToAtomicNumber.contains(symbol())) {
                return std::nullopt;
            }
            return SymbolToAtomicNumber.at(symbol());
        }

        std::optional<Real> mass() const {
            return atomicNumber().transform(
                [](const size_t Z) -> double { return Masses[Z]; }
                );
        }

        bool operator<(const Species& other) const {
            return id < other.id;
        }

        bool operator==(const Species& other) const {
            return id == other.id;
        }

    private:
        inline static std::mutex Mtx{};
        inline static std::map<std::string, uint16_t> SymbolIds = {};
        inline static std::map<uint16_t, std::string> IdSymbols = {};
        // TODO: add properties per ID ??

        uint16_t id;
    };
}

#endif