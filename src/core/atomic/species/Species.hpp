#ifndef JGAP_SPECIES_HPP
#define JGAP_SPECIES_HPP

#include <array>
#include <atomic>
#include <map>
#include <mutex>
#include <optional>
#include <string>

#include "../../Real.hpp"
#include "io/log/CurrentLogger.hpp"

#define NUMBER_OF_ELEMENTS 118

namespace jgap {
    class Species {
    public:
        static std::map<std::string, int16_t> SYMBOL_TO_ATOMIC_NUMBER;
        static std::array<std::string, NUMBER_OF_ELEMENTS+1> ATOMIC_NUMBER_TO_SYMBOL;
        static std::array<Real, NUMBER_OF_ELEMENTS+1> MASSES;

        Species(const Species& other) = default;

        Species(const std::string &symbol) {
            std::lock_guard lock(mtx_);

            if (symbol_ids_.contains(symbol)) {
                id_ = symbol_ids_[symbol];
            } else {
                id_ = static_cast<int16_t>(symbol_ids_.size());
                symbol_ids_[symbol] = id_;
                id_symbols_[id_] = symbol;

                JGAP_LOG_DEBUG("Internally mapped {} with id = {}", symbol, id_);
            }
        }

        std::string symbol() const {
            return id_symbols_.at(id_);
        }

        int16_t id() const {
            return id_;
        }

        std::optional<size_t> atomicNumber() const {
            if (!SYMBOL_TO_ATOMIC_NUMBER.contains(symbol())) {
                return std::nullopt;
            }
            return SYMBOL_TO_ATOMIC_NUMBER.at(symbol());
        }

        std::optional<Real> mass() const {
            return atomicNumber().transform(
                [](const size_t Z) -> double { return MASSES[Z]; }
                );
        }

        bool operator<(const Species& other) const {
            return id_ < other.id_;
        }

        bool operator==(const Species& other) const {
            return id_ == other.id_;
        }

    private:
        inline static std::mutex mtx_{};
        inline static std::map<std::string, int16_t> symbol_ids_ = {};
        inline static std::map<int16_t, std::string> id_symbols_ = {};
        // TODO: add properties per ID ??

        int16_t id_;
    };
}

#endif