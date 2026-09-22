#ifndef JGAP_REGULARIZATIONRULES_HPP
#define JGAP_REGULARIZATIONRULES_HPP

#include "Regularization.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include <memory>
#include <vector>

namespace jgap {
    class RegularizationRules {
    public:
        virtual ~RegularizationRules() = default;

        virtual Regularization determine(const Atoms& atoms) const = 0;

        virtual std::vector<Regularization> determineForAll(const std::vector<Atoms>& structures) const {
            std::vector<Regularization> result;
            result.reserve(structures.size());
            for (const auto& atoms: structures) {
                result.push_back(determine(atoms));
            }
            return result;
        }

        virtual RegularizationRules* clone() const = 0;
    };
}

#endif