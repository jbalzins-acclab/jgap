#ifndef JGAP_MULTIAGGREGATOR_HPP
#define JGAP_MULTIAGGREGATOR_HPP
#include "TransformationAggregator.hpp"

namespace jgap {
    template<size_t Dim>
    class MultiAggregator final : public TransformationAggregator<Dim> {
    public:
        MultiAggregator(std::vector<std::unique_ptr<TransformationAggregator<Dim>>> aggregators)
            : aggregators(std::move(aggregators)) {}

        std::map<size_t, ManyBodyDescriptor<Dim>> aggregate(const NeighbourList& nl) const override {
            std::map<size_t, ManyBodyDescriptor<Dim>> final_descriptors;

            for (const auto& aggregator : aggregators) {
                auto partial_descriptors = aggregator->aggregate(nl);
                for (auto& [atom_idx, descriptor] : partial_descriptors) {
                    if (!final_descriptors.contains(atom_idx)) {
                        final_descriptors[atom_idx] = ManyBodyDescriptor<Dim>(nl.nAtoms());
                    }
                    final_descriptors.at(atom_idx) += descriptor;
                }
            }
            return final_descriptors;
        }

        Cutoffs getCutoffs() const override {
            Cutoffs combined;
            for (const auto& aggregator : aggregators) {
                combined += aggregator->getCutoffs();
            }
            return combined;
        }

    private:
        std::vector<std::unique_ptr<TransformationAggregator<Dim>>> aggregators;
    };
}

#endif
