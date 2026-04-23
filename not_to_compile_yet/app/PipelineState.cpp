#include "PipelineState.hpp"

namespace jgap {
    void PipelineState::addData(std::string label, DataSetWithMeta data) {
        if (data_sets.contains(label) || potentials.contains(label)) {
            JGAP_LOG_AND_THROW("Label '{}' already exists in data sets or potentials.", label);
        }
        data_sets[label] = std::move(data);
    }

    void PipelineState::changeData(std::string label, DataSetWithMeta data) {
        if (!data_sets.contains(label)) {
            JGAP_LOG_AND_THROW("Data set with label '{}' not found for changing.", label);
        }
        data_sets[label] = std::move(data);
    }

    void PipelineState::addPotential(std::string label, PotentialWithMeta potential) {
        if (data_sets.contains(label) || potentials.contains(label)) {
            JGAP_LOG_AND_THROW("Label '{}' already exists in data sets or potentials.", label);
        }
        potentials[label] = std::move(potential);
    }

    DataSetWithMeta & PipelineState::getDataSet(const std::string &label) {
        auto it = data_sets.find(label);
        if (it == data_sets.end()) {
            JGAP_LOG_AND_THROW("Data set with label '{}' not found.", label);
        }
        return it->second;
    }

    PotentialWithMeta & PipelineState::getPotential(const std::string &label) {
        auto it = potentials.find(label);
        if (it == potentials.end()) {
            JGAP_LOG_AND_THROW("Potential with label '{}' not found.", label);
        }
        return it->second;
    }
}
