#ifndef JGAP_PIPELINESTATE_HPP
#define JGAP_PIPELINESTATE_HPP

#include <map>
#include <string>
#include <memory>
#include <stdexcept> 

#include "core/fit/Fit.hpp"
#include "data/DataSetWithMeta.hpp"
#include "data/PotentialWithMeta.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    class PipelineState {
    public:
        PipelineState() : data_sets({}), potentials({}), fits({}) {}
        PipelineState(const PipelineState& other) = default;

        void addData(std::string label, DataSetWithMeta data);
        void changeData(std::string label, DataSetWithMeta data);
        DataSetWithMeta& getDataSet(const std::string& label);

        void addPotential(std::string label, PotentialWithMeta potential);
        PotentialWithMeta& getPotential(const std::string& label);

    private:
        std::map<std::string, DataSetWithMeta> data_sets;
        std::map<std::string, PotentialWithMeta> potentials;
        std::map<std::string, std::shared_ptr<Fit>> fits;
    };
}

#endif
