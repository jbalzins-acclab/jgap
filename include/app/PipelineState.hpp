#ifndef JGAP_STATE_HPP
#define JGAP_STATE_HPP

#include <map>
#include <string>
#include <memory>

#include "core/fit/Fit.hpp"
#include "data/DataSetWithMeta.hpp"
#include "data/PotentialWithMeta.hpp"

namespace jgap {
    struct PipelineState {
        std::map<std::string, DataSetWithMeta> data_sets;
        std::map<std::string, PotentialWithMeta> potentials;
    };
}

#endif
