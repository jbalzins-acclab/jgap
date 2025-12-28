#ifndef JGAP_DATAWITHMETA_HPP
#define JGAP_DATAWITHMETA_HPP

#include "core/fit/Fit.hpp"

namespace jgap {
    struct DataSetWithMeta {
        DataSet data;
        bool saved;
        bool to_be_saved;
    };
}

#endif