#ifndef JGAP_SERIALIZABLE_HPP
#define JGAP_SERIALIZABLE_HPP

#include <string>

#include "../core/DataNode.hpp"

namespace jgap {
    class Serializable {
    public:
        virtual ~Serializable() = default;
        virtual std::string getTypeId() = 0;
        virtual DataNode serialize() = 0;
    };
}

#endif
