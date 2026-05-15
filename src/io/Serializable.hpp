#ifndef JGAP_SERIALIZABLE_HPP
#define JGAP_SERIALIZABLE_HPP

#include <string>

//#include "../core/DataNode.hpp"

namespace jgap {
    class Serializable {
    public:
        static constexpr std::string TYPE_FIELD = "type";

        virtual ~Serializable() = default;

        virtual std::string getTypeId() = 0;

        /*
        virtual DataNode serialize() {
            DataNode result = serializeWithoutType();

            if (result.contains(TYPE_FIELD)) {
                JGAP_LOG_AND_THROW("\"type\" incorrectly specified in serializeWithoutType for {}", getTypeId());
            }

            result[TYPE_FIELD] = getTypeId();

            return result;
        }
        virtual void deserialize(const DataNode& serialized) {
            if (!serialized.contains(TYPE_FIELD)) {
                JGAP_LOG_AND_THROW("\"type\" not specified in {}", serialized.toString());
            }

            if (serialized[TYPE_FIELD].asString() != getTypeId()) {
                JGAP_LOG_AND_THROW("Attempted to deserialize \"{}\" as {}",
                    serialized[TYPE_FIELD].asString(), getTypeId());
            }

            // Otherwise ok
            deserializeNoTypeCheck(serialized);
        }*/

    protected:
        //virtual DataNode serializeWithoutType() = 0;
        //virtual void deserializeNoTypeCheck(const DataNode& serialized) = 0;
    };
}

#endif
