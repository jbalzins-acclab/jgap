#include "SerializationRegistry.hpp"

namespace jgap::detail {
    std::map<std::type_index, std::shared_ptr<void>>& serializationRegistryStorage() {
        static std::map<std::type_index, std::shared_ptr<void>> storage;
        return storage;
    }
}
