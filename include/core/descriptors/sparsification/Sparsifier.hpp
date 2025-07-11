#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>
#include <map>

#include "data/BasicDataTypes.hpp"
#include "io/log/StdoutLogger.hpp"

using namespace std;

namespace jgap {

    template<class TSparsePoints, class TFromData>
    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        //virtual TSparsePoints sparsifyFromData(const vector<AtomicStructure> &fromData) = 0;
        virtual TSparsePoints sparsifyFromData(const TFromData &fromData) = 0;
    };
}

#endif //SPARSIFIER_HPP
